#include "spatial/GodunovKolganRodionovSpatialOperator.hpp"

#include <algorithm>
#include <cctype>
#include <string>

#include "data/DataLayer.hpp"
#include "data/Variables.hpp"
#include "reconstruction/ENOReconstruction.hpp"
#include "reconstruction/P1Reconstruction.hpp"
#include "reconstruction/WENOReconstruction.hpp"
#include "riemann/ExactIdealGasRiemannSolver.hpp"
#include "riemann/HLLCRiemannSolver.hpp"
#include "riemann/HLLRiemannSolver.hpp"
#include "solver/EOS.hpp"

GodunovKolganRodionovSpatialOperator::GodunovKolganRodionovSpatialOperator(
    const Settings& settings)
    : dt_local_(0.0) {
    InitializeReconstruction(settings);
    InitializeRiemannSolver(settings);
}

void GodunovKolganRodionovSpatialOperator::InitializeReconstruction(
    const Settings& settings) {
    std::string name = settings.reconstruction;

    if (name.find("p0") != std::string::npos) {
        std::cout
            << "Godunov-Kolgan-Rodionov method requires at least P1 reconstruction.\n"
            << "Setting reconstruction P1." << '\n';
        name = "p1";
    }

    if (name.find("p1") != std::string::npos) {
        reconstruction_ = std::make_shared<P1Reconstruction>();
        return;
    }

    if (name.starts_with("eno")) {
        int order = 3;
        try {
            order = std::stoi(name.substr(3, std::string::npos));
        } catch (...) {
            std::cout << "Order of ENO not found. Set order to 3" << '\n';
        }
        reconstruction_ = std::make_shared<ENOReconstruction>(order);
        return;
    }

    if (name.starts_with("weno")) {
        int order = 5;
        try {
            order = std::stoi(name.substr(4, std::string::npos));
        } catch (...) {
            std::cout << "Order of WENO not found. Set order to 5" << '\n';
        }
        if (order != 3 && order != 5) {
            std::cout << "WENO supports only orders 3 or 5 for now. Set order to 5"
                << '\n';
        }
        reconstruction_ = std::make_shared<WENOReconstruction>(order);
        return;
    }

    std::cout << "Reconstruction type is unknown.\n"
        << "Setting reconstruction P1." << '\n';
    reconstruction_ = std::make_shared<P1Reconstruction>();
}

void GodunovKolganRodionovSpatialOperator::InitializeRiemannSolver(
    const Settings& settings) {
    std::string name = settings.riemann_solver;
    std::string lower(name.size(), '\0');

    std::transform(name.begin(), name.end(), lower.begin(),
                   [](unsigned char c) -> char {
                       return static_cast<char>(std::tolower(c));
                   });

    if (lower.find("hllc") != std::string::npos) {
        riemann_solver_ = std::make_shared<HLLCRiemannSolver>();
    } else if (lower.find("hll") != std::string::npos) {
        riemann_solver_ = std::make_shared<HLLRiemannSolver>();
    } else if (lower.find("exact") != std::string::npos) {
        riemann_solver_ =
            std::make_shared<ExactIdealGasRiemannSolver>(0.0, settings.Q_user);
    } else {
        riemann_solver_ = std::make_shared<HLLRiemannSolver>();
    }
}

void GodunovKolganRodionovSpatialOperator::SetLocalTimeStep(double dt) {
    dt_local_ = dt;
}

void GodunovKolganRodionovSpatialOperator::ComputeRHS(const DataLayer& layer,
                                                      double dx,
                                                      double gamma,
                                                      xt::xarray<Conservative>& rhs) const {
    const int total_size = layer.GetTotalSize();
    const int core_start = layer.GetCoreStart(0);
    const int core_end = layer.GetCoreEndExclusive(0);
    const int n_core = core_end - core_start;

    if (total_size < 3 || n_core < 2 || dx <= 0.0 || dt_local_ < 0.0) {
        // Degenerate grid, invalid dx, or dt not set: zero RHS
        rhs = xt::xarray<Conservative>::from_shape(
            {static_cast<std::size_t>(std::max(total_size, 0))});
        for (int j = 0; j < total_size; ++j) {
            rhs(j) = Conservative{};
        }
        return;
    }

    const int n_interfaces = total_size - 1;

    rhs = xt::xarray<Conservative>::from_shape(
        {static_cast<std::size_t>(total_size)});

    const double dt_over_dx = dt_local_ / dx;
    const double half_dt_over_dx = 0.5 * dt_over_dx;

    // Interface primitives from reconstruction
    xt::xarray<Primitive> WL_interface;
    xt::xarray<Primitive> WR_interface;
    reconstruction_->ReconstructStates(layer, WL_interface, WR_interface);

    // Cell-edge states and conservative copies
    xt::xarray<Primitive> W_L_cell =
        xt::xarray<Primitive>::from_shape({static_cast<std::size_t>(total_size)});
    xt::xarray<Primitive> W_R_cell =
        xt::xarray<Primitive>::from_shape({static_cast<std::size_t>(total_size)});

    xt::xarray<Conservative> U_L =
        xt::xarray<Conservative>::from_shape({static_cast<std::size_t>(total_size)});
    xt::xarray<Conservative> U_R =
        xt::xarray<Conservative>::from_shape({static_cast<std::size_t>(total_size)});
    xt::xarray<Conservative> U_L_star =
        xt::xarray<Conservative>::from_shape({static_cast<std::size_t>(total_size)});
    xt::xarray<Conservative> U_R_star =
        xt::xarray<Conservative>::from_shape({static_cast<std::size_t>(total_size)});

    xt::xarray<Flux> F_L =
        xt::xarray<Flux>::from_shape({static_cast<std::size_t>(total_size)});
    xt::xarray<Flux> F_R =
        xt::xarray<Flux>::from_shape({static_cast<std::size_t>(total_size)});

    xt::xarray<Flux> fluxes =
        xt::xarray<Flux>::from_shape({static_cast<std::size_t>(n_interfaces)});

    auto clamp_interface = [n_interfaces](int idx) {
        if (idx < 0) {
            return 0;
        }
        if (idx >= n_interfaces) {
            return n_interfaces - 1;
        }
        return idx;
    };

    // Build left/right cell-edge states & corresponding fluxes
    for (int j = 0; j < total_size; ++j) {
        const int left_interface_index = clamp_interface(j - 1);
        const int right_interface_index = clamp_interface(j);

        W_L_cell(j) = WR_interface(left_interface_index);
        W_R_cell(j) = WL_interface(right_interface_index);

        U_L(j) = EOS::PrimToCons(W_L_cell(j), gamma);
        U_R(j) = EOS::PrimToCons(W_R_cell(j), gamma);

        F_L(j) = EulerFlux(W_L_cell(j), gamma);
        F_R(j) = EulerFlux(W_R_cell(j), gamma);
    }

    for (int j = 0; j < total_size; ++j) {
        const Flux dF = Flux::Diff(F_R(j), F_L(j));

        U_L_star(j) = U_L(j);
        U_R_star(j) = U_R(j);

        U_L_star(j) -= half_dt_over_dx * dF;
        U_R_star(j) -= half_dt_over_dx * dF;
    }

    for (int i = 0; i < n_interfaces; ++i) {
        const Primitive WL_star = EOS::ConsToPrim(U_R_star(i), gamma);
        const Primitive WR_star = EOS::ConsToPrim(U_L_star(i + 1), gamma);
        fluxes(i) = riemann_solver_->ComputeFlux(WL_star, WR_star, gamma);
    }

    for (int j = 0; j < total_size; ++j) {
        rhs(j) = Conservative{};
    }

    for (int j = core_start; j < core_end; ++j) {
        const Flux& Fminus = fluxes(j - 1); // F_{j-1/2}
        const Flux& Fplus = fluxes(j); // F_{j+1/2}

        const Flux dF = Flux::Diff(Fplus, Fminus);

        Conservative Lj;
        Lj.rho = -dF.mass / dx;
        Lj.rhoU = -dF.momentum / dx;
        Lj.E = -dF.energy / dx;

        rhs(j) = Lj;
    }
}