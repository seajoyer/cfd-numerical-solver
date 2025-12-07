#include "spatial/GodunovSpatialOperator.hpp"

#include <algorithm>
#include <cctype>
#include <string>

#include "data/DataLayer.hpp"
#include "reconstruction/ENOReconstruction.hpp"
#include "reconstruction/P0Reconstruction.hpp"
#include "reconstruction/P1Reconstruction.hpp"
#include "reconstruction/WENOReconstruction.hpp"
#include "riemann/AcousticRiemannSolver.hpp"
#include "riemann/ExactIdealGasRiemannSolver.hpp"
#include "riemann/HLLCRiemannSolver.hpp"
#include "riemann/HLLRiemannSolver.hpp"
#include "viscosity/VNRArtificialViscosity.hpp"

GodunovSpatialOperator::GodunovSpatialOperator(const Settings& settings) {
    InitializeReconstruction(settings);
    InitializeRiemannSolver(settings);
    viscosity_ = std::make_unique<VNRArtificialViscosity>(settings);
}

void GodunovSpatialOperator::InitializeReconstruction(const Settings& settings) {
    std::string name = settings.reconstruction;

    if (name.find("p1") != std::string::npos) {
        reconstruction_ = std::make_shared<P1Reconstruction>();
    } else if (name.find("p0") != std::string::npos) {
        reconstruction_ = std::make_shared<P0Reconstruction>();
    } else if (name.starts_with("eno")) {
        int order = 3;
        try {
            order = std::stoi(name.substr(3, std::string::npos));
        } catch (...) {
            std::cout << "Order of ENO not found. Set order to 3" << std::endl;
        }
        reconstruction_ = std::make_shared<ENOReconstruction>(order);
    } else if (name.starts_with("weno")) {
        int order = 5;
        try {
            order = std::stoi(name.substr(4, std::string::npos));
        } catch (...) {
            std::cout << "Order of WENO not found. Set order to 5" << std::endl;
        }
        if (order != 3 && order != 5) {
            std::cout << "WENO supports only orders 3 or 5 for now. Set order to 5"
                << std::endl;
        }
        reconstruction_ = std::make_shared<WENOReconstruction>(order);
    } else {
        reconstruction_ = std::make_shared<P0Reconstruction>();
    }
}

void GodunovSpatialOperator::InitializeRiemannSolver(const Settings& settings) {
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
    } else if (lower.find("acoustic") != std::string::npos) {
        riemann_solver_ = std::make_shared<AcousticRiemannSolver>();
    } else {
        riemann_solver_ = std::make_shared<HLLRiemannSolver>();
    }
}

void GodunovSpatialOperator::ComputeRHS(const DataLayer& layer,
                                        double dx,
                                        double gamma,
                                        xt::xarray<Conservative>& rhs) const {

    const int total_size = layer.GetTotalSize();
    if (total_size < 3 || dx <= 0.0) {
        rhs = xt::xarray<Conservative>::from_shape(
            {static_cast<std::size_t>(std::max(total_size, 0))});
        for (int j = 0; j < total_size; ++j) {
            rhs(j) = Conservative{};
        }
        return;
    }

    const int core_start = layer.GetCoreStart(0);
    const int core_end = layer.GetCoreEndExclusive(0);
    const int n_core = core_end - core_start;
    if (n_core < 2) {
        rhs = xt::xarray<Conservative>::from_shape(
            {static_cast<std::size_t>(total_size)});
        for (int j = 0; j < total_size; ++j) {
            rhs(j) = Conservative{};
        }
        return;
    }

    const int n_interfaces = total_size - 1;

    rhs = xt::xarray<Conservative>::from_shape(
        {static_cast<std::size_t>(total_size)});

    xt::xarray<Primitive> left_states;
    xt::xarray<Primitive> right_states;
    reconstruction_->ReconstructStates(layer, left_states, right_states);

    xt::xarray<Flux> fluxes =
        xt::xarray<Flux>::from_shape({static_cast<std::size_t>(n_interfaces)});

    xt::xarray<double> q;
    if (viscosity_) {
        viscosity_->ComputeInterfaceQ(layer, dx, q);
    } else {
        q = xt::zeros<double>({static_cast<std::size_t>(n_interfaces)});
    }

    for (int i = 0; i < n_interfaces; ++i) {
        Primitive WL = left_states(i);
        Primitive WR = right_states(i);

        const double qi = q(i);
        if (qi > 0.0) {
            WL.P += qi;
            WR.P += qi;
        }

        fluxes(i) = riemann_solver_->ComputeFlux(WL, WR, gamma);
    }

    for (int j = 0; j < total_size; ++j) {
        rhs(j) = Conservative{};
    }

    for (int j = core_start; j < core_end; ++j) {

        const Flux& Fminus = fluxes(j - 1); // F_{j-1/2}
        const Flux& Fplus = fluxes(j); // F_{j+1/2}

        const Flux dF = Flux::Diff(Fplus, Fminus);

        rhs(j) = dF * (-1 / dx);
    }
}