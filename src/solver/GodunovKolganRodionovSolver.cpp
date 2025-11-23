#include "solver/GodunovKolganRodionovSolver.hpp"

#include <cmath>

#include "reconstruction/P1Reconstruction.hpp"
#include "riemann/ExactIdealGasRiemannSolver.hpp"
#include "riemann/HLLCRiemannSolver.hpp"
#include "riemann/HLLRiemannSolver.hpp"
#include "solver/EOS.hpp"
#include "solver/PositivityLimiter.hpp"
#include "solver/TimeStepCalculator.hpp"

GodunovKolganRodionovSolver::GodunovKolganRodionovSolver(const Settings& settings)
    : settings_(settings),
      boundary_manager_(settings.dim),
      rho_min_(1e-10),
      p_min_(1e-10) {
    cfl_ = settings_.cfl;
    InitializeReconstruction();
    InitializeRiemannSolver();
}

void GodunovKolganRodionovSolver::AddBoundary(
    int axis, std::shared_ptr<BoundaryCondition> left_bc,
    std::shared_ptr<BoundaryCondition> right_bc) {
    boundary_manager_.Set(axis, std::move(left_bc), std::move(right_bc));
}

void GodunovKolganRodionovSolver::InitializeReconstruction() {
    reconstruction_ = std::make_shared<P1Reconstruction>();
}

void GodunovKolganRodionovSolver::InitializeRiemannSolver() {
    std::string name = settings_.riemann_solver;
    std::string lower(name.size(), '\0');

    std::transform(name.begin(), name.end(), lower.begin(), [](unsigned char c) -> char {
        return static_cast<char>(std::tolower(c));
    });

    if (lower.find("hllc") != std::string::npos) {
        riemann_solver_ = std::make_shared<HLLCRiemannSolver>();
    } else if (lower.find("hll") != std::string::npos) {
        riemann_solver_ = std::make_shared<HLLRiemannSolver>();
    } else if (lower.find("exact") != std::string::npos) {
        riemann_solver_ =
            std::make_shared<ExactIdealGasRiemannSolver>(0., settings_.Q_user);
    } else {
        riemann_solver_ = std::make_shared<HLLRiemannSolver>();
    }
}

auto GodunovKolganRodionovSolver::ComputeDx(const DataLayer& layer) const -> double {
    if (settings_.N > 0 && settings_.L_x > 0.0) {
        return settings_.L_x / static_cast<double>(settings_.N);
    }

    const int core_start = layer.GetCoreStart(0);
    const int core_end = layer.GetCoreEndExclusive(0);

    if (core_end - core_start > 1) {
        return layer.xc(core_start + 1) - layer.xc(core_start);
    }

    return 1.0;
}

auto GodunovKolganRodionovSolver::Step(DataLayer& layer, double& t_cur) -> double {
    boundary_manager_.ApplyAll(layer);

    const double dx = ComputeDx(layer);
    const int total_size = layer.GetTotalSize();
    const int core_start = layer.GetCoreStart(0);
    const int core_end = layer.GetCoreEndExclusive(0);
    const int n_core = core_end - core_start;

    if (n_core < 2 || total_size < 3) {
        return 0.0;
    }

    double dt = TimeStepCalculator::ComputeDt(layer, dx, cfl_, settings_.gamma);
    if (dt <= 0.0) {
        return 0.0;
    }

    if (t_cur + dt > settings_.t_end) {
        dt = settings_.t_end - t_cur;
        if (dt <= 0.0) {
            return 0.0;
        }
    }

    const double dt_over_dx = dt / dx;
    const double half_dt_over_dx = 0.5 * dt_over_dx;
    const double gamma = settings_.gamma;
    const int n_interfaces = total_size - 1;

    xt::xarray<Primitive> WL_interface;
    xt::xarray<Primitive> WR_interface;
    reconstruction_->ReconstructStates(layer, WL_interface, WR_interface);

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

    for (int j = core_start; j < core_end; ++j) {
        const Flux& Fminus = fluxes(j - 1);
        const Flux& Fplus = fluxes(j);

        Primitive w = layer.GetPrimitive(j);
        Conservative Uj = EOS::PrimToCons(w, gamma);

        Uj -= dt_over_dx * Flux::Diff(Fplus, Fminus);

        PositivityLimiter::Apply(Uj, gamma, rho_min_, p_min_);
        StoreConservativeCell(Uj, j, dx, layer);
    }

    t_cur += dt;
    return dt;
}

void GodunovKolganRodionovSolver::StoreConservativeCell(const Conservative& uc,
                                                        const int i, const double dx,
                                                        DataLayer& layer) const {
    const double rho = uc.rho;
    const double rhoU = uc.rhoU;

    const double uvel = (rho > 0.0) ? (rhoU / rho) : 0.0;
    const double P = EOS::Pressure(uc, settings_.gamma);

    layer.rho(i) = rho;
    layer.u(i) = uvel;
    layer.P(i) = P;

    layer.p(i) = rhoU;
    layer.V(i) = (rho > 0.0) ? (1.0 / rho) : 0.0;

    const double kinetic = 0.5 * rho * uvel * uvel;
    const double Eint = uc.E - kinetic;
    const double eint = (rho > 0.0) ? (Eint / rho) : 0.0;
    const double Etot = (rho > 0.0) ? uc.E : 0.0;

    layer.U(i) = eint;
    layer.e(i) = Etot;
    layer.m(i) = rho * dx;
}