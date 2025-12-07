#include "solver/MacCormackSolver.hpp"

#include <algorithm>

#include "bc/BoundaryCondition.hpp"
#include "data/DataLayer.hpp"
#include "data/Variables.hpp"
#include "solver/EOS.hpp"
#include "solver/PositivityLimiter.hpp"
#include "solver/TimeStepCalculator.hpp"

MacCormackSolver::MacCormackSolver(const Settings& settings)
    : settings_(settings),
      boundary_manager_(settings.dim),
      rho_min_(1e-10),
      p_min_(1e-10) {
    cfl_ = settings_.cfl;
}

void MacCormackSolver::AddBoundary(int axis,
                                   std::shared_ptr<BoundaryCondition> left_bc,
                                   std::shared_ptr<BoundaryCondition> right_bc) {
    boundary_manager_.Set(axis, std::move(left_bc), std::move(right_bc));
}

auto MacCormackSolver::ComputeDx(const DataLayer& layer) const -> double {
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

auto MacCormackSolver::Step(DataLayer& layer, double& t_cur) -> double {
    boundary_manager_.ApplyAll(layer);

    const double dx = ComputeDx(layer);
    const int total_size = layer.GetTotalSize();
    const int core_start = std::max(0, layer.GetCoreStart(0));
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
    const double gamma = settings_.gamma;

    xt::xarray<Conservative> U =
        xt::xarray<Conservative>::from_shape({static_cast<std::size_t>(total_size)});
    xt::xarray<Conservative> U_pred =
        xt::xarray<Conservative>::from_shape({static_cast<std::size_t>(total_size)});

    xt::xarray<Flux> F =
        xt::xarray<Flux>::from_shape({static_cast<std::size_t>(total_size)});
    xt::xarray<Flux> F_pred =
        xt::xarray<Flux>::from_shape({static_cast<std::size_t>(total_size)});

    for (int j = 0; j < total_size; ++j) {
        Primitive w = layer.GetPrimitive(j);
        U(j) = EOS::PrimToCons(w, gamma);
        F(j) = EulerFlux(w, gamma);
    }

    for (int j = core_start; j < core_end; ++j) {
        const int jp1 = std::min(j + 1, total_size - 1);
        U_pred(j) = U(j) - dt_over_dx * (F(jp1) - F(j));
    }

    for (int j = 0; j < core_start; ++j) {
        U_pred(j) = U(j);
    }
    for (int j = core_end; j < total_size; ++j) {
        U_pred(j) = U(j);
    }

    for (int j = 0; j < total_size; ++j) {
        Primitive w_pred = EOS::ConsToPrim(U_pred(j), gamma);
        F_pred(j) = EulerFlux(w_pred, gamma);
    }

    for (int j = core_start; j < core_end; ++j) {
        const int jm1 = std::max(j - 1, 0);

        Conservative U_new = U(j);
        U_new += U_pred(j);
        U_new -= dt_over_dx * (F_pred(j) - F_pred(jm1));
        U_new *= 0.5;

        PositivityLimiter::Apply(U_new, gamma, rho_min_, p_min_);

        StoreConservativeCell(U_new, j, dx, layer);
    }

    t_cur += dt;
    return dt;
}

void MacCormackSolver::StoreConservativeCell(const Conservative& uc,
                                             const int i,
                                             const double dx,
                                             DataLayer& layer) const {
    const double rho = uc.rho;
    const double rhoU = uc.rhoU;
    const double uvel = rho > 0.0 ? rhoU / rho : 0.0;
    const double P = EOS::Pressure(uc, settings_.gamma);

    layer.rho(i) = rho;
    layer.u(i) = uvel;
    layer.P(i) = P;

    layer.p(i) = rhoU;
    layer.V(i) = rho > 0.0 ? 1.0 / rho : 0.0;

    const double kinetic = 0.5 * rho * uvel * uvel;
    const double Eint = uc.E - kinetic;
    const double eint = rho > 0.0 ? Eint / rho : 0.0;
    const double Etot = rho > 0.0 ? uc.E : 0.0;

    layer.U(i) = eint;
    layer.e(i) = Etot;
    layer.m(i) = rho * dx;
}
