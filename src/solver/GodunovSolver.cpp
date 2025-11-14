#include "solver/GodunovSolver.hpp"

#include <algorithm>
#include <cctype>
#include <string>
#include <vector>

#include "reconstruction/P0Reconstruction.hpp"
#include "riemann/ExactIdealGasRiemannSolver.hpp"
#include "riemann/HLLCRiemannSolver.hpp"
#include "riemann/HLLRiemannSolver.hpp"
#include "solver/EOS.hpp"
#include "solver/PositivityLimiter.hpp"
#include "solver/TimeStepCalculator.hpp"
#include "solver/Variables.hpp"

GodunovSolver::GodunovSolver(const Settings& settings)
    : settings_(settings),
      boundary_manager_(settings.dim),
      rho_min_(1e-10),
      p_min_(1e-10) {
    cfl_ = settings_.cfl;
    InitializeReconstruction();
    InitializeRiemannSolver();
}

auto GodunovSolver::Step(DataLayer& layer, double& t_cur) -> double {
    boundary_manager_.ApplyAll(layer);

    const double dx = ComputeDx(layer);
    const int total_size = layer.GetTotalSize();
    const int core_start = layer.GetCoreStart(0);
    const int core_end = layer.GetCoreEndExclusive(0);

    if (core_end - core_start < 2 || total_size < 3) {
        return 0.0;
    }

    std::vector<Primitive> prim(static_cast<std::size_t>(total_size));
    for (int i = 0; i < total_size; ++i) {
        Primitive w;
        w.rho = layer.rho(i);
        w.u = layer.u(i);
        w.P = layer.P(i);
        prim[static_cast<std::size_t>(i)] = w;
    }

    std::vector<Primitive> prim_core;
    prim_core.reserve(static_cast<std::size_t>(core_end - core_start));
    for (int i = core_start; i < core_end; ++i) {
        prim_core.push_back(prim[static_cast<std::size_t>(i)]);
    }

    double dt = TimeStepCalculator::ComputeDt(prim_core, dx, cfl_, settings_.gamma);
    if (dt <= 0.0) {
        return 0.0;
    }

    if (t_cur + dt > settings_.t_end) {
        dt = settings_.t_end - t_cur;
        if (dt <= 0.0) {
            return 0.0;
        }
    }

    const int n_interfaces = total_size - 1;
    std::vector<Primitive> left_states(static_cast<std::size_t>(n_interfaces));
    std::vector<Primitive> right_states(static_cast<std::size_t>(n_interfaces));
    reconstruction_->Reconstruct(prim, left_states, right_states);

    std::vector<Flux> fluxes(static_cast<std::size_t>(n_interfaces));
    for (int i = 0; i < n_interfaces; ++i) {
        fluxes[static_cast<std::size_t>(i)] = riemann_solver_->ComputeFlux(
            left_states[static_cast<std::size_t>(i)],
            right_states[static_cast<std::size_t>(i)], settings_.gamma, settings_.Q_user);
    }

    std::vector<Conservative> cons(static_cast<std::size_t>(total_size));
    for (int i = 0; i < total_size; ++i) {
        const Primitive& w = prim[static_cast<std::size_t>(i)];

        Conservative u;
        u.rho = w.rho;
        u.rhoU = w.rho * w.u;
        u.E = w.P / (settings_.gamma - 1.0) + 0.5 * w.rho * w.u * w.u;

        cons[static_cast<std::size_t>(i)] = u;
    }

    std::vector<Conservative> updated = cons;
    const double dt_over_dx = dt / dx;

    for (int j = core_start; j < core_end; ++j) {
        const Flux& f_minus = fluxes[static_cast<std::size_t>(j - 1)];
        const Flux& f_plus = fluxes[static_cast<std::size_t>(j)];

        Conservative u = cons[static_cast<std::size_t>(j)];

        u.rho -= dt_over_dx * (f_plus.mass - f_minus.mass);
        u.rhoU -= dt_over_dx * (f_plus.momentum - f_minus.momentum);
        u.E -= dt_over_dx * (f_plus.energy - f_minus.energy);

        PositivityLimiter::Apply(u, settings_.gamma, rho_min_, p_min_);
        updated[static_cast<std::size_t>(j)] = u;
    }

    StoreConservativeArray(updated, layer);

    t_cur += dt;
    return dt;
}

void GodunovSolver::SetCfl(double cfl) { cfl_ = cfl; }

void GodunovSolver::AddBoundary(int axis, std::shared_ptr<BoundaryCondition> left_bc,
                                std::shared_ptr<BoundaryCondition> right_bc) {
    boundary_manager_.Set(axis, std::move(left_bc), std::move(right_bc));
}

void GodunovSolver::InitializeReconstruction() {
    std::string name = settings_.reconstruction;
    std::string lower(name.size(), '\0');
    std::transform(name.begin(), name.end(), lower.begin(), [](unsigned char c) -> char {
        return static_cast<char>(std::tolower(c));
    });

    reconstruction_ = std::make_shared<P0Reconstruction>();
}

void GodunovSolver::InitializeRiemannSolver() {
    std::string name = settings_.riemann_solver;
    std::string lower(name.size(), '\0');
    std::transform(name.begin(), name.end(), lower.begin(), [](unsigned char c) -> char {
        return static_cast<char>(std::tolower(c));
    });

    if (lower.find("hll") != std::string::npos) {
        riemann_solver_ = std::make_shared<HLLRiemannSolver>();
    } else if (lower.find("exact") != std::string::npos) {
        riemann_solver_ = std::make_shared<ExactIdealGasRiemannSolver>();
    } else {
        riemann_solver_ = std::make_shared<HLLCRiemannSolver>();
    }
}

auto GodunovSolver::ComputeDx(const DataLayer& layer) const -> double {
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

void GodunovSolver::BuildPrimitiveArray(const DataLayer& layer,
                                        std::vector<Primitive>& primitives) const {
    const int total_size = layer.GetTotalSize();
    primitives.resize(static_cast<std::size_t>(total_size));

    for (int i = 0; i < total_size; ++i) {
        Primitive w;
        w.rho = layer.rho(i);
        w.u = layer.u(i);
        w.P = layer.P(i);
        primitives[static_cast<std::size_t>(i)] = w;
    }
}

void GodunovSolver::StoreConservativeArray(
    const std::vector<Conservative>& updated_conservative, DataLayer& layer) const {
    const int core_start = layer.GetCoreStart(0);
    const int core_end = layer.GetCoreEndExclusive(0);
    const int total_size = layer.GetTotalSize();

    if (static_cast<int>(updated_conservative.size()) < total_size) {
        return;
    }

    const double dx = ComputeDx(layer);

    for (int i = core_start; i < core_end; ++i) {
        const Conservative& uc = updated_conservative[static_cast<std::size_t>(i)];

        const double rho = uc.rho;
        const double u_vel = uc.rhoU / rho;
        const double P = EOS::Pressure(uc, settings_.gamma);

        layer.rho(i) = rho;
        layer.u(i) = u_vel;
        layer.P(i) = P;

        layer.p(i) = uc.rhoU;

        layer.V(i) = rho > 0.0 ? 1.0 / rho : 0.0;

        const double kinetic = 0.5 * rho * u_vel * u_vel;
        const double E_int = uc.E - kinetic;
        const double e_int = rho > 0.0 ? E_int / rho : 0.0;
        const double e_tot = rho > 0.0 ? uc.E : 0.0;

        layer.U(i) = e_int;
        layer.e(i) = e_tot;
        layer.m(i) = rho * dx;
    }
}
