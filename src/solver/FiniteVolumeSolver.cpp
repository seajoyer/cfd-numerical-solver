#include "solver/FiniteVolumeSolver.hpp"

#include <utility>

#include "bc/BoundaryCondition.hpp"
#include "solver/TimeStepCalculator.hpp"
#include "solver/EOS.hpp"
#include "solver/PositivityLimiter.hpp"
#include "time/ForwardEulerTimeIntegrator.hpp"
#include "time/MacCormackTimeIntegrator.hpp"
#include "time/SSPRK2TimeIntegrator.hpp"
#include "time/SSPRK3TimeIntegrator.hpp"

FiniteVolumeSolver::FiniteVolumeSolver(const Settings& settings,
                                       std::shared_ptr<SpatialOperator> spatial_operator)
    : settings_(settings),
      spatial_operator_(std::move(spatial_operator)),
      rho_min_(1e-10),
      p_min_(1e-10) {
    cfl_ = settings_.cfl;

    boundary_manager_ = std::make_shared<BoundaryManager>(settings.dim);
    global_limiter_ = nullptr;
    if (settings.global_limiter) {
        global_limiter_ = std::make_unique<GlobalLimiter>();
    }
    diffusion_ = nullptr;
    if (settings.diffusion) {
        diffusion_ = std::make_unique<SolutionFilter>();
    }
    time_integrator_ = nullptr;
    InitializeTimeIntegrator();
    if (time_integrator_) {
        time_integrator_->SetPositivityThresholds(rho_min_, p_min_);
    }
    vacuum_fix_limiter_ = std::make_unique<VacuumFixLimiter>();
}

void FiniteVolumeSolver::InitializeTimeIntegrator() {
    std::string t_integrator = settings_.time_integrator;

    if (t_integrator == "ssprk2") {
        time_integrator_ = std::make_shared<SSPRK2TimeIntegrator>(boundary_manager_);
    } else if (t_integrator == "ssprk3") {
        time_integrator_ = std::make_shared<SSPRK3TimeIntegrator>(boundary_manager_);
    } else if (t_integrator == "euler") {
        time_integrator_ = std::make_shared<ForwardEulerTimeIntegrator>();
    } else if (settings_.solver == "maccormack") {
        time_integrator_ = std::make_shared<MacCormackTimeIntegrator>(boundary_manager_, settings_);
    }

    if (time_integrator_) return;

    std::cerr << "\nError: Unknown time integrator: " << settings_.time_integrator << "\n\n";
    throw std::runtime_error("Unknown time integrator type: " + settings_.time_integrator);
}

void FiniteVolumeSolver::SetCfl(double cfl) { cfl_ = cfl; }

void FiniteVolumeSolver::AddBoundary(int axis,
                                     std::shared_ptr<BoundaryCondition> left_bc,
                                     std::shared_ptr<BoundaryCondition> right_bc) {
    boundary_manager_->Set(axis, std::move(left_bc), std::move(right_bc));
}

auto FiniteVolumeSolver::ComputeDx(const DataLayer& layer) const -> double {
    const int nx = (settings_.GetNx() > 0) ? settings_.GetNx() : settings_.N;
    if (nx > 0 && settings_.L_x > 0.0) {
        return settings_.L_x / static_cast<double>(nx);
    }
    // Fallback: deduce from DataLayer
    const int core_start = layer.GetCoreStart(0);
    const int core_end = layer.GetCoreEndExclusive(0);
    if (layer.GetDim() <= 1 && core_end - core_start > 1) {
        return layer.xc(core_start + 1) - layer.xc(core_start);
    }
    return 1.0;
}

auto FiniteVolumeSolver::ComputeDy(const DataLayer& layer) const -> double {
    const int ny = (settings_.GetNy() > 0) ? settings_.GetNy() : settings_.N;
    if (ny > 0 && settings_.L_y > 0.0) {
        return settings_.L_y / static_cast<double>(ny);
    }
    return 1.0;
}

auto FiniteVolumeSolver::Step(DataLayer& layer, double& t_cur) -> double {
    if (!time_integrator_) return 0.0;

    boundary_manager_->ApplyAll(layer);

    const double dx = ComputeDx(layer);
    double dt;

    if (layer.GetDim() >= 2) {
        const double dy = ComputeDy(layer);
        dt = TimeStepCalculator::ComputeDt2D(layer, dx, dy, cfl_, settings_.gamma);
    } else {
        const int total_size = layer.GetTotalSize();
        const int core_start = layer.GetCoreStart(0);
        const int core_end = layer.GetCoreEndExclusive(0);
        const int n_core = core_end - core_start;
        if (n_core < 2 || total_size < 3) return 0.0;
        dt = TimeStepCalculator::ComputeDt(layer, dx, cfl_, settings_.gamma);
    }

    if (dt <= 0.0) return 0.0;

    if (settings_.t_end > 0.0 && t_cur + dt > settings_.t_end) {
        dt = settings_.t_end - t_cur;
        if (dt <= 0.0) return 0.0;
    }

    if (layer.GetDim() >= 2) {
        // 2D time advance
        Step2D(layer, dt, dx, ComputeDy(layer));
    } else {
        // 1D time advance (existing path)
        time_integrator_->Advance(layer, dt, dx, settings_, *spatial_operator_);
        if (diffusion_) {
            solution_filter_.Apply(layer);
        }
        if (global_limiter_) {
            global_limiter_->Apply(layer, dx, settings_);
        }
    }

    t_cur += dt;
    return dt;
}

void FiniteVolumeSolver::Step2D(DataLayer& layer, double dt, double dx, double dy) {
    // Forward Euler for 2D: U^{n+1} = U^n + dt * RHS(U^n)
    const int tx = layer.GetTotalSize(0);
    const int ty = layer.GetTotalSize(1);
    const int cs_x = layer.GetCoreStart(0);
    const int ce_x = layer.GetCoreEndExclusive(0);
    const int cs_y = layer.GetCoreStart(1);
    const int ce_y = layer.GetCoreEndExclusive(1);
    const double gamma = settings_.gamma;

    // Extract conservative state
    xt::xarray<Conservative> U = xt::xarray<Conservative>::from_shape(
        {static_cast<std::size_t>(tx), static_cast<std::size_t>(ty)});

    for (int i = 0; i < tx; ++i) {
        for (int j = 0; j < ty; ++j) {
            U(i, j) = layer.GetConservative2D(i, j, gamma);
        }
    }

    // Compute RHS using the 2D spatial operator
    xt::xarray<Conservative> rhs;

    // Use the 2D-specific ComputeRHS2D if available, otherwise fall back
    // We need to cast to GodunovSpatialOperator2D to access ComputeRHS2D
    // For now, use the overloaded ComputeRHS which detects dim internally
    spatial_operator_->ComputeRHS(layer, dx, gamma, rhs);

    // Update core cells
    for (int i = cs_x; i < ce_x; ++i) {
        for (int j = cs_y; j < ce_y; ++j) {
            Conservative U_new = U(i, j);
            U_new.rho  += dt * rhs(i, j).rho;
            U_new.rhoU += dt * rhs(i, j).rhoU;
            U_new.rhoV += dt * rhs(i, j).rhoV;
            U_new.E    += dt * rhs(i, j).E;

            // Positivity limiting
            PositivityLimiter::Apply(U_new, gamma, rho_min_, p_min_);

            // Store back to DataLayer
            layer.SetConservative2D(i, j, U_new, gamma, dx, dy);
        }
    }
}
