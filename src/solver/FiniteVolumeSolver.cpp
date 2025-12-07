#include "solver/FiniteVolumeSolver.hpp"

#include <utility>

#include "bc/BoundaryCondition.hpp"
#include "solver/TimeStepCalculator.hpp"
#include "time/ForwardEulerTimeIntegrator.hpp"
#include "time/MacCormackTimeIntegrator.hpp"
#include "time/SSPRK2TimeIntegrator.hpp"
#include "time/SSPRK3TimeIntegrator.hpp"

FiniteVolumeSolver::FiniteVolumeSolver(const Settings& settings,
                                       std::shared_ptr<SpatialOperator> spatial_operator)
    : settings_(settings),
      boundary_manager_(settings.dim),
      spatial_operator_(std::move(spatial_operator)),
      rho_min_(1e-10),
      p_min_(1e-10) {
    cfl_ = settings_.cfl;

    time_integrator_ = nullptr;
    InitializeTimeIntegrator();

    if (time_integrator_) {
        time_integrator_->SetPositivityThresholds(rho_min_, p_min_);
    }
}

void FiniteVolumeSolver::InitializeTimeIntegrator() {
    std::string t_integrator = settings_.time_integrator;

    if (t_integrator == "ssprk2") {
        time_integrator_ = std::make_shared<SSPRK2TimeIntegrator>();
    } else if (t_integrator == "ssprk3") {
        time_integrator_ = std::make_shared<SSPRK3TimeIntegrator>();
    } else if (t_integrator == "euler") {
        time_integrator_ = std::make_shared<ForwardEulerTimeIntegrator>();
    } else if (settings_.solver == "maccormack") {
        time_integrator_ = std::make_shared<MacCormackTimeIntegrator>();
    }

    if (time_integrator_) {
        return;
    }
    std::cerr << "\nError: Unknown time integrator: " << settings_.time_integrator <<
        "\n\n";
    std::cerr << "Supported time integrators:\n";
    std::cerr << "  - Euler\n";
    std::cerr << "  - SSPRK2\n";
    std::cerr << "  - SSPRK3\n\n";
    std::cerr << "To see all time integrator options and their compatiblity, use:\n";
    std::cerr << "  ./cfd_numerical_solver --list-time-integrators\n\n";

    throw std::runtime_error("Unknown time integrator type: " + settings_.time_integrator);
}


void FiniteVolumeSolver::SetCfl(double cfl) { cfl_ = cfl; }

void FiniteVolumeSolver::AddBoundary(int axis,
                                     std::shared_ptr<BoundaryCondition> left_bc,
                                     std::shared_ptr<BoundaryCondition> right_bc) {
    boundary_manager_.Set(axis, std::move(left_bc), std::move(right_bc));
}

auto FiniteVolumeSolver::ComputeDx(const DataLayer& layer) const -> double {
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

auto FiniteVolumeSolver::Step(DataLayer& layer, double& t_cur) -> double {
    if (!spatial_operator_ || !time_integrator_) {
        return 0.0;
    }

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

    time_integrator_->Advance(layer, dt, dx, settings_, *spatial_operator_);

    t_cur += dt;
    return dt;
}