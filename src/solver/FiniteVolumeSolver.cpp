#include "solver/FiniteVolumeSolver.hpp"

#include "bc/BoundaryCondition.hpp"

FiniteVolumeSolver::FiniteVolumeSolver(const Settings& settings,
                                       std::shared_ptr<SpatialOperator> spatial_operator,
                                       std::shared_ptr<TimeIntegrator> time_integrator)
    : settings_(settings),
      spatial_operator_(std::move(spatial_operator)),
      time_integrator_(std::move(time_integrator)) {
    cfl_ = settings_.cfl;

    if (!spatial_operator_) {
        throw std::runtime_error("FiniteVolumeSolver: spatial_operator is null");
    }
    if (!time_integrator_) {
        throw std::runtime_error("FiniteVolumeSolver: time_integrator is null");
    }

    // Важно: boundary_manager НЕ передаем каждый шаг — он хранится внутри SpatialOperator.
    // Нужен один из вариантов:
    //  (A) SpatialOperator имеет SetBoundaryManager(shared_ptr<BoundaryManager>)
    //  (B) BoundaryManager передается в ctor SpatialOperator
    // Ниже — вариант (A). Если у тебя его нет, скажи — адаптирую под твой SpatialOperator.
    // spatial_operator_->SetBoundaryManager(boundary_manager_);

    time_integrator_->SetPositivityThresholds(rho_min_, p_min_);
}

void FiniteVolumeSolver::SetCfl(const double cfl) {
    cfl_ = cfl;
}

void FiniteVolumeSolver::EnsureWorkspaceSized(const DataLayer& layer) {
    workspace_.ResizeFrom(layer);
}

auto FiniteVolumeSolver::Step(DataLayer& layer, double& t_cur) -> double {
    EnsureWorkspaceSized(layer);

    // dt из conservative U + метрик DataLayer (inv_dx/...)
    double dt = TimeStepCalculator::ComputeDt(layer, settings_.gamma, cfl_);
    if (dt <= 0.0) {
        return 0.0;
    }

    if (t_cur + dt > settings_.t_end) {
        dt = settings_.t_end - t_cur;
        if (dt <= 0.0) {
            return 0.0;
        }
    }

    // TimeIntegrator НЕ трогает BC/halo и не знает про оси.
    // Он вызывает SpatialOperator::ComputeRhs, который внутри:
    // UpdateHalo -> ApplyPhysicalBc -> U->W -> divergence per axis -> rhs
    time_integrator_->Advance(layer, workspace_, dt, settings_.gamma, *spatial_operator_);

    t_cur += dt;
    return dt;
}