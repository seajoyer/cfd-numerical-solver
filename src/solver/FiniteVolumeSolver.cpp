#include "solver/FiniteVolumeSolver.hpp"

FiniteVolumeSolver::FiniteVolumeSolver(const Settings& settings,
                                       Mesh mesh,
                                       std::shared_ptr<SpatialOperator> spatial_operator,
                                       std::shared_ptr<TimeIntegrator> time_integrator,
                                       const MPIContext* mpi_context) : settings_(settings),
      mesh_(std::move(mesh)),
      spatial_operator_(std::move(spatial_operator)),
      time_integrator_(std::move(time_integrator)),
      mpi_context_(mpi_context) {
    cfl_ = settings_.cfl;

    if (!spatial_operator_) {
        throw std::runtime_error("FiniteVolumeSolver: spatial_operator is null");
    }
    if (!time_integrator_) {
        throw std::runtime_error("FiniteVolumeSolver: time_integrator is null");
    }

    diffusion_ = nullptr;
    if (settings.diffusion) {
        diffusion_ = std::make_unique<SolutionFilter>(settings_);
    }

    time_integrator_->SetPositivityThresholds(rho_min_, p_min_);
}

void FiniteVolumeSolver::SetCfl(const double cfl) {
    cfl_ = cfl;
}

const Mesh& FiniteVolumeSolver::GetMesh() const {
    return mesh_;
}

Mesh& FiniteVolumeSolver::GetMesh() {
    return mesh_;
}

void FiniteVolumeSolver::EnsureWorkspaceSized() {
    workspace_.ResizeFrom(mesh_);
}

auto FiniteVolumeSolver::Step(DataLayer& layer, double& t_cur) -> double {
    EnsureWorkspaceSized();

    double dt_local = TimeStepCalculator::ComputeDt(layer, mesh_, settings_.gamma, cfl_);
    double dt = dt_local;

    if (mpi_context_) {
        const double large_dt = 1e300;
        if (dt_local <= 0.0) {
            dt_local = large_dt;
        }
        dt = mpi_context_->GlobalMin(dt_local);
        if (dt >= large_dt) {
            dt = 0.0;
        }
    }

    if (dt <= 0.0) {
        return 0.0;
    }

    if (t_cur + dt > settings_.t_end) {
        dt = settings_.t_end - t_cur;
        if (dt <= 0.0) {
            return 0.0;
        }
    }

    time_integrator_->Advance(layer, mesh_, workspace_, dt, settings_.gamma, *spatial_operator_);

    if (diffusion_) {
        diffusion_->Apply(layer, mesh_, settings_.gamma);
    }

    t_cur += dt;
    return dt;
}
