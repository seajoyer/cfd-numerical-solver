#include "pressure_velocity/PimpleSolver.hpp"

#include <algorithm>
#include <stdexcept>

PimpleSolver::PimpleSolver(const Settings& settings,
                           Mesh mesh,
                           std::shared_ptr<BoundaryManager> boundary_manager,
                           const MPIContext* mpi_context)
    : PressureVelocitySolver(settings,
                             std::move(mesh),
                             std::move(boundary_manager),
                             mpi_context,
                             /*steady=*/false) {}

auto PimpleSolver::Step(DataLayer& layer, double& t_cur) -> double {
    (void)layer;

    EnsureStorageSized();

    if (mesh_.GetDim() != 2) {
        throw std::runtime_error("PimpleSolver: only dim=2 is supported for now");
    }

    state_.CopyCurrentToOld();

    const double dt = ComputeDt(t_cur);
    if (dt <= 0.0) {
        return 0.0;
    }

    const double nu = settings_.kinematic_viscosity;

    const int n_outer = std::max(1, settings_.n_outer_correctors);
    const int n_corr = std::max(1, settings_.n_pressure_correctors);

    const double alpha_u_relaxed = settings_.alpha_u;
    const double alpha_p_relaxed = settings_.alpha_p;

    if (alpha_u_relaxed <= 0.0 || alpha_u_relaxed > 1.0) {
        throw std::runtime_error("PimpleSolver: alpha_u must be in (0, 1]");
    }
    if (alpha_p_relaxed <= 0.0 || alpha_p_relaxed > 1.0) {
        throw std::runtime_error("PimpleSolver: alpha_p must be in (0, 1]");
    }

    for (int outer = 0; outer < n_outer; ++outer) {
        const bool last_outer = (outer == n_outer - 1);

        const double alpha_u = last_outer ? 1.0 : alpha_u_relaxed;
        const double alpha_p = last_outer ? 1.0 : alpha_p_relaxed;

        BuildMomentumCoefficients(dt, nu);
        SolveMomentumPredictor(dt, alpha_u);

        for (int corr = 0; corr < n_corr; ++corr) {
            BuildPressureCorrectionEquation(dt, nu);
            SolvePressureCorrection();
            ApplyPressureCorrection(alpha_p);

            workspace_.UxStar() = state_.Ux();
            workspace_.VyStar() = state_.Vy();
        }

        state_.PressureOld() = state_.Pressure();
        state_.UxOld() = state_.Ux();
        state_.VyOld() = state_.Vy();
    }

    t_cur += dt;
    return dt;
}
