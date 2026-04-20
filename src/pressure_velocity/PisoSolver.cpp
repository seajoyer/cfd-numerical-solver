#include "pressure_velocity/PisoSolver.hpp"

#include <algorithm>
#include <stdexcept>

PisoSolver::PisoSolver(const Settings& settings,
                       Mesh mesh,
                       std::shared_ptr<BoundaryManager> boundary_manager,
                       const MPIContext* mpi_context)
    : PressureVelocitySolver(settings,
                             std::move(mesh),
                             std::move(boundary_manager),
                             mpi_context,
                             /*steady=*/false) {}

auto PisoSolver::Step(DataLayer& layer, double& t_cur) -> double {
    (void)layer;

    EnsureStorageSized();

    if (mesh_.GetDim() != 2) {
        throw std::runtime_error("PisoSolver: only dim=2 is supported for now");
    }

    state_.CopyCurrentToOld();

    const double dt = ComputeDt(t_cur);
    if (dt <= 0.0) {
        return 0.0;
    }

    const double nu = settings_.kinematic_viscosity;

    BuildMomentumCoefficients(dt, nu);
    SolveMomentumPredictor(dt, 1.0);

    const int n_corr = std::max(1, settings_.n_pressure_correctors);

    for (int corr = 0; corr < n_corr; ++corr) {
        BuildPressureCorrectionEquation(dt, nu);
        SolvePressureCorrection();
        ApplyPressureCorrection(1.0);

        workspace_.UxStar() = state_.Ux();
        workspace_.VyStar() = state_.Vy();
    }

    t_cur += dt;
    return dt;
}
