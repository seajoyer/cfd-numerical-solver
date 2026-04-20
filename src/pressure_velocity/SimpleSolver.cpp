#include "pressure_velocity/SimpleSolver.hpp"

#include <stdexcept>

SimpleSolver::SimpleSolver(const Settings& settings,
                           Mesh mesh,
                           std::shared_ptr<BoundaryManager> boundary_manager,
                           const MPIContext* mpi_context)
    : PressureVelocitySolver(settings,
                             std::move(mesh),
                             std::move(boundary_manager),
                             mpi_context,
                             /*steady=*/true) {}

auto SimpleSolver::Step(DataLayer& layer, double& t_cur) -> double {
    (void)layer;

    EnsureStorageSized();

    if (mesh_.GetDim() != 2) {
        throw std::runtime_error("SimpleSolver: only dim=2 is supported for now");
    }

    state_.CopyCurrentToOld();

    const double dt = 1.0;
    const double nu = settings_.kinematic_viscosity;
    const double alpha_u = settings_.alpha_u;
    const double alpha_p = settings_.alpha_p;

    if (alpha_u <= 0.0 || alpha_u > 1.0) {
        throw std::runtime_error("SimpleSolver: alpha_u must be in (0, 1]");
    }
    if (alpha_p <= 0.0 || alpha_p > 1.0) {
        throw std::runtime_error("SimpleSolver: alpha_p must be in (0, 1]");
    }

    BuildMomentumCoefficients(dt, nu);
    SolveMomentumPredictor(dt, alpha_u);

    BuildPressureCorrectionEquation(dt, nu);
    SolvePressureCorrection();
    ApplyPressureCorrection(alpha_p);

    t_cur += 1.0;
    return dt;
}
