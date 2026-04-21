#ifndef PISOSOLVER_HPP
#define PISOSOLVER_HPP

#include "solver/PressureVelocitySolver.hpp"

class PisoSolver final : public PressureVelocitySolver {
public:
    PisoSolver(const Settings& settings,
               Mesh mesh,
               std::shared_ptr<BoundaryManager> boundary_manager,
               const MPIContext* mpi_context);
    // CHECK: PISO_LOOP
    auto Step(DataLayer& layer, double& t_cur) -> double override;
};

#endif
