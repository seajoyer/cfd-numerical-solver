#ifndef SIMPLESOLVER_HPP
#define SIMPLESOLVER_HPP

#include "solver/PressureVelocitySolver.hpp"

class SimpleSolver final : public PressureVelocitySolver {
public:
    SimpleSolver(const Settings& settings,
                 Mesh mesh,
                 std::shared_ptr<BoundaryManager> boundary_manager,
                 const MPIContext* mpi_context);
    // CHECK: SIMPLE_LOOP
    auto Step(DataLayer& layer, double& t_cur) -> double override;
};

#endif
