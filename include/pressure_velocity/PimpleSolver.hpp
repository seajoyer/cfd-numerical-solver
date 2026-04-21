#ifndef PIMPLESOLVER_HPP
#define PIMPLESOLVER_HPP

#include "solver/PressureVelocitySolver.hpp"

class PimpleSolver final : public PressureVelocitySolver {
public:
    PimpleSolver(const Settings& settings,
                 Mesh mesh,
                 std::shared_ptr<BoundaryManager> boundary_manager,
                 const MPIContext* mpi_context);
    // CHECK: PIMPLE_LOOP
    auto Step(DataLayer& layer, double& t_cur) -> double override;
};

#endif
