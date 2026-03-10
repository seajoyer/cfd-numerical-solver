#ifndef SOLVERFACTORY_HPP
#define SOLVERFACTORY_HPP

#include <memory>

#include "bc/BoundaryManager.hpp"
#include "bc/BoundaryFactory.hpp"
#include "config/Settings.hpp"
#include "data/Mesh.hpp"
#include "data/Variables.hpp"
#include "solver/Solver.hpp"

/**
 * @class SolverFactory
 * @brief Factory class for creating solver instances based on string identifiers.
 *
 * This class provides a centralized way to create different solver types
 * without coupling the Simulation class to specific solver implementations.
 */
class SolverFactory {
public:
    SolverFactory() : boundary_manager_(std::make_shared<BoundaryManager>()) {}

    /**
     * @brief Creates a solver instance based on the solver type string.
     *
     * @param settings Settings for solver construction.
     * @param mesh Structured mesh for solver construction.
     * @param boundary_manager Boundary manager for spatial operator.
     * @return Unique pointer to the created solver.
     * @throws std::runtime_error if solver type is not recognized.
     */
    static auto Create(const Settings& settings,
                       Mesh mesh,
                       const std::shared_ptr<BoundaryManager>& boundary_manager,
                       const MPIContext* mpi_context) -> std::unique_ptr<Solver>;

    /** @brief Set boundary conditions for an axis. */
    void AddBoundary(Axis axis,
                     std::shared_ptr<BoundaryCondition> left_bc,
                     std::shared_ptr<BoundaryCondition> right_bc);

private:
    std::shared_ptr<BoundaryManager> boundary_manager_;
};

#endif  // SOLVERFACTORY_HPP
