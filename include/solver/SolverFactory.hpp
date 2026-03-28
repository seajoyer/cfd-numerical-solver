#ifndef SOLVERFACTORY_HPP
#define SOLVERFACTORY_HPP

#include <memory>

#include "bc/BoundaryManager.hpp"
#include "config/Settings.hpp"
#include "geometry/Mesh.hpp"
#include "solver/Solver.hpp"

class MPIContext;

/**
 * @class SolverFactory
 * @brief Factory for constructing solver pipeline from runtime settings.
 *
 * Assembles:
 * - spatial operator
 * - time integrator
 * - finite-volume solver
 *
 * Boundary conditions are assumed to be already created and registered in
 * BoundaryManager by boundary tag.
 */
class SolverFactory {
public:
    /**
     * @brief Create solver instance from runtime settings.
     *
     * @param settings Runtime settings.
     * @param mesh Generic mesh.
     * @param boundary_manager Boundary manager indexed by boundary tag.
     * @param mpi_context Optional MPI context placeholder for future use.
     * @return Constructed solver instance.
     */
    static std::unique_ptr<Solver> Create(const Settings& settings,
                                          Mesh mesh,
                                          const std::shared_ptr<BoundaryManager>& boundary_manager,
                                          const MPIContext* mpi_context);
};

#endif  // SOLVERFACTORY_HPP
