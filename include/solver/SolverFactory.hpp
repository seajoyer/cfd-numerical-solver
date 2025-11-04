#ifndef SOLVERFACTORY_HPP
#define SOLVERFACTORY_HPP

#include <memory>
#include <string>

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
    /**
     * @brief Creates a solver instance based on the solver type string.
     *
     * @param solver_type String identifier for the solver (e.g., "godunov").
     * @param dim Spatial dimension (1, 2, or 3).
     * @return Unique pointer to the created solver.
     * @throws std::runtime_error if solver type is not recognized.
     */
    static auto Create(const std::string& solver_type, int dim)
        -> std::unique_ptr<Solver>;
};

#endif  // SOLVERFACTORY_HPP
