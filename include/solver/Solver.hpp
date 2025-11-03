#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <memory>
#include "bc/BoundaryManager.hpp"

struct DataLayer;
class BoundaryCondition;


/**
 * @class Solver
 * @brief Abstract base class for all numerical solvers.
 *
 * Provides a unified interface for advancing physical quantities in time.
 * Concrete solvers (e.g., GodunovSolver) implement numerical schemes
 * for different orders of accuracy or flux computation methods.
 *
 * This class follows the Dependency Inversion Principle (DIP) —
 * higher-level modules depend on abstractions rather than specific implementations.
 *
 * @note Solver does not store configuration of the grid or time controller directly;
 *       it operates on an external DataLayer that contains all variable arrays.
 */
class Solver {
public:
    virtual ~Solver() = default;


    virtual void Step(DataLayer &layer, double &time, double t_end) = 0;


    virtual void SetCfl(double value) = 0;


    virtual void AddBoundary(int axis,
                             std::shared_ptr<BoundaryCondition> min,
                             std::shared_ptr<BoundaryCondition> max) = 0;
};

#endif // SOLVER_HPP
