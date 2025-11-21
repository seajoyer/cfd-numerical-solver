#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <memory>

struct DataLayer;
class BoundaryCondition;

/**
 * @class Solver
 * @brief Abstract base class for numerical solvers.
 *
 * This class defines the interface for time-stepping algorithms
 * used to advance the simulation state. Concrete implementations
 * (e.g., GodunovSolver) provide specific numerical schemes.
 */
class Solver {
   public:
    virtual ~Solver() = default;

    /**
     * @brief Advances the simulation by one time step.
     * @param layer Data layer containing current state (modified in-place)
     * @param t_cur Current simulation time (will be updated)
     * @return Actual time step taken (dt)
     */
    virtual auto Step(DataLayer& layer, double& t_cur) -> double = 0;

    /**
     * @brief Sets the CFL number for time step calculation.
     * @param cfl CFL number (typically 0.4-0.9)
     */
    virtual void SetCfl(double cfl) = 0;

    /**
     * @brief Adds boundary conditions for a specific axis.
     * @param axis Spatial axis index (0=x, 1=y, 2=z)
     * @param left_bc Boundary condition for left/lower boundary
     * @param right_bc Boundary condition for right/upper boundary
     */
    virtual void AddBoundary(int axis, std::shared_ptr<BoundaryCondition> left_bc,
                             std::shared_ptr<BoundaryCondition> right_bc) = 0;

   protected:
    double cfl_ = 0.5;
};

#endif  // SOLVER_HPP
