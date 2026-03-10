#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <memory>

class DataLayer;
class BoundaryCondition;
enum class Axis : std::uint8_t;

/**
 * @class Solver
 * @brief Abstract base class for explicit time-stepping solvers.
 */
class Solver {
public:
    virtual ~Solver() = default;

    /**
     * @brief Advances simulation by one step.
     * @param layer DataLayer (updated in-place).
     * @param t_cur Current time (incremented by dt).
     * @return dt actually taken, or 0.0 if step is skipped.
     */
    virtual auto Step(DataLayer& layer, double& t_cur) -> double = 0;

    /** @brief Set CFL number. */
    virtual void SetCfl(double cfl) = 0;

protected:
    double cfl_ = 0.5;
};

#endif  // SOLVER_HPP