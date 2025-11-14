#ifndef TIMESTEPCALCULATOR_HPP
#define TIMESTEPCALCULATOR_HPP

#include <vector>
#include "Variables.hpp"

/**
 * @class TimeStepCalculator
 * @brief Utility class for CFL-based time step computation.
 *
 * Provides a static method to compute the stable time step size
 * for explicit finite-volume schemes solving the 1D Euler equations.
 */
class TimeStepCalculator {
public:
    /**
     * @brief Computes the time step from CFL condition.
     *
     * Uses the maximal signal speed over all cells:
     * \f$ s_{\max} = \max_j ( |u_j| + a_j ) \f$
     * and returns:
     * \f$ \Delta t = \mathrm{CFL} \cdot \Delta x / s_{\max} \f$.
     *
     * @param cellStates Cell-centered primitive variables.
     * @param dx Cell size (assumed uniform).
     * @param cfl CFL number (0 < cfl <= 1).
     * @param gamma Ratio of specific heats.
     * @return Stable time step size.
     */
    static auto ComputeDt(const std::vector<Primitive> &cell_states,
                            double dx,
                            double cfl,
                            double gamma) -> double;
};

#endif  // TIMESTEPCALCULATOR_HPP
