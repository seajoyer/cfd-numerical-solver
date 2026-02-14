#ifndef TIMESTEPCALCULATOR_HPP
#define TIMESTEPCALCULATOR_HPP

#include "data/DataLayer.hpp"

/**
 * @class TimeStepCalculator
 * @brief Utility for CFL-based timestep selection (1D and 2D).
 *
 * 1D: dt = cfl * dx / max(|u| + c)
 * 2D: dt = cfl / max( (|u|+c)/dx + (|v|+c)/dy )
 */
class TimeStepCalculator {
   public:
    /**
     * @brief Computes stable timestep (1D or 2D with dx=dy).
     */
    static auto ComputeDt(const DataLayer& layer, double dx, double cfl, double gamma)
        -> double;

    /**
     * @brief Computes stable timestep for 2D with separate dx, dy.
     */
    static auto ComputeDt2D(const DataLayer& layer, double dx, double dy,
                            double cfl, double gamma) -> double;
};

#endif  // TIMESTEPCALCULATOR_HPP
