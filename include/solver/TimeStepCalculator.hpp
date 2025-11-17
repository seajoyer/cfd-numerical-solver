#ifndef TIMESTEPCALCULATOR_HPP
#define TIMESTEPCALCULATOR_HPP

#include "data/DataLayer.hpp"

/**
 * @class TimeStepCalculator
 * @brief Utility for CFL-based timestep selection.
 *
 * This helper class provides a static function to compute a stable
 * explicit timestep from the CFL condition, scanning the primitive
 * fields stored in a DataLayer instance.
 *
 * The maximal signal speed is estimated as
 *
 *  |u| + c,   where c = sqrt(gamma * P / rho)
 *
 * over all physical cells, and the timestep is chosen as
 *
 *  dt = cfl * dx / maxSpeed.
 */
class TimeStepCalculator {
public:
    /**
     * @brief Computes a stable timestep from the CFL condition.
     *
     * The function iterates over the physical cells (core region) of the
     * DataLayer and estimates the maximal signal speed |u| + c. Cells
     * with non-positive rho or P are skipped to avoid invalid sound speeds.
     *
     * If no valid cell is found or the maximal speed is non-positive,
     * the function returns dt = 0.0.
     *
     * @param layer DataLayer with current primitive fields (rho, u, P).
     * @param dx    Cell size in x-direction.
     * @param cfl   CFL number (0 < cfl <= 1).
     * @param gamma Ratio of specific heats.
     * @return Stable timestep dt; 0.0 if no stable timestep can be computed.
     */
    static double ComputeDt(const DataLayer &layer,
                            double dx,
                            double cfl,
                            double gamma);
};

#endif  // TIMESTEPCALCULATOR_HPP
