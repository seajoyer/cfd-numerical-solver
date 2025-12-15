#ifndef VACUUMFIX_HPP
#define VACUUMFIX_HPP

#include "data/DataLayer.hpp"
#include "config/Settings.hpp"

/**
 * @class VacuumFixLimiter
 * @brief Extremely simple post-processing fix for near-vacuum regions.
 *
 * This fix does NOT modify conservative variables.
 * It only suppresses computation of sensitive derived quantities
 * (specific internal energy, specific volume, optionally velocity)
 * when density and pressure are very small.
 *
 * Purpose:
 * - avoid Eint/rho blow-up
 * - avoid 1/rho blow-up
 * - stabilize plots and post-processing
 */
class VacuumFixLimiter {
public:
    VacuumFixLimiter() = default;

    /**
     * @brief Applies vacuum fix directly to DataLayer.
     *
     * @param layer    DataLayer to modify.
     * @param dx       Cell size.
     * @param settings Global settings (gamma used for consistency).
     */
    void Apply(DataLayer& layer, double dx, const Settings& settings) const;
};

#endif  // VACUUMFIX_HPP
