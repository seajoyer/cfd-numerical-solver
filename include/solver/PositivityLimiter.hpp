#ifndef POSITIVITYLIMITER_HPP
#define POSITIVITYLIMITER_HPP

#include "Variables.hpp"

/**
 * @class PositivityLimiter
 * @brief Ensures physically admissible conservative states.
 *
 * Provides a simple correction mechanism to enforce positivity of
 * density and pressure in conservative variables for the 1D Euler
 * equations. Intended as a safety net after numerical updates.
 */
class PositivityLimiter {
public:
    /**
     * @brief Applies positivity corrections to a conservative state.
     *
     * If density or pressure (computed from conservative variables)
     * fall below specified thresholds, they are corrected to restore
     * a physically meaningful state.
     *
     * @param u Conservative variables to be adjusted in-place.
     * @param gamma Ratio of specific heats.
     * @param rhoMin Minimal allowed density.
     * @param pMin Minimal allowed pressure.
     */
    static void Apply(Conservative &u,
                      double gamma,
                      double rho_min,
                      double p_min);
};

#endif  // POSITIVITYLIMITER_HPP
