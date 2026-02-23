// PositivityLimiter.hpp
#ifndef POSITIVITYLIMITER_HPP
#define POSITIVITYLIMITER_HPP

#include "data/DataLayer.hpp"

/**
 * @class PositivityLimiter
 * @brief Ensures physically admissible conservative Euler state U (rho and P floors).
 *
 * Applies in-place corrections on conservative U(var,i,j,k):
 *  - rho is clamped to rho_min (preserving velocities when possible)
 *  - pressure P(U) is clamped to p_min by adjusting total energy E
 *
 * Intended as a safety net after explicit updates.
 */
class PositivityLimiter final {
public:
    /**
     * @brief Apply positivity corrections in-place on core cells.
     * @param layer DataLayer containing U.
     * @param gamma Ratio of specific heats.
     * @param rho_min Minimum density.
     * @param p_min Minimum pressure.
     */
    static void Apply(DataLayer& layer, double gamma, double rho_min, double p_min);
};

#endif  // POSITIVITYLIMITER_HPP