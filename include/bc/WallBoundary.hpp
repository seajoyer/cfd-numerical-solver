#ifndef WALLBOUNDARY_HPP
#define WALLBOUNDARY_HPP

#include "bc/BoundaryCondition.hpp"

/**
 * @class WallBoundary
 * @brief Stationary wall boundary for conservative Euler state U.
 *
 * Implementation:
 * - copies rho and E from the nearest interior core layer
 * - sets all momentum components (rhoU, rhoV, rhoW) to zero in ghost cells
 *
 * Note:
 * For inviscid Euler, reflective/slip-wall is usually more physical.
 * This boundary enforces zero velocity in ghost cells.
 */
class WallBoundary final : public BoundaryCondition {
public:
    void Apply(DataLayer& layer, const Mesh& mesh, Axis axis, Side side) const override;
};

#endif  // WALLBOUNDARY_HPP