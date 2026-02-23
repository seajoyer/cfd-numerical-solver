#ifndef WALLBOUNDARY_HPP
#define WALLBOUNDARY_HPP

#include "bc/BoundaryCondition.hpp"
#include "data/Variables.hpp"

/**
 * @class WallBoundary
 * @brief Stationary wall boundary for conservative Euler state U.
 *
 * Implementation (Euler-compatible kinematic wall):
 *  - Copy rho and E from nearest interior core layer (zero-gradient for scalars).
 *  - Set all momentum components (rhoU, rhoV, rhoW) to zero in ghost cells.
 *
 * Note: For inviscid Euler, the physically standard wall is slip (ReflectiveBoundary).
 * This class enforces a strict zero-velocity condition in ghost cells.
 */
class WallBoundary final : public BoundaryCondition {
public:
    void Apply(DataLayer& layer, Axis axis, Side side) const override;
};

#endif  // WALLBOUNDARY_HPP