#ifndef WALLBOUNDARY_HPP
#define WALLBOUNDARY_HPP

#include "bc/BoundaryCondition.hpp"

/**
 * @class WallBoundary
 * @brief Implements impermeable wall boundary conditions.
 *
 * This boundary condition models a solid, stationary, and impermeable wall.
 * It enforces zero normal velocity at the wall surface (no mass flux through
 * the boundary), while scalar variables such as density and pressure are
 * mirrored to maintain continuity.
 *
 * In 1D, this reduces to reversing the sign of the velocity component normal
 * to the wall while copying other variables from adjacent interior cells.
 *
 * Physically, it represents a rigid wall where the fluid cannot penetrate
 * but may exert pressure. It is suitable for inviscid flow simulations
 * (Euler equations) and can serve as a simplified alternative to
 * no-slip walls in viscous models.
 *
 * @see ReflectiveBoundary
 * @see SymmetryBoundary
 */
class WallBoundary : public BoundaryCondition {
   public:
    /**
     * @brief Applies wall boundary condition along a specified axis.
     *
     * Enforces zero normal velocity at the wall by inverting the normal
     * component of velocity in ghost cells, while scalar quantities are
     * copied from interior cells adjacent to the wall.
     *
     * @param layer Reference to the DataLayer being modified.
     * @param axis Axis index (0 for 1D problems).
     * @param side Boundary side to which the condition is applied:
     *             Side::Min (left) or Side::Max (right).
     *
     * @note For 1D problems, this acts as a simple reflection.
     *       For multidimensional cases, only the velocity component normal
     *       to the wall is reversed.
     */
    void Apply(DataLayer& layer, int axis, Side side) const override;
};

#endif  // WALLBOUNDARY_HPP
