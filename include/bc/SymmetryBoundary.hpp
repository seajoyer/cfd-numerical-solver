#ifndef SYMMETRYBOUNDARY_HPP
#define SYMMETRYBOUNDARY_HPP

#include "bc/BoundaryCondition.hpp"

/**
 * @class SymmetryBoundary
 * @brief Implements symmetry boundary conditions.
 *
 * This boundary condition enforces mirror symmetry of all physical variables
 * relative to a plane (or line in 1D) of symmetry. It ensures that the normal
 * component of velocity and fluxes vanish across the symmetry plane, while all
 * scalar fields are mirrored with zero normal gradients.
 *
 * In 1D, the symmetry boundary is equivalent to a reflective wall:
 *  - The normal velocity component (u) is reversed in sign.
 *  - Scalar quantities (ρ, P, e, etc.) are copied symmetrically.
 *
 * Physically, this condition is used to exploit geometric or flow symmetry,
 * reducing computational cost by modeling only half (or a fraction) of a domain.
 *
 *
 * @see ReflectiveBoundary
 * @see WallBoundary
 */
class SymmetryBoundary : public BoundaryCondition {
public:
    /**
     * @brief Applies symmetry boundary condition along a specified axis.
     *
     * Mirrors scalar fields and reverses the normal velocity component
     * to enforce symmetry relative to the domain boundary.
     *
     * @param layer Reference to the DataLayer being modified.
     * @param axis Axis index (0 for 1D problems).
     * @param side Which side to apply: Side::Min (left) or Side::Max (right).
     *
     * @note For 1D cases, this operation is equivalent to a reflective boundary.
     *       For higher dimensions, the logic extends to the axis normal direction.
     */
    void Apply(DataLayer &layer, int axis, Side side) const override;
};

#endif // SYMMETRYBOUNDARY_HPP