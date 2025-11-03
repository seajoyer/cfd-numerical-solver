#ifndef REFLECTIVEBOUNDARY_HPP
#define REFLECTIVEBOUNDARY_HPP

#include "bc/BoundaryCondition.hpp"

/**
 * @class ReflectiveBoundary
 * @brief Implements reflective (wall) boundary conditions.
 *
 * This boundary condition models a perfectly reflecting solid wall,
 * where the normal component of velocity (or momentum) is reversed in sign,
 * while all scalar quantities (density, pressure, energy, etc.) are mirrored.
 *
 * Physically, it enforces the condition of zero normal mass flux through
 * the wall — meaning the wall is impermeable but may reflect the flow.
 *
 * In 1D:
 *  - Left ghost cells copy values from the leftmost interior region
 *    but with the velocity component inverted.
 *  - Right ghost cells copy from the rightmost interior region
 *    and also invert the velocity component.
 *
 *
 * @see PeriodicBoundary
 * @see SymmetryBoundary
 * @see WallBoundary
 */
class ReflectiveBoundary : public BoundaryCondition {
public:
    /**
     * @brief Applies reflective boundary condition along the specified axis.
     *
     * Copies scalar variables into ghost cells as mirror images of the inner cells,
     * while inverting the sign of the normal velocity component.
     *
     * @param layer Reference to the DataLayer being updated.
     * @param axis Axis along which the reflection is applied (0 for 1D).
     * @param side Side of the domain (Side::Min or Side::Max).
     *
     * @note Currently implemented for 1D domains; extension to higher dimensions
     *       would involve reversing the normal component of velocity along
     *       the corresponding axis.
     */
    void Apply(DataLayer &layer, int axis, Side side) const override;
};

#endif // REFLECTIVEBOUNDARY_HPP
