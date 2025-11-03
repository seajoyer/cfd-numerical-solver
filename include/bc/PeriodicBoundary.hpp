#ifndef PERIODICBOUNDARY_HPP
#define PERIODICBOUNDARY_HPP

#include "bc/BoundaryCondition.hpp"

/**
 * @class PeriodicBoundary
 * @brief Implements periodic boundary conditions for the computational domain.
 *
 * This boundary condition enforces periodicity by copying data from the
 * opposite side of the domain into ghost cells. Physically, it represents
 * an infinitely repeating medium where flow exiting one boundary re-enters
 * from the other side.
 *
 * In 1D, this means:
 *  - Left ghost cells are filled with values from the right edge of the domain.
 *  - Right ghost cells are filled with values from the left edge of the domain.
 *
 *
 * @see ReflectiveBoundary
 * @see InletBoundary
 * @see OutletBoundary
 */
class PeriodicBoundary : public BoundaryCondition {
public:
    /**
     * @brief Applies the periodic boundary condition along a specified axis.
     *
     * Copies data from the opposite side of the computational domain into the
     * corresponding ghost cells to ensure continuity across boundaries.
     *
     * @param layer Reference to the DataLayer whose ghost cells are updated.
     * @param axis Axis along which periodicity is applied (0 for 1D).
     * @param side Side of the domain being processed (Side::Min or Side::Max).
     *
     * @note Currently implemented for 1D grids. For higher dimensions, the same logic
     *       will be applied independently along each coordinate direction.
     */
    void Apply(DataLayer &layer, int axis, Side side) const override;
};

#endif // PERIODICBOUNDARY_HPP
