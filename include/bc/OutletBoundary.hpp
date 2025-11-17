#ifndef OUTLETBOUNDARY_HPP
#define OUTLETBOUNDARY_HPP

#include "bc/BoundaryCondition.hpp"

/**
 * @class OutletBoundary
 * @brief Implements outlet (outflow) boundary conditions.
 *
 * This boundary condition allows flow to freely exit the computational
 * domain without reflection or artificial pressure buildup. It enforces
 * a zero-gradient (Neumann) condition by copying the nearest interior
 * cell values into the ghost cells.
 *
 * Physically, it assumes that disturbances leaving the domain do not
 * affect the interior flow — the state at the outlet is determined
 * solely by upstream dynamics.
 *
 * In 1D:
 *  - Left side (Side::Min): ghost cells copy from the first physical cell.
 *  - Right side (Side::Max): ghost cells copy from the last physical cell.
 *
 *
 * @see InletBoundary
 * @see FreeStreamBoundary
 * @see NonReflectingBoundary
 */
class OutletBoundary : public BoundaryCondition {
   public:
    /**
     * @brief Applies outlet boundary condition along the specified axis.
     *
     * Fills ghost cells by copying values from adjacent interior cells,
     * thereby enforcing a zero spatial gradient at the boundary.
     *
     * @param layer Reference to the DataLayer being modified.
     * @param axis Axis index (0 for 1D problems).
     * @param side Which side to apply: Side::Min (left) or Side::Max (right).
     *
     * @note For multidimensional extensions, the same logic is applied
     *       independently for each axis.
     */
    void Apply(DataLayer& layer, int axis, Side side) const override;
};

#endif  // OUTLETBOUNDARY_HPP
