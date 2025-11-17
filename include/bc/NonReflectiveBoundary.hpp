#ifndef NONREFLECTINGBOUNDARY_HPP
#define NONREFLECTINGBOUNDARY_HPP

#include "bc/BoundaryCondition.hpp"

/**
 * @class NonReflectingBoundary
 * @brief Implements non-reflecting (radiation-type) boundary conditions.
 *
 * This boundary condition minimizes artificial reflections of acoustic
 * or pressure waves at the computational boundaries. It allows disturbances
 * to freely exit the domain while preventing incoming reflected waves
 * from re-entering.
 *
 * Physically, it approximates an "open" boundary through characteristic-based
 * relations that decompose the flow into incoming and outgoing wave components.
 * In simple implementations, it may fall back to zero-gradient behavior when
 * characteristic information is insufficient.
 *
 * In 1D:
 *  - On outflow, extrapolates variables to preserve outgoing wave behavior.
 *  - On inflow, may combine fixed external data with characteristic relations.
 *
 * @see OutletBoundary
 * @see FreeStreamBoundary
 * @see InletBoundary
 */
class NonReflectiveBoundary : public BoundaryCondition {
   public:
    /**
     * @brief Applies non-reflective boundary condition along the specified axis.
     *
     * Updates ghost cells based on outgoing wave information to suppress
     * artificial reflections. Typically, uses characteristic decomposition
     * or simplified extrapolation of flow variables.
     *
     * @param layer Reference to the DataLayer being modified.
     * @param axis Axis index (0 for 1D problems).
     * @param side Boundary side where condition is applied (Side::Min or Side::Max).
     *
     * @note In the current implementation, this serves as a placeholder.
     *       Future extensions may compute characteristic wave propagation
     *       using local sound speed and flow direction.
     */
    void Apply(DataLayer& layer, int axis, Side side) const override;
};

#endif  // NONREFLECTINGBOUNDARY_HPP
