#ifndef OUTLETBOUNDARY_HPP
#define OUTLETBOUNDARY_HPP

#include "bc/BoundaryCondition.hpp"

/**
 * @class OutletBoundary
 * @brief Zero-gradient outlet boundary.
 *
 * Exterior state equals interior state.
 */
class OutletBoundary final : public BoundaryCondition {
public:
    [[nodiscard]] PrimitiveCell BuildExteriorState(const DataLayer& layer,
                                                   const Mesh& mesh,
                                                   const Face& face,
                                                   const PrimitiveCell& interior_state) const override;
};

#endif  // OUTLETBOUNDARY_HPP