#ifndef PERIODICBOUNDARY_HPP
#define PERIODICBOUNDARY_HPP

#include "bc/BoundaryCondition.hpp"

/**
 * @class PeriodicBoundary
 * @brief Placeholder periodic boundary.
 *
 * Proper periodic support should be implemented through explicit face pairing
 * in the mesh or a dedicated periodic mapping layer.
 */
class PeriodicBoundary final : public BoundaryCondition {
public:
    [[nodiscard]] PrimitiveCell BuildExteriorState(const DataLayer& layer,
                                                   const Mesh& mesh,
                                                   const Face& face,
                                                   const PrimitiveCell& interior_state) const override;
};

#endif  // PERIODICBOUNDARY_HPP
