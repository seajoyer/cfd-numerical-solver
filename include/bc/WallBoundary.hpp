#ifndef WALLBOUNDARY_HPP
#define WALLBOUNDARY_HPP

#include "bc/BoundaryCondition.hpp"

/**
 * @class WallBoundary
 * @brief Solid slip wall boundary for inviscid Euler equations.
 */
class WallBoundary final : public BoundaryCondition {
public:
    [[nodiscard]] PrimitiveCell BuildExteriorState(const DataLayer& layer,
                                                   const Mesh& mesh,
                                                   const Face& face,
                                                   const PrimitiveCell& interior_state) const override;
};

#endif  // WALLBOUNDARY_HPP
