#ifndef SLIPWALLINTERNALBOUNDARY_HPP
#define SLIPWALLINTERNALBOUNDARY_HPP

#include "bc/InternalBoundaryCondition.hpp"

/**
 * @class SlipWallInternalBoundary
 * @brief Slip-wall treatment for immersed internal fluid-solid interfaces.
 *
 * The wall is impermeable and frictionless:
 *  - normal velocity component is reflected
 *  - tangential velocity components are preserved
 *  - density and pressure are preserved
 *
 * If the face normal is not normalized, it is normalized internally.
 * If the normal magnitude is degenerate, the fluid state is copied unchanged.
 */
class SlipWallInternalBoundary final : public InternalBoundaryCondition {
public:
    SlipWallInternalBoundary() = default;
    ~SlipWallInternalBoundary() override = default;

    void BuildBoundaryState(const PrimitiveCell& fluid_state,
                            const ImmersedFaceInfo& face_info,
                            PrimitiveCell& boundary_state) const override;
};

#endif  // SLIPWALLINTERNALBOUNDARY_HPP