#ifndef WALLINTERNALBOUNDARY_HPP
#define WALLINTERNALBOUNDARY_HPP

#include "bc/InternalBoundaryCondition.hpp"

/**
 * @class WallInternalBoundary
 * @brief No-slip impermeable wall treatment for immersed internal fluid-solid interfaces.
 *
 * Compressible branch:
 *  - velocity is mirrored in all components for a no-slip wall surrogate
 *  - density and pressure are preserved
 *
 * Pressure-velocity branch:
 *  - normal staggered face velocities on fluid-solid interfaces are forced to zero
 *  - momentum coefficients near solid are corrected as no-slip half-cell walls
 */
class WallInternalBoundary final : public InternalBoundaryCondition {
public:
    WallInternalBoundary() = default;
    ~WallInternalBoundary() override = default;

    void BuildBoundaryState(const PrimitiveCell& fluid_state,
                            const ImmersedFaceInfo& face_info,
                            PrimitiveCell& boundary_state) const override;

    void ApplyPressureVelocityVelocityConstraints(PressureVelocityState& state,
                                                  PressureVelocityWorkspace& workspace,
                                                  const Mesh& mesh) const override;

    void ApplyPressureVelocityMomentumCorrections(PressureVelocityState& state,
                                                  PressureVelocityWorkspace& workspace,
                                                  const Mesh& mesh,
                                                  bool steady,
                                                  double dt,
                                                  double nu) const override;

private:
    double k_eps_ = 1e-14;
};

#endif  // WALLINTERNALBOUNDARY_HPP
