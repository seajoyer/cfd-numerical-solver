#ifndef WALLBOUNDARY_HPP
#define WALLBOUNDARY_HPP

#include "bc/BoundaryCondition.hpp"

/**
 * @class WallBoundary
 * @brief Stationary no-slip wall boundary.
 *
 * Conservative branch:
 *  - copies scalar quantities from the nearest interior layer
 *  - sets momentum components in ghost cells to zero
 *
 * Pressure-velocity branch:
 *  - pressure uses zero normal gradient
 *  - normal face velocity at the wall is zero
 *  - tangential face velocities use odd reflection in ghost layers
 *
 * Pressure-velocity assembly:
 *  - modifies near-wall momentum coefficients using half-cell diffusion distance
 *  - pressure-correction stage keeps zero-normal-gradient behavior
 */
class WallBoundary final : public BoundaryCondition {
public:
    void Apply(DataLayer& layer, const Mesh& mesh, Axis axis, Side side) const override;

    void Apply(PressureVelocityState& state,
               const Mesh& mesh,
               Axis axis,
               Side side) const override;

    void ApplyPressureVelocityBoundary(PressureVelocityState& state,
                                       PressureVelocityWorkspace& workspace,
                                       const Mesh& mesh,
                                       Axis axis,
                                       Side side,
                                       PvAssemblyStage stage,
                                       bool steady,
                                       double dt,
                                       double nu) const override;
};

#endif  // WALLBOUNDARY_HPP
