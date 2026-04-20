#ifndef OUTLETBOUNDARY_HPP
#define OUTLETBOUNDARY_HPP

#include "bc/BoundaryCondition.hpp"

/**
 * @class OutletBoundary
 * @brief Zero-gradient outlet boundary.
 *
 * Conservative branch:
 *  - copies the nearest interior layer into ghost layers
 *
 * Pressure-velocity branch:
 *  - applies zero normal gradient to pressure and velocity components
 *
 * Assembly branch:
 *  - currently leaves momentum and pressure-correction operators unchanged
 *    and relies on natural zero-gradient behavior of the discrete stencil
 *
 * @note
 * This class represents a Neumann-type outlet.
 * A fixed-pressure outlet should be implemented as a separate boundary class.
 */
class OutletBoundary final : public BoundaryCondition {
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

#endif  // OUTLETBOUNDARY_HPP
