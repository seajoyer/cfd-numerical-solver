#ifndef INLETBOUNDARY_HPP
#define INLETBOUNDARY_HPP

#include "bc/BoundaryCondition.hpp"
#include "bc/BoundaryFactory.hpp"

/**
 * @class InletBoundary
 * @brief Velocity inlet boundary with prescribed inflow state.
 *
 * Conservative branch:
 *  - fills ghost cells with prescribed conservative inflow state
 *
 * Pressure-velocity branch:
 *  - prescribes inlet velocity components
 *  - uses zero normal gradient for pressure
 *
 * Assembly branch:
 *  - modifies near-boundary momentum coefficients using half-cell diffusion distance
 */
class InletBoundary final : public BoundaryCondition {
public:
    explicit InletBoundary(const FarfieldConservative& inflow_U);
    explicit InletBoundary(const BoundaryStateSettings& primitive_state);

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

private:
    FarfieldConservative inflow_U_;
    BoundaryStateSettings primitive_state_;
};

#endif  // INLETBOUNDARY_HPP
