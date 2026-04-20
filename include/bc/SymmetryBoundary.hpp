#ifndef SYMMETRYBOUNDARY_HPP
#define SYMMETRYBOUNDARY_HPP

#include "bc/BoundaryCondition.hpp"

/**
 * @class SymmetryBoundary
 * @brief Symmetry boundary.
 *
 * Conservative branch:
 *  - ghost cells mirror interior values
 *  - normal momentum changes sign
 *  - tangential momenta and scalars are copied
 *
 * Pressure-velocity branch:
 *  - pressure uses zero normal gradient
 *  - normal face velocity is zero
 *  - tangential face velocities use zero normal gradient
 *
 * Assembly branch:
 *  - tangential momentum equations use zero diffusive flux through symmetry plane
 *  - pressure-correction equation keeps natural zero-normal-gradient behavior
 */
class SymmetryBoundary final : public BoundaryCondition {
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

#endif  // SYMMETRYBOUNDARY_HPP
