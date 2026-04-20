#ifndef INTERNALBOUNDARYCONDITION_HPP
#define INTERNALBOUNDARYCONDITION_HPP

#include "data/Mesh.hpp"
#include "data/PressureVelocityState.hpp"
#include "data/PressureVelocityWorkspace.hpp"
#include "data/Variables.hpp"

/**
 * @class InternalBoundaryCondition
 * @brief Abstract base class for immersed internal boundary treatment.
 *
 * Internal boundary conditions are applied on fluid-solid interfaces inside the domain.
 *
 * Two supported usages:
 *  - compressible / Riemann-based branch via BuildBoundaryState(...)
 *  - pressure-velocity branch via ApplyPressureVelocity... hooks
 */
class InternalBoundaryCondition {
public:
    virtual ~InternalBoundaryCondition() = default;

    /**
     * @brief Build boundary-compatible primitive state at an immersed face.
     */
    virtual void BuildBoundaryState(const PrimitiveCell& fluid_state,
                                    const ImmersedFaceInfo& face_info,
                                    PrimitiveCell& boundary_state) const = 0;

    /**
     * @brief Apply immersed-face velocity constraints for pressure-velocity solvers.
     *
     * Default implementation does nothing.
     */
    virtual void ApplyPressureVelocityVelocityConstraints(PressureVelocityState& state,
                                                          PressureVelocityWorkspace& workspace,
                                                          const Mesh& mesh) const {
        (void)state;
        (void)workspace;
        (void)mesh;
    }

    /**
     * @brief Apply immersed-wall momentum coefficient corrections for pressure-velocity solvers.
     *
     * Default implementation does nothing.
     */
    virtual void ApplyPressureVelocityMomentumCorrections(PressureVelocityState& state,
                                                          PressureVelocityWorkspace& workspace,
                                                          const Mesh& mesh,
                                                          bool steady,
                                                          double dt,
                                                          double nu) const {
        (void)state;
        (void)workspace;
        (void)mesh;
        (void)steady;
        (void)dt;
        (void)nu;
    }
};

#endif  // INTERNALBOUNDARYCONDITION_HPP
