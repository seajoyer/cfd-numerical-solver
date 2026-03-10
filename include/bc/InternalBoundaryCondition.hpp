#ifndef INTERNALBOUNDARYCONDITION_HPP
#define INTERNALBOUNDARYCONDITION_HPP

#include "data/Mesh.hpp"
#include "data/Variables.hpp"

/**
 * @class InternalBoundaryCondition
 * @brief Abstract base class for immersed internal boundary treatment.
 *
 * Internal boundary conditions are applied on fluid-solid interfaces inside the domain.
 * Unlike external boundary conditions, they are not tied to global domain sides and do
 * not fill outer ghost layers. Instead, they construct a boundary-compatible state for
 * an immersed face adjacent to a fluid cell.
 *
 * Typical usage:
 *  - Spatial operator detects an immersed face for a fluid cell.
 *  - Fluid-side primitive state is known.
 *  - InternalBoundaryCondition builds a mirrored / wall-compatible state on the solid side.
 *  - Riemann solver uses (fluid_state, boundary_state) at the immersed face.
 */
class InternalBoundaryCondition {
public:
    /** @brief Virtual destructor for safe polymorphic deletion. */
    virtual ~InternalBoundaryCondition() = default;

    /**
     * @brief Build boundary-compatible primitive state at an immersed face.
     *
     * @param fluid_state Primitive state in the adjacent fluid cell.
     * @param face_info Immersed-face geometry data.
     * @param boundary_state Output primitive state representing the solid-side boundary treatment.
     */
    virtual void BuildBoundaryState(const PrimitiveCell& fluid_state,
                                    const ImmersedFaceInfo& face_info,
                                    PrimitiveCell& boundary_state) const = 0;
};

#endif  // INTERNALBOUNDARYCONDITION_HPP