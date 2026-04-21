#ifndef BOUNDARYCONDITION_HPP
#define BOUNDARYCONDITION_HPP

#include "data/Variables.hpp"

class DataLayer;
class Mesh;
class Face;
// CHECK: GHOST_STATE
/**
 * @class BoundaryCondition
 * @brief Abstract physical boundary condition for one boundary face.
 *
 * Contract:
 * - Works on one boundary face at a time.
 * - Does not modify mesh or storage.
 * - Returns the external primitive state used by the Riemann solver.
 */
class BoundaryCondition {
public:
    virtual ~BoundaryCondition() = default;

    /**
     * @brief Build external primitive state for a boundary face.
     *
     * @param layer Current numerical solution.
     * @param mesh Mesh with geometry and connectivity.
     * @param face Boundary face.
     * @param interior_state Primitive state in owner cell adjacent to the face.
     * @return External primitive state at the boundary.
     */
    [[nodiscard]] virtual PrimitiveCell BuildExteriorState(const DataLayer& layer,
                                                           const Mesh& mesh,
                                                           const Face& face,
                                                           const PrimitiveCell& interior_state) const = 0;
};

#endif  // BOUNDARYCONDITION_HPP
