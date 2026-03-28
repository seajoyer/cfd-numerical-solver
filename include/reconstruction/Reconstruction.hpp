#ifndef RECONSTRUCTION_HPP
#define RECONSTRUCTION_HPP

#include "data/Variables.hpp"

class Mesh;
class Workspace;
class Face;

/**
 * @class Reconstruction
 * @brief Interface for face-based reconstruction on generic cell-centered meshes.
 *
 * Reconstruction works on one mesh face at a time and uses cell-centered
 * primitive cache stored in Workspace::W(cell,var).
 *
 * Contract:
 * - No allocations in hot path.
 * - Does not modify mesh or workspace.
 * - For internal faces returns left/right states adjacent to the face.
 * - For boundary faces returns only the interior state on the owner side.
 *
 * State orientation:
 * - "left" means owner-side state
 * - "right" means neighbor-side state for internal faces
 *
 * Boundary handling:
 * - Reconstruction does not build exterior boundary state.
 * - Exterior state is built later by BoundaryCondition from the reconstructed
 *   interior boundary state.
 */
class Reconstruction {
public:
    virtual ~Reconstruction() = default;

    /**
     * @brief Reconstruct owner-side and neighbor-side states on one internal face.
     *
     * @param mesh Mesh with geometry and connectivity.
     * @param workspace Workspace containing primitive cache W(cell,var).
     * @param face Internal face.
     * @param owner_state Reconstructed state on owner side of the face.
     * @param neighbor_state Reconstructed state on neighbor side of the face.
     */
    virtual void ReconstructInteriorFace(const Mesh& mesh,
                                         const Workspace& workspace,
                                         const Face& face,
                                         PrimitiveCell& owner_state,
                                         PrimitiveCell& neighbor_state) const = 0;

    /**
     * @brief Reconstruct owner-side interior state on one boundary face.
     *
     * @param mesh Mesh with geometry and connectivity.
     * @param workspace Workspace containing primitive cache W(cell,var).
     * @param face Boundary face.
     * @param interior_state Reconstructed owner-side state adjacent to the face.
     */
    virtual void ReconstructBoundaryFaceInterior(const Mesh& mesh,
                                                 const Workspace& workspace,
                                                 const Face& face,
                                                 PrimitiveCell& interior_state) const = 0;
};

#endif  // RECONSTRUCTION_HPP
