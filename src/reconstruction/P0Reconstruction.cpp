#include "reconstruction/P0Reconstruction.hpp"

#include <stdexcept>

#include "data/Workspace.hpp"
#include "geometry/Face.hpp"
#include "geometry/Mesh.hpp"

PrimitiveCell P0Reconstruction::LoadCellPrimitive(const Workspace& workspace,
                                                  const std::size_t cell_id) const {
    const auto& W = workspace.W();

    PrimitiveCell state;
    state.rho = W(cell_id, Workspace::k_rho);
    state.u = W(cell_id, Workspace::k_u);
    state.v = W(cell_id, Workspace::k_v);
    state.w = W(cell_id, Workspace::k_w);
    state.P = W(cell_id, Workspace::k_p);

    return state;
}

void P0Reconstruction::ReconstructInteriorFace(const Mesh& mesh,
                                               const Workspace& workspace,
                                               const Face& face,
                                               PrimitiveCell& owner_state,
                                               PrimitiveCell& neighbor_state) const {
    (void)mesh;

    if (!(face.IsInternal() || face.IsMPIBoundary())) {
        throw std::runtime_error(
            "P0Reconstruction::ReconstructInteriorFace: face is not internal or MPI boundary"
        );
    }

    owner_state = LoadCellPrimitive(workspace, face.owner_cell_id);
    neighbor_state = LoadCellPrimitive(workspace, face.neighbor_cell_id);
}

void P0Reconstruction::ReconstructBoundaryFaceInterior(const Mesh& mesh,
                                                       const Workspace& workspace,
                                                       const Face& face,
                                                       PrimitiveCell& interior_state) const {
    (void)mesh;

    if (!face.IsPhysicalBoundary()) {
        throw std::runtime_error(
            "P0Reconstruction::ReconstructBoundaryFaceInterior: face is not physical boundary"
        );
    }

    interior_state = LoadCellPrimitive(workspace, face.owner_cell_id);
}
