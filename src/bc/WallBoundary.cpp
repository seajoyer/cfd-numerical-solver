#include "bc/WallBoundary.hpp"

#include "data/DataLayer.hpp"
#include "geometry/Face.hpp"
#include "geometry/Mesh.hpp"

PrimitiveCell WallBoundary::BuildExteriorState(const DataLayer& layer,
                                               const Mesh& mesh,
                                               const Face& face,
                                               const PrimitiveCell& interior_state) const {
    (void)layer;
    (void)mesh;

    const FaceNormal normal{face.normal_x, face.normal_y, face.normal_z};
    const double vn = NormalVelocity(interior_state, normal);

    PrimitiveCell exterior = interior_state;
    exterior.u = interior_state.u - 2.0 * vn * normal.x;
    exterior.v = interior_state.v - 2.0 * vn * normal.y;
    exterior.w = interior_state.w - 2.0 * vn * normal.z;

    return exterior;
}
