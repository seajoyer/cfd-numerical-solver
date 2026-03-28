#include "bc/OutletBoundary.hpp"

#include "data/DataLayer.hpp"
#include "geometry/Face.hpp"
#include "geometry/Mesh.hpp"

PrimitiveCell OutletBoundary::BuildExteriorState(const DataLayer& layer,
                                                 const Mesh& mesh,
                                                 const Face& face,
                                                 const PrimitiveCell& interior_state) const {
    (void)layer;
    (void)mesh;
    (void)face;
    return interior_state;
}
