#include "bc/FreeStreamBoundary.hpp"

#include <utility>

#include "data/DataLayer.hpp"
#include "geometry/Face.hpp"
#include "geometry/Mesh.hpp"

FreeStreamBoundary::FreeStreamBoundary(PrimitiveCell farfield_state)
    : farfield_state_(std::move(farfield_state)) {}

PrimitiveCell FreeStreamBoundary::BuildExteriorState(const DataLayer& layer,
                                                     const Mesh& mesh,
                                                     const Face& face,
                                                     const PrimitiveCell& interior_state) const {
    (void)layer;
    (void)mesh;
    (void)face;
    (void)interior_state;
    return farfield_state_;
}
