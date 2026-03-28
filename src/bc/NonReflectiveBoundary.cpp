#include "bc/NonReflectiveBoundary.hpp"

#include <utility>

#include "data/DataLayer.hpp"
#include "geometry/Face.hpp"
#include "geometry/Mesh.hpp"

NonReflectiveBoundary::NonReflectiveBoundary(PrimitiveCell farfield_state,
                                             const double gamma)
    : farfield_state_(std::move(farfield_state)),
      gamma_(gamma) {}

PrimitiveCell NonReflectiveBoundary::BuildExteriorState(const DataLayer& layer,
                                                        const Mesh& mesh,
                                                        const Face& face,
                                                        const PrimitiveCell& interior_state) const {
    (void)layer;
    (void)mesh;
    (void)face;
    (void)interior_state;
    (void)gamma_;

    return farfield_state_;
}