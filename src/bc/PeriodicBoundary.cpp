#include "bc/PeriodicBoundary.hpp"

#include <stdexcept>

#include "data/DataLayer.hpp"
#include "geometry/Face.hpp"
#include "geometry/Mesh.hpp"

PrimitiveCell PeriodicBoundary::BuildExteriorState(const DataLayer& layer,
                                                   const Mesh& mesh,
                                                   const Face& face,
                                                   const PrimitiveCell& interior_state) const {
    (void)layer;
    (void)mesh;
    (void)face;
    (void)interior_state;

    throw std::runtime_error(
                             "PeriodicBoundary: periodic boundary requires explicit periodic face mapping and is not implemented yet"
                            );
}
