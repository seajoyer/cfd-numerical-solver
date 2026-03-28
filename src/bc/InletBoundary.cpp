#include "bc/InletBoundary.hpp"

#include <utility>

#include "data/DataLayer.hpp"
#include "geometry/Face.hpp"
#include "geometry/Mesh.hpp"

InletBoundary::InletBoundary(PrimitiveCell inflow_state)
    : inflow_state_(std::move(inflow_state)) {}

PrimitiveCell InletBoundary::BuildExteriorState(const DataLayer& layer,
                                                const Mesh& mesh,
                                                const Face& face,
                                                const PrimitiveCell& interior_state) const {
    (void)layer;
    (void)mesh;
    (void)face;
    (void)interior_state;
    return inflow_state_;
}
