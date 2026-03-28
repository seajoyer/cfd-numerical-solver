#include "bc/BoundaryManager.hpp"

#include <stdexcept>

#include "bc/BoundaryCondition.hpp"
#include "data/DataLayer.hpp"
#include "geometry/Face.hpp"
#include "geometry/Mesh.hpp"

void BoundaryManager::Register(const int boundary_tag,
                               std::shared_ptr<BoundaryCondition> boundary_condition) {
    if (!boundary_condition) {
        throw std::invalid_argument("BoundaryManager::Register: boundary_condition is null");
    }

    by_tag_[boundary_tag] = std::move(boundary_condition);
}

bool BoundaryManager::Has(const int boundary_tag) const {
    return by_tag_.find(boundary_tag) != by_tag_.end();
}

const BoundaryCondition& BoundaryManager::Get(const int boundary_tag) const {
    const auto it = by_tag_.find(boundary_tag);
    if (it == by_tag_.end()) {
        throw std::runtime_error(
                                 "BoundaryManager::Get: no boundary condition registered for tag " +
                                 std::to_string(boundary_tag)
                                );
    }

    return *it->second;
}

PrimitiveCell BoundaryManager::BuildExteriorState(const DataLayer& layer,
                                                  const Mesh& mesh,
                                                  const Face& face,
                                                  const PrimitiveCell& interior_state) const {
    if (!face.IsBoundary()) {
        throw std::runtime_error(
                                 "BoundaryManager::BuildExteriorState: face is not a boundary face"
                                );
    }

    return Get(face.boundary_tag).BuildExteriorState(layer, mesh, face, interior_state);
}
