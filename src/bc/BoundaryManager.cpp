#include "bc/BoundaryManager.hpp"

#include "bc/BoundaryCondition.hpp"

BoundaryManager::BoundaryManager(int dim) : axes_(static_cast<std::size_t>(dim)) {}

void BoundaryManager::Set(int axis, std::shared_ptr<BoundaryCondition> left_bc,
                          std::shared_ptr<BoundaryCondition> right_bc) {
    auto idx = static_cast<std::size_t>(axis);

    axes_.at(idx).left_bc = std::move(left_bc);
    axes_.at(idx).right_bc = std::move(right_bc);
}

void BoundaryManager::ApplyAll(DataLayer& layer) const {
    for (int axis = 0; axis < static_cast<int>(axes_.size()); ++axis) {
        const auto& bc = axes_[static_cast<std::size_t>(axis)];

        if (bc.left_bc) bc.left_bc->Apply(layer, axis, Side::kLeft);
        if (bc.right_bc) bc.right_bc->Apply(layer, axis, Side::kRight);
    }
}
