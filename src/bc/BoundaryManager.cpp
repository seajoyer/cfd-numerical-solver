#include "bc/BoundaryManager.hpp"
#include "bc/BoundaryCondition.hpp"

// TODO: implement support of multiple dimensions

BoundaryManager::BoundaryManager(int dim) : axes_(static_cast<std::size_t>(dim)) {
}

void BoundaryManager::Set(int axis,
                          std::shared_ptr<BoundaryCondition> left,
                          std::shared_ptr<BoundaryCondition> right) {
    auto idx = static_cast<std::size_t>(axis);
    axes_.at(idx).left = std::move(left);
    axes_.at(idx).right = std::move(right);
}

void BoundaryManager::ApplyAll(DataLayer &layer) const {
    for (int axis = 0; axis < static_cast<int>(axes_.size()); ++axis) {
        const auto &bc = axes_[static_cast<std::size_t>(axis)];
        if (bc.left) bc.left->Apply(layer, axis, Side::kLeft);
        if (bc.right) bc.right->Apply(layer, axis, Side::kRight);
    }
}
