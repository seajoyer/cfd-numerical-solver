#include "bc/BoundaryManager.hpp"
#include "bc/BoundaryCondition.hpp"

BoundaryManager::BoundaryManager(int dim) : axes_(static_cast<std::size_t>(dim)) {
}

void BoundaryManager::Set(int axis,
                          std::shared_ptr<BoundaryCondition> min,
                          std::shared_ptr<BoundaryCondition> max) {
    auto idx = static_cast<std::size_t>(axis);
    axes_.at(idx).min = std::move(min);
    axes_.at(idx).max = std::move(max);
}

void BoundaryManager::ApplyAll(DataLayer &layer) const {
    for (int axis = 0; axis < static_cast<int>(axes_.size()); ++axis) {
        const auto &bc = axes_[static_cast<std::size_t>(axis)];
        if (bc.min) bc.min->Apply(layer, axis, Side::kMin);
        if (bc.max) bc.max->Apply(layer, axis, Side::kMax);
    }
}
