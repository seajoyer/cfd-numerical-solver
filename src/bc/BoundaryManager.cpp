#include "bc/BoundaryManager.hpp"

#include "bc/BoundaryCondition.hpp"

BoundaryManager::BoundaryManager() : axes_(3) {}

void BoundaryManager::Set(const Axis axis,
                          std::shared_ptr<BoundaryCondition> left_bc,
                          std::shared_ptr<BoundaryCondition> right_bc) {
    axes_.at(static_cast<std::size_t>(axis)).left_bc = std::move(left_bc);
    axes_.at(static_cast<std::size_t>(axis)).right_bc = std::move(right_bc);
}

void BoundaryManager::UpdateHalo(DataLayer& /*layer*/) const {
    // MPI halo exchange will live here later.
    // For single-domain smoke-test: no-op.
}

void BoundaryManager::ApplyPhysicalBc(DataLayer& layer) const {
    const int dim = layer.GetDim();

    // X axis always active
    {
        const AxisBc& bc = axes_[static_cast<std::size_t>(Axis::X)];
        if (bc.left_bc && layer.IsGlobalBoundary(Axis::X, Side::Left)) {
            bc.left_bc->Apply(layer, Axis::X, Side::Left);
        }
        if (bc.right_bc && layer.IsGlobalBoundary(Axis::X, Side::Right)) {
            bc.right_bc->Apply(layer, Axis::X, Side::Right);
        }
    }

    // Y axis only if dim >= 2
    if (dim >= 2) {
        const AxisBc& bc = axes_[static_cast<std::size_t>(Axis::Y)];
        if (bc.left_bc && layer.IsGlobalBoundary(Axis::Y, Side::Left)) {
            bc.left_bc->Apply(layer, Axis::Y, Side::Left);
        }
        if (bc.right_bc && layer.IsGlobalBoundary(Axis::Y, Side::Right)) {
            bc.right_bc->Apply(layer, Axis::Y, Side::Right);
        }
    }

    // Z axis only if dim >= 3
    if (dim >= 3) {
        const AxisBc& bc = axes_[static_cast<std::size_t>(Axis::Z)];
        if (bc.left_bc && layer.IsGlobalBoundary(Axis::Z, Side::Left)) {
            bc.left_bc->Apply(layer, Axis::Z, Side::Left);
        }
        if (bc.right_bc && layer.IsGlobalBoundary(Axis::Z, Side::Right)) {
            bc.right_bc->Apply(layer, Axis::Z, Side::Right);
        }
    }
}

const AxisBc& BoundaryManager::Get(const Axis axis) const {
    return axes_.at(static_cast<std::size_t>(axis));
}