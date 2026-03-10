#include "bc/BoundaryManager.hpp"

#include "bc/BoundaryCondition.hpp"
#include "data/DataLayer.hpp"
#include "data/Mesh.hpp"

BoundaryManager::BoundaryManager(std::shared_ptr<HaloExchange> halo_exchange)
    : axes_(3), halo_exchange_(std::move(halo_exchange)) {}

void BoundaryManager::Set(const Axis axis,
                          std::shared_ptr<BoundaryCondition> left_bc,
                          std::shared_ptr<BoundaryCondition> right_bc) {
    axes_.at(static_cast<std::size_t>(axis)).left_bc = std::move(left_bc);
    axes_.at(static_cast<std::size_t>(axis)).right_bc = std::move(right_bc);
}

void BoundaryManager::UpdateHalo(DataLayer& layer, const Mesh& mesh) const {
    if (!halo_exchange_) {
        return;
    }

    halo_exchange_->Exchange(layer, mesh);
}

void BoundaryManager::ApplyPhysicalBc(DataLayer& layer, const Mesh& mesh) const {
    const int dim = mesh.GetDim();

    {
        const AxisBc& bc = axes_[static_cast<std::size_t>(Axis::X)];
        if (bc.left_bc && mesh.IsGlobalBoundary(Axis::X, Side::Left)) {
            bc.left_bc->Apply(layer, mesh, Axis::X, Side::Left);
        }
        if (bc.right_bc && mesh.IsGlobalBoundary(Axis::X, Side::Right)) {
            bc.right_bc->Apply(layer, mesh, Axis::X, Side::Right);
        }
    }

    if (dim >= 2) {
        const AxisBc& bc = axes_[static_cast<std::size_t>(Axis::Y)];
        if (bc.left_bc && mesh.IsGlobalBoundary(Axis::Y, Side::Left)) {
            bc.left_bc->Apply(layer, mesh, Axis::Y, Side::Left);
        }
        if (bc.right_bc && mesh.IsGlobalBoundary(Axis::Y, Side::Right)) {
            bc.right_bc->Apply(layer, mesh, Axis::Y, Side::Right);
        }
    }

    if (dim >= 3) {
        const AxisBc& bc = axes_[static_cast<std::size_t>(Axis::Z)];
        if (bc.left_bc && mesh.IsGlobalBoundary(Axis::Z, Side::Left)) {
            bc.left_bc->Apply(layer, mesh, Axis::Z, Side::Left);
        }
        if (bc.right_bc && mesh.IsGlobalBoundary(Axis::Z, Side::Right)) {
            bc.right_bc->Apply(layer, mesh, Axis::Z, Side::Right);
        }
    }
}

const AxisBc& BoundaryManager::Get(const Axis axis) const {
    return axes_.at(static_cast<std::size_t>(axis));
}
