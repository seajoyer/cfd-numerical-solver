#include "bc/BoundaryManager.hpp"

#include <utility>

#include "data/DataLayer.hpp"
#include "data/Mesh.hpp"
#include "data/PressureVelocityState.hpp"
#include "data/PressureVelocityWorkspace.hpp"

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

void BoundaryManager::UpdateHalo(PressureVelocityState& state, const Mesh& mesh) const {
    if (!halo_exchange_) {
        return;
    }

    halo_exchange_->Exchange(state, mesh);
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

void BoundaryManager::ApplyPhysicalBc(PressureVelocityState& state, const Mesh& mesh) const {
    const int dim = mesh.GetDim();

    {
        const AxisBc& bc = axes_[static_cast<std::size_t>(Axis::X)];
        if (bc.left_bc && mesh.IsGlobalBoundary(Axis::X, Side::Left)) {
            bc.left_bc->Apply(state, mesh, Axis::X, Side::Left);
        }
        if (bc.right_bc && mesh.IsGlobalBoundary(Axis::X, Side::Right)) {
            bc.right_bc->Apply(state, mesh, Axis::X, Side::Right);
        }
    }

    if (dim >= 2) {
        const AxisBc& bc = axes_[static_cast<std::size_t>(Axis::Y)];
        if (bc.left_bc && mesh.IsGlobalBoundary(Axis::Y, Side::Left)) {
            bc.left_bc->Apply(state, mesh, Axis::Y, Side::Left);
        }
        if (bc.right_bc && mesh.IsGlobalBoundary(Axis::Y, Side::Right)) {
            bc.right_bc->Apply(state, mesh, Axis::Y, Side::Right);
        }
    }

    if (dim >= 3) {
        const AxisBc& bc = axes_[static_cast<std::size_t>(Axis::Z)];
        if (bc.left_bc && mesh.IsGlobalBoundary(Axis::Z, Side::Left)) {
            bc.left_bc->Apply(state, mesh, Axis::Z, Side::Left);
        }
        if (bc.right_bc && mesh.IsGlobalBoundary(Axis::Z, Side::Right)) {
            bc.right_bc->Apply(state, mesh, Axis::Z, Side::Right);
        }
    }
}

void BoundaryManager::ApplyPressureVelocityBoundary(PressureVelocityState& state,
                                                    PressureVelocityWorkspace& workspace,
                                                    const Mesh& mesh,
                                                    const PvAssemblyStage stage,
                                                    const bool steady,
                                                    const double dt,
                                                    const double nu) const {
    const int dim = mesh.GetDim();

    {
        const AxisBc& bc = axes_[static_cast<std::size_t>(Axis::X)];
        if (bc.left_bc && mesh.IsGlobalBoundary(Axis::X, Side::Left)) {
            bc.left_bc->ApplyPressureVelocityBoundary(
                state, workspace, mesh, Axis::X, Side::Left, stage, steady, dt, nu);
        }
        if (bc.right_bc && mesh.IsGlobalBoundary(Axis::X, Side::Right)) {
            bc.right_bc->ApplyPressureVelocityBoundary(
                state, workspace, mesh, Axis::X, Side::Right, stage, steady, dt, nu);
        }
    }

    if (dim >= 2) {
        const AxisBc& bc = axes_[static_cast<std::size_t>(Axis::Y)];
        if (bc.left_bc && mesh.IsGlobalBoundary(Axis::Y, Side::Left)) {
            bc.left_bc->ApplyPressureVelocityBoundary(
                state, workspace, mesh, Axis::Y, Side::Left, stage, steady, dt, nu);
        }
        if (bc.right_bc && mesh.IsGlobalBoundary(Axis::Y, Side::Right)) {
            bc.right_bc->ApplyPressureVelocityBoundary(
                state, workspace, mesh, Axis::Y, Side::Right, stage, steady, dt, nu);
        }
    }

    if (dim >= 3) {
        const AxisBc& bc = axes_[static_cast<std::size_t>(Axis::Z)];
        if (bc.left_bc && mesh.IsGlobalBoundary(Axis::Z, Side::Left)) {
            bc.left_bc->ApplyPressureVelocityBoundary(
                state, workspace, mesh, Axis::Z, Side::Left, stage, steady, dt, nu);
        }
        if (bc.right_bc && mesh.IsGlobalBoundary(Axis::Z, Side::Right)) {
            bc.right_bc->ApplyPressureVelocityBoundary(
                state, workspace, mesh, Axis::Z, Side::Right, stage, steady, dt, nu);
        }
    }
}

const AxisBc& BoundaryManager::Get(const Axis axis) const {
    return axes_.at(static_cast<std::size_t>(axis));
}

const BoundaryCondition* BoundaryManager::GetCondition(const Axis axis, const Side side) const {
    const AxisBc& bc = axes_.at(static_cast<std::size_t>(axis));
    const auto& ptr = side == Side::Left ? bc.left_bc : bc.right_bc;
    return ptr ? ptr.get() : nullptr;
}

bool BoundaryManager::IsPeriodic(const Axis axis, const Side side) const {
    const BoundaryCondition* bc = GetCondition(axis, side);
    return bc && bc->IsPeriodic();
}
