#ifndef BOUNDARYMANAGER_HPP
#define BOUNDARYMANAGER_HPP

#include <cstdint>
#include <memory>
#include <vector>

#include "bc/BoundaryCondition.hpp"
#include "data/Variables.hpp"
#include "parallel/HaloExchange.hpp"

class DataLayer;
class Mesh;
class PressureVelocityState;
class PressureVelocityWorkspace;

/**
 * @struct AxisBc
 * @brief Boundary condition pair for one axis (left/right side).
 */
struct AxisBc final {
    std::shared_ptr<BoundaryCondition> left_bc;
    std::shared_ptr<BoundaryCondition> right_bc;
};

/**
 * @class BoundaryManager
 * @brief Manages halo update and physical boundary conditions.
 */
class BoundaryManager final {
public:
    explicit BoundaryManager(std::shared_ptr<HaloExchange> halo_exchange = nullptr);

    void Set(Axis axis,
             std::shared_ptr<BoundaryCondition> left_bc,
             std::shared_ptr<BoundaryCondition> right_bc);

    void UpdateHalo(DataLayer& layer, const Mesh& mesh) const;
    void UpdateHalo(PressureVelocityState& state, const Mesh& mesh) const;

    void ApplyPhysicalBc(DataLayer& layer, const Mesh& mesh) const;
    void ApplyPhysicalBc(PressureVelocityState& state, const Mesh& mesh) const;

    void ApplyPressureVelocityBoundary(PressureVelocityState& state,
                                       PressureVelocityWorkspace& workspace,
                                       const Mesh& mesh,
                                       PvAssemblyStage stage,
                                       bool steady,
                                       double dt,
                                       double nu) const;

    [[nodiscard]] const AxisBc& Get(Axis axis) const;

    /**
     * @brief Get BC object for one axis/side.
     * @return Raw pointer or nullptr if unset.
     */
    [[nodiscard]] const BoundaryCondition* GetCondition(Axis axis, Side side) const;

    /**
     * @brief Check whether a boundary side is periodic.
     */
    [[nodiscard]] bool IsPeriodic(Axis axis, Side side) const;

private:
    std::vector<AxisBc> axes_;
    std::shared_ptr<HaloExchange> halo_exchange_;
};

#endif  // BOUNDARYMANAGER_HPP
