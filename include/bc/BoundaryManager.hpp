#ifndef BOUNDARYMANAGER_HPP
#define BOUNDARYMANAGER_HPP

#include <cstdint>
#include <memory>
#include <vector>
#include "parallel/HaloExchange.hpp"
#include "data/Variables.hpp"



class DataLayer;
class Mesh;
class BoundaryCondition;

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
 * @brief Manages halo update (stub for MPI) and physical boundary conditions.
 *
 * Contract:
 *  - UpdateHalo() handles internal subdomain interfaces (MPI later). For now: no-op.
 *  - ApplyPhysicalBc() applies only physical BC on global external boundaries,
 *    as indicated by Mesh boundary flags.
 */
class BoundaryManager final {
public:
    /** @brief Constructs boundary manager for three axes. */
    explicit BoundaryManager(std::shared_ptr<HaloExchange> halo_exchange = nullptr);

    /**
     * @brief Assign boundary conditions for an axis.
     * @param axis Axis (X/Y/Z).
     * @param left_bc Boundary at lower/min side.
     * @param right_bc Boundary at upper/max side.
     */
    void Set(Axis axis,
             std::shared_ptr<BoundaryCondition> left_bc,
             std::shared_ptr<BoundaryCondition> right_bc);

    /**
     * @brief Halo exchange/update for internal interfaces (MPI later).
     * @details For now: no-op.
     */
    void UpdateHalo(DataLayer& layer, const Mesh& mesh) const;

    /**
     * @brief Apply physical boundary conditions on global external boundaries.
     * @details Uses Mesh::IsGlobalBoundary(axis, side) to decide whether to apply.
     */
    void ApplyPhysicalBc(DataLayer& layer, const Mesh& mesh) const;

    /** @brief Get boundary pair for an axis. */
    [[nodiscard]] const AxisBc& Get(Axis axis) const;

private:
    std::vector<AxisBc> axes_;
    std::shared_ptr<HaloExchange> halo_exchange_;
};

#endif  // BOUNDARYMANAGER_HPP