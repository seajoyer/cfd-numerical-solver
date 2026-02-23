#ifndef SPATIALOPERATOR_HPP
#define SPATIALOPERATOR_HPP

#include <memory>

#include "data/DataLayer.hpp"
#include "data/Workspace.hpp"

class BoundaryManager;

/**
 * @class SpatialOperator
 * @brief Abstract semi-discrete FV operator in conservative form: dU/dt = L(U).
 *
 * Contract:
 *  - Must apply halo + physical BC internally (BoundaryManager).
 *  - Must write RHS into workspace.Rhs() (no allocations in hot path).
 *  - DataLayer stores only conservative U.
 */
class SpatialOperator {
public:
    explicit SpatialOperator(std::shared_ptr<BoundaryManager> boundary_manager) : boundary_manager_(boundary_manager) {}
    virtual ~SpatialOperator() = default;

    /**
     * @brief Computes RHS dU/dt = L(U) on padded grid (only core is used for update).
     * @param layer Conservative state owner (ghosts will be written by BC).
     * @param workspace Scratch buffers (W and rhs).
     * @param gamma Ratio of specific heats.
     * @param dt Local timestep.
     */
    virtual void ComputeRHS(DataLayer& layer, Workspace& workspace, double gamma, double dt) const = 0;

protected:
    std::shared_ptr<BoundaryManager> boundary_manager_;
};

#endif  // SPATIALOPERATOR_HPP
