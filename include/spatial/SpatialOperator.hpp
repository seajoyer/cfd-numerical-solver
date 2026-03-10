#ifndef SPATIALOPERATOR_HPP
#define SPATIALOPERATOR_HPP

#include <memory>

class DataLayer;
class Mesh;
class Workspace;
class BoundaryManager;

/**
 * @class SpatialOperator
 * @brief Abstract semi-discrete FV operator in conservative form: dU/dt = L(U).
 *
 * Contract:
 *  - Must apply halo + physical BC internally through BoundaryManager.
 *  - Must write RHS into workspace.Rhs() (no allocations in hot path).
 *  - DataLayer stores only conservative U.
 *  - Mesh stores geometry, metrics, domain ranges, and cell classification.
 */
class SpatialOperator {
public:
    explicit SpatialOperator(std::shared_ptr<BoundaryManager> boundary_manager) : boundary_manager_(std::move(boundary_manager)) {}
    virtual ~SpatialOperator() = default;

    /**
     * @brief Computes RHS dU/dt = L(U) on padded grid.
     * @param layer Conservative state owner (ghosts may be written by BC).
     * @param mesh Structured mesh with geometry, metrics, and cell types.
     * @param workspace Scratch buffers (W and rhs).
     * @param gamma Ratio of specific heats.
     * @param dt Local timestep.
     */
    virtual void ComputeRHS(DataLayer& layer,
                            const Mesh& mesh,
                            Workspace& workspace,
                            double gamma,
                            double dt) const = 0;

protected:
    std::shared_ptr<BoundaryManager> boundary_manager_;
};

#endif  // SPATIALOPERATOR_HPP