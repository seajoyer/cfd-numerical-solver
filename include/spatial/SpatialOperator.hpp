#ifndef SPATIALOPERATOR_HPP
#define SPATIALOPERATOR_HPP

#include <memory>

class DataLayer;
class Mesh;
class Workspace;
class BoundaryManager;

/**
 * @class SpatialOperator
 * @brief Abstract semi-discrete finite-volume operator in conservative form.
 *
 * Computes:
 *   dU/dt = L(U)
 *
 * Contract:
 * - Works on generic face-based meshes.
 * - Does not use ghost cells or structured indexing.
 * - Does not modify mesh.
 * - Writes conservative RHS into Workspace::Rhs().
 * - Uses BoundaryManager to build exterior states on boundary faces.
 */
class SpatialOperator {
public:
    explicit SpatialOperator(std::shared_ptr<BoundaryManager> boundary_manager)
        : boundary_manager_(std::move(boundary_manager)) {}

    virtual ~SpatialOperator() = default;

    /**
     * @brief Compute conservative RHS for the current state.
     *
     * @param layer Conservative solution storage.
     * @param mesh Mesh with cells, faces, geometry, and connectivity.
     * @param workspace Reusable scratch buffers and output RHS storage.
     * @param gamma Ratio of specific heats.
     * @param dt Current timestep size, available for optional model terms.
     */
    virtual void ComputeRHS(const DataLayer& layer,
                            const Mesh& mesh,
                            Workspace& workspace,
                            double gamma,
                            double dt) const = 0;

protected:
    std::shared_ptr<BoundaryManager> boundary_manager_;
};

#endif  // SPATIALOPERATOR_HPP
