#ifndef GODUNOVSPATIALOPERATOR_HPP
#define GODUNOVSPATIALOPERATOR_HPP

#include <memory>

#include "config/Settings.hpp"
#include "spatial/SpatialOperator.hpp"
#include "data/Variables.hpp"

class Reconstruction;
class RiemannSolver;
class Face;

/**
 * @class GodunovSpatialOperator
 * @brief Face-based Godunov finite-volume operator for generic meshes.
 *
 * Workflow of one RHS evaluation:
 * - convert conservative cell state U(cell,var) to primitive cache W(cell,var)
 * - loop over mesh faces
 * - reconstruct owner/neighbor face states
 * - for boundary face build exterior state through BoundaryManager
 * - compute numerical flux through face normal
 * - accumulate flux contribution into cell-centered conservative RHS
 *
 * Notes:
 * - No ghost cells.
 * - No structured indexing.
 * - Boundary handling is performed face-by-face.
 */
class GodunovSpatialOperator final : public SpatialOperator {
public:
    GodunovSpatialOperator(const Settings& settings,
                           std::shared_ptr<BoundaryManager> boundary_manager,
                           const StateSynchronizer* synchronizer = nullptr);

    void ComputeRHS(const DataLayer& layer,
                    const Mesh& mesh,
                    Workspace& workspace,
                    double gamma,
                    double dt) const override;

private:
    std::shared_ptr<Reconstruction> reconstruction_;
    std::shared_ptr<RiemannSolver> riemann_solver_;

    void InitializeReconstruction(const Settings& settings);
    void InitializeRiemannSolver(const Settings& settings);

    void FillPrimitiveCache(const DataLayer& layer,
                            const Mesh& mesh,
                            Workspace& workspace,
                            double gamma) const;

    [[nodiscard]] FaceNormal BuildFaceNormal(const Face& face) const;

    void AccumulateInternalFace(const DataLayer& layer,
                                const Mesh& mesh,
                                const Face& face,
                                Workspace& workspace,
                                double gamma) const;

    void AccumulateBoundaryFace(const DataLayer& layer,
                                const Mesh& mesh,
                                const Face& face,
                                Workspace& workspace,
                                double gamma) const;

    void AccumulateFluxToOwner(const Mesh& mesh,
                               const Face& face,
                               const ConservativeCell& flux,
                               Workspace& workspace) const;

    void AccumulateFluxToNeighbor(const Mesh& mesh,
                                  const Face& face,
                                  const ConservativeCell& flux,
                                  Workspace& workspace) const;
};

#endif  // GODUNOVSPATIALOPERATOR_HPP
