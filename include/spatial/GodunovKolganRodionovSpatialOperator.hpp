#ifndef GODUNOVKOLGANRODIONOVSPATIALOPERATOR_HPP
#define GODUNOVKOLGANRODIONOVSPATIALOPERATOR_HPP

#include <memory>
#include <xtensor.hpp>

#include "config/Settings.hpp"
#include "data/Variables.hpp"
#include "spatial/SpatialOperator.hpp"

class Reconstruction;
class RiemannSolver;
class Face;

/**
 * @class GodunovKolganRodionovSpatialOperator
 * @brief Face-based second-order Godunov-Kolgan-Rodionov type operator.
 *
 * Predictor-corrector workflow:
 * - build primitive cache from current conservative state
 * - compute first-order predictor RHS with piecewise-constant face states
 * - build half-step conservative state U_half
 * - rebuild primitive cache from U_half
 * - compute final RHS using selected reconstruction (typically P1)
 *
 * Notes:
 * - Works on generic face-based meshes.
 * - No ghost cells.
 * - Boundary handling is done face-by-face through BoundaryManager.
 */
class GodunovKolganRodionovSpatialOperator final : public SpatialOperator {
public:
    GodunovKolganRodionovSpatialOperator(const Settings& settings,
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

    void FillPrimitiveCacheFromConservative(const xt::xtensor<double, 2>& U,
                                            const Mesh& mesh,
                                            Workspace& workspace,
                                            double gamma) const;

    [[nodiscard]] PrimitiveCell LoadCellPrimitive
    (
        const Workspace& workspace,
        std::size_t cell_id
    )
    const;

    [[nodiscard]] FaceNormal BuildFaceNormal(const Face& face) const;

    void ComputePredictorRhs(const DataLayer& layer,
                             const Mesh& mesh,
                             Workspace& workspace,
                             double gamma) const;

    void AccumulatePredictorInternalFace(const DataLayer& layer,
                                         const Mesh& mesh,
                                         const Face& face,
                                         Workspace& workspace,
                                         double gamma) const;

    void AccumulatePredictorBoundaryFace(const DataLayer& layer,
                                         const Mesh& mesh,
                                         const Face& face,
                                         Workspace& workspace,
                                         double gamma) const;

    void ComputeFinalRhs(const DataLayer& layer,
                         const Mesh& mesh,
                         Workspace& workspace,
                         double gamma) const;

    void AccumulateFinalInternalFace(const DataLayer& layer,
                                     const Mesh& mesh,
                                     const Face& face,
                                     Workspace& workspace,
                                     double gamma) const;

    void AccumulateFinalBoundaryFace(const DataLayer& layer,
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

#endif  // GODUNOVKOLGANRODIONOVSPATIALOPERATOR_HPP
