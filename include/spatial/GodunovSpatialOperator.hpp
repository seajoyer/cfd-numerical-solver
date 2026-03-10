#ifndef GODUNOVSPATIALOPERATOR_HPP
#define GODUNOVSPATIALOPERATOR_HPP

#include <memory>

#include "config/Settings.hpp"
#include "data/DataLayer.hpp"
#include "data/Mesh.hpp"
#include "data/Workspace.hpp"
#include "spatial/SpatialOperator.hpp"

class Reconstruction;
class RiemannSolver;
class ArtificialViscosity;
class InternalBoundaryCondition;

/**
 * @class GodunovSpatialOperator
 * @brief Godunov FV operator: reconstruction + Riemann solver, summed per axis.
 *
 * Uses Workspace:
 *  - Converts U -> W once per RHS evaluation.
 *  - Computes flux divergence per axis and accumulates into workspace.Rhs().
 *
 * Embedded internal boundaries are handled on fluid-solid faces via
 * InternalBoundaryCondition.
 */
class GodunovSpatialOperator final : public SpatialOperator {
public:
    GodunovSpatialOperator(const Settings& settings,
                           std::shared_ptr<BoundaryManager> boundary_manager);

    void ComputeRHS(DataLayer& layer,
                    const Mesh& mesh,
                    Workspace& workspace,
                    double gamma,
                    double dt) const override;

private:
    std::shared_ptr<Reconstruction> reconstruction_;
    std::shared_ptr<RiemannSolver> riemann_solver_;
    std::shared_ptr<ArtificialViscosity> viscosity_;
    std::shared_ptr<InternalBoundaryCondition> internal_boundary_condition_;

    void InitializeReconstruction(const Settings& settings);
    void InitializeRiemannSolver(const Settings& settings);

    void AccumulateAxisFluxDivergence(const DataLayer& layer,
                                      const Mesh& mesh,
                                      const xt::xtensor<double, 4>& W,
                                      xt::xtensor<double, 4>& rhs,
                                      double gamma,
                                      Axis axis) const;

    [[nodiscard]] const xt::xtensor<double, 1>& InvMetric(const Mesh& mesh, Axis axis) const;
};

#endif  // GODUNOVSPATIALOPERATOR_HPP