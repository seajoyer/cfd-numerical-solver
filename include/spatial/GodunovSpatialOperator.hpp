// GodunovSpatialOperator.hpp
#ifndef GODUNOVSPATIALOPERATOR_HPP
#define GODUNOVSPATIALOPERATOR_HPP

#include <memory>

#include "spatial/SpatialOperator.hpp"
#include "config/Settings.hpp"

class Reconstruction;
class RiemannSolver;
class ArtificialViscosity;

/**
 * @class GodunovSpatialOperator
 * @brief Godunov FV operator: reconstruction + Riemann solver, summed per axis.
 *
 * Uses Workspace:
 *  - Converts U -> W once per RHS evaluation.
 *  - Computes flux divergence per axis and accumulates into workspace.Rhs().
 */
class GodunovSpatialOperator final : public SpatialOperator {
public:
    GodunovSpatialOperator(const Settings& settings,
                           std::shared_ptr<BoundaryManager> boundary_manager);

    void ComputeRHS(DataLayer& layer, Workspace& workspace, double gamma, double dt) const override;

private:
    std::shared_ptr<Reconstruction> reconstruction_;
    std::shared_ptr<RiemannSolver> riemann_solver_;
    std::shared_ptr<ArtificialViscosity> viscosity_;

    void InitializeReconstruction(const Settings& settings);
    void InitializeRiemannSolver(const Settings& settings);

    void AccumulateAxisFluxDivergence(const DataLayer& layer,
                                  const xt::xtensor<double, 4>& W,
                                  xt::xtensor<double, 4>& rhs,
                                  double gamma,
                                  Axis axis) const;

    [[nodiscard]] const xt::xtensor<double, 1>& InvMetric(const DataLayer& layer, Axis axis) const;
};

#endif  // GODUNOVSPATIALOPERATOR_HPP