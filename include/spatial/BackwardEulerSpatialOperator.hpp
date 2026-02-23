#ifndef BACKWARDEULERSPATIALOPERATOR_HPP
#define BACKWARDEULERSPATIALOPERATOR_HPP

#include <memory>

#include "spatial/SpatialOperator.hpp"
#include "config/Settings.hpp"
#include "data/Variables.hpp"

class ArtificialViscosity;

/**
 * @class BackwardEulerSpatialOperator
 * @brief Backward-difference flux-divergence operator for MacCormack corrector.
 *
 * Computes per axis:
 *   L(U) = - (F_i - F_{i-1}) / dx
 * where F_i is the physical Euler flux evaluated at cell centers (from primitives).
 *
 * Contract:
 *  - Applies halo + physical BC internally (BoundaryManager).
 *  - Converts U -> W once per call using Workspace.
 *  - Writes RHS into workspace.Rhs().
 *
 * Artificial viscosity is not wired yet (viscosity_ stays null).
 */
class BackwardEulerSpatialOperator final : public SpatialOperator {
public:
    explicit BackwardEulerSpatialOperator(const Settings& settings,
                                          std::shared_ptr<BoundaryManager> boundary_manager);

    void ComputeRHS(DataLayer& layer, Workspace& workspace, double gamma, double dt) const override;

private:
    std::shared_ptr<ArtificialViscosity> viscosity_; // TODO: wire later

    void AccumulateAxis(const DataLayer& layer,
                        const xt::xtensor<double, 4>& W,
                        xt::xtensor<double, 4>& rhs,
                        double gamma,
                        Axis axis) const;

    [[nodiscard]] const xt::xtensor<double, 1>& InvMetric(const DataLayer& layer, Axis axis) const;

    static PrimitiveCell LoadPrimitive(const xt::xtensor<double, 4>& W, int i, int j, int k);
};

#endif  // BACKWARDEULERSPATIALOPERATOR_HPP