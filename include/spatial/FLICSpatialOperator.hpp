#ifndef FLICSPATIALOPERATOR_HPP
#define FLICSPATIALOPERATOR_HPP

#include <memory>

#include "config/Settings.hpp"
#include "spatial/SpatialOperator.hpp"
#include "data/DataLayer.hpp"
#include "data/Mesh.hpp"
#include "data/Workspace.hpp"

class ArtificialViscosity;

/**
 * @class FLICSpatialOperator
 * @brief Fluid-in-Cell (FLIC) style split operator: Lagrangian (pressure) + Eulerian remap (advection).
 *
 * Algorithm (per ComputeRHS):
 *  1) Apply halo + physical BC on U.
 *  2) Convert U -> W on full padded domain.
 *  3) Compute Lagrangian part:
 *       d(rho)/dt = 0
 *       d(rho u)/dt = - dP/dx
 *       d(rho v)/dt = - dP/dy
 *       d(rho w)/dt = - dP/dz
 *       dE/dt       = - div(P * v)
 *  4) Build intermediate conservative state U* = U + dt * RHS_lag.
 *  5) Compute Eulerian remap (donor-cell upwind advection) on U*:
 *       dU/dt += - div( v * U )
 *  6) Optionally add artificial viscosity contribution to rhs.
 *
 * Notes:
 *  - This is a simple, robust baseline FLIC-like method, 3D-ready.
 *  - Requires ng >= 1.
 */
class FLICSpatialOperator final : public SpatialOperator {
public:
    FLICSpatialOperator(const Settings& settings,
                        std::shared_ptr<BoundaryManager> boundary_manager);

    void ComputeRHS(DataLayer& layer,
                    const Mesh& mesh,
                    Workspace& workspace,
                    double gamma,
                    double dt) const override;

private:
    std::shared_ptr<ArtificialViscosity> viscosity_;

    void ComputeLagrangianRhs(const Mesh& mesh,
                              const xt::xtensor<double, 4>& W,
                              xt::xtensor<double, 4>& rhs) const;

    void BuildUStar(const DataLayer& layer,
                    const Mesh& mesh,
                    const xt::xtensor<double, 4>& U,
                    const xt::xtensor<double, 4>& rhs_lag,
                    double dt,
                    xt::xtensor<double, 4>& U_star) const;

    void ComputeEulerianAdvectionRhs(const Mesh& mesh,
                                     const xt::xtensor<double, 4>& U_star,
                                     xt::xtensor<double, 4>& rhs) const;

    void AccumulateAxisAdvection(const Mesh& mesh,
                                 const xt::xtensor<double, 4>& U_star,
                                 xt::xtensor<double, 4>& rhs,
                                 Axis axis) const;

    [[nodiscard]] const xt::xtensor<double, 1>& InvMetric(const Mesh& mesh, Axis axis) const;

    static double CellVelocityComponent(const xt::xtensor<double, 4>& U_star,
                                        int i, int j, int k,
                                        Axis axis,
                                        double rho_floor = 1e-14);
};

#endif  // FLICSPATIALOPERATOR_HPP