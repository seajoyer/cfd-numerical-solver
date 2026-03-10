#include "time/SSPRK3TimeIntegrator.hpp"

#include "data/DataLayer.hpp"
#include "data/Mesh.hpp"
#include "data/Workspace.hpp"
#include "solver/PositivityLimiter.hpp"
#include "spatial/SpatialOperator.hpp"

void SSPRK3TimeIntegrator::Advance(DataLayer& layer,
                                   const Mesh& mesh,
                                   Workspace& workspace,
                                   const double dt,
                                   const double gamma,
                                   const SpatialOperator& op) const {
    if (dt <= 0.0) {
        return;
    }

    workspace.ResizeFrom(mesh);

    auto& U = layer.U();
    auto& rhs = workspace.Rhs();

    const int i0 = mesh.GetCoreStartX();
    const int i1 = mesh.GetCoreEndExclusiveX();
    const int j0 = mesh.GetCoreStartY();
    const int j1 = mesh.GetCoreEndExclusiveY();
    const int k0 = mesh.GetCoreStartZ();
    const int k1 = mesh.GetCoreEndExclusiveZ();

    xt::xtensor<double, 4> U0 = xt::eval(U);

    op.ComputeRHS(layer, mesh, workspace, gamma, dt);

    for (int k = k0; k < k1; ++k) {
        for (int j = j0; j < j1; ++j) {
            for (int i = i0; i < i1; ++i) {
                if (!mesh.IsFluidCell(i, j, k)) {
                    continue;
                }

                for (std::size_t v = 0; v < DataLayer::k_nvar; ++v) {
                    U(v, i, j, k) = U0(v, i, j, k) + dt * rhs(v, i, j, k);
                }
            }
        }
    }

    op.ComputeRHS(layer, mesh, workspace, gamma, dt);

    for (int k = k0; k < k1; ++k) {
        for (int j = j0; j < j1; ++j) {
            for (int i = i0; i < i1; ++i) {
                if (!mesh.IsFluidCell(i, j, k)) {
                    continue;
                }

                for (std::size_t v = 0; v < DataLayer::k_nvar; ++v) {
                    U(v, i, j, k) =
                        3.0 / 4.0 * U0(v, i, j, k) +
                        1.0 / 4.0 * (U(v, i, j, k) + dt * rhs(v, i, j, k));
                }
            }
        }
    }

    op.ComputeRHS(layer, mesh, workspace, gamma, dt);

    for (int k = k0; k < k1; ++k) {
        for (int j = j0; j < j1; ++j) {
            for (int i = i0; i < i1; ++i) {
                if (!mesh.IsFluidCell(i, j, k)) {
                    continue;
                }

                for (std::size_t v = 0; v < DataLayer::k_nvar; ++v) {
                    U(v, i, j, k) =
                        1.0 / 3.0 * U0(v, i, j, k) +
                        2.0 / 3.0 * (U(v, i, j, k) + dt * rhs(v, i, j, k));
                }
            }
        }
    }

    PositivityLimiter::Apply(layer, mesh, gamma, rho_min_, p_min_);
}