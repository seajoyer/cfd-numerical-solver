#include "time/SSPRK3TimeIntegrator.hpp"

#include "spatial/SpatialOperator.hpp"
#include "solver/PositivityLimiter.hpp"

void SSPRK3TimeIntegrator::Advance(DataLayer& layer,
                                   Workspace& workspace,
                                   const double dt,
                                   const double gamma,
                                   const SpatialOperator& op) const {
    if (dt <= 0.0) {
        return;
    }

    workspace.ResizeFrom(layer);

    auto& U = layer.U();
    auto& rhs = workspace.Rhs();

    const int i0 = layer.GetCoreStartX();
    const int i1 = layer.GetCoreEndExclusiveX();
    const int j0 = layer.GetCoreStartY();
    const int j1 = layer.GetCoreEndExclusiveY();
    const int k0 = layer.GetCoreStartZ();
    const int k1 = layer.GetCoreEndExclusiveZ();

    // Save U^n (simple + safe).
    xt::xtensor<double, 4> U0 = xt::eval(U);

    // ---------- Stage 1: U = U0 + dt * L(U0) ----------
    op.ComputeRHS(layer, workspace, gamma, dt);

    xt::view(U, xt::all(),
             xt::range(i0, i1),
             xt::range(j0, j1),
             xt::range(k0, k1)) =
        xt::view(U0, xt::all(),
                 xt::range(i0, i1),
                 xt::range(j0, j1),
                 xt::range(k0, k1)) +
        dt * xt::view(rhs, xt::all(),
                      xt::range(i0, i1),
                      xt::range(j0, j1),
                      xt::range(k0, k1));

    // ---------- Stage 2: U = 3/4 U0 + 1/4 (U + dt*L(U)) ----------
    op.ComputeRHS(layer, workspace, gamma, dt);

    xt::view(U, xt::all(),
             xt::range(i0, i1),
             xt::range(j0, j1),
             xt::range(k0, k1)) =
        3.0 / 4.0 * xt::view(U0, xt::all(),
                             xt::range(i0, i1),
                             xt::range(j0, j1),
                             xt::range(k0, k1)) +
        1.0 / 4.0 * (xt::view(U, xt::all(),
                              xt::range(i0, i1),
                              xt::range(j0, j1),
                              xt::range(k0, k1)) +
            dt * xt::view(rhs, xt::all(),
                          xt::range(i0, i1),
                          xt::range(j0, j1),
                          xt::range(k0, k1)));

    // ---------- Stage 3: U = 1/3 U0 + 2/3 (U + dt*L(U)) ----------
    op.ComputeRHS(layer, workspace, gamma, dt);

    xt::view(U, xt::all(),
             xt::range(i0, i1),
             xt::range(j0, j1),
             xt::range(k0, k1)) =
        1.0 / 3.0 * xt::view(U0, xt::all(),
                             xt::range(i0, i1),
                             xt::range(j0, j1),
                             xt::range(k0, k1)) +
        2.0 / 3.0 * (xt::view(U, xt::all(),
                              xt::range(i0, i1),
                              xt::range(j0, j1),
                              xt::range(k0, k1)) +
            dt * xt::view(rhs, xt::all(),
                          xt::range(i0, i1),
                          xt::range(j0, j1),
                          xt::range(k0, k1)));

    PositivityLimiter::Apply(layer, gamma, rho_min_, p_min_);
}
