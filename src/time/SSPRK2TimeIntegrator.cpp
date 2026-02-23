#include "time/SSPRK2TimeIntegrator.hpp"

#include <stdexcept>

#include "spatial/SpatialOperator.hpp"
#include "solver/PositivityLimiter.hpp"

void SSPRK2TimeIntegrator::Advance(DataLayer& layer,
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

    // Save U^n on core (we keep a full tensor for simplicity; could also store only core view).
    xt::xtensor<double, 4> U0 = xt::eval(U); // includes ghosts; safe and simple

    // -------- Stage 1: U1 = U0 + dt * L(U0) --------
    op.ComputeRHS(layer, workspace, gamma, dt);

    // U(core) = U0(core) + dt * rhs(core)
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

    // -------- Stage 2: U^{n+1} = 0.5*U0 + 0.5*(U + dt*L(U)) --------
    op.ComputeRHS(layer, workspace, gamma, dt);

    xt::view(U, xt::all(),
             xt::range(i0, i1),
             xt::range(j0, j1),
             xt::range(k0, k1)) =
        0.5 * xt::view(U0, xt::all(),
                       xt::range(i0, i1),
                       xt::range(j0, j1),
                       xt::range(k0, k1)) +
        0.5 * (xt::view(U, xt::all(),
                        xt::range(i0, i1),
                        xt::range(j0, j1),
                        xt::range(k0, k1)) +
               dt * xt::view(rhs, xt::all(),
                             xt::range(i0, i1),
                             xt::range(j0, j1),
                             xt::range(k0, k1)));

    PositivityLimiter::Apply(layer, gamma, rho_min_, p_min_);
}