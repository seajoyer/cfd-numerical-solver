#include "time/MacCormackTimeIntegrator.hpp"

#include "solver/PositivityLimiter.hpp"

MacCormackTimeIntegrator::MacCormackTimeIntegrator(const Settings& settings,
                                                   std::shared_ptr<BoundaryManager> boundary_manager)
    : forward_op_(settings, boundary_manager),
      backward_op_(settings, boundary_manager) {}

void MacCormackTimeIntegrator::Advance(DataLayer& layer,
                                       Workspace& workspace,
                                       const double dt,
                                       const double gamma,
                                       const SpatialOperator& op) const {
    (void)op;

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

    // Save U^n (simple and safe).
    xt::xtensor<double, 4> U0 = xt::eval(U);

    // ---------- Predictor: U = U0 + dt * L_fwd(U0) ----------
    forward_op_.ComputeRHS(layer, workspace, gamma, dt);

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

    // Save U* (predictor result) for the corrector formula.
    xt::xtensor<double, 4> Upred = xt::eval(U);

    // ---------- Corrector RHS: L_bwd(U*) computed at current layer.U() = U* ----------
    backward_op_.ComputeRHS(layer, workspace, gamma, dt);

    // U^{n+1} = 0.5 * (U0 + Upred + dt*rhs)
    xt::view(U, xt::all(),
             xt::range(i0, i1),
             xt::range(j0, j1),
             xt::range(k0, k1)) =
        0.5 * (xt::view(U0, xt::all(),
                        xt::range(i0, i1),
                        xt::range(j0, j1),
                        xt::range(k0, k1)) +
            xt::view(Upred, xt::all(),
                     xt::range(i0, i1),
                     xt::range(j0, j1),
                     xt::range(k0, k1)) +
            dt * xt::view(rhs, xt::all(),
                          xt::range(i0, i1),
                          xt::range(j0, j1),
                          xt::range(k0, k1)));

    PositivityLimiter::Apply(layer, gamma, rho_min_, p_min_);
}
