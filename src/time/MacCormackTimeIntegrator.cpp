#include "time/MacCormackTimeIntegrator.hpp"
#include "data/DataLayer.hpp"
#include "data/Variables.hpp"
#include "solver/EOS.hpp"
#include "solver/PositivityLimiter.hpp"


MacCormackTimeIntegrator::MacCormackTimeIntegrator() {
    rho_min_ = 1e-10;
    p_min_ = 1e-10;

    forward_operator_ = std::make_unique<ForwardEulerSpatialOperator>();
    backward_operator_ = std::make_unique<BackwardEulerSpatialOperator>();
}

void MacCormackTimeIntegrator::Advance(DataLayer& layer,
                                       double dt,
                                       double dx,
                                       const Settings& settings,
                                       const SpatialOperator& op) const {
    (void)op;

    const int total_size = layer.GetTotalSize();
    if (total_size <= 0 || dt <= 0.0 || dx <= 0.0) {
        return;
    }

    if (!forward_operator_ || !backward_operator_) {
        return;
    }

    const int core_start = layer.GetCoreStart(0);
    const int core_end = layer.GetCoreEndExclusive(0);
    const int n_core = core_end - core_start;

    if (n_core < 2 || total_size < 3) {
        return;
    }

    const double gamma = settings.gamma;

    xt::xarray<Conservative> U0 =
        xt::xarray<Conservative>::from_shape({static_cast<std::size_t>(total_size)});
    xt::xarray<Conservative> U_pred =
        xt::xarray<Conservative>::from_shape({static_cast<std::size_t>(total_size)});
    xt::xarray<Conservative> rhs =
        xt::xarray<Conservative>::from_shape({static_cast<std::size_t>(total_size)});

    for (int j = 0; j < total_size; ++j) {
        Primitive w = layer.GetPrimitive(j);
        U0(j) = EOS::PrimToCons(w, gamma);
    }

    forward_operator_->ComputeRHS(layer, dx, gamma, rhs);

    for (int j = 0; j < total_size; ++j) {
        U_pred(j) = U0(j);
        U_pred(j) += dt * rhs(j);
    }

    for (int j = core_start; j < core_end; ++j) {
        StoreConservativeCell(U_pred(j), j, dx, settings, layer);
    }

    backward_operator_->ComputeRHS(layer, dx, gamma, rhs);

    for (int j = core_start; j < core_end; ++j) {
        Conservative U_new = U0(j);
        U_new += U_pred(j);
        U_new += dt * rhs(j);
        U_new *= 0.5;

        PositivityLimiter::Apply(U_new, gamma, rho_min_, p_min_);
        StoreConservativeCell(U_new, j, dx, settings, layer);
    }
}