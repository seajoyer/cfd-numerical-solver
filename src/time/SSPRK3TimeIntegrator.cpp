#include "time/SSPRK3TimeIntegrator.hpp"

#include "solver/EOS.hpp"
#include "solver/PositivityLimiter.hpp"

SSPRK3TimeIntegrator::SSPRK3TimeIntegrator() {
    rho_min_ = 1e-10;
    p_min_ = 1e-10;
}

void SSPRK3TimeIntegrator::Advance(DataLayer& layer,
                                   double dt,
                                   double dx,
                                   const Settings& settings,
                                   const SpatialOperator& op) const {
    const int total_size = layer.GetTotalSize();
    if (total_size <= 0 || dt <= 0.0 || dx <= 0.0) {
        return;
    }

    const int core_start = layer.GetCoreStart(0);
    const int core_end = layer.GetCoreEndExclusive(0);
    const double gamma = settings.gamma;

    xt::xarray<Conservative> U0 =
        xt::xarray<Conservative>::from_shape({static_cast<std::size_t>(total_size)});
    xt::xarray<Conservative> U1 =
        xt::xarray<Conservative>::from_shape({static_cast<std::size_t>(total_size)});
    xt::xarray<Conservative> U2 =
        xt::xarray<Conservative>::from_shape({static_cast<std::size_t>(total_size)});
    xt::xarray<Conservative> rhs =
        xt::xarray<Conservative>::from_shape({static_cast<std::size_t>(total_size)});

    for (int j = 0; j < total_size; ++j) {
        Primitive w = layer.GetPrimitive(j);
        U0(j) = EOS::PrimToCons(w, gamma);
    }

    op.ComputeRHS(layer, dx, gamma, rhs);

    for (int j = 0; j < total_size; ++j) {
        Primitive w = layer.GetPrimitive(j);
        U0(j) = EOS::PrimToCons(w, gamma);
        U1(j) = U0(j);
        U1(j) += dt * rhs(j);
    }

    for (int j = core_start; j < core_end; ++j) {
        StoreConservativeCell(U1(j), j, dx, settings, layer);
    }

    op.ComputeRHS(layer, dx, gamma, rhs);

    for (int j = 0; j < total_size; ++j) {
        Conservative tmp = U1(j);
        tmp += dt * rhs(j);

        U2(j) = 3.0 / 4.0 * U0(j) + 1.0 / 4.0 * tmp;
    }

    for (int j = core_start; j < core_end; ++j) {
        StoreConservativeCell(U2(j), j, dx, settings, layer);
    }

    op.ComputeRHS(layer, dx, gamma, rhs);

    for (int j = core_start; j < core_end; ++j) {
        Conservative tmp = U2(j);
        tmp += dt * rhs(j);

        Conservative U_new = 1.0 / 3.0 * U0(j) + 2.0 / 3.0 * tmp;

        PositivityLimiter::Apply(U_new, gamma, rho_min_, p_min_);
        StoreConservativeCell(U_new, j, dx, settings, layer);
    }
}