#include "time/ForwardEulerTimeIntegrator.hpp"

#include "solver/EOS.hpp"
#include "solver/PositivityLimiter.hpp"

ForwardEulerTimeIntegrator::ForwardEulerTimeIntegrator() {
    rho_min_ = 1e-10;
    p_min_ = 1e-10;
}

void ForwardEulerTimeIntegrator::Advance(DataLayer& layer,
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

    xt::xarray<Conservative> U =
        xt::xarray<Conservative>::from_shape({static_cast<std::size_t>(total_size)});
    xt::xarray<Conservative> rhs =
        xt::xarray<Conservative>::from_shape({static_cast<std::size_t>(total_size)});

    for (int j = 0; j < total_size; ++j) {
        Primitive w = layer.GetPrimitive(j);
        U(j) = EOS::PrimToCons(w, gamma);
    }

    op.ComputeRHS(layer, dx, gamma, rhs);

    for (int j = core_start; j < core_end; ++j) {
        Conservative U_new = U(j);
        U_new += dt * rhs(j);

        PositivityLimiter::Apply(U_new, gamma, rho_min_, p_min_);
        StoreConservativeCell(U_new, j, dx, settings, layer);
    }
}
