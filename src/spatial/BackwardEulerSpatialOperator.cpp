#include "spatial/BackwardEulerSpatialOperator.hpp"

#include "data/DataLayer.hpp"
#include "data/Variables.hpp"

BackwardEulerSpatialOperator::BackwardEulerSpatialOperator(const Settings& settings) {
    viscosity_ = nullptr;
    if (settings.viscosity) {
        viscosity_ = std::make_shared<VNRArtificialViscosity>(settings);
    }
}

void BackwardEulerSpatialOperator::ComputeRHS(const DataLayer& layer,
                                              double dx,
                                              double gamma,
                                              xt::xarray<Conservative>& rhs) const {
    const int total_size = layer.GetTotalSize();
    const int core_start = layer.GetCoreStart(0);
    const int core_end = layer.GetCoreEndExclusive(0);
    const int n_core = core_end - core_start;

    if (total_size < 3 || n_core < 2 || dx <= 0.0) {
        rhs = xt::xarray<Conservative>::from_shape(
            {static_cast<std::size_t>(std::max(total_size, 0))});
        for (int j = 0; j < total_size; ++j) {
            rhs(j) = Conservative{};
        }
        return;
    }

    rhs = xt::xarray<Conservative>::from_shape(
        {static_cast<std::size_t>(total_size)});
    for (int j = 0; j < total_size; ++j) {
        rhs(j) = Conservative{};
    }

    xt::xarray<Flux> fluxes =
        xt::xarray<Flux>::from_shape({static_cast<std::size_t>(total_size)});

    xt::xarray<double> q;
    if (viscosity_) {
        viscosity_->ComputeInterfaceQ(layer, dx, q);
    } else {
        q = xt::zeros<double>({static_cast<std::size_t>(total_size - 1)});
    }

    for (int j = 0; j < total_size; ++j) {
        Primitive w = layer.GetPrimitive(j);
        w.P += q(j - 1);
        fluxes(j) = EulerFlux(w, gamma);
    }

    const double inv_dx = 1.0 / dx;

    for (int j = core_start; j < core_end; ++j) {
        const int jm1 = j - 1;

        const Flux dF = Flux::Diff(fluxes(j), fluxes(jm1));

        Conservative lj;
        lj -= dF * inv_dx;

        rhs(j) = lj;
    }
}
