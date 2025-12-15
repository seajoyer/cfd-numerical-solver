#include "spatial/ForwardEulerSpatialOperator.hpp"

#include <algorithm>

#include "data/DataLayer.hpp"
#include "data/Variables.hpp"
#include "viscosity/VNRArtificialViscosity.hpp"

ForwardEulerSpatialOperator::ForwardEulerSpatialOperator(const Settings& settings) {
    viscosity_ = nullptr;
    if (settings.viscosity) {
        viscosity_ = std::make_shared<VNRArtificialViscosity>(settings);
    }
}

void ForwardEulerSpatialOperator::ComputeRHS(const DataLayer& layer,
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
        for (int j = 0; j < total_size; ++j) rhs(j) = Conservative{};
        return;
    }

    rhs = xt::xarray<Conservative>::from_shape({static_cast<std::size_t>(total_size)});
    for (int j = 0; j < total_size; ++j) rhs(j) = Conservative{};

    xt::xarray<double> q_face;
    if (viscosity_) {
        viscosity_->ComputeInterfaceQ(layer, dx, q_face);
    } else {
        q_face = xt::zeros<double>({static_cast<std::size_t>(total_size - 1)});
    }

    xt::xarray<double> q_cell = xt::zeros<double>({static_cast<std::size_t>(total_size)});
    for (int j = 1; j < total_size - 1; ++j) {
        q_cell(j) = 0.5 * (q_face(j - 1) + q_face(j));
    }

    xt::xarray<Flux> fluxes =
        xt::xarray<Flux>::from_shape({static_cast<std::size_t>(total_size)});

    for (int j = 0; j < total_size; ++j) {
        Primitive w = layer.GetPrimitive(j);
        w.P += q_cell(j);
        fluxes(j) = EulerFlux(w, gamma);
    }

    const double inv_dx = 1.0 / dx;

    for (int j = core_start; j < core_end; ++j) {
        const Flux dF = Flux::Diff(fluxes(j + 1), fluxes(j));

        Conservative lj;
        lj = -1 * dF * inv_dx;

        rhs(j) = lj;
    }
}