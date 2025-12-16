#include "spatial/ForwardEulerSpatialOperator.hpp"

#include <algorithm>
#include <xtensor.hpp>

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
    const int N = layer.GetTotalSize();
    const int cs = layer.GetCoreStart(0);
    const int ce = layer.GetCoreEndExclusive(0);
    const int n_core = ce - cs;

    if (N < 3 || n_core < 2 || dx <= 0.0) {
        rhs = xt::xarray<Conservative>::from_shape( {static_cast<std::size_t>(std::max(N, 0))});
        for (int j = 0; j < N; ++j)
            rhs(j) = Conservative{};
        return;
    }

    xt::xarray<double> q;
    if (viscosity_) {
        viscosity_->ComputeInterfaceQ(layer, dx, q);
    } else {
        q = xt::zeros<double>({(size_t)(N - 1)});
    }

    xt::xarray<Flux> F =
        xt::xarray<Flux>::from_shape({(size_t)(N - 1)});

    for (int i = 0; i < N - 1; ++i)
    {
        Primitive wL = layer.GetPrimitive(i);
        if (q(i) > 0.0) {
            wL.P += q(i);
        }
        F(i) = EulerFlux(wL, gamma);
    }

    const double inv_dx = 1.0 / dx;
    for (int j = cs; j < ce; ++j) {
        rhs(j) = (F(j - 1) - F(j)) * inv_dx;
    }
}
