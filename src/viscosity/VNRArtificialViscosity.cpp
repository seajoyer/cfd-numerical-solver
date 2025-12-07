#include "viscosity/VNRArtificialViscosity.hpp"

#include <algorithm>
#include <cmath>

#include "data/DataLayer.hpp"
#include "data/Variables.hpp"
#include "solver/EOS.hpp"

void VNRArtificialViscosity::ComputeInterfaceQ(const DataLayer& layer,
                                               double dx,
                                               xt::xarray<double>& q) const {
    (void)dx;

    const int total_size = layer.GetTotalSize();
    if (total_size <= 1) {
        q = xt::xarray<double>::from_shape({std::size_t(0)});
        return;
    }

    const int n_interfaces = total_size - 1;
    q = xt::zeros<double>({static_cast<std::size_t>(n_interfaces)});

    const double gamma = settings_.gamma;

    const int core_start = layer.GetCoreStart(0);
    const int core_end = layer.GetCoreEndExclusive(0);

    const int i_min = std::max(0, core_start - 1);
    const int i_max = std::min(n_interfaces - 1, core_end - 1);

    for (int i = i_min; i <= i_max; ++i) {
        const int jl = i;
        const int jr = i + 1;

        const double ul = layer.u(jl);
        const double ur = layer.u(jr);
        const double du = ur - ul;

        // Only act in compression: Δu < 0
        if (du >= 0.0) {
            continue;
        }

        const double rho_l = layer.rho(jl);
        const double rho_r = layer.rho(jr);
        const double rho_bar = 0.5 * (rho_l + rho_r);

        if (rho_bar <= 0.0) {
            continue;
        }

        const Primitive wl = layer.GetPrimitive(jl);
        const Primitive wr = layer.GetPrimitive(jr);

        const double c_l = EOS::SoundSpeed(wl, gamma);
        const double c_r = EOS::SoundSpeed(wr, gamma);
        const double c_bar = 0.5 * (c_l + c_r);

        const double du_abs = std::fabs(du);

        double q_i = 0.0;
        q_i += C2_ * rho_bar * du * du;
        q_i += C1_ * rho_bar * c_bar * du_abs;

        if (q_i > 0.0) {
            q(i) = q_i;
        }
    }
}