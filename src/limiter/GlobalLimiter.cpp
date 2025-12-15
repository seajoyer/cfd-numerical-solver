#include "limiter/GlobalLimiter.hpp"

#include <algorithm>
#include <cmath>

#include "solver/EOS.hpp"
#include "solver/PositivityLimiter.hpp"

void GlobalLimiter::Apply(DataLayer& layer, const double dx, const Settings& settings) const {
    const int total_size = layer.GetTotalSize();
    if (total_size < 5 || dx <= 0.0) {
        return;
    }

    const int core_start = layer.GetCoreStart(0);
    const int core_end = layer.GetCoreEndExclusive(0);
    const int n_core = core_end - core_start;
    if (n_core < 3) {
        return;
    }

    const double gamma = settings.gamma;

    const double rho_min = 1e-10;
    const double p_min = 1e-10;

    for (int i = core_start + 1; i < core_end - 1; ++i) {
        const double em = layer.e(i - 1);
        const double ec = layer.e(i);
        const double ep = layer.e(i + 1);

        const double dm = ec - em;
        const double dp = ep - ec;

        const bool is_extremum = dm * dp < 0.0;
        if (!is_extremum) {
            continue;
        }

        const double target = 0.5 * (em + ep);
        const double e_new_raw = ec + strength_ * (target - ec);

        Conservative uc;
        uc.rho = layer.rho(i);
        uc.rhoU = layer.p(i);
        uc.E = e_new_raw;

        PositivityLimiter::Apply(uc, gamma, rho_min, p_min);

        const double rho = uc.rho;
        const double uvel = rho > 0.0 ? uc.rhoU / rho : 0.0;

        layer.e(i) = uc.E;
        layer.p(i) = uc.rhoU;
        layer.u(i) = uvel;
        layer.rho(i) = rho;

        const double P = EOS::Pressure(uc, gamma);
        layer.P(i) = P;

        const double kinetic = 0.5 * rho * uvel * uvel;
        const double Eint = uc.E - kinetic;
        layer.U(i) = rho > 0.0 ? Eint / rho : 0.0;

        layer.V(i) = rho > 0.0 ? 1.0 / rho : 0.0;
        layer.m(i) = rho * dx;
    }
}
