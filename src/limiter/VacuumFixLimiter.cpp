#include "limiter/VacuumFixLimiter.hpp"

#include <cmath>

void VacuumFixLimiter::Apply(DataLayer& layer, double dx, const Settings& settings) const {
    const int total_size = layer.GetTotalSize();
    if (total_size <= 0) {
        return;
    }

    const int core_start = layer.GetCoreStart(0);
    const int core_end = layer.GetCoreEndExclusive(0);

    constexpr double rho_vac = 1.0e-2;
    constexpr double p_vac   = 1.0e-1;
    constexpr double u_cap   = 50.0;

    for (int i = core_start; i < core_end; ++i) {
        const double rho = layer.rho(i);
        const double P   = layer.P(i);

        if (rho <= 0.0) {
            continue;
        }

        const bool near_vacuum = (rho < rho_vac) && (P < p_vac);
        if (!near_vacuum) {
            continue;
        }

        // -------------------------------
        // 1) Kill specific volume
        // -------------------------------
        layer.V(i) = 0.0;

        // -------------------------------
        // 2) Kill specific internal energy
        // -------------------------------
        layer.U(i) = 0.0;

        // -------------------------------
        // 3) Optionally damp velocity
        // (keeps momentum untouched)
        // -------------------------------
        double u = layer.u(i);
        if (std::abs(u) > u_cap) {
            u = (u > 0.0 ? u_cap : -u_cap);
            layer.u(i) = u;
        }

        // -------------------------------
        // 4) Recompute mass safely
        // -------------------------------
        layer.m(i) = rho * dx;
    }
}
