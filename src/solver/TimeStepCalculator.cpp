#include "solver/TimeStepCalculator.hpp"

#include <cmath>
#include <algorithm>

auto TimeStepCalculator::ComputeDt(const DataLayer& layer, const double dx,
                                   const double cfl, const double gamma) -> double {
    // 1D case (backward compatible)
    if (layer.GetDim() <= 1) {
        const int core_start = layer.GetCoreStart(0);
        const int core_end = layer.GetCoreEndExclusive(0);

        double max_speed = 0.0;
        for (int i = core_start; i < core_end; ++i) {
            const double rho = layer.rho(i);
            const double u = layer.u(i);
            const double P = layer.P(i);

            if (rho <= 0.0 || P <= 0.0) continue;

            const double c = std::sqrt(gamma * P / rho);
            const double s = std::abs(u) + c;
            if (s > max_speed) max_speed = s;
        }

        if (max_speed <= 0.0) return 0.0;
        return cfl * dx / max_speed;
    }

    // For 2D, fall back to the method that takes dx and dy explicitly
    return ComputeDt2D(layer, dx, dx, cfl, gamma);
}

auto TimeStepCalculator::ComputeDt2D(const DataLayer& layer,
                                      double dx, double dy,
                                      double cfl, double gamma) -> double {
    const int core_start_x = layer.GetCoreStart(0);
    const int core_end_x = layer.GetCoreEndExclusive(0);
    const int core_start_y = layer.GetCoreStart(1);
    const int core_end_y = layer.GetCoreEndExclusive(1);

    double max_spectral_radius = 0.0;

    for (int i = core_start_x; i < core_end_x; ++i) {
        for (int j = core_start_y; j < core_end_y; ++j) {
            const double rho = layer.rho(i, j);
            const double u_val = layer.u(i, j);
            const double v_val = layer.v(i, j);
            const double P_val = layer.P(i, j);

            if (rho <= 0.0 || P_val <= 0.0) continue;

            const double c = std::sqrt(gamma * P_val / rho);

            // 2D CFL: spectral radius = |u|/dx + |v|/dy + c*(1/dx + 1/dy)
            // But the standard formulation for unsplit methods is:
            //   dt = CFL / ( (|u|+c)/dx + (|v|+c)/dy )
            const double sr = (std::abs(u_val) + c) / dx + (std::abs(v_val) + c) / dy;

            if (sr > max_spectral_radius) {
                max_spectral_radius = sr;
            }
        }
    }

    if (max_spectral_radius <= 0.0) return 0.0;

    return cfl / max_spectral_radius;
}
