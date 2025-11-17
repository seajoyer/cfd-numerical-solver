#include "solver/TimeStepCalculator.hpp"

auto TimeStepCalculator::ComputeDt(const DataLayer& layer, const double dx,
                                   const double cfl, const double gamma) -> double {
    const int core_start = layer.GetCoreStart(0);
    const int core_end = layer.GetCoreEndExclusive(0);

    double max_speed = 0.0;

    for (int i = core_start; i < core_end; ++i) {
        const double rho = layer.rho(i);
        const double u = layer.u(i);
        const double P = layer.P(i);

        if (rho <= 0.0 || P <= 0.0) {
            continue;
        }

        const double c = std::sqrt(gamma * P / rho);
        const double s = std::abs(u) + c;
        if (s > max_speed) {
            max_speed = s;
        }
    }

    if (max_speed <= 0.0) {
        return 0.0;
    }

    return cfl * dx / max_speed;
}
