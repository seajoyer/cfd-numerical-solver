#include "filter/SolutionFilter.hpp"

#include <vector>
#include <algorithm>

void SolutionFilter::Apply(DataLayer& layer) const
{
    const int total_size = layer.GetTotalSize();
    if (total_size < 5) {
        return;
    }

    const int core_start = layer.GetCoreStart(0);
    const int core_end = layer.GetCoreEndExclusive(0);

    constexpr double eps_diff = 0.05;
    constexpr double eps_anti = 0.03;

    std::vector<double> rho_f(total_size);
    std::vector<double> P_f(total_size);

    for (int j = core_start; j < core_end; ++j) {
        rho_f[j] = layer.rho(j)
            + eps_diff * (layer.rho(j + 1) - 2.0 * layer.rho(j) + layer.rho(j - 1));

        P_f[j] = layer.P(j)
            + eps_diff * (layer.P(j + 1) - 2.0 * layer.P(j) + layer.P(j - 1));
    }

    for (int j = 0; j < core_start; ++j) {
        rho_f[j] = layer.rho(j);
        P_f[j] = layer.P(j);
    }
    for (int j = core_end; j < total_size; ++j)
    {
        rho_f[j] = layer.rho(j);
        P_f[j] = layer.P(j);
    }

    for (int j = core_start + 1; j < core_end - 1; ++j) {
        layer.rho(j) = rho_f[j]
            - eps_anti * (rho_f[j + 1] - 2.0 * rho_f[j] + rho_f[j - 1]);

        layer.P(j) = P_f[j]
            - eps_anti * (P_f[j + 1] - 2.0 * P_f[j] + P_f[j - 1]);
    }

    constexpr double rho_floor = 1e-10;
    constexpr double p_floor = 1e-10;

    for (int j = core_start; j < core_end; ++j) {
        layer.rho(j) = std::max(layer.rho(j), rho_floor);
        layer.P(j) = std::max(layer.P(j), p_floor);
    }
}
