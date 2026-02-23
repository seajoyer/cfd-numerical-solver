#include "solver/TimeStepCalculator.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

auto TimeStepCalculator::ComputeDt(const DataLayer& layer, const double gamma, const double cfl) -> double {
    if (cfl <= 0.0) {
        return 0.0;
    }

    const int dim = layer.GetDim();

    const int i0 = layer.GetCoreStartX();
    const int i1 = layer.GetCoreEndExclusiveX();
    const int j0 = layer.GetCoreStartY();
    const int j1 = layer.GetCoreEndExclusiveY();
    const int k0 = layer.GetCoreStartZ();
    const int k1 = layer.GetCoreEndExclusiveZ();

    const bool active_x = (i1 - i0) >= 2;
    const bool active_y = (dim >= 2) && ((j1 - j0) >= 2);
    const bool active_z = (dim >= 3) && ((k1 - k0) >= 2);

    if (!active_x && !active_y && !active_z) {
        return 0.0;
    }

    const auto& U = layer.U();
    const auto& dx = layer.Dx();
    const auto& dy = layer.Dy();
    const auto& dz = layer.Dz();

    double dt_min = std::numeric_limits<double>::infinity();
    bool has_dt = false;

    for (int k = k0; k < k1; ++k) {
        for (int j = j0; j < j1; ++j) {
            for (int i = i0; i < i1; ++i) {
                const double rho = U(DataLayer::k_rho, i, j, k);
                if (rho <= 0.0) {
                    continue;
                }

                const double inv_rho = 1.0 / rho;
                const double u = U(DataLayer::k_rhoU, i, j, k) * inv_rho;
                const double v = U(DataLayer::k_rhoV, i, j, k) * inv_rho;
                const double w = U(DataLayer::k_rhoW, i, j, k) * inv_rho;

                const double kinetic = 0.5 * rho * (u * u + v * v + w * w);
                const double E = U(DataLayer::k_E, i, j, k);
                const double eint = E - kinetic;
                const double P = (gamma - 1.0) * eint;

                if (P <= 0.0) {
                    continue;
                }

                const double c2 = gamma * P * inv_rho;
                if (c2 <= 0.0) {
                    continue;
                }
                const double c = std::sqrt(c2);

                if (active_x) {
                    const double s = std::abs(u) + c;
                    if (s > 0.0) {
                        const double dxi = dx(static_cast<std::size_t>(i));
                        if (dxi > 0.0) {
                            dt_min = std::min(dt_min, cfl * (dxi / s));
                            has_dt = true;
                        }
                    }
                }

                if (active_y) {
                    const double s = std::abs(v) + c;
                    if (s > 0.0) {
                        const double dyj = dy(static_cast<std::size_t>(j));
                        if (dyj > 0.0) {
                            dt_min = std::min(dt_min, cfl * (dyj / s));
                            has_dt = true;
                        }
                    }
                }

                if (active_z) {
                    const double s = std::abs(w) + c;
                    if (s > 0.0) {
                        const double dzk = dz(static_cast<std::size_t>(k));
                        if (dzk > 0.0) {
                            dt_min = std::min(dt_min, cfl * (dzk / s));
                            has_dt = true;
                        }
                    }
                }
            }
        }
    }

    if (!has_dt || !std::isfinite(dt_min) || dt_min <= 0.0) {
        return 0.0;
    }
    return dt_min;
}
