#include "solver/TimeStepCalculator.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "config/Settings.hpp"
#include "data/DataLayer.hpp"
#include "data/Mesh.hpp"
#include "solver/EOS.hpp"

auto TimeStepCalculator::ComputeDt(const DataLayer& layer,
                                   const Mesh& mesh,
                                   const Settings& settings,
                                   std::shared_ptr<EOS> eos) -> double {
    const double cfl = settings.cfl;
    if (cfl <= 0.0) {
        return 0.0;
    }

    const int dim = mesh.GetDim();

    const int i0 = mesh.GetCoreStartX();
    const int i1 = mesh.GetCoreEndExclusiveX();
    const int j0 = mesh.GetCoreStartY();
    const int j1 = mesh.GetCoreEndExclusiveY();
    const int k0 = mesh.GetCoreStartZ();
    const int k1 = mesh.GetCoreEndExclusiveZ();

    const bool active_x = (i1 - i0) >= 2;
    const bool active_y = (dim >= 2) && ((j1 - j0) >= 2);
    const bool active_z = (dim >= 3) && ((k1 - k0) >= 2);

    if (!active_x && !active_y && !active_z) {
        return 0.0;
    }

    const auto& U = layer.U();
    const auto& dx = mesh.Dx();
    const auto& dy = mesh.Dy();
    const auto& dz = mesh.Dz();

    const bool has_chemistry = settings.chemistry_enabled;

    double dt_min = std::numeric_limits<double>::infinity();
    bool has_dt = false;

    for (int k = k0; k < k1; ++k) {
        for (int j = j0; j < j1; ++j) {
            for (int i = i0; i < i1; ++i) {
                if (!mesh.IsFluidCell(i, j, k)) {
                    continue;
                }

                const double rho = U(DataLayer::k_rho, i, j, k);
                if (rho <= 0.0) {
                    continue;
                }

                const double inv_rho = 1.0 / rho;
                const double u = U(DataLayer::k_rhoU, i, j, k) * inv_rho;
                const double v = U(DataLayer::k_rhoV, i, j, k) * inv_rho;
                const double w = U(DataLayer::k_rhoW, i, j, k) * inv_rho;

                const double kinetic = 0.5 * (u * u + v * v + w * w);
                const double E_in = U(DataLayer::k_E, i, j, k);

                const double I_cell = std::max(E_in * inv_rho - kinetic, 0.0);

                double lambda = 0.0;
                if (has_chemistry) {
                    lambda = layer.ReactantMassFraction()(i, j, k);
                }

                EosCellInput eos_in{rho, I_cell, lambda};
                EosCellOutput eos_out = eos->Evaluate(eos_in);

                if (eos_out.P <= 0.0 || eos_out.c <= 0.0) {
                    continue;
                }

                const double c = eos_out.c;

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