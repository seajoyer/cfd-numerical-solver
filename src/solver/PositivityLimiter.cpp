// PositivityLimiter.cpp
#include "solver/PositivityLimiter.hpp"

#include <algorithm>

void PositivityLimiter::Apply(DataLayer& layer, const double gamma, const double rho_min, const double p_min) {
    if (rho_min <= 0.0 && p_min <= 0.0) {
        return;
    }

    auto& U = layer.U();

    const int i0 = layer.GetCoreStartX();
    const int i1 = layer.GetCoreEndExclusiveX();
    const int j0 = layer.GetCoreStartY();
    const int j1 = layer.GetCoreEndExclusiveY();
    const int k0 = layer.GetCoreStartZ();
    const int k1 = layer.GetCoreEndExclusiveZ();

    for (int k = k0; k < k1; ++k) {
        for (int j = j0; j < j1; ++j) {
            for (int i = i0; i < i1; ++i) {
                double rho = U(DataLayer::k_rho, i, j, k);
                double rhoU = U(DataLayer::k_rhoU, i, j, k);
                double rhoV = U(DataLayer::k_rhoV, i, j, k);
                double rhoW = U(DataLayer::k_rhoW, i, j, k);
                double E    = U(DataLayer::k_E, i, j, k);

                // Compute velocities from current state if possible
                double u = 0.0, v = 0.0, w = 0.0;
                if (rho > 0.0) {
                    const double inv_rho = 1.0 / rho;
                    u = rhoU * inv_rho;
                    v = rhoV * inv_rho;
                    w = rhoW * inv_rho;
                }

                // 1) Clamp density, preserving velocities
                if (rho_min > 0.0 && rho < rho_min) {
                    rho = rho_min;
                    rhoU = rho * u;
                    rhoV = rho * v;
                    rhoW = rho * w;
                }

                // 2) Clamp pressure by adjusting energy E
                if (p_min > 0.0) {
                    const double kinetic = 0.5 * rho * (u * u + v * v + w * w);
                    const double eint = E - kinetic;
                    const double P = (gamma - 1.0) * eint;

                    if (P < p_min) {
                        const double eint_min = p_min / (gamma - 1.0);
                        E = eint_min + kinetic;
                    }
                }

                U(DataLayer::k_rho,  i, j, k) = rho;
                U(DataLayer::k_rhoU, i, j, k) = rhoU;
                U(DataLayer::k_rhoV, i, j, k) = rhoV;
                U(DataLayer::k_rhoW, i, j, k) = rhoW;
                U(DataLayer::k_E,    i, j, k) = E;
            }
        }
    }
}