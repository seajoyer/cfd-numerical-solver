#include "solver/PositivityLimiter.hpp"

#include "data/DataLayer.hpp"
#include "data/Mesh.hpp"

void PositivityLimiter::Apply(DataLayer& layer,
                              const Mesh& mesh,
                              const double gamma,
                              const double rho_min,
                              const double p_min) {
    if (rho_min <= 0.0 && p_min <= 0.0) {
        return;
    }

    auto& U = layer.U();

    const int i0 = mesh.GetCoreStartX();
    const int i1 = mesh.GetCoreEndExclusiveX();
    const int j0 = mesh.GetCoreStartY();
    const int j1 = mesh.GetCoreEndExclusiveY();
    const int k0 = mesh.GetCoreStartZ();
    const int k1 = mesh.GetCoreEndExclusiveZ();

    for (int k = k0; k < k1; ++k) {
        for (int j = j0; j < j1; ++j) {
            for (int i = i0; i < i1; ++i) {
                if (!mesh.IsFluidCell(i, j, k)) {
                    continue;
                }

                double rho = U(DataLayer::k_rho, i, j, k);
                double rhoU = U(DataLayer::k_rhoU, i, j, k);
                double rhoV = U(DataLayer::k_rhoV, i, j, k);
                double rhoW = U(DataLayer::k_rhoW, i, j, k);
                double E = U(DataLayer::k_E, i, j, k);

                double u = 0.0;
                double v = 0.0;
                double w = 0.0;

                if (rho > 0.0) {
                    const double inv_rho = 1.0 / rho;
                    u = rhoU * inv_rho;
                    v = rhoV * inv_rho;
                    w = rhoW * inv_rho;
                }

                if (rho_min > 0.0 && rho < rho_min) {
                    rho = rho_min;
                    rhoU = rho * u;
                    rhoV = rho * v;
                    rhoW = rho * w;
                }

                if (p_min > 0.0) {
                    const double kinetic = 0.5 * rho * (u * u + v * v + w * w);
                    const double eint = E - kinetic;
                    const double P = (gamma - 1.0) * eint;

                    if (P < p_min) {
                        const double eint_min = p_min / (gamma - 1.0);
                        E = eint_min + kinetic;
                    }
                }

                U(DataLayer::k_rho, i, j, k) = rho;
                U(DataLayer::k_rhoU, i, j, k) = rhoU;
                U(DataLayer::k_rhoV, i, j, k) = rhoV;
                U(DataLayer::k_rhoW, i, j, k) = rhoW;
                U(DataLayer::k_E, i, j, k) = E;
            }
        }
    }
}