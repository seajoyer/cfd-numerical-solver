#include "solver/PositivityLimiter.hpp"

#include "data/DataLayer.hpp"
#include "geometry/Mesh.hpp"

void PositivityLimiter::Apply(DataLayer& layer,
                              const Mesh& mesh,
                              const double gamma,
                              const double rho_min,
                              const double p_min) {
    if (rho_min <= 0.0 && p_min <= 0.0) {
        return;
    }

    auto& U = layer.U();

    for (std::size_t cell_id = 0; cell_id < mesh.GetOwnedCellCount(); ++cell_id) {
        double rho = U(cell_id, DataLayer::k_rho);
        double rhoU = U(cell_id, DataLayer::k_rhoU);
        double rhoV = U(cell_id, DataLayer::k_rhoV);
        double rhoW = U(cell_id, DataLayer::k_rhoW);
        double E = U(cell_id, DataLayer::k_E);

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
            const double internal_energy_density = E - kinetic;
            const double P = (gamma - 1.0) * internal_energy_density;

            if (P < p_min) {
                const double internal_energy_density_min = p_min / (gamma - 1.0);
                E = internal_energy_density_min + kinetic;
            }
        }

        U(cell_id, DataLayer::k_rho) = rho;
        U(cell_id, DataLayer::k_rhoU) = rhoU;
        U(cell_id, DataLayer::k_rhoV) = rhoV;
        U(cell_id, DataLayer::k_rhoW) = rhoW;
        U(cell_id, DataLayer::k_E) = E;
    }
}
