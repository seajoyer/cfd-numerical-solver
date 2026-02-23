#include "bc/InletBoundary.hpp"

#include "data/DataLayer.hpp"

#include <algorithm>

namespace {
    inline double NormalVelocityFromU(
        const xt::xtensor<double, 4>& U, Axis axis,
        int i, int j, int k,
        double rho_floor = 1e-14
    ) {
        const double rho_in = U(DataLayer::k_rho, i, j, k);
        const double rho = (rho_in > rho_floor) ? rho_in : rho_floor;
        const double inv_rho = 1.0 / rho;

        if (axis == Axis::X) return U(DataLayer::k_rhoU, i, j, k) * inv_rho;
        if (axis == Axis::Y) return U(DataLayer::k_rhoV, i, j, k) * inv_rho;
        return U(DataLayer::k_rhoW, i, j, k) * inv_rho; // Axis::Z
    }

    inline void SetPlaneToInflowX(
        xt::xtensor<double, 4>& U,
        int dst_i,
        const FarfieldConservative& inflow
    ) {
        xt::view(U, DataLayer::k_rho, dst_i, xt::all(), xt::all()) = inflow.rho;
        xt::view(U, DataLayer::k_rhoU, dst_i, xt::all(), xt::all()) = inflow.rhoU;
        xt::view(U, DataLayer::k_rhoV, dst_i, xt::all(), xt::all()) = inflow.rhoV;
        xt::view(U, DataLayer::k_rhoW, dst_i, xt::all(), xt::all()) = inflow.rhoW;
        xt::view(U, DataLayer::k_E, dst_i, xt::all(), xt::all()) = inflow.E;
    }

    inline void SetPlaneToInflowY(
        xt::xtensor<double, 4>& U,
        int dst_j,
        const FarfieldConservative& inflow
    ) {
        xt::view(U, DataLayer::k_rho, xt::all(), dst_j, xt::all()) = inflow.rho;
        xt::view(U, DataLayer::k_rhoU, xt::all(), dst_j, xt::all()) = inflow.rhoU;
        xt::view(U, DataLayer::k_rhoV, xt::all(), dst_j, xt::all()) = inflow.rhoV;
        xt::view(U, DataLayer::k_rhoW, xt::all(), dst_j, xt::all()) = inflow.rhoW;
        xt::view(U, DataLayer::k_E, xt::all(), dst_j, xt::all()) = inflow.E;
    }

    inline void SetPlaneToInflowZ(
        xt::xtensor<double, 4>& U,
        int dst_k,
        const FarfieldConservative& inflow
    ) {
        xt::view(U, DataLayer::k_rho, xt::all(), xt::all(), dst_k) = inflow.rho;
        xt::view(U, DataLayer::k_rhoU, xt::all(), xt::all(), dst_k) = inflow.rhoU;
        xt::view(U, DataLayer::k_rhoV, xt::all(), xt::all(), dst_k) = inflow.rhoV;
        xt::view(U, DataLayer::k_rhoW, xt::all(), xt::all(), dst_k) = inflow.rhoW;
        xt::view(U, DataLayer::k_E, xt::all(), xt::all(), dst_k) = inflow.E;
    }
} // namespace

InletBoundary::InletBoundary(const FarfieldConservative& inflow_U)
    : inflow_U_(inflow_U) {}

void InletBoundary::Apply(DataLayer& layer, const Axis axis, const Side side) const {
    const int ng = layer.GetPadding();
    if (ng == 0) return;

    const int dim = layer.GetDim();
    if (axis == Axis::Y && dim < 2) return;
    if (axis == Axis::Z && dim < 3) return;

    auto& U = layer.U();

    const int i0 = layer.GetCoreStartX();
    const int i1 = layer.GetCoreEndExclusiveX();
    const int j0 = layer.GetCoreStartY();
    const int j1 = layer.GetCoreEndExclusiveY();
    const int k0 = layer.GetCoreStartZ();
    const int k1 = layer.GetCoreEndExclusiveZ();

    // Determine inward/outward based on normal velocity at the nearest interior core cell.
    bool inward = false;

    if (axis == Axis::X) {
        const int ii = (side == Side::Left) ? i0 : (i1 - 1);
        const double un = NormalVelocityFromU(U, axis, ii, j0, k0);
        inward = (side == Side::Left) ? (un > 0.0) : (un < 0.0);

        for (int g = 0; g < ng; ++g) {
            const int dst_i = (side == Side::Left) ? (i0 - 1 - g) : (i1 + g);
            const int src_i = (side == Side::Left) ? i0 : (i1 - 1);

            if (inward) {
                SetPlaneToInflowX(U, dst_i, inflow_U_);
            }
            else {
                xt::view(U, xt::all(), dst_i, xt::all(), xt::all()) =
                    xt::view(U, xt::all(), src_i, xt::all(), xt::all());
            }
        }
        return;
    }

    if (axis == Axis::Y) {
        const int jj = (side == Side::Left) ? j0 : (j1 - 1);
        const double un = NormalVelocityFromU(U, axis, i0, jj, k0);
        inward = (side == Side::Left) ? (un > 0.0) : (un < 0.0);

        for (int g = 0; g < ng; ++g) {
            const int dst_j = (side == Side::Left) ? (j0 - 1 - g) : (j1 + g);
            const int src_j = (side == Side::Left) ? j0 : (j1 - 1);

            if (inward) {
                SetPlaneToInflowY(U, dst_j, inflow_U_);
            }
            else {
                xt::view(U, xt::all(), xt::all(), dst_j, xt::all()) =
                    xt::view(U, xt::all(), xt::all(), src_j, xt::all());
            }
        }
        return;
    }

    // Axis::Z
    const int kk = (side == Side::Left) ? k0 : (k1 - 1);
    const double un = NormalVelocityFromU(U, axis, i0, j0, kk);
    inward = (side == Side::Left) ? (un > 0.0) : (un < 0.0);

    for (int g = 0; g < ng; ++g) {
        const int dst_k = (side == Side::Left) ? (k0 - 1 - g) : (k1 + g);
        const int src_k = (side == Side::Left) ? k0 : (k1 - 1);

        if (inward) {
            SetPlaneToInflowZ(U, dst_k, inflow_U_);
        }
        else {
            xt::view(U, xt::all(), xt::all(), xt::all(), dst_k) =
                xt::view(U, xt::all(), xt::all(), xt::all(), src_k);
        }
    }
}
