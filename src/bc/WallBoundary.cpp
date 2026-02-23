#include "bc/WallBoundary.hpp"

#include "data/DataLayer.hpp"


namespace {

inline void CopyScalarsAndZeroMomentaX(
    xt::xtensor<double, 4>& U,
    int dst_i, int src_i
) {
    // Copy scalars: rho, E
    xt::view(U, DataLayer::k_rho, dst_i, xt::all(), xt::all()) =
        xt::view(U, DataLayer::k_rho, src_i, xt::all(), xt::all());
    xt::view(U, DataLayer::k_E, dst_i, xt::all(), xt::all()) =
        xt::view(U, DataLayer::k_E, src_i, xt::all(), xt::all());

    // Zero momenta: rhoU, rhoV, rhoW
    xt::view(U, DataLayer::k_rhoU, dst_i, xt::all(), xt::all()) = 0.0;
    xt::view(U, DataLayer::k_rhoV, dst_i, xt::all(), xt::all()) = 0.0;
    xt::view(U, DataLayer::k_rhoW, dst_i, xt::all(), xt::all()) = 0.0;
}

inline void CopyScalarsAndZeroMomentaY(
    xt::xtensor<double, 4>& U,
    int dst_j, int src_j
) {
    xt::view(U, DataLayer::k_rho, xt::all(), dst_j, xt::all()) =
        xt::view(U, DataLayer::k_rho, xt::all(), src_j, xt::all());
    xt::view(U, DataLayer::k_E, xt::all(), dst_j, xt::all()) =
        xt::view(U, DataLayer::k_E, xt::all(), src_j, xt::all());

    xt::view(U, DataLayer::k_rhoU, xt::all(), dst_j, xt::all()) = 0.0;
    xt::view(U, DataLayer::k_rhoV, xt::all(), dst_j, xt::all()) = 0.0;
    xt::view(U, DataLayer::k_rhoW, xt::all(), dst_j, xt::all()) = 0.0;
}

inline void CopyScalarsAndZeroMomentaZ(
    xt::xtensor<double, 4>& U,
    int dst_k, int src_k
) {
    xt::view(U, DataLayer::k_rho, xt::all(), xt::all(), dst_k) =
        xt::view(U, DataLayer::k_rho, xt::all(), xt::all(), src_k);
    xt::view(U, DataLayer::k_E, xt::all(), xt::all(), dst_k) =
        xt::view(U, DataLayer::k_E, xt::all(), xt::all(), src_k);

    xt::view(U, DataLayer::k_rhoU, xt::all(), xt::all(), dst_k) = 0.0;
    xt::view(U, DataLayer::k_rhoV, xt::all(), xt::all(), dst_k) = 0.0;
    xt::view(U, DataLayer::k_rhoW, xt::all(), xt::all(), dst_k) = 0.0;
}

}  // namespace

void WallBoundary::Apply(DataLayer& layer, const Axis axis, const Side side) const {
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

    if (axis == Axis::X) {
        const int src_i = (side == Side::Left) ? i0 : (i1 - 1);
        for (int g = 0; g < ng; ++g) {
            const int dst_i = (side == Side::Left) ? (i0 - 1 - g) : (i1 + g);
            CopyScalarsAndZeroMomentaX(U, dst_i, src_i);
        }
        return;
    }

    if (axis == Axis::Y) {
        const int src_j = (side == Side::Left) ? j0 : (j1 - 1);
        for (int g = 0; g < ng; ++g) {
            const int dst_j = (side == Side::Left) ? (j0 - 1 - g) : (j1 + g);
            CopyScalarsAndZeroMomentaY(U, dst_j, src_j);
        }
        return;
    }

    // Axis::Z
    const int src_k = (side == Side::Left) ? k0 : (k1 - 1);
    for (int g = 0; g < ng; ++g) {
        const int dst_k = (side == Side::Left) ? (k0 - 1 - g) : (k1 + g);
        CopyScalarsAndZeroMomentaZ(U, dst_k, src_k);
    }
}