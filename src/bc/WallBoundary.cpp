#include "bc/WallBoundary.hpp"

#include "data/DataLayer.hpp"
#include "data/Mesh.hpp"
#include "data/Variables.hpp"

namespace {
    void CopyScalarsAndZeroMomentaX(xt::xtensor<double, 4>& U, const int dst_i, const int src_i) {
        xt::view(U, DataLayer::k_rho, dst_i, xt::all(), xt::all()) =
            xt::view(U, DataLayer::k_rho, src_i, xt::all(), xt::all());
        xt::view(U, DataLayer::k_E, dst_i, xt::all(), xt::all()) =
            xt::view(U, DataLayer::k_E, src_i, xt::all(), xt::all());

        xt::view(U, DataLayer::k_rhoU, dst_i, xt::all(), xt::all()) = 0.0;
        xt::view(U, DataLayer::k_rhoV, dst_i, xt::all(), xt::all()) = 0.0;
        xt::view(U, DataLayer::k_rhoW, dst_i, xt::all(), xt::all()) = 0.0;
    }

    void CopyScalarsAndZeroMomentaY(xt::xtensor<double, 4>& U, const int dst_j, const int src_j) {
        xt::view(U, DataLayer::k_rho, xt::all(), dst_j, xt::all()) =
            xt::view(U, DataLayer::k_rho, xt::all(), src_j, xt::all());
        xt::view(U, DataLayer::k_E, xt::all(), dst_j, xt::all()) =
            xt::view(U, DataLayer::k_E, xt::all(), src_j, xt::all());

        xt::view(U, DataLayer::k_rhoU, xt::all(), dst_j, xt::all()) = 0.0;
        xt::view(U, DataLayer::k_rhoV, xt::all(), dst_j, xt::all()) = 0.0;
        xt::view(U, DataLayer::k_rhoW, xt::all(), dst_j, xt::all()) = 0.0;
    }

    void CopyScalarsAndZeroMomentaZ(xt::xtensor<double, 4>& U, const int dst_k, const int src_k) {
        xt::view(U, DataLayer::k_rho, xt::all(), xt::all(), dst_k) =
            xt::view(U, DataLayer::k_rho, xt::all(), xt::all(), src_k);
        xt::view(U, DataLayer::k_E, xt::all(), xt::all(), dst_k) =
            xt::view(U, DataLayer::k_E, xt::all(), xt::all(), src_k);

        xt::view(U, DataLayer::k_rhoU, xt::all(), xt::all(), dst_k) = 0.0;
        xt::view(U, DataLayer::k_rhoV, xt::all(), xt::all(), dst_k) = 0.0;
        xt::view(U, DataLayer::k_rhoW, xt::all(), xt::all(), dst_k) = 0.0;
    }
} // namespace

void WallBoundary::Apply(DataLayer& layer, const Mesh& mesh, const Axis axis, const Side side) const {
    const int ng = mesh.GetPadding();
    if (ng == 0) {
        return;
    }

    const int dim = mesh.GetDim();
    if (axis == Axis::Y && dim < 2) {
        return;
    }
    if (axis == Axis::Z && dim < 3) {
        return;
    }

    auto& U = layer.U();

    const int i0 = mesh.GetCoreStartX();
    const int i1 = mesh.GetCoreEndExclusiveX();
    const int j0 = mesh.GetCoreStartY();
    const int j1 = mesh.GetCoreEndExclusiveY();
    const int k0 = mesh.GetCoreStartZ();
    const int k1 = mesh.GetCoreEndExclusiveZ();

    if (axis == Axis::X) {
        const int src_i = side == Side::Left ? i0 : (i1 - 1);
        for (int g = 0; g < ng; ++g) {
            const int dst_i = side == Side::Left ? (i0 - 1 - g) : (i1 + g);
            CopyScalarsAndZeroMomentaX(U, dst_i, src_i);
        }
        return;
    }

    if (axis == Axis::Y) {
        const int src_j = side == Side::Left ? j0 : (j1 - 1);
        for (int g = 0; g < ng; ++g) {
            const int dst_j = side == Side::Left ? (j0 - 1 - g) : (j1 + g);
            CopyScalarsAndZeroMomentaY(U, dst_j, src_j);
        }
        return;
    }

    const int src_k = side == Side::Left ? k0 : (k1 - 1);
    for (int g = 0; g < ng; ++g) {
        const int dst_k = side == Side::Left ? (k0 - 1 - g) : (k1 + g);
        CopyScalarsAndZeroMomentaZ(U, dst_k, src_k);
    }
}
