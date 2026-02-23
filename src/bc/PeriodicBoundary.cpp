#include "bc/PeriodicBoundary.hpp"

#include "data/DataLayer.hpp"


void PeriodicBoundary::Apply(DataLayer& layer, const Axis axis, const Side side) const {
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

    for (int g = 0; g < ng; ++g) {
        if (axis == Axis::X) {
            const int dst_i = (side == Side::Left) ? (i0 - 1 - g) : (i1 + g);
            const int src_i = (side == Side::Left) ? (i1 - 1 - g) : (i0 + g);

            xt::view(U, xt::all(), dst_i, xt::all(), xt::all()) =
                xt::view(U, xt::all(), src_i, xt::all(), xt::all());
            continue;
        }

        if (axis == Axis::Y) {
            const int dst_j = (side == Side::Left) ? (j0 - 1 - g) : (j1 + g);
            const int src_j = (side == Side::Left) ? (j1 - 1 - g) : (j0 + g);

            xt::view(U, xt::all(), xt::all(), dst_j, xt::all()) =
                xt::view(U, xt::all(), xt::all(), src_j, xt::all());
            continue;
        }

        // Axis::Z
        const int dst_k = (side == Side::Left) ? (k0 - 1 - g) : (k1 + g);
        const int src_k = (side == Side::Left) ? (k1 - 1 - g) : (k0 + g);

        xt::view(U, xt::all(), xt::all(), xt::all(), dst_k) =
            xt::view(U, xt::all(), xt::all(), xt::all(), src_k);
    }
}