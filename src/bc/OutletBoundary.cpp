#include "bc/OutletBoundary.hpp"

#include "data/DataLayer.hpp"
#include "data/Mesh.hpp"
#include "data/Variables.hpp"

void OutletBoundary::Apply(DataLayer& layer, const Mesh& mesh, const Axis axis, const Side side) const {
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

    for (int g = 0; g < ng; ++g) {
        if (axis == Axis::X) {
            const int dst_i = side == Side::Left ? (i0 - 1 - g) : (i1 + g);
            const int src_i = side == Side::Left ? i0 : (i1 - 1);

            xt::view(U, xt::all(), dst_i, xt::all(), xt::all()) =
                xt::view(U, xt::all(), src_i, xt::all(), xt::all());
            continue;
        }

        if (axis == Axis::Y) {
            const int dst_j = side == Side::Left ? (j0 - 1 - g) : (j1 + g);
            const int src_j = side == Side::Left ? j0 : (j1 - 1);

            xt::view(U, xt::all(), xt::all(), dst_j, xt::all()) =
                xt::view(U, xt::all(), xt::all(), src_j, xt::all());
            continue;
        }

        const int dst_k = side == Side::Left ? (k0 - 1 - g) : (k1 + g);
        const int src_k = side == Side::Left ? k0 : (k1 - 1);

        xt::view(U, xt::all(), xt::all(), xt::all(), dst_k) =
            xt::view(U, xt::all(), xt::all(), xt::all(), src_k);
    }
}