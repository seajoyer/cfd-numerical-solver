#include "bc/ReflectiveBoundary.hpp"

#include "data/DataLayer.hpp"
#include "data/Mesh.hpp"
#include "data/PressureVelocityState.hpp"

void ReflectiveBoundary::Apply(DataLayer& layer, const Mesh& mesh, const Axis axis, const Side side) const {
    const int ng = mesh.GetPadding();
    if (ng == 0) {
        return;
    }

    auto& U = layer.U();
    auto& reactant = layer.ReactantMassFraction();

    const int i0 = mesh.GetCoreStartX();
    const int i1 = mesh.GetCoreEndExclusiveX();
    const int j0 = mesh.GetCoreStartY();
    const int j1 = mesh.GetCoreEndExclusiveY();
    const int k0 = mesh.GetCoreStartZ();
    const int k1 = mesh.GetCoreEndExclusiveZ();

    for (int g = 0; g < ng; ++g) {
        if (axis == Axis::X) {
            const int dst_i = side == Side::Left ? (i0 - 1 - g) : (i1 + g);
            const int src_i = side == Side::Left ? (i0 + g) : (i1 - 1 - g);

            xt::view(U, xt::all(), dst_i, xt::all(), xt::all()) =
                xt::view(U, xt::all(), src_i, xt::all(), xt::all());

            xt::view(reactant, dst_i, xt::all(), xt::all()) =
                xt::view(reactant, src_i, xt::all(), xt::all());

            xt::view(U, DataLayer::k_rhoU, dst_i, xt::all(), xt::all()) *= -1.0;
            continue;
        }

        if (axis == Axis::Y) {
            const int dst_j = side == Side::Left ? (j0 - 1 - g) : (j1 + g);
            const int src_j = side == Side::Left ? (j0 + g) : (j1 - 1 - g);

            xt::view(U, xt::all(), xt::all(), dst_j, xt::all()) =
                xt::view(U, xt::all(), xt::all(), src_j, xt::all());

            xt::view(reactant, xt::all(), dst_j, xt::all()) =
                xt::view(reactant, xt::all(), src_j, xt::all());

            xt::view(U, DataLayer::k_rhoV, xt::all(), dst_j, xt::all()) *= -1.0;
            continue;
        }

        const int dst_k = side == Side::Left ? (k0 - 1 - g) : (k1 + g);
        const int src_k = side == Side::Left ? (k0 + g) : (k1 - 1 - g);

        xt::view(U, xt::all(), xt::all(), xt::all(), dst_k) =
            xt::view(U, xt::all(), xt::all(), xt::all(), src_k);

        xt::view(reactant, xt::all(), xt::all(), dst_k) =
            xt::view(reactant, xt::all(), xt::all(), src_k);

        xt::view(U, DataLayer::k_rhoW, xt::all(), xt::all(), dst_k) *= -1.0;
    }
}


void ReflectiveBoundary::Apply(PressureVelocityState& state,
                               const Mesh& mesh,
                               const Axis axis,
                               const Side side) const {
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

    auto& p = state.Pressure();
    auto& ux = state.Ux();
    auto& vy = state.Vy();
    auto& wz = state.Wz();

    const int i0 = mesh.GetCoreStartX();
    const int i1 = mesh.GetCoreEndExclusiveX();
    const int j0 = mesh.GetCoreStartY();
    const int j1 = mesh.GetCoreEndExclusiveY();
    const int k0 = mesh.GetCoreStartZ();
    const int k1 = mesh.GetCoreEndExclusiveZ();

    for (int g = 0; g < ng; ++g) {
        if (axis == Axis::X) {
            const int dst_i_p = side == Side::Left ? (i0 - 1 - g) : (i1 + g);
            const int src_i_p = side == Side::Left ? (i0 + g) : (i1 - 1 - g);

            xt::view(p, dst_i_p, xt::all(), xt::all()) =
                xt::view(p, src_i_p, xt::all(), xt::all());

            const int dst_i_ux = side == Side::Left ? (i0 - g) : (i1 + g);
            const int src_i_ux = side == Side::Left ? (i0 + 1 + g) : (i1 - 1 - g);

            xt::view(ux, dst_i_ux, xt::all(), xt::all()) =
                xt::view(ux, src_i_ux, xt::all(), xt::all());
            xt::view(ux, dst_i_ux, xt::all(), xt::all()) *= -1.0;

            xt::view(vy, dst_i_p, xt::all(), xt::all()) =
                xt::view(vy, src_i_p, xt::all(), xt::all());

            xt::view(wz, dst_i_p, xt::all(), xt::all()) =
                xt::view(wz, src_i_p, xt::all(), xt::all());

            continue;
        }

        if (axis == Axis::Y) {
            const int dst_j_p = side == Side::Left ? (j0 - 1 - g) : (j1 + g);
            const int src_j_p = side == Side::Left ? (j0 + g) : (j1 - 1 - g);

            xt::view(p, xt::all(), dst_j_p, xt::all()) =
                xt::view(p, xt::all(), src_j_p, xt::all());

            xt::view(ux, xt::all(), dst_j_p, xt::all()) =
                xt::view(ux, xt::all(), src_j_p, xt::all());

            const int dst_j_vy = side == Side::Left ? (j0 - g) : (j1 + g);
            const int src_j_vy = side == Side::Left ? (j0 + 1 + g) : (j1 - 1 - g);

            xt::view(vy, xt::all(), dst_j_vy, xt::all()) =
                xt::view(vy, xt::all(), src_j_vy, xt::all());
            xt::view(vy, xt::all(), dst_j_vy, xt::all()) *= -1.0;

            xt::view(wz, xt::all(), dst_j_p, xt::all()) =
                xt::view(wz, xt::all(), src_j_p, xt::all());

            continue;
        }

        const int dst_k_p = side == Side::Left ? (k0 - 1 - g) : (k1 + g);
        const int src_k_p = side == Side::Left ? (k0 + g) : (k1 - 1 - g);

        xt::view(p, xt::all(), xt::all(), dst_k_p) =
            xt::view(p, xt::all(), xt::all(), src_k_p);

        xt::view(ux, xt::all(), xt::all(), dst_k_p) =
            xt::view(ux, xt::all(), xt::all(), src_k_p);

        xt::view(vy, xt::all(), xt::all(), dst_k_p) =
            xt::view(vy, xt::all(), xt::all(), src_k_p);

        const int dst_k_wz = side == Side::Left ? (k0 - g) : (k1 + g);
        const int src_k_wz = side == Side::Left ? (k0 + 1 + g) : (k1 - 1 - g);

        xt::view(wz, xt::all(), xt::all(), dst_k_wz) =
            xt::view(wz, xt::all(), xt::all(), src_k_wz);
        xt::view(wz, xt::all(), xt::all(), dst_k_wz) *= -1.0;
    }
}
