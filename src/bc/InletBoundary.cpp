#include "bc/InletBoundary.hpp"

#include "data/DataLayer.hpp"
#include "data/Mesh.hpp"
#include "data/Variables.hpp"

namespace {
void SetPlaneToStateX(xt::xtensor<double, 4>& U, const int dst_i, const FarfieldConservative& s) {
    xt::view(U, DataLayer::k_rho,  dst_i, xt::all(), xt::all()) = s.rho;
    xt::view(U, DataLayer::k_rhoU, dst_i, xt::all(), xt::all()) = s.rhoU;
    xt::view(U, DataLayer::k_rhoV, dst_i, xt::all(), xt::all()) = s.rhoV;
    xt::view(U, DataLayer::k_rhoW, dst_i, xt::all(), xt::all()) = s.rhoW;
    xt::view(U, DataLayer::k_E,    dst_i, xt::all(), xt::all()) = s.E;
}

void SetPlaneToStateY(xt::xtensor<double, 4>& U, const int dst_j, const FarfieldConservative& s) {
    xt::view(U, DataLayer::k_rho,  xt::all(), dst_j, xt::all()) = s.rho;
    xt::view(U, DataLayer::k_rhoU, xt::all(), dst_j, xt::all()) = s.rhoU;
    xt::view(U, DataLayer::k_rhoV, xt::all(), dst_j, xt::all()) = s.rhoV;
    xt::view(U, DataLayer::k_rhoW, xt::all(), dst_j, xt::all()) = s.rhoW;
    xt::view(U, DataLayer::k_E,    xt::all(), dst_j, xt::all()) = s.E;
}

void SetPlaneToStateZ(xt::xtensor<double, 4>& U, const int dst_k, const FarfieldConservative& s) {
    xt::view(U, DataLayer::k_rho,  xt::all(), xt::all(), dst_k) = s.rho;
    xt::view(U, DataLayer::k_rhoU, xt::all(), xt::all(), dst_k) = s.rhoU;
    xt::view(U, DataLayer::k_rhoV, xt::all(), xt::all(), dst_k) = s.rhoV;
    xt::view(U, DataLayer::k_rhoW, xt::all(), xt::all(), dst_k) = s.rhoW;
    xt::view(U, DataLayer::k_E,    xt::all(), xt::all(), dst_k) = s.E;
}
}  // namespace

InletBoundary::InletBoundary(const FarfieldConservative& inflow_U)
    : inflow_U_(inflow_U) {}

void InletBoundary::Apply(DataLayer& layer, const Mesh& mesh, const Axis axis, const Side side) const {
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
            SetPlaneToStateX(U, dst_i, inflow_U_);
            continue;
        }

        if (axis == Axis::Y) {
            const int dst_j = side == Side::Left ? (j0 - 1 - g) : (j1 + g);
            SetPlaneToStateY(U, dst_j, inflow_U_);
            continue;
        }

        const int dst_k = side == Side::Left ? (k0 - 1 - g) : (k1 + g);
        SetPlaneToStateZ(U, dst_k, inflow_U_);
    }
}