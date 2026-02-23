#include "reconstruction/P0Reconstruction.hpp"

namespace {
    PrimitiveCell LoadCellW(const xt::xtensor<double, 4>& W, const int i, const int j, const int k) {
        PrimitiveCell w;
        w.rho = W(var::u_rho, i, j, k);
        w.u = W(var::u_u, i, j, k);
        w.v = W(var::u_v, i, j, k);
        w.w = W(var::u_w, i, j, k);
        w.P = W(var::u_P, i, j, k);
        return w;
    }
} // namespace

void P0Reconstruction::ReconstructFace(const xt::xtensor<double, 4>& W,
                                       const Axis axis,
                                       const int i, const int j, const int k,
                                       PrimitiveCell& WL,
                                       PrimitiveCell& WR) const {
    WL = LoadCellW(W, i, j, k);

    if (axis == Axis::X) {
        WR = LoadCellW(W, i + 1, j, k);
        return;
    }
    if (axis == Axis::Y) {
        WR = LoadCellW(W, i, j + 1, k);
        return;
    }
    // Axis::Z
    WR = LoadCellW(W, i, j, k + 1);
}
