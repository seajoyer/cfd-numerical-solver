#include "bc/OutletBoundary.hpp"

#include "data/DataLayer.hpp"

void OutletBoundary::Apply(DataLayer& layer, int axis, Side side) const {
    if (layer.GetDim() >= 2) {
        // 2D: zero-gradient extrapolation along the given axis
        const int pad = layer.GetPadding();
        const int cs_x = layer.GetCoreStart(0);
        const int ce_x = layer.GetCoreEndExclusive(0);
        const int cs_y = layer.GetCoreStart(1);
        const int ce_y = layer.GetCoreEndExclusive(1);
        const int tx = layer.GetTotalSize(0);
        const int ty = layer.GetTotalSize(1);

        if (axis == 0) {
            for (int j = 0; j < ty; ++j) {
                int src = (side == Side::kLeft) ? cs_x : (ce_x - 1);
                for (int g = 0; g < pad; ++g) {
                    int dst = (side == Side::kLeft) ? g : (ce_x + g);
                    layer.rho(dst, j) = layer.rho(src, j);
                    layer.u(dst, j)   = layer.u(src, j);
                    layer.v(dst, j)   = layer.v(src, j);
                    layer.P(dst, j)   = layer.P(src, j);
                    layer.p(dst, j)   = layer.p(src, j);
                    layer.q(dst, j)   = layer.q(src, j);
                    layer.e(dst, j)   = layer.e(src, j);
                    layer.U(dst, j)   = layer.U(src, j);
                    layer.V(dst, j)   = layer.V(src, j);
                    layer.m(dst, j)   = layer.m(src, j);
                }
            }
        } else {
            for (int i = 0; i < tx; ++i) {
                int src = (side == Side::kLeft) ? cs_y : (ce_y - 1);
                for (int g = 0; g < pad; ++g) {
                    int dst = (side == Side::kLeft) ? g : (ce_y + g);
                    layer.rho(i, dst) = layer.rho(i, src);
                    layer.u(i, dst)   = layer.u(i, src);
                    layer.v(i, dst)   = layer.v(i, src);
                    layer.P(i, dst)   = layer.P(i, src);
                    layer.p(i, dst)   = layer.p(i, src);
                    layer.q(i, dst)   = layer.q(i, src);
                    layer.e(i, dst)   = layer.e(i, src);
                    layer.U(i, dst)   = layer.U(i, src);
                    layer.V(i, dst)   = layer.V(i, src);
                    layer.m(i, dst)   = layer.m(i, src);
                }
            }
        }
        return;
    }

    // --- Original 1D logic ---
    (void)axis;

    const int pad = layer.GetPadding();
    const int core_start = layer.GetCoreStart();
    const int core_end = layer.GetCoreEndExclusive();

    if (side == Side::kLeft) {
        int src = core_start;
        for (int g = 0; g < pad; ++g) {
            int dst = g;
            layer.rho(dst) = layer.rho(src);
            layer.u(dst) = layer.u(src);
            layer.P(dst) = layer.P(src);
            layer.p(dst) = layer.p(src);
            layer.e(dst) = layer.e(src);
            layer.U(dst) = layer.U(src);
            layer.V(dst) = layer.V(src);
            layer.m(dst) = layer.m(src);
            layer.xb(dst) = layer.xb(src);
            layer.xc(dst) = layer.xc(src);
        }
    } else {
        int src = core_end - 1;
        for (int g = 0; g < pad; ++g) {
            int dst = core_end + g;
            layer.rho(dst) = layer.rho(src);
            layer.u(dst) = layer.u(src);
            layer.P(dst) = layer.P(src);
            layer.p(dst) = layer.p(src);
            layer.e(dst) = layer.e(src);
            layer.U(dst) = layer.U(src);
            layer.V(dst) = layer.V(src);
            layer.m(dst) = layer.m(src);
            layer.xb(dst) = layer.xb(src);
            layer.xc(dst) = layer.xc(src);
        }
    }
}
