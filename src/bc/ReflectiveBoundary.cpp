#include "bc/ReflectiveBoundary.hpp"

#include "data/DataLayer.hpp"

void ReflectiveBoundary::Apply(DataLayer& layer, int axis, Side side) const {
    if (layer.GetDim() >= 2) {
        Apply2D(layer, axis, side);
        return;
    }

    // --- Original 1D logic ---
    const int pad = layer.GetPadding();
    const int core_start = layer.GetCoreStart();
    const int core_end = layer.GetCoreEndExclusive();

    if (side == Side::kLeft) {
        for (int g = 0; g < pad; ++g) {
            int dst = pad - 1 - g;
            int src = core_start + g;
            layer.rho(dst) = layer.rho(src);
            layer.u(dst)   = -layer.u(src);  // reflect normal velocity
            layer.P(dst)   = layer.P(src);
            layer.p(dst)   = -layer.p(src);
            layer.e(dst)   = layer.e(src);
            layer.U(dst)   = layer.U(src);
            layer.V(dst)   = layer.V(src);
            layer.m(dst)   = layer.m(src);
        }
    } else {
        for (int g = 0; g < pad; ++g) {
            int dst = core_end + g;
            int src = core_end - 1 - g;
            layer.rho(dst) = layer.rho(src);
            layer.u(dst)   = -layer.u(src);
            layer.P(dst)   = layer.P(src);
            layer.p(dst)   = -layer.p(src);
            layer.e(dst)   = layer.e(src);
            layer.U(dst)   = layer.U(src);
            layer.V(dst)   = layer.V(src);
            layer.m(dst)   = layer.m(src);
        }
    }
}

void ReflectiveBoundary::Apply2D(DataLayer& layer, int axis, Side side) const {
    const int pad = layer.GetPadding();
    const int cs_x = layer.GetCoreStart(0);
    const int ce_x = layer.GetCoreEndExclusive(0);
    const int cs_y = layer.GetCoreStart(1);
    const int ce_y = layer.GetCoreEndExclusive(1);
    const int tx = layer.GetTotalSize(0);
    const int ty = layer.GetTotalSize(1);

    if (axis == 0) {
        // X-axis reflective: mirror across x-boundary, reflect u, keep v
        for (int j = 0; j < ty; ++j) {
            for (int g = 0; g < pad; ++g) {
                int dst, src;
                if (side == Side::kLeft) {
                    dst = pad - 1 - g;
                    src = cs_x + g;
                } else {
                    dst = ce_x + g;
                    src = ce_x - 1 - g;
                }
                layer.rho(dst, j) = layer.rho(src, j);
                layer.u(dst, j)   = -layer.u(src, j);   // reflect x-velocity
                layer.v(dst, j)   = layer.v(src, j);     // keep y-velocity
                layer.P(dst, j)   = layer.P(src, j);
                layer.p(dst, j)   = -layer.p(src, j);
                layer.q(dst, j)   = layer.q(src, j);
                layer.e(dst, j)   = layer.e(src, j);
                layer.U(dst, j)   = layer.U(src, j);
                layer.V(dst, j)   = layer.V(src, j);
                layer.m(dst, j)   = layer.m(src, j);
            }
        }
    } else {
        // Y-axis reflective: mirror across y-boundary, reflect v, keep u
        for (int i = 0; i < tx; ++i) {
            for (int g = 0; g < pad; ++g) {
                int dst, src;
                if (side == Side::kLeft) {
                    dst = pad - 1 - g;
                    src = cs_y + g;
                } else {
                    dst = ce_y + g;
                    src = ce_y - 1 - g;
                }
                layer.rho(i, dst) = layer.rho(i, src);
                layer.u(i, dst)   = layer.u(i, src);     // keep x-velocity
                layer.v(i, dst)   = -layer.v(i, src);    // reflect y-velocity
                layer.P(i, dst)   = layer.P(i, src);
                layer.p(i, dst)   = layer.p(i, src);
                layer.q(i, dst)   = -layer.q(i, src);
                layer.e(i, dst)   = layer.e(i, src);
                layer.U(i, dst)   = layer.U(i, src);
                layer.V(i, dst)   = layer.V(i, src);
                layer.m(i, dst)   = layer.m(i, src);
            }
        }
    }
}
