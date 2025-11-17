#include "bc/SymmetryBoundary.hpp"

#include "data/DataLayer.hpp"

void SymmetryBoundary::Apply(DataLayer& layer, int axis, Side side) const {
    (void)axis;

    const int padding = layer.GetPadding();
    const int core_start = layer.GetCoreStart();
    const int core_end = layer.GetCoreEndExclusive();

    if (side == Side::kLeft) {
        for (int g = 0; g < padding; ++g) {
            int src = core_start + g;
            int dst = padding - 1 - g;

            layer.rho(dst) = layer.rho(src);
            layer.P(dst) = layer.P(src);
            layer.p(dst) = layer.p(src);
            layer.e(dst) = layer.e(src);
            layer.U(dst) = layer.U(src);
            layer.V(dst) = layer.V(src);
            layer.m(dst) = layer.m(src);
            layer.xb(dst) = layer.xb(src);
            layer.xc(dst) = layer.xc(src);

            layer.u(dst) = -layer.u(src);
        }
    } else {
        for (int g = 0; g < padding; ++g) {
            int src = core_end - 1 - g;
            int dst = core_end + g;

            layer.rho(dst) = layer.rho(src);
            layer.P(dst) = layer.P(src);
            layer.p(dst) = layer.p(src);
            layer.e(dst) = layer.e(src);
            layer.U(dst) = layer.U(src);
            layer.V(dst) = layer.V(src);
            layer.m(dst) = layer.m(src);
            layer.xb(dst) = layer.xb(src);
            layer.xc(dst) = layer.xc(src);

            layer.u(dst) = -layer.u(src);
        }
    }
}
