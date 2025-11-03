#include "bc/OutletBoundary.hpp"
#include "data/DataLayer.hpp"


void OutletBoundary::Apply(DataLayer &layer, int axis, Side side) const {
    (void) axis;

    const int pad = layer.GetPadding();
    const int coreStart = layer.GetCoreStart();
    const int coreEnd = layer.GetCoreEndExclusive();

    if (side == Side::Min) {
        int src = coreStart;
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
        int src = coreEnd - 1;
        for (int g = 0; g < pad; ++g) {
            int dst = coreEnd + g;
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