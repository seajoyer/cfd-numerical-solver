#include "bc/SymmetryBoundary.hpp"
#include "data/DataLayer.hpp"


void SymmetryBoundary::Apply(DataLayer &layer, int axis, Side side) const {
    (void) axis;

    const int pad = layer.GetPadding();
    const int coreStart = layer.GetCoreStart();
    const int coreEnd = layer.GetCoreEndExclusive();

    if (side == Side::Min) {
        for (int g = 0; g < pad; ++g) {
            int src = coreStart + g;
            int dst = pad - 1 - g;

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
        for (int g = 0; g < pad; ++g) {
            int src = coreEnd - 1 - g;
            int dst = coreEnd + g;

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
