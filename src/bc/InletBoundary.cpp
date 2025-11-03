#include "bc/InletBoundary.hpp"
#include "data/DataLayer.hpp"


void InletBoundary::Apply(DataLayer &layer, int axis, Side side) const {
    (void) axis;

    const int pad = layer.GetPadding();
    const int coreStart = layer.GetCoreStart();
    const int coreEnd = layer.GetCoreEndExclusive();

    bool inward = false;
    if (side == Side::Min) {
        double ui = layer.u(coreStart);
        inward = (ui > 0.0);
    } else {
        double ui = layer.u(coreEnd - 1);
        inward = (ui < 0.0);
    }

    if (side == Side::Min) {
        if (inward) {
            for (int g = 0; g < pad; ++g) {
                int dst = g;
                layer.rho(dst) = rhoIn;
                layer.u(dst) = uIn;
                layer.P(dst) = pIn;


                int src = coreStart;
                layer.p(dst) = layer.p(src);
                layer.e(dst) = layer.e(src);
                layer.U(dst) = layer.U(src);
                layer.V(dst) = layer.V(src);
                layer.m(dst) = layer.m(src);
                layer.xb(dst) = layer.xb(src);
                layer.xc(dst) = layer.xc(src);
            }
        } else {
            for (int g = 0; g < pad; ++g) {
                int dst = g;
                int src = coreStart;
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
    } else {
        if (inward) {
            for (int g = 0; g < pad; ++g) {
                int dst = coreEnd + g;
                layer.rho(dst) = rhoIn;
                layer.u(dst) = uIn;
                layer.P(dst) = pIn;

                int src = coreEnd - 1;
                layer.p(dst) = layer.p(src);
                layer.e(dst) = layer.e(src);
                layer.U(dst) = layer.U(src);
                layer.V(dst) = layer.V(src);
                layer.m(dst) = layer.m(src);
                layer.xb(dst) = layer.xb(src);
                layer.xc(dst) = layer.xc(src);
            }
        } else {
            // как Outlet
            for (int g = 0; g < pad; ++g) {
                int dst = coreEnd + g;
                int src = coreEnd - 1;
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
}