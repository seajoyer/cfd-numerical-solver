#include "bc/ReflectiveBoundary.hpp"
#include "data/DataLayer.hpp"


void ReflectiveBoundary::Apply(DataLayer &layer, int axis, Side side) const {
    (void) axis;

    const int pad = layer.GetPadding();
    const int core_start = layer.GetCoreStart();
    const int core_end = layer.GetCoreEndExclusive();

    if (side == Side::kLeft) {
        for (int g = 0; g < pad; ++g) {
            const int mirror = core_start + g;
            const int dst = pad - 1 - g;

            layer.rho(dst) = layer.rho(mirror);
            layer.P(dst) = layer.P(mirror);
            layer.p(dst) = layer.p(mirror);
            layer.e(dst) = layer.e(mirror);
            layer.U(dst) = layer.U(mirror);
            layer.V(dst) = layer.V(mirror);
            layer.m(dst) = layer.m(mirror);
            layer.xb(dst) = layer.xb(mirror);
            layer.xc(dst) = layer.xc(mirror);

            layer.u(dst) = -layer.u(mirror);
        }
    } else {
        for (int g = 0; g < pad; ++g) {
            const int mirror = core_end - 1 - g;
            const int dst = core_end + g;

            layer.rho(dst) = layer.rho(mirror);
            layer.P(dst) = layer.P(mirror);
            layer.p(dst) = layer.p(mirror);
            layer.e(dst) = layer.e(mirror);
            layer.U(dst) = layer.U(mirror);
            layer.V(dst) = layer.V(mirror);
            layer.m(dst) = layer.m(mirror);
            layer.xb(dst) = layer.xb(mirror);
            layer.xc(dst) = layer.xc(mirror);

            layer.u(dst) = -layer.u(mirror);
        }
    }
}
