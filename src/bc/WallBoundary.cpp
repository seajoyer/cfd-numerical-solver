#include "bc/WallBoundary.hpp"
#include "data/DataLayer.hpp"

void WallBoundary::Apply(DataLayer &layer, int axis, Side side) const {
    (void) axis;

    const int padding = layer.GetPadding();
    const int core_start = layer.GetCoreStart();
    const int core_end = layer.GetCoreEndExclusive();

    if (side == Side::kLeft) {
        const int src = core_start;
        for (int g = 0; g < padding; ++g) {
            const int dst = padding - 1 - g;

            const double rho = layer.rho(src);
            const double P = layer.P(src);

            layer.rho(dst) = rho;
            layer.u(dst) = 0.0;
            layer.P(dst) = P;

            layer.p(dst) = 0.0;

            layer.U(dst) = layer.U(src);
            layer.e(dst) = layer.e(src);
            layer.V(dst) = layer.V(src);
            layer.m(dst) = layer.m(src);

            layer.xb(dst) = layer.xb(src);
            layer.xc(dst) = layer.xc(src);
        }
    } else {
        const int src = core_end - 1;
        for (int g = 0; g < padding; ++g) {
            const int dst = core_end + g;

            const double rho = layer.rho(src);
            const double P = layer.P(src);

            layer.rho(dst) = rho;
            layer.u(dst) = 0.0;
            layer.P(dst) = P;

            layer.p(dst) = 0.0;

            layer.U(dst) = layer.U(src);
            layer.e(dst) = layer.e(src);
            layer.V(dst) = layer.V(src);
            layer.m(dst) = layer.m(src);

            layer.xb(dst) = layer.xb(src);
            layer.xc(dst) = layer.xc(src);
        }
    }
}
