#include "bc/FreeStreamBoundary.hpp"

#include "data/DataLayer.hpp"

FreeStreamBoundary::FreeStreamBoundary(double rho_inf, double u_inf, double p_inf)
    : rho_inf_(rho_inf), u_inf_(u_inf), p_inf_(p_inf) {}

void FreeStreamBoundary::Apply(DataLayer& layer, int axis, Side side) const {
    (void)axis;

    const int pad = layer.GetPadding();
    const int core_start = layer.GetCoreStart();
    const int core_end = layer.GetCoreEndExclusive();

    bool inward = false;
    if (side == Side::kLeft) {
        double ui = layer.u(core_start);
        inward = (ui > 0.0);
    } else {
        double ui = layer.u(core_end - 1);
        inward = (ui < 0.0);
    }

    if (side == Side::kLeft) {
        if (inward) {
            for (int g = 0; g < pad; ++g) {
                int dst = g;
                layer.rho(dst) = rho_inf_;
                layer.u(dst) = u_inf_;
                layer.P(dst) = p_inf_;

                int src = core_start;
                layer.p(dst) = layer.p(src);
                layer.e(dst) = layer.e(src);
                layer.U(dst) = layer.U(src);
                layer.V(dst) = layer.V(src);
                layer.m(dst) = layer.m(src);
                layer.xb(dst) = layer.xb(src);
                layer.xc(dst) = layer.xc(src);
            }
        } else {
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
        }
    } else {
        if (inward) {
            for (int g = 0; g < pad; ++g) {
                int dst = core_end + g;
                layer.rho(dst) = rho_inf_;
                layer.u(dst) = u_inf_;
                layer.P(dst) = p_inf_;

                int src = core_end - 1;
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
}
