#include "bc/InletBoundary.hpp"

#include "data/DataLayer.hpp"

InletBoundary::InletBoundary(double rho_in, double u_in, double p_in)
    : rho_in_(rho_in), u_in_(u_in), v_in_(0.0), p_in_(p_in) {}

InletBoundary::InletBoundary(double rho_in, double u_in, double v_in, double p_in)
    : rho_in_(rho_in), u_in_(u_in), v_in_(v_in), p_in_(p_in) {}

void InletBoundary::Apply(DataLayer& layer, int axis, Side side) const {
    if (layer.GetDim() >= 2) {
        // ===== 2D inlet =====
        const int pad = layer.GetPadding();
        const int cs_x = layer.GetCoreStart(0);
        const int ce_x = layer.GetCoreEndExclusive(0);
        const int cs_y = layer.GetCoreStart(1);
        const int ce_y = layer.GetCoreEndExclusive(1);
        const int tx = layer.GetTotalSize(0);
        const int ty = layer.GetTotalSize(1);

        if (axis == 0) {
            // X-axis boundary: normal velocity is u
            for (int j = 0; j < ty; ++j) {
                // Check normal velocity at nearest core cell
                int ref;
                if (side == Side::kLeft) {
                    ref = cs_x;
                } else {
                    ref = ce_x - 1;
                }
                double u_ref = layer.u(ref, j);
                bool inward = (side == Side::kLeft) ? (u_ref > 0.0) : (u_ref < 0.0);

                for (int g = 0; g < pad; ++g) {
                    int dst = (side == Side::kLeft) ? g : (ce_x + g);
                    int src = ref;

                    if (inward) {
                        // Impose prescribed inlet values
                        layer.rho(dst, j) = rho_in_;
                        layer.u(dst, j)   = u_in_;
                        layer.v(dst, j)   = v_in_;
                        layer.P(dst, j)   = p_in_;
                        // Copy derived quantities from nearest core cell
                        layer.p(dst, j) = layer.p(src, j);
                        layer.q(dst, j) = layer.q(src, j);
                        layer.e(dst, j) = layer.e(src, j);
                        layer.U(dst, j) = layer.U(src, j);
                        layer.V(dst, j) = layer.V(src, j);
                        layer.m(dst, j) = layer.m(src, j);
                    } else {
                        // Zero-gradient (outlet behavior)
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
            }
        } else {
            // Y-axis boundary: normal velocity is v
            for (int i = 0; i < tx; ++i) {
                int ref;
                if (side == Side::kLeft) {
                    ref = cs_y;
                } else {
                    ref = ce_y - 1;
                }
                double v_ref = layer.v(i, ref);
                bool inward = (side == Side::kLeft) ? (v_ref > 0.0) : (v_ref < 0.0);

                for (int g = 0; g < pad; ++g) {
                    int dst = (side == Side::kLeft) ? g : (ce_y + g);
                    int src = ref;

                    if (inward) {
                        layer.rho(i, dst) = rho_in_;
                        layer.u(i, dst)   = u_in_;
                        layer.v(i, dst)   = v_in_;
                        layer.P(i, dst)   = p_in_;
                        layer.p(i, dst) = layer.p(i, src);
                        layer.q(i, dst) = layer.q(i, src);
                        layer.e(i, dst) = layer.e(i, src);
                        layer.U(i, dst) = layer.U(i, src);
                        layer.V(i, dst) = layer.V(i, src);
                        layer.m(i, dst) = layer.m(i, src);
                    } else {
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
        }
        return;
    }

    // ===== Original 1D logic =====
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
                layer.rho(dst) = rho_in_;
                layer.u(dst) = u_in_;
                layer.P(dst) = p_in_;

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
            for (int g = 0; g < pad; ++g) {
                int dst = g;
                int src = core_start;
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
                layer.rho(dst) = rho_in_;
                layer.u(dst) = u_in_;
                layer.P(dst) = p_in_;

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
            for (int g = 0; g < pad; ++g) {
                int dst = core_end + g;
                int src = core_end - 1;
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
