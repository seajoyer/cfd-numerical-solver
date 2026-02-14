#include "bc/FreeStreamBoundary.hpp"

#include "data/DataLayer.hpp"

FreeStreamBoundary::FreeStreamBoundary(double rho_inf, double u_inf, double p_inf)
    : rho_inf_(rho_inf), u_inf_(u_inf), p_inf_(p_inf), v_inf_(0.0) {}

FreeStreamBoundary::FreeStreamBoundary(double rho_inf, double u_inf, double v_inf, double p_inf)
    : rho_inf_(rho_inf), u_inf_(u_inf), p_inf_(p_inf), v_inf_(v_inf) {}

void FreeStreamBoundary::Apply(DataLayer& layer, int axis, Side side) const {
    if (layer.GetDim() >= 2) {
        Apply2D(layer, axis, side);
        return;
    }

    // --- Original 1D logic (unchanged) ---
    const int pad = layer.GetPadding();
    const int core_start = layer.GetCoreStart();
    const int core_end = layer.GetCoreEndExclusive();

    bool inward = false;
    if (side == Side::kLeft) {
        inward = (layer.u(core_start) > 0.0);
    } else {
        inward = (layer.u(core_end - 1) < 0.0);
    }

    if (side == Side::kLeft) {
        for (int g = 0; g < pad; ++g) {
            int dst = g;
            int src = core_start;
            if (inward) {
                layer.rho(dst) = rho_inf_;
                layer.u(dst)   = u_inf_;
                layer.P(dst)   = p_inf_;
            } else {
                layer.rho(dst) = layer.rho(src);
                layer.u(dst)   = layer.u(src);
                layer.P(dst)   = layer.P(src);
            }
            layer.p(dst)  = layer.p(src);
            layer.e(dst)  = layer.e(src);
            layer.U(dst)  = layer.U(src);
            layer.V(dst)  = layer.V(src);
            layer.m(dst)  = layer.m(src);
            layer.xb(dst) = layer.xb(src);
            layer.xc(dst) = layer.xc(src);
        }
    } else {
        for (int g = 0; g < pad; ++g) {
            int dst = core_end + g;
            int src = core_end - 1;
            if (inward) {
                layer.rho(dst) = rho_inf_;
                layer.u(dst)   = u_inf_;
                layer.P(dst)   = p_inf_;
            } else {
                layer.rho(dst) = layer.rho(src);
                layer.u(dst)   = layer.u(src);
                layer.P(dst)   = layer.P(src);
            }
            layer.p(dst)  = layer.p(src);
            layer.e(dst)  = layer.e(src);
            layer.U(dst)  = layer.U(src);
            layer.V(dst)  = layer.V(src);
            layer.m(dst)  = layer.m(src);
            layer.xb(dst) = layer.xb(src);
            layer.xc(dst) = layer.xc(src);
        }
    }
}

void FreeStreamBoundary::Apply2D(DataLayer& layer, int axis, Side side) const {
    const int pad = layer.GetPadding();
    const int cs_x = layer.GetCoreStart(0);
    const int ce_x = layer.GetCoreEndExclusive(0);
    const int cs_y = layer.GetCoreStart(1);
    const int ce_y = layer.GetCoreEndExclusive(1);
    const int tx = layer.GetTotalSize(0);
    const int ty = layer.GetTotalSize(1);

    if (axis == 0) {
        // X-axis: fill ghost columns (left or right)
        for (int j = 0; j < ty; ++j) {
            int src_i = (side == Side::kLeft) ? cs_x : (ce_x - 1);
            double normal_vel = layer.u(src_i, j);
            bool inward = (side == Side::kLeft) ? (normal_vel > 0.0) : (normal_vel < 0.0);

            for (int g = 0; g < pad; ++g) {
                int dst_i = (side == Side::kLeft) ? g : (ce_x + g);

                if (inward) {
                    layer.rho(dst_i, j) = rho_inf_;
                    layer.u(dst_i, j)   = u_inf_;
                    layer.v(dst_i, j)   = v_inf_;
                    layer.P(dst_i, j)   = p_inf_;
                } else {
                    layer.rho(dst_i, j) = layer.rho(src_i, j);
                    layer.u(dst_i, j)   = layer.u(src_i, j);
                    layer.v(dst_i, j)   = layer.v(src_i, j);
                    layer.P(dst_i, j)   = layer.P(src_i, j);
                }
                // Copy auxiliary fields
                layer.p(dst_i, j) = layer.p(src_i, j);
                layer.q(dst_i, j) = layer.q(src_i, j);
                layer.e(dst_i, j) = layer.e(src_i, j);
                layer.U(dst_i, j) = layer.U(src_i, j);
                layer.V(dst_i, j) = layer.V(src_i, j);
                layer.m(dst_i, j) = layer.m(src_i, j);
            }
        }
    } else {
        // Y-axis: fill ghost rows (bottom or top)
        for (int i = 0; i < tx; ++i) {
            int src_j = (side == Side::kLeft) ? cs_y : (ce_y - 1);
            double normal_vel = layer.v(i, src_j);
            bool inward = (side == Side::kLeft) ? (normal_vel > 0.0) : (normal_vel < 0.0);

            for (int g = 0; g < pad; ++g) {
                int dst_j = (side == Side::kLeft) ? g : (ce_y + g);

                if (inward) {
                    layer.rho(i, dst_j) = rho_inf_;
                    layer.u(i, dst_j)   = u_inf_;
                    layer.v(i, dst_j)   = v_inf_;
                    layer.P(i, dst_j)   = p_inf_;
                } else {
                    layer.rho(i, dst_j) = layer.rho(i, src_j);
                    layer.u(i, dst_j)   = layer.u(i, src_j);
                    layer.v(i, dst_j)   = layer.v(i, src_j);
                    layer.P(i, dst_j)   = layer.P(i, src_j);
                }
                layer.p(i, dst_j) = layer.p(i, src_j);
                layer.q(i, dst_j) = layer.q(i, src_j);
                layer.e(i, dst_j) = layer.e(i, src_j);
                layer.U(i, dst_j) = layer.U(i, src_j);
                layer.V(i, dst_j) = layer.V(i, src_j);
                layer.m(i, dst_j) = layer.m(i, src_j);
            }
        }
    }
}
