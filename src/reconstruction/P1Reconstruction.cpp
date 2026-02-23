#include "reconstruction/P1Reconstruction.hpp"

#include <algorithm>
#include <cmath>

namespace {
    inline PrimitiveCell LoadCellW(const xt::xtensor<double, 4>& W, int i, int j, int k) {
        PrimitiveCell w;
        w.rho = W(var::u_rho, i, j, k);
        w.u = W(var::u_u, i, j, k);
        w.v = W(var::u_v, i, j, k);
        w.w = W(var::u_w, i, j, k);
        w.P = W(var::u_P, i, j, k);
        return w;
    }

    inline void AddScaled(PrimitiveCell& out, const PrimitiveCell& a, const double s) {
        out.rho += s * a.rho;
        out.u += s * a.u;
        out.v += s * a.v;
        out.w += s * a.w;
        out.P += s * a.P;
    }
} // namespace

auto P1Reconstruction::Minmod(const double a, const double b) -> double {
    if (a * b <= 0.0) return 0.0;
    return (std::abs(a) < std::abs(b)) ? a : b;
}

auto P1Reconstruction::Mc(const double a, const double b) -> double {
    if (a * b <= 0.0) return 0.0;
    const double sign = (a > 0.0) ? 1.0 : -1.0;
    const double aa = std::abs(a);
    const double bb = std::abs(b);
    return sign * std::min({2.0 * aa, 2.0 * bb, 0.5 * (aa + bb)});
}

auto P1Reconstruction::Superbee(const double a, const double b) -> double {
    if (a * b <= 0.0) return 0.0;
    const double sign = (a > 0.0) ? 1.0 : -1.0;
    const double aa = std::abs(a);
    const double bb = std::abs(b);
    const double t1 = std::min(2.0 * aa, bb);
    const double t2 = std::min(aa, 2.0 * bb);
    return sign * std::max(t1, t2);
}

auto P1Reconstruction::ApplyLimiter(const double a, const double b) const -> double {
    switch (limiter_type_) {
    case LimiterType::kMc: return Mc(a, b);
    case LimiterType::kSuperbee: return Superbee(a, b);
    case LimiterType::kMinmod:
    default: return Minmod(a, b);
    }
}

void P1Reconstruction::ReconstructFace(const xt::xtensor<double, 4>& W,
                                       const Axis axis,
                                       const int i, const int j, const int k,
                                       PrimitiveCell& WL,
                                       PrimitiveCell& WR) const {
    const AxisStride s = AxisStride::FromAxis(axis);

    // Left cell index = (i,j,k); right cell = (i+ix, j+iy, k+iz)
    const PrimitiveCell Wm1 = LoadCellW(W, i - s.di, j - s.dj, k - s.dk);
    const PrimitiveCell W0 = LoadCellW(W, i, j, k);
    const PrimitiveCell Wp1 = LoadCellW(W, i + s.di, j + s.dj, k + s.dk);
    const PrimitiveCell Wp2 = LoadCellW(W, i + 2 * s.di, j + 2 * s.dj, k + 2 * s.dk);

    // slope in left cell (centered around W0): limiter(W0 - Wm1, Wp1 - W0)
    PrimitiveCell slopeL;
    slopeL.rho = ApplyLimiter(W0.rho - Wm1.rho, Wp1.rho - W0.rho);
    slopeL.u = ApplyLimiter(W0.u - Wm1.u, Wp1.u - W0.u);
    slopeL.v = ApplyLimiter(W0.v - Wm1.v, Wp1.v - W0.v);
    slopeL.w = ApplyLimiter(W0.w - Wm1.w, Wp1.w - W0.w);
    slopeL.P = ApplyLimiter(W0.P - Wm1.P, Wp1.P - W0.P);

    // slope in right cell (centered around Wp1): limiter(Wp1 - W0, Wp2 - Wp1)
    PrimitiveCell slopeR;
    slopeR.rho = ApplyLimiter(Wp1.rho - W0.rho, Wp2.rho - Wp1.rho);
    slopeR.u = ApplyLimiter(Wp1.u - W0.u, Wp2.u - Wp1.u);
    slopeR.v = ApplyLimiter(Wp1.v - W0.v, Wp2.v - Wp1.v);
    slopeR.w = ApplyLimiter(Wp1.w - W0.w, Wp2.w - Wp1.w);
    slopeR.P = ApplyLimiter(Wp1.P - W0.P, Wp2.P - Wp1.P);

    // Reconstruct
    WL = W0;
    WR = Wp1;

    AddScaled(WL, slopeL, 0.5);
    AddScaled(WR, slopeR, -0.5);
}
