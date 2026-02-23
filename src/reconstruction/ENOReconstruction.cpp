#include "reconstruction/ENOReconstruction.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "data/Variables.hpp"  // AxisStride, var::u_*

namespace {

// Notes:
// - We reconstruct in uniform index space along the chosen axis.
// - For face with left cell index (i,j,k):
//     WL: base = left cell,  x_eval = base + 0.5
//     WR: base = right cell, x_eval = base - 0.5
//
// Required padding:
// - For ENO order r, you generally need ng >= r to safely read W along the stencil near core boundaries.

struct LineBuffers5 final {
    std::vector<double> rho;
    std::vector<double> u;
    std::vector<double> v;
    std::vector<double> w;
    std::vector<double> p;

    void Resize(const int N) {
        const std::size_t n = static_cast<std::size_t>(N);
        rho.resize(n);
        u.resize(n);
        v.resize(n);
        w.resize(n);
        p.resize(n);
    }
};

// ---------------------- small math helpers ----------------------

static inline bool IsFinite(const double x) {
    return std::isfinite(x);
}

static double LagrangeInterpolate(const std::vector<double>& f,
                                  const std::vector<int>& stencil,
                                  const double x_eval) {
    const std::size_t m = stencil.size();
    if (m == 0) return 0.0;
    if (m == 1) return f[static_cast<std::size_t>(stencil[0])];

    double result = 0.0;

    for (std::size_t k = 0; k < m; ++k) {
        const double xk = static_cast<double>(stencil[k]);
        double Lk = 1.0;
        for (std::size_t j = 0; j < m; ++j) {
            if (j == k) continue;
            const double xj = static_cast<double>(stencil[j]);
            const double denom = xk - xj;
            if (denom != 0.0) {
                Lk *= (x_eval - xj) / denom;
            } else {
                // degenerate stencil; should not happen for distinct indices
                Lk = 0.0;
            }
        }
        result += f[static_cast<std::size_t>(stencil[k])] * Lk;
    }

    return result;
}

// Compute divided difference of order (m-1) on indices idx[0..m-1] in uniform index space.
static double DividedDifference(const std::vector<double>& f,
                                const std::vector<int>& idx,
                                const int m) {
    if (m <= 0) return 0.0;
    if (m == 1) return f[static_cast<std::size_t>(idx[0])];

    thread_local std::vector<double> g;
    if (static_cast<int>(g.size()) < m) g.resize(static_cast<std::size_t>(m));

    for (int t = 0; t < m; ++t) {
        g[static_cast<std::size_t>(t)] = f[static_cast<std::size_t>(idx[static_cast<std::size_t>(t)])];
    }

    for (int kk = 1; kk < m; ++kk) {
        for (int t = 0; t < m - kk; ++t) {
            const double x0 = static_cast<double>(idx[static_cast<std::size_t>(t)]);
            const double x1 = static_cast<double>(idx[static_cast<std::size_t>(t + kk)]);
            const double denom = x1 - x0;
            g[static_cast<std::size_t>(t)] =
                (denom != 0.0) ? (g[static_cast<std::size_t>(t + 1)] - g[static_cast<std::size_t>(t)]) / denom
                               : 0.0;
        }
    }

    return g[0];
}

// Build ENO stencil of size r around base_index inside [0..N-1].
// The algorithm grows the stencil by choosing left/right extension based on smoothness criterion.
static void BuildENOStencil(const std::vector<double>& f,
                            const int base_index,
                            const int r,
                            const int N,
                            std::vector<int>& stencil_out) {
    stencil_out.clear();

    if (N <= 0 || r <= 0) return;
    if (base_index < 0 || base_index >= N) return;

    stencil_out.push_back(base_index);

    int left = base_index;
    int right = base_index;

    const int max_r = std::min(r, N);

    thread_local std::vector<int> cand_left;
    thread_local std::vector<int> cand_right;

    for (int level = 1; level < max_r; ++level) {
        const int try_left = left - 1;
        const int try_right = right + 1;

        const bool has_left = (try_left >= 0);
        const bool has_right = (try_right < N);

        if (!has_left && !has_right) break;

        int choose = 0; // -1 left, +1 right

        if (has_left && !has_right) {
            choose = -1;
        } else if (!has_left && has_right) {
            choose = +1;
        } else {
            // Compare |divided difference| of order (level) (i.e., m = level+1 points)
            cand_left.clear();
            cand_left.reserve(static_cast<std::size_t>(level + 1));
            cand_left.push_back(try_left);
            cand_left.insert(cand_left.end(), stencil_out.begin(), stencil_out.end());
            const double dd_left = std::abs(DividedDifference(f, cand_left, level + 1));

            cand_right = stencil_out;
            cand_right.push_back(try_right);
            const double dd_right = std::abs(DividedDifference(f, cand_right, level + 1));

            if (dd_left < dd_right) choose = -1;
            else if (dd_right < dd_left) choose = +1;
            else choose = -1; // deterministic tie-break
        }

        if (choose < 0) {
            left = try_left;
            stencil_out.insert(stencil_out.begin(), left);
        } else {
            right = try_right;
            stencil_out.push_back(right);
        }
    }
}

static double ENOReconstructScalar(const std::vector<double>& f,
                                  const int base_index,
                                  const double x_eval,
                                  const int r,
                                  const int N) {
    thread_local std::vector<int> stencil;
    BuildENOStencil(f, base_index, r, N, stencil);
    return LagrangeInterpolate(f, stencil, x_eval);
}

// Fill 5 primitive line buffers along an axis, starting at (start_i,start_j,start_k),
// reading N points with stride st.
static void FillLineBuffers5(const xt::xtensor<double, 4>& W,
                             const AxisStride& st,
                             const int start_i,
                             const int start_j,
                             const int start_k,
                             const int N,
                             LineBuffers5& out) {
    out.Resize(N);

    for (int s = 0; s < N; ++s) {
        const int ii = start_i + s * st.di;
        const int jj = start_j + s * st.dj;
        const int kk = start_k + s * st.dk;

        out.rho[static_cast<std::size_t>(s)] = W(var::u_rho, ii, jj, kk);
        out.u  [static_cast<std::size_t>(s)] = W(var::u_u,   ii, jj, kk);
        out.v  [static_cast<std::size_t>(s)] = W(var::u_v,   ii, jj, kk);
        out.w  [static_cast<std::size_t>(s)] = W(var::u_w,   ii, jj, kk);
        out.p  [static_cast<std::size_t>(s)] = W(var::u_P,   ii, jj, kk);
    }
}

static void SanitizePrimitive(PrimitiveCell& w) {
    // ENO may overshoot; keep states finite. Floors are handled later in Riemann/Variables,
    // but we avoid NaNs here.
    if (!IsFinite(w.rho)) w.rho = 0.0;
    if (!IsFinite(w.u))   w.u   = 0.0;
    if (!IsFinite(w.v))   w.v   = 0.0;
    if (!IsFinite(w.w))   w.w   = 0.0;
    if (!IsFinite(w.P))   w.P   = 0.0;
}

} // namespace

ENOReconstruction::ENOReconstruction(const int order) : order_(order) {
    if (order_ < 1) {
        throw std::invalid_argument("ENOReconstruction: order must be >= 1");
    }
}

void ENOReconstruction::SetOrder(const int order) {
    if (order < 1) {
        throw std::invalid_argument("ENOReconstruction: order must be >= 1");
    }
    order_ = order;
}

void ENOReconstruction::ReconstructFace(const xt::xtensor<double, 4>& W,
                                        const Axis axis,
                                        const int i, const int j, const int k,
                                        PrimitiveCell& WL,
                                        PrimitiveCell& WR) const {
    const int r = order_;
    const AxisStride st = AxisStride::FromAxis(axis);

    // Build a local uniform-index buffer around the face:
    // local s = 0 corresponds to global index (i - (r-1)) along the axis.
    //
    // We take N = 2*r + 1 points:
    //   covers global indices [i-(r-1), ..., i+(r+1)] along the axis.
    // This is safe for the stencil growth and gives a little extra room.
    const int N = 2 * r + 1;

    const int start_i = i - (r - 1) * st.di;
    const int start_j = j - (r - 1) * st.dj;
    const int start_k = k - (r - 1) * st.dk;

    const int base_left = (r - 1);     // maps to global left cell (i,j,k)
    const int base_right = base_left + 1; // maps to global right cell

    const double xL = static_cast<double>(base_left) + 0.5;
    const double xR = static_cast<double>(base_right) - 0.5;

    thread_local LineBuffers5 line;

    // NOTE: This assumes caller ensured sufficient padding (ng >= r).
    // If you want a defensive check here, you'd need layer sizes; we only have W.
    FillLineBuffers5(W, st, start_i, start_j, start_k, N, line);

    WL.rho = ENOReconstructScalar(line.rho, base_left,  xL, r, N);
    WL.u   = ENOReconstructScalar(line.u,   base_left,  xL, r, N);
    WL.v   = ENOReconstructScalar(line.v,   base_left,  xL, r, N);
    WL.w   = ENOReconstructScalar(line.w,   base_left,  xL, r, N);
    WL.P   = ENOReconstructScalar(line.p,   base_left,  xL, r, N);

    WR.rho = ENOReconstructScalar(line.rho, base_right, xR, r, N);
    WR.u   = ENOReconstructScalar(line.u,   base_right, xR, r, N);
    WR.v   = ENOReconstructScalar(line.v,   base_right, xR, r, N);
    WR.w   = ENOReconstructScalar(line.w,   base_right, xR, r, N);
    WR.P   = ENOReconstructScalar(line.p,   base_right, xR, r, N);

    SanitizePrimitive(WL);
    SanitizePrimitive(WR);
}