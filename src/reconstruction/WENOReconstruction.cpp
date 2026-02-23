#include "reconstruction/WENOReconstruction.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "data/Variables.hpp"  // AxisStride, var::u_*

WENOReconstruction::WENOReconstruction(const int order) : order_(order) {
    SetOrder(order);
}

void WENOReconstruction::SetOrder(const int order) {
    if (order != 3 && order != 5) {
        throw std::invalid_argument("WENOReconstruction: supported orders are 3 and 5");
    }
    order_ = order;
}

void WENOReconstruction::SetEpsilon(const double epsilon) {
    if (!(epsilon > 0.0) || !std::isfinite(epsilon)) {
        throw std::invalid_argument("WENOReconstruction: epsilon must be positive");
    }
    epsilon_ = epsilon;
}

void WENOReconstruction::SetNonlinearWeightPower(const int p) {
    if (p < 1) {
        throw std::invalid_argument("WENOReconstruction: p must be >= 1");
    }
    p_ = p;
}

// ----------------------- private helpers (methods) -----------------------

double WENOReconstruction::PowInt(const double x, const int p) {
    if (p == 1) return x;
    if (p == 2) return x * x;
    return std::pow(x, static_cast<double>(p));
}

void WENOReconstruction::ComputeNonlinearWeights(const double* beta,
                                                 const double* d,
                                                 const int r,
                                                 const double eps,
                                                 const int p,
                                                 double* omega) {
    double alpha_sum = 0.0;
    for (int k = 0; k < r; ++k) {
        const double denom = PowInt(eps + beta[k], p);
        const double alpha = d[k] / denom;
        omega[k] = alpha;
        alpha_sum += alpha;
    }

    if (!(alpha_sum > 0.0) || !std::isfinite(alpha_sum)) {
        for (int k = 0; k < r; ++k) omega[k] = d[k];
        return;
    }

    const double inv = 1.0 / alpha_sum;
    for (int k = 0; k < r; ++k) omega[k] *= inv;
}

double WENOReconstruction::LoadAt(const xt::xtensor<double, 4>& W,
                                  const std::size_t v,
                                  const int i, const int j, const int k,
                                  const AxisStride& st,
                                  const int offset) {
    return W(v,
             i + offset * st.di,
             j + offset * st.dj,
             k + offset * st.dk);
}

void WENOReconstruction::StorePrimitiveComponent(PrimitiveCell& w, const std::size_t v, const double val) {
    if (v == var::u_rho) w.rho = val;
    else if (v == var::u_u) w.u = val;
    else if (v == var::u_v) w.v = val;
    else if (v == var::u_w) w.w = val;
    else w.P = val;
}

// ----- WENO3 (JS) -----

double WENOReconstruction::Weno3Left(const double fim1, const double fi, const double fip1,
                                     const double eps, const int p) {
    const double d[2] = {1.0 / 3.0, 2.0 / 3.0};

    double beta[2];
    beta[0] = (fi - fim1) * (fi - fim1);
    beta[1] = (fip1 - fi) * (fip1 - fi);

    const double q0 = -0.5 * fim1 + 1.5 * fi;
    const double q1 =  0.5 * fi   + 0.5 * fip1;

    double omega[2];
    ComputeNonlinearWeights(beta, d, 2, eps, p, omega);

    return omega[0] * q0 + omega[1] * q1;
}

double WENOReconstruction::Weno3Right(const double fi, const double fip1, const double fip2,
                                      const double eps, const int p) {
    const double d[2] = {2.0 / 3.0, 1.0 / 3.0};

    double beta[2];
    beta[0] = (fip2 - fip1) * (fip2 - fip1);
    beta[1] = (fip1 - fi)   * (fip1 - fi);

    const double q0 = 0.5 * fip1 + 0.5 * fip2;
    const double q1 = 1.5 * fip1 - 0.5 * fi;

    double omega[2];
    ComputeNonlinearWeights(beta, d, 2, eps, p, omega);

    return omega[0] * q0 + omega[1] * q1;
}

// ----- WENO5 (JS) -----

double WENOReconstruction::Weno5Left(const double fim2, const double fim1, const double fi,
                                     const double fip1, const double fip2,
                                     const double eps, const int p) {
    const double d[3] = {0.1, 0.6, 0.3};

    double beta[3];
    {
        const double t1 = fim2 - 2.0 * fim1 + fi;
        const double t2 = fim2 - 4.0 * fim1 + 3.0 * fi;
        beta[0] = 13.0 / 12.0 * t1 * t1 + 0.25 * t2 * t2;
    }
    {
        const double t1 = fim1 - 2.0 * fi + fip1;
        const double t2 = fim1 - fip1;
        beta[1] = 13.0 / 12.0 * t1 * t1 + 0.25 * t2 * t2;
    }
    {
        const double t1 = fi - 2.0 * fip1 + fip2;
        const double t2 = 3.0 * fi - 4.0 * fip1 + fip2;
        beta[2] = 13.0 / 12.0 * t1 * t1 + 0.25 * t2 * t2;
    }

    const double q0 = (1.0 / 3.0) * fim2 + (-7.0 / 6.0) * fim1 + (11.0 / 6.0) * fi;
    const double q1 = (-1.0 / 6.0) * fim1 + (5.0 / 6.0) * fi + (1.0 / 3.0) * fip1;
    const double q2 = (1.0 / 3.0) * fi + (5.0 / 6.0) * fip1 + (-1.0 / 6.0) * fip2;

    double omega[3];
    ComputeNonlinearWeights(beta, d, 3, eps, p, omega);

    return omega[0] * q0 + omega[1] * q1 + omega[2] * q2;
}

double WENOReconstruction::Weno5Right(const double fim1, const double fi, const double fip1,
                                      const double fip2, const double fip3,
                                      const double eps, const int p) {
    const double d[3] = {0.3, 0.6, 0.1};

    double beta[3];
    {
        const double t1 = fip3 - 2.0 * fip2 + fip1;
        const double t2 = fip3 - 4.0 * fip2 + 3.0 * fip1;
        beta[0] = 13.0 / 12.0 * t1 * t1 + 0.25 * t2 * t2;
    }
    {
        const double t1 = fip2 - 2.0 * fip1 + fi;
        const double t2 = fip2 - fi;
        beta[1] = 13.0 / 12.0 * t1 * t1 + 0.25 * t2 * t2;
    }
    {
        const double t1 = fip1 - 2.0 * fi + fim1;
        const double t2 = 3.0 * fip1 - 4.0 * fi + fim1;
        beta[2] = 13.0 / 12.0 * t1 * t1 + 0.25 * t2 * t2;
    }

    const double q0 = (-1.0 / 6.0) * fip3 + (5.0 / 6.0) * fip2 + (1.0 / 3.0) * fip1;
    const double q1 = (1.0 / 3.0) * fip2 + (5.0 / 6.0) * fip1 + (-1.0 / 6.0) * fi;
    const double q2 = (11.0 / 6.0) * fip1 + (-7.0 / 6.0) * fi + (1.0 / 3.0) * fim1;

    double omega[3];
    ComputeNonlinearWeights(beta, d, 3, eps, p, omega);

    return omega[0] * q0 + omega[1] * q1 + omega[2] * q2;
}

// ----------------------- main API -----------------------

void WENOReconstruction::ReconstructFace(const xt::xtensor<double, 4>& W,
                                         const Axis axis,
                                         const int i, const int j, const int k,
                                         PrimitiveCell& WL,
                                         PrimitiveCell& WR) const {
    const AxisStride st = AxisStride::FromAxis(axis);

    for (std::size_t v = 0; v < var::nvar; ++v) {
        double wl = 0.0;
        double wr = 0.0;

        if (order_ == 3) {
            // WL at face i+1/2 from left cell i: uses offsets {-1,0,+1}
            const double fim1 = LoadAt(W, v, i, j, k, st, -1);
            const double fi   = LoadAt(W, v, i, j, k, st,  0);
            const double fip1 = LoadAt(W, v, i, j, k, st, +1);

            // WR at same face from right cell (i+1): needs {0,+1,+2} in left-cell indexing
            const double fip2 = LoadAt(W, v, i, j, k, st, +2);

            wl = Weno3Left (fim1, fi,   fip1, epsilon_, p_);
            wr = Weno3Right(fi,   fip1, fip2, epsilon_, p_);
        } else {
            // order_ == 5
            // WL uses {-2,-1,0,+1,+2}
            const double fim2 = LoadAt(W, v, i, j, k, st, -2);
            const double fim1 = LoadAt(W, v, i, j, k, st, -1);
            const double fi   = LoadAt(W, v, i, j, k, st,  0);
            const double fip1 = LoadAt(W, v, i, j, k, st, +1);
            const double fip2 = LoadAt(W, v, i, j, k, st, +2);

            // WR uses {-1,0,+1,+2,+3}
            const double fip3 = LoadAt(W, v, i, j, k, st, +3);

            wl = Weno5Left (fim2, fim1, fi,   fip1, fip2, epsilon_, p_);
            wr = Weno5Right(fim1, fi,   fip1, fip2, fip3, epsilon_, p_);
        }

        StorePrimitiveComponent(WL, v, wl);
        StorePrimitiveComponent(WR, v, wr);
    }
}