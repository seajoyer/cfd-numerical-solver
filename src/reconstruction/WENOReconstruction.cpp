#include "reconstruction/WENOReconstruction.hpp"

#include "data/DataLayer.hpp"
#include "data/Variables.hpp"

#include <cmath>
#include <stdexcept>

namespace {
// Compute nonlinear WENO-JS weights from smoothness indicators and linear weights.
void ComputeNonlinearWeights(const double* beta,
                             const double* d,
                             int r,
                             double epsilon,
                             int p,
                             double* omega) {
    double alpha_sum = 0.0;

    for (int k = 0; k < r; ++k) {
        const double denom = std::pow(epsilon + beta[k], p);
        const double alpha = d[k] / denom;
        omega[k] = alpha;
        alpha_sum += alpha;
    }

    if (alpha_sum <= 0.0) {
        // Fallback to linear weights if something goes wrong.
        for (int k = 0; k < r; ++k) {
            omega[k] = d[k];
        }
    } else {
        const double inv = 1.0 / alpha_sum;
        for (int k = 0; k < r; ++k) {
            omega[k] *= inv;
        }
    }
}

// ----- WENO3: left-biased scalar reconstruction at interface i+1/2 -----

double Weno3LeftScalar(const xt::xarray<double>& f,
                       int i,
                       int total_size,
                       double epsilon,
                       int p) {
    // Need cells: i-1, i, i+1
    const int n_interfaces = total_size > 0 ? total_size - 1 : 0;
    const int half = 1; // (3-1)/2

    if (i < half || i > n_interfaces - 1 - half) {
        // Near boundaries: fall back to piecewise-constant.
        return f(i);
    }

    const double fim1 = f(i - 1);
    const double fi = f(i);
    const double fip1 = f(i + 1);

    // Linear weights for WENO3 (JS)
    const int r = 2;
    const double d[r] = {1.0 / 3.0, 2.0 / 3.0};

    // Smoothness indicators
    double beta[r];
    beta[0] = (fi - fim1) * (fi - fim1);
    beta[1] = (fip1 - fi) * (fip1 - fi);

    // Candidate polynomials at i+1/2
    //
    // q0 = -0.5*f_{i-1} + 1.5*f_{i}
    // q1 =  0.5*f_{i}   + 0.5*f_{i+1}
    static const double C[2][3] = {
        {-0.5, 1.5, 0.0}, // q0: f_{i-1}, f_i, f_{i+1}
        {0.0, 0.5, 0.5} // q1: f_{i-1}, f_i, f_{i+1}
    };

    const double q0 = C[0][0] * fim1 + C[0][1] * fi + C[0][2] * fip1;
    const double q1 = C[1][0] * fim1 + C[1][1] * fi + C[1][2] * fip1;

    double omega[r];
    ComputeNonlinearWeights(beta, d, r, epsilon, p, omega);

    return omega[0] * q0 + omega[1] * q1;
}

// ----- WENO5: left-biased scalar reconstruction at interface i+1/2 -----

double Weno5LeftScalar(const xt::xarray<double>& f,
                       int i,
                       int total_size,
                       double epsilon,
                       int p) {
    // Need cells: i-2, i-1, i, i+1, i+2
    const int n_interfaces = total_size > 0 ? total_size - 1 : 0;
    const int half = 2;

    if (i < half || i > n_interfaces - 1 - half) {
        // Near boundaries: fall back to piecewise-constant.
        return f(i);
    }

    const double fim2 = f(i - 2);
    const double fim1 = f(i - 1);
    const double fi = f(i);
    const double fip1 = f(i + 1);
    const double fip2 = f(i + 2);

    // Linear weights for WENO5 (JS)
    const int r = 3;
    const double d[r] = {0.1, 0.6, 0.3}; // {1/10, 6/10, 3/10}

    // Smoothness indicators (Jiang–Shu)
    double beta[r];

    // β0
    {
        const double t1 = fim2 - 2.0 * fim1 + fi;
        const double t2 = fim2 - 4.0 * fim1 + 3.0 * fi;
        beta[0] = 13.0 / 12.0 * t1 * t1 + 0.25 * t2 * t2;
    }

    // β1
    {
        const double t1 = fim1 - 2.0 * fi + fip1;
        const double t2 = fim1 - fip1;
        beta[1] = 13.0 / 12.0 * t1 * t1 + 0.25 * t2 * t2;
    }

    // β2
    {
        const double t1 = fi - 2.0 * fip1 + fip2;
        const double t2 = 3.0 * fi - 4.0 * fip1 + fip2;
        beta[2] = 13.0 / 12.0 * t1 * t1 + 0.25 * t2 * t2;
    }

    // Candidate polynomials at i+1/2 written as
    // qk = sum_{m=-2}^{2} C[k][m+2] * f_{i+m}
    static const double C[3][5] = {
        // q0 uses stencil {i-2, i-1, i}
        {1.0 / 3.0, -7.0 / 6.0, 11.0 / 6.0, 0.0, 0.0},
        // q1 uses stencil {i-1, i, i+1}
        {0.0, -1.0 / 6.0, 5.0 / 6.0, 1.0 / 3.0, 0.0},
        // q2 uses stencil {i, i+1, i+2}
        {0.0, 0.0, 1.0 / 3.0, 5.0 / 6.0, -1.0 / 6.0}
    };

    const double q0 = C[0][0] * fim2 + C[0][1] * fim1 + C[0][2] * fi +
                      C[0][3] * fip1 + C[0][4] * fip2;

    const double q1 = C[1][0] * fim2 + C[1][1] * fim1 + C[1][2] * fi +
                      C[1][3] * fip1 + C[1][4] * fip2;

    const double q2 = C[2][0] * fim2 + C[2][1] * fim1 + C[2][2] * fi +
                      C[2][3] * fip1 + C[2][4] * fip2;

    double omega[r];
    ComputeNonlinearWeights(beta, d, r, epsilon, p, omega);

    return omega[0] * q0 + omega[1] * q1 + omega[2] * q2;
}

// ----- Generic left-biased dispatcher (scalar) -----

double WenoLeftScalar(const xt::xarray<double>& f,
                      int i,
                      int total_size,
                      int order,
                      double epsilon,
                      int p) {
    switch (order) {
        case 3:
            return Weno3LeftScalar(f, i, total_size, epsilon, p);
        case 5:
            return Weno5LeftScalar(f, i, total_size, epsilon, p);
        default:
            // Should not happen if user guarantees {3,5}.
            throw std::runtime_error("WENOReconstruction: unsupported order.");
    }
}

// Reconstruct left states on a given triple of scalar fields.
void ReconstructLeftStatesOnScalars(const xt::xarray<double>& rho,
                                    const xt::xarray<double>& u,
                                    const xt::xarray<double>& P,
                                    int order,
                                    double epsilon,
                                    int p,
                                    xt::xarray<Primitive>& left_states) {
    const int total_size = static_cast<int>(rho.shape()[0]);
    const int n_interfaces = total_size > 0 ? total_size - 1 : 0;

    for (int i = 0; i < n_interfaces; ++i) {
        Primitive w;
        w.rho = WenoLeftScalar(rho, i, total_size, order, epsilon, p);
        w.u = WenoLeftScalar(u, i, total_size, order, epsilon, p);
        w.P = WenoLeftScalar(P, i, total_size, order, epsilon, p);
        left_states(i) = w;
    }
}
} // namespace

// ----- WENOReconstruction methods -----

void WENOReconstruction::SetOrder(const int order) {
    order_ = order;
}

void WENOReconstruction::SetEpsilon(const double eps) {
    epsilon_ = eps;
}

void WENOReconstruction::SetNonlinearWeightPower(const int p) {
    p_ = p;
}

void WENOReconstruction::ReconstructStates(const DataLayer& layer,
                                           xt::xarray<Primitive>& left_states,
                                           xt::xarray<Primitive>& right_states) const {
    const int total_size = layer.GetTotalSize();
    const int n_interfaces = total_size > 0 ? total_size - 1 : 0;
    const std::size_t n_int = static_cast<std::size_t>(n_interfaces);

    left_states = xt::xarray<Primitive>::from_shape({n_int});
    right_states = xt::xarray<Primitive>::from_shape({n_int});

    if (n_interfaces <= 0) {
        return;
    }

    // Left states: WENO on original arrays.
    ReconstructLeftStatesOnScalars(layer.rho, layer.u, layer.P, order_, epsilon_, p_,
                                   left_states);

    // Right states: apply the same left-biased reconstruction on a reversed grid
    // and then mirror the interface indices back.
    //
    // Reverse primitive fields: f_rev(k) = f(N-1-k).
    xt::xarray<double> rho_rev = xt::xarray<double>::from_shape(
        {static_cast<std::size_t>(total_size)});
    xt::xarray<double> u_rev = xt::xarray<double>::from_shape(
        {static_cast<std::size_t>(total_size)});
    xt::xarray<double> P_rev = xt::xarray<double>::from_shape(
        {static_cast<std::size_t>(total_size)});

    for (int k = 0; k < total_size; ++k) {
        const int kk = total_size - 1 - k;
        rho_rev(k) = layer.rho(kk);
        u_rev(k) = layer.u(kk);
        P_rev(k) = layer.P(kk);
    }

    // Temporary left states on reversed grid.
    xt::xarray<Primitive> left_rev =
        xt::xarray<Primitive>::from_shape({n_int});

    ReconstructLeftStatesOnScalars(rho_rev, u_rev, P_rev, order_, epsilon_, p_, left_rev);

    // Map them back: interface i in original grid corresponds to
    // interface (n_interfaces - 1 - i) in the reversed grid.
    for (int i = 0; i < n_interfaces; ++i) {
        const int j = n_interfaces - 1 - i;
        right_states(i) = left_rev(j);
    }
}