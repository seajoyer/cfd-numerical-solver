#include "riemann/OsherRiemannSolver.hpp"

#include <cmath>
#include <algorithm>
#include <iostream>
#include <ostream>

#include "solver/EOS.hpp"

namespace {
auto Solve3x3(const Conservative& r1,
              const Conservative& r2,
              const Conservative& r3,
              const Conservative& b,
              double& x1, double& x2, double& x3) -> bool {
    double a11 = r1.rho;
    double a12 = r2.rho;
    double a13 = r3.rho;
    double a21 = r1.rhoU;
    double a22 = r2.rhoU;
    double a23 = r3.rhoU;
    double a31 = r1.E;
    double a32 = r2.E;
    double a33 = r3.E;

    double b1 = b.rho;
    double b2 = b.rhoU;
    double b3 = b.E;

    const double eps = 1e-14;

    int pivot_row = 0;
    double max_abs = std::fabs(a11);
    if (std::fabs(a21) > max_abs) {
        max_abs = std::fabs(a21);
        pivot_row = 1;
    }
    if (std::fabs(a31) > max_abs) {
        max_abs = std::fabs(a31);
        pivot_row = 2;
    }
    if (max_abs < eps) {
        return false;
    }

    auto swap_rows = [&](int r) {
        if (r == 0) {
            return;
        }
        if (r == 1) {
            std::swap(a11, a21);
            std::swap(a12, a22);
            std::swap(a13, a23);
            std::swap(b1, b2);
        } else {
            std::swap(a11, a31);
            std::swap(a12, a32);
            std::swap(a13, a33);
            std::swap(b1, b3);
        }
    };
    swap_rows(pivot_row);

    double m21 = a21 / a11;
    double m31 = a31 / a11;

    a21 -= m21 * a11;
    a22 -= m21 * a12;
    a23 -= m21 * a13;
    b2 -= m21 * b1;

    a31 -= m31 * a11;
    a32 -= m31 * a12;
    a33 -= m31 * a13;
    b3 -= m31 * b1;

    pivot_row = 1;
    max_abs = std::fabs(a22);
    if (std::fabs(a32) > max_abs) {
        max_abs = std::fabs(a32);
        pivot_row = 2;
    }
    if (max_abs < eps) {
        return false;
    }
    if (pivot_row == 2) {
        std::swap(a22, a32);
        std::swap(a23, a33);
        std::swap(b2, b3);
    }

    double m32 = a32 / a22;
    a32 -= m32 * a22;
    a33 -= m32 * a23;
    b3 -= m32 * b2;

    if (std::fabs(a33) < eps) {
        return false;
    }

    x3 = b3 / a33;
    x2 = (b2 - a23 * x3) / a22;
    x1 = (b1 - a12 * x2 - a13 * x3) / a11;

    return std::isfinite(x1) && std::isfinite(x2) && std::isfinite(x3);
}
} // namespace

auto OsherRiemannSolver::ApplyAbsJacobian(const Primitive& w,
                                          const double gamma,
                                          const Conservative& dU,
                                          Conservative& result) -> bool {
    const double rho = w.rho;
    const double u = w.u;
    const double p = w.P;

    if (rho <= 0.0 || p <= 0.0) {
        return false;
    }

    const double a = std::sqrt(gamma * p / rho);
    const double H = (p / (gamma - 1.0) + 0.5 * rho * u * u + p) / rho;

    const double lambda1 = u - a;
    const double lambda2 = u;
    const double lambda3 = u + a;

    Conservative r1;
    r1.rho = 1.0;
    r1.rhoU = u - a;
    r1.E = H - u * a;

    Conservative r2;
    r2.rho = 1.0;
    r2.rhoU = u;
    r2.E = 0.5 * u * u;

    Conservative r3;
    r3.rho = 1.0;
    r3.rhoU = u + a;
    r3.E = H + u * a;

    double alpha1 = 0.0;
    double alpha2 = 0.0;
    double alpha3 = 0.0;
    if (!Solve3x3(r1, r2, r3, dU, alpha1, alpha2, alpha3)) {
        return false;
    }

    const double s1 = std::fabs(lambda1) * alpha1;
    const double s2 = std::fabs(lambda2) * alpha2;
    const double s3 = std::fabs(lambda3) * alpha3;

    result.rho = s1 * r1.rho + s2 * r2.rho + s3 * r3.rho;
    result.rhoU = s1 * r1.rhoU + s2 * r2.rhoU + s3 * r3.rhoU;
    result.E = s1 * r1.E + s2 * r2.E + s3 * r3.E;

    return std::isfinite(result.rho) &&
           std::isfinite(result.rhoU) &&
           std::isfinite(result.E);
}

auto OsherRiemannSolver::RusanovFlux(const Primitive& left,
                                     const Primitive& right,
                                     const double gamma) -> Flux {
    const Flux f_l = EulerFlux(left, gamma);
    const Flux f_r = EulerFlux(right, gamma);

    const Conservative u_l = EOS::PrimToCons(left, gamma);
    const Conservative u_r = EOS::PrimToCons(right, gamma);

    const double a_l = std::sqrt(gamma * left.P / left.rho);
    const double a_r = std::sqrt(gamma * right.P / right.rho);

    const double s_l = std::fabs(left.u) + a_l;
    const double s_r = std::fabs(right.u) + a_r;
    const double a_max = std::max(s_l, s_r);

    Conservative dU;
    dU.rho = u_r.rho - u_l.rho;
    dU.rhoU = u_r.rhoU - u_l.rhoU;
    dU.E = u_r.E - u_l.E;

    Flux flux;
    flux.mass = 0.5 * (f_l.mass + f_r.mass) - 0.5 * a_max * dU.rho;
    flux.momentum = 0.5 * (f_l.momentum + f_r.momentum) - 0.5 * a_max * dU.rhoU;
    flux.energy = 0.5 * (f_l.energy + f_r.energy) - 0.5 * a_max * dU.E;

    return flux;
}

auto OsherRiemannSolver::ComputeFlux(const Primitive& left,
                                     const Primitive& right,
                                     const double gamma) const -> Flux {
    const Conservative U_L = EOS::PrimToCons(left, gamma);
    const Conservative U_R = EOS::PrimToCons(right, gamma);

    const Flux F_L = EulerFlux(left, gamma);
    const Flux F_R = EulerFlux(right, gamma);

    Conservative dU = U_R - U_L;

    const double norm_dU = std::fabs(dU.rho) + std::fabs(dU.rhoU) + std::fabs(dU.E);
    if (norm_dU < 1e-14) {
        return F_L;
    }

    const double s[3] = {0.5 * (1.0 - std::sqrt(3.0 / 5.0)),
                         0.5,
                         0.5 * (1.0 + std::sqrt(3.0 / 5.0))};
    const double w[3] = {5.0 / 18.0,
                         4.0 / 9.0,
                         5.0 / 18.0};

    Conservative integral;

    for (int m = 0; m < 3; ++m) {
        const double theta = s[m];

        Primitive w_theta = (1.0 - theta) * left + theta * right;

        Conservative absA_dU;
        if (!ApplyAbsJacobian(w_theta, gamma, dU, absA_dU)) {
            return RusanovFlux(left, right, gamma);
        }
        integral += w[m] * absA_dU;
    }

    Flux flux = 0.5 * (F_L + F_R) - 0.5 * integral;

    return flux;
}