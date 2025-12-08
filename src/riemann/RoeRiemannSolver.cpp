#include "riemann/RoeRiemannSolver.hpp"

#include <algorithm>
#include <cmath>

#include "solver/EOS.hpp"

namespace {
auto RusanovFallback(const Primitive& left,
                     const Primitive& right,
                     double gamma) -> Flux {
    const double rho_l = left.rho;
    const double u_l = left.u;
    const double p_l = left.P;

    const double rho_r = right.rho;
    const double u_r = right.u;
    const double p_r = right.P;

    const Conservative ul = EOS::PrimToCons(left, gamma);
    const Conservative ur = EOS::PrimToCons(right, gamma);

    const Flux fl = EulerFlux(left, gamma);
    const Flux fr = EulerFlux(right, gamma);

    const double a_l = std::sqrt(gamma * p_l / rho_l);
    const double a_r = std::sqrt(gamma * p_r / rho_r);

    const double s_l = std::abs(u_l) + a_l;
    const double s_r = std::abs(u_r) + a_r;
    const double a_max = std::max(s_l, s_r);

    return 0.5 * (fl + fr + a_max * (ul - ur));
}

auto EntropyFix(double lambda, double lambda_L, double lambda_R) -> double {
    const double delta = std::max(0.0, lambda_R - lambda_L);
    if (delta <= 0.0) {
        return std::abs(lambda);
    }
    const double abs_lambda = std::abs(lambda);
    if (abs_lambda >= delta) {
        return abs_lambda;
    }
    return 0.5 * (lambda * lambda / delta + delta);
}
} // namespace

auto RoeRiemannSolver::ComputeFlux(const Primitive& left,
                                   const Primitive& right,
                                   double gamma) const -> Flux {
    const double rho_l = left.rho;
    const double u_l = left.u;
    const double p_l = left.P;

    const double rho_r = right.rho;
    const double u_r = right.u;
    const double p_r = right.P;

    const Conservative ul = EOS::PrimToCons(left, gamma);
    const Conservative ur = EOS::PrimToCons(right, gamma);

    const Flux fl = EulerFlux(left, gamma);
    const Flux fr = EulerFlux(right, gamma);

    if (rho_l <= 0.0 || rho_r <= 0.0) {
        return RusanovFallback(left, right, gamma);
    }

    const double sqrt_rho_l = std::sqrt(rho_l);
    const double sqrt_rho_r = std::sqrt(rho_r);
    const double denom = sqrt_rho_l + sqrt_rho_r;

    if (denom <= 0.0) {
        return RusanovFallback(left, right, gamma);
    }

    const double el = ul.E;
    const double er = ur.E;

    const double H_l = (el + p_l) / rho_l;
    const double H_r = (er + p_r) / rho_r;

    const double u_tilde = (sqrt_rho_l * u_l + sqrt_rho_r * u_r) / denom;
    const double H_tilde = (sqrt_rho_l * H_l + sqrt_rho_r * H_r) / denom;
    const double rho_tilde = sqrt_rho_l * sqrt_rho_r;

    const double a2_tilde =
        (gamma - 1.0) * (H_tilde - 0.5 * u_tilde * u_tilde);

    if (a2_tilde <= 0.0) {
        return RusanovFallback(left, right, gamma);
    }

    const double a_tilde = std::sqrt(a2_tilde);

    double lambda1 = u_tilde - a_tilde;
    double lambda2 = u_tilde;
    double lambda3 = u_tilde + a_tilde;

    const double a_l = std::sqrt(gamma * p_l / rho_l);
    const double a_r = std::sqrt(gamma * p_r / rho_r);

    const double lambda1_L = u_l - a_l;
    const double lambda1_R = u_r - a_r;
    const double lambda3_L = u_l + a_l;
    const double lambda3_R = u_r + a_r;

    const double a1 = EntropyFix(lambda1, lambda1_L, lambda1_R);
    const double a2 = std::abs(lambda2);
    const double a3 = EntropyFix(lambda3, lambda3_L, lambda3_R);

    const double drho = rho_r - rho_l;
    const double du = u_r - u_l;
    const double dp = p_r - p_l;

    const double alpha2 = drho - dp / a2_tilde;
    const double alpha1 = 0.5 / a2_tilde * (dp - rho_tilde * a_tilde * du);
    const double alpha3 = 0.5 / a2_tilde * (dp + rho_tilde * a_tilde * du);

    const double r1_0 = 1.0;
    const double r1_1 = u_tilde - a_tilde;
    const double r1_2 = H_tilde - u_tilde * a_tilde;

    const double r2_0 = 1.0;
    const double r2_1 = u_tilde;
    const double r2_2 = 0.5 * u_tilde * u_tilde;

    const double r3_0 = 1.0;
    const double r3_1 = u_tilde + a_tilde;
    const double r3_2 = H_tilde + u_tilde * a_tilde;

    const double dF0 = a1 * alpha1 * r1_0 +
                       a2 * alpha2 * r2_0 +
                       a3 * alpha3 * r3_0;

    const double dF1 = a1 * alpha1 * r1_1 +
                       a2 * alpha2 * r2_1 +
                       a3 * alpha3 * r3_1;

    const double dF2 = a1 * alpha1 * r1_2 +
                       a2 * alpha2 * r2_2 +
                       a3 * alpha3 * r3_2;

    Flux f;
    f.mass = 0.5 * (fl.mass + fr.mass) - 0.5 * dF0;
    f.momentum = 0.5 * (fl.momentum + fr.momentum) - 0.5 * dF1;
    f.energy = 0.5 * (fl.energy + fr.energy) - 0.5 * dF2;

    return f;
}