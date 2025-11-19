#include "riemann/ExactIdealGasRiemannSolver.hpp"

#include <algorithm>
#include <cmath>

namespace {
struct State {
    double rho;
    double u;
    double p;
    double a;
};

auto MakeState(const Primitive& w, double gamma) -> State {
    State s;
    s.rho = w.rho;
    s.u = w.u;
    s.p = w.P;
    s.a = std::sqrt(gamma * s.p / s.rho);
    return s;
}

// ================= Phi_k(p) and derivatives (Toro) =================

auto PhiRarefaction(double p, const State& s, double gamma) -> double {
    const double pr = p / s.p;
    const double exp = (gamma - 1.0) / (2.0 * gamma);
    return 2.0 * s.a / (gamma - 1.0) * (std::pow(pr, exp) - 1.0);
}

auto PhiRarefactionDerivative(double p, const State& s, double gamma) -> double {
    const double pr = p / s.p;
    const double exp = (gamma - 1.0) / (2.0 * gamma);
    const double coeff = (2.0 * s.a / (gamma - 1.0)) * (exp / s.p);
    return coeff * std::pow(pr, exp - 1.0);
}

auto PhiShock(double p, const State& s, double gamma) -> double {
    const double a = 2.0 / ((gamma + 1.0) * s.rho);
    const double b = (gamma - 1.0) / (gamma + 1.0) * s.p;
    return (p - s.p) * std::sqrt(a / (p + b));
}

auto PhiShockDerivative(double p, const State& s, double gamma) -> double {
    const double a = 2.0 / ((gamma + 1.0) * s.rho);
    const double b = (gamma - 1.0) / (gamma + 1.0) * s.p;
    const double sqrt_term = std::sqrt(a / (p + b));
    const double term = (p - s.p) / (2.0 * (p + b));
    return sqrt_term * (1.0 - term);
}

auto SidePhi(double p, const State& s, double gamma) -> double {
    return p <= s.p ? PhiRarefaction(p, s, gamma) : PhiShock(p, s, gamma);
}

auto SidePhiDerivative(double p, const State& s, double gamma) -> double {
    return p <= s.p ? PhiRarefactionDerivative(p, s, gamma)
                    : PhiShockDerivative(p, s, gamma);
}

// ================= Initial guess (PVRS / TRRS / TSRS) =================

auto InitialGuess(const State& l, const State& r, double Q_user, double gamma) -> double {
    const double p_pvrs =
        0.5 * (l.p + r.p) - 0.125 * (r.u - l.u) * (l.rho + r.rho) * (l.a + r.a);

    const double p_min = std::min(l.p, r.p);
    const double p_max = std::max(l.p, r.p);
    const double Q = p_max / p_min;

    double p_star = std::max(p_pvrs, 1e-16);

    // PVRS region
    if (p_min < p_pvrs && p_pvrs < p_max && Q < Q_user) {
        return p_star;
    }

    // Two-rarefaction
    if (p_pvrs <= p_min) {
        const double z = (gamma - 1.0) / (2.0 * gamma);
        const double num = l.a + r.a - 0.5 * (gamma - 1.0) * (r.u - l.u);
        const double den = l.a / std::pow(l.p, z) + r.a / std::pow(r.p, z);
        p_star = std::pow(num / den, 1.0 / z);
        return std::max(p_star, 1e-16);
    }

    // Two-shock
    const auto a = [gamma](const State& s) -> double {
        return 2.0 / ((gamma + 1.0) * s.rho);
    };
    const auto b = [gamma](const State& s) -> double {
        return (gamma - 1.0) / (gamma + 1.0) * s.p;
    };
    const auto g = [&](const State& s) -> double {
        return std::sqrt(a(s) / (p_pvrs + b(s)));
    };

    const double g_l = g(l);
    const double g_r = g(r);

    p_star = (g_l * l.p + g_r * r.p - (r.u - l.u)) / (g_l + g_r);

    return std::max(p_star, 1e-16);
}

// ================= Solve p* (with correct vacuum check) =================

auto SolveStarPressure(const State& l, const State& r, double gamma, double q_user)
    -> double {
    // Correct vacuum condition:
    // u_R - u_L >= 2 (a_L + a_R) / (gamma - 1)
    const double du = r.u - l.u;
    const double crit = 2.0 * (l.a + r.a) / (gamma - 1.0);
    if (du >= crit) {
        return 0.0;  // vacuum (double rarefaction)
    }

    double p = InitialGuess(l, r, q_user, gamma);
    const double p_min = 1e-16;

    for (int iter = 0; iter < 40; ++iter) {
        const double f_l = SidePhi(p, l, gamma);
        const double f_r = SidePhi(p, r, gamma);
        const double f = f_l + f_r + (r.u - l.u);

        if (std::fabs(f) < 1e-10) {
            break;
        }

        const double d_l = SidePhiDerivative(p, l, gamma);
        const double d_r = SidePhiDerivative(p, r, gamma);
        const double df = d_l + d_r;

        if (df == 0.0 || !std::isfinite(df)) {
            break;
        }

        double p_new = p - f / df;

        if (!std::isfinite(p_new) || p_new < p_min) {
            p_new = p_min;
        }

        if (std::fabs(p_new - p) <= 1e-8 * (p + p_min)) {
            p = p_new;
            break;
        }

        p = p_new;
    }

    if (!std::isfinite(p) || p <= 0.0) {
        return 0.0;  // treat as vacuum/failure
    }

    return p;
}

// ================= Sampling helpers =================

// --- Vacuum / double rarefaction sampling
auto SampleVacuum(double xi, const State& l, const State& r, double gamma) -> Primitive {
    const double a_l = l.a;
    const double a_r = r.a;

    const double shl = l.u - a_l;
    const double svl = l.u + 2.0 * a_l / (gamma - 1.0);

    const double shr = r.u + a_r;
    const double svr = r.u - 2.0 * a_r / (gamma - 1.0);

    if (xi <= shl) {
        return Primitive{l.rho, l.u, l.p};
    }
    if (xi >= shr) {
        return Primitive{r.rho, r.u, r.p};
    }

    // Left fan
    if (xi > shl && xi < svl) {
        const double u = 2.0 / (gamma + 1.0) * (a_l + 0.5 * (gamma - 1.0) * l.u + xi);
        const double a = 2.0 / (gamma + 1.0) * (a_l + 0.5 * (gamma - 1.0) * (l.u - xi));
        const double rho = l.rho * std::pow(a / a_l, 2.0 / (gamma - 1.0));
        const double p = l.p * std::pow(a / a_l, 2.0 * gamma / (gamma - 1.0));
        return Primitive{rho, u, p};
    }

    // Right fan
    if (xi > svr && xi < shr) {
        const double u = 2.0 / (gamma + 1.0) * (-a_r + 0.5 * (gamma - 1.0) * r.u + xi);
        const double a = 2.0 / (gamma + 1.0) *
                         (a_r - 0.5 * (gamma - 1.0) * r.u + 0.5 * (gamma - 1.0) * xi);
        const double rho = r.rho * std::pow(a / a_r, 2.0 / (gamma - 1.0));
        const double p = r.p * std::pow(a / a_r, 2.0 * gamma / (gamma - 1.0));
        return Primitive{rho, u, p};
    }

    // True vacuum between SVL and SVR
    return Primitive{0.0, 0.0, 0.0};
}

// --- Non-vacuum standard Riemann fan
auto SampleNonVacuum(double xi, double p_star, double u_star, const State& l,
                     const State& r, double gamma) -> Primitive {
    // LEFT of contact
    if (xi <= u_star) {
        if (p_star > l.p) {
            // Left shock
            const double q = p_star / l.p;
            const double sl = l.u - l.a * std::sqrt(0.5 * ((gamma + 1.0) / gamma * q +
                                                           (gamma - 1.0) / gamma));
            if (xi <= sl) {
                return Primitive{l.rho, l.u, l.p};
            }
            const double factor = (q + (gamma - 1.0) / (gamma + 1.0)) /
                                  ((gamma - 1.0) / (gamma + 1.0) * q + 1.0);
            const double rho_star_l = l.rho * factor;
            return Primitive{rho_star_l, u_star, p_star};
        }
        // Left rarefaction
        const double a_l = l.a;
        const double shl = l.u - a_l;
        const double rho_star_l = l.rho * std::pow(p_star / l.p, 1.0 / gamma);
        const double a_star_l = std::sqrt(gamma * p_star / rho_star_l);
        const double stl = u_star - a_star_l;

        if (xi <= shl) {
            return Primitive{l.rho, l.u, l.p};
        }
        if (xi >= stl) {
            return Primitive{rho_star_l, u_star, p_star};
        }

        const double u = 2.0 / (gamma + 1.0) * (a_l + 0.5 * (gamma - 1.0) * l.u + xi);
        const double a = 2.0 / (gamma + 1.0) * (a_l + 0.5 * (gamma - 1.0) * (l.u - xi));
        const double rho = l.rho * std::pow(a / a_l, 2.0 / (gamma - 1.0));
        const double p = l.p * std::pow(a / a_l, 2.0 * gamma / (gamma - 1.0));
        return Primitive{rho, u, p};
    }

    // RIGHT of contact (xi > uStar)
    if (p_star > r.p) {
        // Right shock
        const double q = p_star / r.p;
        const double sr =
            r.u +
            r.a * std::sqrt(0.5 * ((gamma + 1.0) / gamma * q + (gamma - 1.0) / gamma));
        if (xi >= sr) {
            return Primitive{r.rho, r.u, r.p};
        }
        const double factor = (q + (gamma - 1.0) / (gamma + 1.0)) /
                              ((gamma - 1.0) / (gamma + 1.0) * q + 1.0);
        const double rho_star_r = r.rho * factor;
        return Primitive{rho_star_r, u_star, p_star};
    }

    // Right rarefaction (pStar <= p_R)
    const double a_r = r.a;
    const double shr = r.u + a_r;

    const double rho_star_r = r.rho * std::pow(p_star / r.p, 1.0 / gamma);
    const double a_star_r = std::sqrt(gamma * p_star / rho_star_r);
    const double str = u_star + a_star_r;

    if (xi >= shr) {
        // Right initial state region
        return Primitive{r.rho, r.u, r.p};
    }
    if (xi <= str) {
        // Right star state
        return Primitive{rho_star_r, u_star, p_star};
    }

    // Self-similar right rarefaction fan between STR and SHR
    const double u = 2.0 / (gamma + 1.0) * (-a_r + 0.5 * (gamma - 1.0) * r.u + xi);
    const double a = 2.0 / (gamma + 1.0) *
                     (a_r - 0.5 * (gamma - 1.0) * r.u + 0.5 * (gamma - 1.0) * xi);
    const double rho = r.rho * std::pow(a / a_r, 2.0 / (gamma - 1.0));
    const double p = r.p * std::pow(a / a_r, 2.0 * gamma / (gamma - 1.0));
    return Primitive{rho, u, p};
}
}  // namespace

// ================= Public methods =================

ExactIdealGasRiemannSolver::ExactIdealGasRiemannSolver() : xi_(0.0) {}

void ExactIdealGasRiemannSolver::SetXi(double xi) { xi_ = xi; }

void ExactIdealGasRiemannSolver::SetQ(double Q) { Q_user_ = Q; }

auto ExactIdealGasRiemannSolver::Sample(const Primitive& left, const Primitive& right,
                                        double gamma, double xi) const -> Primitive {
    const State l = MakeState(left, gamma);
    const State r = MakeState(right, gamma);

    const double p_star = SolveStarPressure(l, r, gamma, Q_user_);

    // Vacuum / double rarefaction handled separately
    if (p_star <= 0.0) {
        return SampleVacuum(xi, l, r, gamma);
    }

    const double f_l = SidePhi(p_star, l, gamma);
    const double f_r = SidePhi(p_star, r, gamma);
    const double u_star = 0.5 * (l.u + r.u + f_r - f_l);

    return SampleNonVacuum(xi, p_star, u_star, l, r, gamma);
}

auto ExactIdealGasRiemannSolver::ComputeFlux(const Primitive& left,
                                             const Primitive& right, double gamma) const
    -> Flux {
    const Primitive sample = Sample(left, right, gamma, xi_);

    return EulerFlux(sample, gamma);
}
