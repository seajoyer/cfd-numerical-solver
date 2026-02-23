#include "riemann/ExactIdealGasRiemannSolver.hpp"

#include <algorithm>
#include <cmath>

#include "riemann/RiemannHelpers.hpp"
#include "data/Variables.hpp"  // SoundSpeed, EulerFlux

namespace {
    // 1D state along interface normal
    struct State1D final {
        double rho = 0.0;
        double un = 0.0;
        double p = 0.0;
        double a = 0.0;
    };

    struct Primitive1D final {
        double rho = 0.0;
        double un = 0.0;
        double p = 0.0;
    };

    State1D MakeState1D(const PrimitiveCell& w, const Axis axis, const double gamma) {
        double un, ut1, ut2;
        riemann::SplitVelocity(w, axis, un, ut1, ut2);

        State1D s;
        s.rho = w.rho;
        s.un = un;
        s.p = w.P;
        s.a = std::sqrt(std::max(gamma * s.p / s.rho, 0.0));
        return s;
    }

    // ================= Phi_k(p) and derivatives (Toro) =================

    double PhiRarefaction(const double p, const State1D& s, const double gamma) {
        const double pr = p / s.p;
        const double exp = (gamma - 1.0) / (2.0 * gamma);
        return 2.0 * s.a / (gamma - 1.0) * (std::pow(pr, exp) - 1.0);
    }

    double PhiRarefactionDerivative(const double p, const State1D& s, const double gamma) {
        const double pr = p / s.p;
        const double exp = (gamma - 1.0) / (2.0 * gamma);
        const double coeff = (2.0 * s.a / (gamma - 1.0)) * (exp / s.p);
        return coeff * std::pow(pr, exp - 1.0);
    }

    double PhiShock(const double p, const State1D& s, const double gamma) {
        const double A = 2.0 / ((gamma + 1.0) * s.rho);
        const double B = (gamma - 1.0) / (gamma + 1.0) * s.p;
        return (p - s.p) * std::sqrt(A / (p + B));
    }

    double PhiShockDerivative(const double p, const State1D& s, const double gamma) {
        const double A = 2.0 / ((gamma + 1.0) * s.rho);
        const double B = (gamma - 1.0) / (gamma + 1.0) * s.p;
        const double sqrt_term = std::sqrt(A / (p + B));
        const double term = (p - s.p) / (2.0 * (p + B));
        return sqrt_term * (1.0 - term);
    }

    double SidePhi(const double p, const State1D& s, const double gamma) {
        return (p <= s.p) ? PhiRarefaction(p, s, gamma) : PhiShock(p, s, gamma);
    }

    double SidePhiDerivative(const double p, const State1D& s, const double gamma) {
        return (p <= s.p)
                   ? PhiRarefactionDerivative(p, s, gamma)
                   : PhiShockDerivative(p, s, gamma);
    }

    // ================= Initial guess (PVRS / TRRS / TSRS) =================

    double InitialGuess(const State1D& l, const State1D& r, const double Q_user, const double gamma) {
        const double p_pvrs =
            0.5 * (l.p + r.p) - 0.125 * (r.un - l.un) * (l.rho + r.rho) * (l.a + r.a);

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
            const double num = l.a + r.a - 0.5 * (gamma - 1.0) * (r.un - l.un);
            const double den = l.a / std::pow(l.p, z) + r.a / std::pow(r.p, z);
            p_star = std::pow(num / den, 1.0 / z);
            return std::max(p_star, 1e-16);
        }

        // Two-shock
        const auto A = [gamma](const State1D& s) -> double {
            return 2.0 / ((gamma + 1.0) * s.rho);
        };
        const auto B = [gamma](const State1D& s) -> double {
            return (gamma - 1.0) / (gamma + 1.0) * s.p;
        };
        const auto g = [&](const State1D& s) -> double {
            return std::sqrt(A(s) / (p_pvrs + B(s)));
        };

        const double g_l = g(l);
        const double g_r = g(r);

        p_star = (g_l * l.p + g_r * r.p - (r.un - l.un)) / (g_l + g_r);
        return std::max(p_star, 1e-16);
    }

    // ================= Solve p* (with vacuum check) =================

    double SolveStarPressure(const State1D& l, const State1D& r, const double gamma, const double q_user) {
        // Correct vacuum condition: u_R - u_L >= 2 (a_L + a_R) / (gamma - 1)
        const double du = r.un - l.un;
        const double crit = 2.0 * (l.a + r.a) / (gamma - 1.0);
        if (du >= crit) {
            return 0.0; // vacuum (double rarefaction)
        }

        double p = InitialGuess(l, r, q_user, gamma);
        const double p_min = 1e-16;

        for (int iter = 0; iter < 40; ++iter) {
            const double f_l = SidePhi(p, l, gamma);
            const double f_r = SidePhi(p, r, gamma);
            const double f = f_l + f_r + (r.un - l.un);

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
            return 0.0; // treat as vacuum/failure
        }

        return p;
    }

    // ================= Sampling helpers =================

    // Vacuum / double rarefaction sampling
    Primitive1D SampleVacuum(const double xi, const State1D& l, const State1D& r, const double gamma) {
        const double a_l = l.a;
        const double a_r = r.a;

        const double shl = l.un - a_l;
        const double svl = l.un + 2.0 * a_l / (gamma - 1.0);

        const double shr = r.un + a_r;
        const double svr = r.un - 2.0 * a_r / (gamma - 1.0);

        if (xi <= shl) return Primitive1D{l.rho, l.un, l.p};
        if (xi >= shr) return Primitive1D{r.rho, r.un, r.p};

        // Left fan
        if (xi > shl && xi < svl) {
            const double un = 2.0 / (gamma + 1.0) * (a_l + 0.5 * (gamma - 1.0) * l.un + xi);
            const double a = 2.0 / (gamma + 1.0) * (a_l + 0.5 * (gamma - 1.0) * (l.un - xi));
            const double rho = l.rho * std::pow(a / a_l, 2.0 / (gamma - 1.0));
            const double p = l.p * std::pow(a / a_l, 2.0 * gamma / (gamma - 1.0));
            return Primitive1D{rho, un, p};
        }

        // Right fan
        if (xi > svr && xi < shr) {
            const double un = 2.0 / (gamma + 1.0) * (-a_r + 0.5 * (gamma - 1.0) * r.un + xi);
            const double a = 2.0 / (gamma + 1.0) *
                (a_r - 0.5 * (gamma - 1.0) * r.un + 0.5 * (gamma - 1.0) * xi);
            const double rho = r.rho * std::pow(a / a_r, 2.0 / (gamma - 1.0));
            const double p = r.p * std::pow(a / a_r, 2.0 * gamma / (gamma - 1.0));
            return Primitive1D{rho, un, p};
        }

        // True vacuum between SVL and SVR
        return Primitive1D{0.0, 0.0, 0.0};
    }

    // Non-vacuum standard sampling
    Primitive1D SampleNonVacuum(const double xi,
                                const double p_star,
                                const double u_star,
                                const State1D& l,
                                const State1D& r,
                                const double gamma) {
        // LEFT of contact
        if (xi <= u_star) {
            if (p_star > l.p) {
                // Left shock
                const double q = p_star / l.p;
                const double sl = l.un - l.a * std::sqrt(0.5 * ((gamma + 1.0) / gamma * q +
                    (gamma - 1.0) / gamma));
                if (xi <= sl) return Primitive1D{l.rho, l.un, l.p};

                const double factor = (q + (gamma - 1.0) / (gamma + 1.0)) /
                    ((gamma - 1.0) / (gamma + 1.0) * q + 1.0);
                const double rho_star_l = l.rho * factor;
                return Primitive1D{rho_star_l, u_star, p_star};
            }

            // Left rarefaction
            const double a_l = l.a;
            const double shl = l.un - a_l;
            const double rho_star_l = l.rho * std::pow(p_star / l.p, 1.0 / gamma);
            const double a_star_l = std::sqrt(std::max(gamma * p_star / rho_star_l, 0.0));
            const double stl = u_star - a_star_l;

            if (xi <= shl) return Primitive1D{l.rho, l.un, l.p};
            if (xi >= stl) return Primitive1D{rho_star_l, u_star, p_star};

            const double un = 2.0 / (gamma + 1.0) * (a_l + 0.5 * (gamma - 1.0) * l.un + xi);
            const double a = 2.0 / (gamma + 1.0) * (a_l + 0.5 * (gamma - 1.0) * (l.un - xi));
            const double rho = l.rho * std::pow(a / a_l, 2.0 / (gamma - 1.0));
            const double p = l.p * std::pow(a / a_l, 2.0 * gamma / (gamma - 1.0));
            return Primitive1D{rho, un, p};
        }

        // RIGHT of contact (xi > uStar)
        if (p_star > r.p) {
            // Right shock
            const double q = p_star / r.p;
            const double sr =
                r.un + r.a * std::sqrt(0.5 * ((gamma + 1.0) / gamma * q + (gamma - 1.0) / gamma));
            if (xi >= sr) return Primitive1D{r.rho, r.un, r.p};

            const double factor = (q + (gamma - 1.0) / (gamma + 1.0)) /
                ((gamma - 1.0) / (gamma + 1.0) * q + 1.0);
            const double rho_star_r = r.rho * factor;
            return Primitive1D{rho_star_r, u_star, p_star};
        }

        // Right rarefaction (pStar <= p_R)
        const double a_r = r.a;
        const double shr = r.un + a_r;

        const double rho_star_r = r.rho * std::pow(p_star / r.p, 1.0 / gamma);
        const double a_star_r = std::sqrt(std::max(gamma * p_star / rho_star_r, 0.0));
        const double str = u_star + a_star_r;

        if (xi >= shr) return Primitive1D{r.rho, r.un, r.p};
        if (xi <= str) return Primitive1D{rho_star_r, u_star, p_star};

        const double un = 2.0 / (gamma + 1.0) * (-a_r + 0.5 * (gamma - 1.0) * r.un + xi);
        const double a = 2.0 / (gamma + 1.0) *
            (a_r - 0.5 * (gamma - 1.0) * r.un + 0.5 * (gamma - 1.0) * xi);
        const double rho = r.rho * std::pow(a / a_r, 2.0 / (gamma - 1.0));
        const double p = r.p * std::pow(a / a_r, 2.0 * gamma / (gamma - 1.0));
        return Primitive1D{rho, un, p};
    }
} // namespace

// ================= Public methods =================

ExactIdealGasRiemannSolver::ExactIdealGasRiemannSolver() {}

void ExactIdealGasRiemannSolver::SetXi(const double xi) {
    xi_ = xi;
}

void ExactIdealGasRiemannSolver::SetQ(const double Q) {
    Q_user_ = Q;
}

PrimitiveCell ExactIdealGasRiemannSolver::Sample(const PrimitiveCell& left,
                                                 const PrimitiveCell& right,
                                                 const double gamma,
                                                 const double xi,
                                                 const Axis axis) const {
    const State1D l = MakeState1D(left, axis, gamma);
    const State1D r = MakeState1D(right, axis, gamma);

    const double p_star = SolveStarPressure(l, r, gamma, Q_user_);

    // Vacuum / double rarefaction
    Primitive1D s1d{};
    double u_star = 0.0;

    if (p_star <= 0.0) {
        s1d = SampleVacuum(xi, l, r, gamma);
        // In vacuum case, choose tangentials by xi position relative to 0 (arbitrary but stable):
        // use left if xi <= 0 else right.
        u_star = 0.0;
    }
    else {
        const double f_l = SidePhi(p_star, l, gamma);
        const double f_r = SidePhi(p_star, r, gamma);
        u_star = 0.5 * (l.un + r.un + f_r - f_l);
        s1d = SampleNonVacuum(xi, p_star, u_star, l, r, gamma);
    }

    // Choose tangential velocities from the side of the contact.
    double unL, ut1L, ut2L;
    double unR, ut1R, ut2R;
    riemann::SplitVelocity(left, axis, unL, ut1L, ut2L);
    riemann::SplitVelocity(right, axis, unR, ut1R, ut2R);

    const bool use_left_tangential = (p_star <= 0.0) ? (xi <= u_star) : (xi <= u_star);

    const double ut1 = use_left_tangential ? ut1L : ut1R;
    const double ut2 = use_left_tangential ? ut2L : ut2R;

    PrimitiveCell sample{};
    sample.rho = s1d.rho;
    sample.P = s1d.p;

    riemann::ComposeVelocity(s1d.un, ut1, ut2, axis, sample.u, sample.v, sample.w);
    return sample;
}

auto ExactIdealGasRiemannSolver::ComputeFlux(const PrimitiveCell& left,
                                             const PrimitiveCell& right,
                                             const double gamma,
                                             const Axis axis) const -> FluxCell {
    const PrimitiveCell sample = Sample(left, right, gamma, xi_, axis);
    return EulerFlux(sample, gamma, axis);
}
