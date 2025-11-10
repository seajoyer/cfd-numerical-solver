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

auto MakeState(const Primitive &w, double gamma) -> State {
    State s;
    s.rho = w.rho;
    s.u = w.u;
    s.p = w.P;
    s.a = std::sqrt(gamma * s.p / s.rho);
    return s;
}

auto PhiRarefaction(double p, const State &s, double gamma) -> double {
    const double pratio = p / s.p;
    const double exponent = (gamma - 1.0) / (2.0 * gamma);
    return (2.0 * s.a / (gamma - 1.0)) * (std::pow(pratio, exponent) - 1.0);
}

auto PhiRarefactionDerivative(double p, const State &s, double gamma) -> double {
    const double pratio = p / s.p;
    const double exponent = (gamma - 1.0) / (2.0 * gamma);
    const double coeff = (2.0 * s.a / (gamma - 1.0)) * (exponent / s.p);
    return coeff * std::pow(pratio, exponent - 1.0);
}

auto PhiShock(double p, const State &s, double gamma) -> double {
    const double A = 2.0 / ((gamma + 1.0) * s.rho);
    const double B = (gamma - 1.0) / (gamma + 1.0) * s.p;
    return (p - s.p) * std::sqrt(A / (p + B));
}

auto PhiShockDerivative(double p, const State &s, double gamma) -> double {
    const double A = 2.0 / ((gamma + 1.0) * s.rho);
    const double B = (gamma - 1.0) / (gamma + 1.0) * s.p;
    const double sqrt_term = std::sqrt(A / (p + B));
    const double term = (p - s.p) / (2.0 * (p + B));
    return sqrt_term * (1.0 - term);
}

auto SidePhi(double p, const State &s, double gamma) -> double {
    if (p <= s.p) {
        return PhiRarefaction(p, s, gamma);
    }
    return PhiShock(p, s, gamma);
}

auto SidePhiDerivative(double p, const State &s, double gamma) -> double {
    if (p <= s.p) {
        return PhiRarefactionDerivative(p, s, gamma);
    }
    return PhiShockDerivative(p, s, gamma);
}

auto InitialGuess(const State &L, const State &R, double Q_user, double gamma) -> double {
    const double p_pvrs = 0.5 * (L.p + R.p)
                       - 0.125 * (R.u - L.u) * (L.rho + R.rho) * (L.a + R.a);
    const double p_min = std::min(L.p, R.p);
    const double p_max = std::max(L.p, R.p);

    const double z = (gamma - 1) / (2 * gamma);
    const double Q = p_max / p_min;

    double p_star = p_pvrs;
    if (p_pvrs < 1e-16) {
        p_star = 1e-16;
        return p_star;
    }

    if (p_min < p_pvrs and p_pvrs < p_max and Q < Q_user) {
        // PVRS
        return p_star;
    }
    if (p_pvrs < p_min) {
        // TRRS
        p_star = pow((L.a + R.a - (gamma - 1) / 2 * (R.u - L.u)) / (L.a / pow(L.p, z) + R.a / pow(R.p, z)), 1 / z);
        return p_star;
    }
    
    // TSRS
    auto A = [gamma](const auto& side) -> auto { return 2 / ((gamma + 1) * side.rho); };
    auto B = [gamma](const auto& side) -> auto { return (gamma - 1) / (gamma + 1) * side.p; };
    auto g = [p_pvrs, A, B](const auto& side) -> auto { return pow(A(side) / (p_pvrs + B(side)), 0.5); };

    p_star = (g(L) * L.p + g(R) * R.p - (R.u - L.u)) / (g(L) + g(R));

    return p_star;
}

auto SolveStarPressure(const State &L, const State &R, double Q_user, double gamma) -> double {
    double p = InitialGuess(L, R, Q_user, gamma);

    for (int iter = 0; iter < 40; ++iter) {
        const double f_L = SidePhi(p, L, gamma);
        const double f_R = SidePhi(p, R, gamma);
        const double f = f_L + f_R + (R.u - L.u);

        if (std::fabs(f) < 1e-16) {
            break;
        }

        const double dL = SidePhiDerivative(p, L, gamma);
        const double dR = SidePhiDerivative(p, R, gamma);
        const double df = dL + dR;

        if (df == 0.0 || std::isnan(df)) {
            break;
        }

        double p_new = p - f / df;
        if (p_new < 1e-16 || std::isnan(p_new) || std::isinf(p_new)) {
            p_new = 1e-16;
        }

        if (std::fabs(p_new - p) <= 1e-16 * (p + 1e-16)) {
            p = p_new;
            break;
        }

        p = p_new;
    }

    if (!(p > 0.0) || std::isnan(p) || std::isinf(p)) {
        p = 1e-16;
    }

    return p;
}

auto SampleAtOrigin(double p_star,
                         double u_star,
                         const State &L,
                         const State &R,
                         double gamma) -> Primitive {
    const double xi = 0.0;

    if (xi <= u_star) {
        if (p_star <= L.p) {
            const double a_L = L.a;
            const double s_head = L.u - a_L;
            const double rho_star_L = L.rho * std::pow(p_star / L.p, 1.0 / gamma);
            const double a_star_L = std::sqrt(gamma * p_star / rho_star_L);
            const double s_tail = u_star - a_star_L;

            if (xi <= s_head) {
                Primitive w = {.rho=L.rho, .u=L.u, .P=L.p};
                return w;
            }
            if (xi >= s_tail) {
                Primitive w = {.rho=rho_star_L, .u=u_star, .P=p_star};
                return w;
            }

            const double u = (2.0 / (gamma + 1.0)) *
                             (a_L + 0.5 * (gamma - 1.0) * L.u + xi);
            const double a = (2.0 / (gamma + 1.0)) *
                             (a_L + 0.5 * (gamma - 1.0) * (L.u - xi));
            const double rho = L.rho * std::pow(a / a_L, 2.0 / (gamma - 1.0));
            const double p = L.p * std::pow(a / a_L, 2.0 * gamma / (gamma - 1.0));
            Primitive w = {.rho=rho, .u=u, .P=p};
            return w;
        }

        const double q = p_star / L.p;
        const double s_L = L.u - L.a * std::sqrt(0.5 * ((gamma + 1.0) * q + (gamma - 1.0)) / gamma);

        if (xi <= s_L) {
            Primitive w = {.rho=L.rho, .u=L.u, .P=L.p};
            return w;
        }

        const double factor = (q + (gamma - 1.0) / (gamma + 1.0))
                              / ((gamma - 1.0) / (gamma + 1.0) * q + 1.0);
        const double rho_star_L = L.rho * factor;
        Primitive w = {.rho=rho_star_L, .u=u_star, .P=p_star};
        return w;
    }

    if (p_star <= R.p) {
        const double a_R = R.a;
        const double s_head = R.u + a_R;
        const double rho_star_R = R.rho * std::pow(p_star / R.p, 1.0 / gamma);
        const double a_star_R = std::sqrt(gamma * p_star / rho_star_R);
        const double s_tail = u_star + a_star_R;

        if (xi >= s_head) {
            Primitive w = {.rho=R.rho, .u=R.u, .P=R.p};
            return w;
        }
        if (xi <= s_tail) {
            Primitive w = {.rho=rho_star_R, .u=u_star, .P=p_star};
            return w;
        }

        const double u = (2.0 / (gamma + 1.0)) *
                         (-a_R + 0.5 * (gamma - 1.0) * R.u + xi);
        const double a = (2.0 / (gamma + 1.0)) *
                         (-a_R + 0.5 * (gamma - 1.0) * (R.u - xi));
        const double rho = R.rho * std::pow(a / a_R, 2.0 / (gamma - 1.0));
        const double p = R.p * std::pow(a / a_R, 2.0 * gamma / (gamma - 1.0));
        Primitive w = {.rho=rho, .u=u, .P=p};
        return w;
    }

    const double q = p_star / R.p;
    const double s_R = R.u + R.a * std::sqrt(0.5 * ((gamma + 1.0) * q + (gamma - 1.0)) / gamma);

    if (xi >= s_R) {
        Primitive w = {.rho=R.rho, .u=R.u, .P=R.p};
        return w;
    }

    const double factor = (q + (gamma - 1.0) / (gamma + 1.0))
                          / ((gamma - 1.0) / (gamma + 1.0) * q + 1.0);
    const double rho_star_R = R.rho * factor;
    Primitive w = {.rho=rho_star_R, .u=u_star, .P=p_star};
    return w;
}

auto HllFallback(const Primitive &left,
                 const Primitive &right,
                 double gamma) -> Flux {
    const double rho_L = left.rho;
    const double u_L = left.u;
    const double p_L = left.P;
    const double a_L = std::sqrt(gamma * p_L / rho_L);

    const double rho_R = right.rho;
    const double u_R = right.u;
    const double p_R = right.P;
    const double a_R = std::sqrt(gamma * p_R / rho_R);

    const Conservative U_L = EOS::PrimToCons(left, gamma);
    const Conservative U_R = EOS::PrimToCons(right, gamma);

    const Flux F_L = EulerFlux(left, gamma);
    const Flux F_R = EulerFlux(right, gamma);

    const double S_L = std::min(u_L - a_L, u_R - a_R);
    const double S_R = std::max(u_L + a_L, u_R + a_R);

    if (S_L >= 0.0) {
        return F_L;
    }
    if (S_R <= 0.0) {
        return F_R;
    }

    const double inv = 1.0 / (S_R - S_L);

    Flux F;
    F.mass = (S_R * F_L.mass - S_L * F_R.mass + S_L * S_R * (U_R.rho - U_L.rho)) * inv;
    F.momentum = (S_R * F_L.momentum - S_L * F_R.momentum + S_L * S_R * (U_R.rhoU - U_L.rhoU)) * inv;
    F.energy = (S_R * F_L.energy - S_L * F_R.energy + S_L * S_R * (U_R.E - U_L.E)) * inv;
    return F;
}
}

auto ExactIdealGasRiemannSolver::ComputeFlux(const Primitive &left,
                                             const Primitive &right,
                                             double Q_user,
                                             double gamma) const -> Flux {
    if (left.rho <= 0.0 || left.P <= 0.0 ||
        right.rho <= 0.0 || right.P <= 0.0) {
        return HllFallback(left, right, gamma);
    }

    const State L = MakeState(left, gamma);
    const State R = MakeState(right, gamma);

    const double p_star = SolveStarPressure(L, R, Q_user, gamma);
    if (!(p_star > 0.0) || std::isnan(p_star) || std::isinf(p_star)) {
        return HllFallback(left, right, gamma);
    }

    const double f_L = SidePhi(p_star, L, gamma);
    const double f_R = SidePhi(p_star, R, gamma);
    const double u_star = 0.5 * (L.u + R.u + f_R - f_L);

    Primitive sample = SampleAtOrigin(p_star, u_star, L, R, gamma);

    if (sample.rho <= 0.0 || sample.P <= 0.0 ||
        std::isnan(sample.rho) || std::isnan(sample.P)) {
        return HllFallback(left, right, gamma);
    }

    return EulerFlux(sample, gamma);
}
