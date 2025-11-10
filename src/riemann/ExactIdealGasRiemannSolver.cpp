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

State MakeState(const Primitive &w, double gamma) {
    State s;
    s.rho = w.rho;
    s.u = w.u;
    s.p = w.P;
    s.a = std::sqrt(gamma * s.p / s.rho);
    return s;
}

double PhiRarefaction(double p, const State &s, double gamma) {
    const double pratio = p / s.p;
    const double exponent = (gamma - 1.0) / (2.0 * gamma);
    return (2.0 * s.a / (gamma - 1.0)) * (std::pow(pratio, exponent) - 1.0);
}

double PhiRarefactionDerivative(double p, const State &s, double gamma) {
    const double pratio = p / s.p;
    const double exponent = (gamma - 1.0) / (2.0 * gamma);
    const double coeff = (2.0 * s.a / (gamma - 1.0)) * (exponent / s.p);
    return coeff * std::pow(pratio, exponent - 1.0);
}

double PhiShock(double p, const State &s, double gamma) {
    const double A = 2.0 / ((gamma + 1.0) * s.rho);
    const double B = (gamma - 1.0) / (gamma + 1.0) * s.p;
    return (p - s.p) * std::sqrt(A / (p + B));
}

double PhiShockDerivative(double p, const State &s, double gamma) {
    const double A = 2.0 / ((gamma + 1.0) * s.rho);
    const double B = (gamma - 1.0) / (gamma + 1.0) * s.p;
    const double sqrtTerm = std::sqrt(A / (p + B));
    const double term = (p - s.p) / (2.0 * (p + B));
    return sqrtTerm * (1.0 - term);
}

double SidePhi(double p, const State &s, double gamma) {
    if (p <= s.p) {
        return PhiRarefaction(p, s, gamma);
    }
    return PhiShock(p, s, gamma);
}

double SidePhiDerivative(double p, const State &s, double gamma) {
    if (p <= s.p) {
        return PhiRarefactionDerivative(p, s, gamma);
    }
    return PhiShockDerivative(p, s, gamma);
}

double InitialGuess(const State &L, const State &R, double gamma) {
    const double pPV = 0.5 * (L.p + R.p)
                       - 0.125 * (R.u - L.u) * (L.rho + R.rho) * (L.a + R.a);
    const double pMin = std::min(L.p, R.p);
    const double pMax = std::max(L.p, R.p);

    double p0 = pPV;
    if (p0 < 1e-16) {
        p0 = 1e-16;
    }
    if (p0 < pMin) {
        p0 = pMin;
    }
    if (p0 > pMax) {
        p0 = pMax;
    }
    return p0;
}

double SolveStarPressure(const State &L, const State &R, double gamma) {
    double p = InitialGuess(L, R, gamma);

    for (int iter = 0; iter < 40; ++iter) {
        const double fL = SidePhi(p, L, gamma);
        const double fR = SidePhi(p, R, gamma);
        const double f = fL + fR + (R.u - L.u);

        if (std::fabs(f) < 1e-16) {
            break;
        }

        const double dL = SidePhiDerivative(p, L, gamma);
        const double dR = SidePhiDerivative(p, R, gamma);
        const double df = dL + dR;

        if (df == 0.0 || std::isnan(df)) {
            break;
        }

        double pNew = p - f / df;
        if (pNew < 1e-16 || std::isnan(pNew) || std::isinf(pNew)) {
            pNew = 1e-16;
        }

        if (std::fabs(pNew - p) <= 1e-16 * (p + 1e-16)) {
            p = pNew;
            break;
        }

        p = pNew;
    }

    if (!(p > 0.0) || std::isnan(p) || std::isinf(p)) {
        p = 1e-16;
    }

    return p;
}

Primitive SampleAtOrigin(double pStar,
                         double uStar,
                         const State &L,
                         const State &R,
                         double gamma) {
    const double xi = 0.0;

    if (xi <= uStar) {
        if (pStar <= L.p) {
            const double aL = L.a;
            const double sHead = L.u - aL;
            const double rhoStarL = L.rho * std::pow(pStar / L.p, 1.0 / gamma);
            const double aStarL = std::sqrt(gamma * pStar / rhoStarL);
            const double sTail = uStar - aStarL;

            if (xi <= sHead) {
                Primitive w = {L.rho, L.u, L.p};
                return w;
            }
            if (xi >= sTail) {
                Primitive w = {rhoStarL, uStar, pStar};
                return w;
            }

            const double u = (2.0 / (gamma + 1.0)) *
                             (aL + 0.5 * (gamma - 1.0) * L.u + xi);
            const double a = (2.0 / (gamma + 1.0)) *
                             (aL + 0.5 * (gamma - 1.0) * (L.u - xi));
            const double rho = L.rho * std::pow(a / aL, 2.0 / (gamma - 1.0));
            const double p = L.p * std::pow(a / aL, 2.0 * gamma / (gamma - 1.0));
            Primitive w = {rho, u, p};
            return w;
        }

        const double q = pStar / L.p;
        const double sL = L.u - L.a * std::sqrt(0.5 * ((gamma + 1.0) * q + (gamma - 1.0)) / gamma);

        if (xi <= sL) {
            Primitive w = {L.rho, L.u, L.p};
            return w;
        }

        const double factor = (q + (gamma - 1.0) / (gamma + 1.0))
                              / ((gamma - 1.0) / (gamma + 1.0) * q + 1.0);
        const double rhoStarL = L.rho * factor;
        Primitive w = {rhoStarL, uStar, pStar};
        return w;
    }

    if (pStar <= R.p) {
        const double aR = R.a;
        const double sHead = R.u + aR;
        const double rhoStarR = R.rho * std::pow(pStar / R.p, 1.0 / gamma);
        const double aStarR = std::sqrt(gamma * pStar / rhoStarR);
        const double sTail = uStar + aStarR;

        if (xi >= sHead) {
            Primitive w = {R.rho, R.u, R.p};
            return w;
        }
        if (xi <= sTail) {
            Primitive w = {rhoStarR, uStar, pStar};
            return w;
        }

        const double u = (2.0 / (gamma + 1.0)) *
                         (-aR + 0.5 * (gamma - 1.0) * R.u + xi);
        const double a = (2.0 / (gamma + 1.0)) *
                         (-aR + 0.5 * (gamma - 1.0) * (R.u - xi));
        const double rho = R.rho * std::pow(a / aR, 2.0 / (gamma - 1.0));
        const double p = R.p * std::pow(a / aR, 2.0 * gamma / (gamma - 1.0));
        Primitive w = {rho, u, p};
        return w;
    }

    const double q = pStar / R.p;
    const double sR = R.u + R.a * std::sqrt(0.5 * ((gamma + 1.0) * q + (gamma - 1.0)) / gamma);

    if (xi >= sR) {
        Primitive w = {R.rho, R.u, R.p};
        return w;
    }

    const double factor = (q + (gamma - 1.0) / (gamma + 1.0))
                          / ((gamma - 1.0) / (gamma + 1.0) * q + 1.0);
    const double rhoStarR = R.rho * factor;
    Primitive w = {rhoStarR, uStar, pStar};
    return w;
}

Flux HllFallback(const Primitive &left,
                 const Primitive &right,
                 double gamma) {
    const double rhoL = left.rho;
    const double uL = left.u;
    const double pL = left.P;
    const double aL = std::sqrt(gamma * pL / rhoL);

    const double rhoR = right.rho;
    const double uR = right.u;
    const double pR = right.P;
    const double aR = std::sqrt(gamma * pR / rhoR);

    const Conservative UL = EOS::PrimToCons(left, gamma);
    const Conservative UR = EOS::PrimToCons(right, gamma);

    const Flux FL = EulerFlux(left, gamma);
    const Flux FR = EulerFlux(right, gamma);

    const double SL = std::min(uL - aL, uR - aR);
    const double SR = std::max(uL + aL, uR + aR);

    if (SL >= 0.0) {
        return FL;
    }
    if (SR <= 0.0) {
        return FR;
    }

    const double inv = 1.0 / (SR - SL);

    Flux F;
    F.mass = (SR * FL.mass - SL * FR.mass + SL * SR * (UR.rho - UL.rho)) * inv;
    F.momentum = (SR * FL.momentum - SL * FR.momentum + SL * SR * (UR.rhoU - UL.rhoU)) * inv;
    F.energy = (SR * FL.energy - SL * FR.energy + SL * SR * (UR.E - UL.E)) * inv;
    return F;
}
}

Flux ExactIdealGasRiemannSolver::ComputeFlux(const Primitive &left,
                                             const Primitive &right,
                                             double gamma) const {
    if (left.rho <= 0.0 || left.P <= 0.0 ||
        right.rho <= 0.0 || right.P <= 0.0) {
        return HllFallback(left, right, gamma);
    }

    const State L = MakeState(left, gamma);
    const State R = MakeState(right, gamma);

    const double pStar = SolveStarPressure(L, R, gamma);
    if (!(pStar > 0.0) || std::isnan(pStar) || std::isinf(pStar)) {
        return HllFallback(left, right, gamma);
    }

    const double fL = SidePhi(pStar, L, gamma);
    const double fR = SidePhi(pStar, R, gamma);
    const double uStar = 0.5 * (L.u + R.u + fR - fL);

    Primitive sample = SampleAtOrigin(pStar, uStar, L, R, gamma);

    if (sample.rho <= 0.0 || sample.P <= 0.0 ||
        std::isnan(sample.rho) || std::isnan(sample.P)) {
        return HllFallback(left, right, gamma);
    }

    return EulerFlux(sample, gamma);
}
