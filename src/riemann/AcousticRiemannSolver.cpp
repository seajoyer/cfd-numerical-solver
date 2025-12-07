#include "riemann/AcousticRiemannSolver.hpp"

#include <algorithm>
#include <cmath>

AcousticRiemannSolver::AcousticRiemannSolver() : rho_min_(1e-10), p_min_(1e-10) {}

static auto CentralFlux(const Primitive& left, const Primitive& right, double gamma)
    -> Flux {
    const Flux FL = EulerFlux(left, gamma);
    const Flux FR = EulerFlux(right, gamma);
    return 0.5 * (FL + FR);
}

auto AcousticRiemannSolver::ComputeFlux(const Primitive& left, const Primitive& right,
                                        const double gamma) const -> Flux {
    const double rhoL = left.rho;
    const double uL = left.u;
    const double pL = left.P;

    const double rhoR = right.rho;
    const double uR = right.u;
    const double pR = right.P;

    if (rhoL <= rho_min_ || rhoR <= rho_min_ || pL <= p_min_ || pR <= p_min_) {
        return CentralFlux(left, right, gamma);
    }

    const double cL = std::sqrt(gamma * pL / rhoL);
    const double cR = std::sqrt(gamma * pR / rhoR);

    if (!std::isfinite(cL) || !std::isfinite(cR) || cL <= 0.0 || cR <= 0.0) {
        return CentralFlux(left, right, gamma);
    }

    const double ZL = rhoL * cL;
    const double ZR = rhoR * cR;
    const double denom = ZL + ZR;

    if (!std::isfinite(denom) || denom <= 0.0) {
        return CentralFlux(left, right, gamma);
    }

    const double u_star = (uL * ZL + uR * ZR + (pL - pR)) / denom;

    const double term1 = ZL * ZR * (uL - uR);
    const double term2 = pR * ZR + pL * ZL;
    const double p_star = (term1 + term2) / denom;

    if (!std::isfinite(u_star) || !std::isfinite(p_star)) {
        return CentralFlux(left, right, gamma);
    }

    double rhoL_star = rhoL + (p_star - pL) / (cL * cL);
    double rhoR_star = rhoR + (p_star - pR) / (cR * cR);

    rhoL_star = std::max(rhoL_star, rho_min_);
    rhoR_star = std::max(rhoR_star, rho_min_);

    double rho_star = 0.5 * (rhoL_star + rhoR_star);
    if (!std::isfinite(rho_star) || rho_star <= 0.0) {
        rho_star = std::max(0.5 * (rhoL + rhoR), rho_min_);
    }

    double p_clamped = std::max(p_star, p_min_);
    if (!std::isfinite(p_clamped)) {
        p_clamped = std::max(0.5 * (pL + pR), p_min_);
    }

    Primitive w_star;
    w_star.rho = rho_star;
    w_star.u = u_star;
    w_star.P = p_clamped;

    return EulerFlux(w_star, gamma);
}
