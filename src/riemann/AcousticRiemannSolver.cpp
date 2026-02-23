#include "riemann/AcousticRiemannSolver.hpp"

#include <algorithm>
#include <cmath>

#include "data/Variables.hpp"
#include "riemann/RiemannHelpers.hpp"

AcousticRiemannSolver::AcousticRiemannSolver() : rho_min_(1e-10), p_min_(1e-10) {}

namespace {
    inline FluxCell CentralFlux(const PrimitiveCell& left, const PrimitiveCell& right, double gamma, Axis axis) {
        const FluxCell FL = EulerFlux(left, gamma, axis);
        const FluxCell FR = EulerFlux(right, gamma, axis);

        FluxCell F;
        F.mass = 0.5 * (FL.mass + FR.mass);
        F.mom_x = 0.5 * (FL.mom_x + FR.mom_x);
        F.mom_y = 0.5 * (FL.mom_y + FR.mom_y);
        F.mom_z = 0.5 * (FL.mom_z + FR.mom_z);
        F.energy = 0.5 * (FL.energy + FR.energy);
        return F;
    }
} // namespace

auto AcousticRiemannSolver::ComputeFlux(const PrimitiveCell& left,
                                        const PrimitiveCell& right,
                                        const double gamma,
                                        const Axis axis) const -> FluxCell {
    // Work in normal/tangential split.
    double unL, ut1L, ut2L;
    double unR, ut1R, ut2R;
    riemann::SplitVelocity(left, axis, unL, ut1L, ut2L);
    riemann::SplitVelocity(right, axis, unR, ut1R, ut2R);

    const double rhoL = left.rho;
    const double pL = left.P;
    const double rhoR = right.rho;
    const double pR = right.P;

    if (rhoL <= rho_min_ || rhoR <= rho_min_ || pL <= p_min_ || pR <= p_min_) {
        return CentralFlux(left, right, gamma, axis);
    }

    const double cL = std::sqrt(gamma * pL / rhoL);
    const double cR = std::sqrt(gamma * pR / rhoR);

    if (!std::isfinite(cL) || !std::isfinite(cR) || cL <= 0.0 || cR <= 0.0) {
        return CentralFlux(left, right, gamma, axis);
    }

    const double ZL = rhoL * cL;
    const double ZR = rhoR * cR;
    const double denom = ZL + ZR;

    if (!std::isfinite(denom) || denom <= 0.0) {
        return CentralFlux(left, right, gamma, axis);
    }

    const double un_star = (unL * ZL + unR * ZR + (pL - pR)) / denom;

    const double term1 = ZL * ZR * (unL - unR);
    const double term2 = pR * ZR + pL * ZL;
    const double p_star = (term1 + term2) / denom;

    if (!std::isfinite(un_star) || !std::isfinite(p_star)) {
        return CentralFlux(left, right, gamma, axis);
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

    // Tangential velocities: for weak waves, keep them approximately continuous.
    const double ut1_star = 0.5 * (ut1L + ut1R);
    const double ut2_star = 0.5 * (ut2L + ut2R);

    PrimitiveCell w_star{};
    w_star.rho = rho_star;
    w_star.P = p_clamped;

    // Map (un*,ut1*,ut2*) back to (u,v,w) in global components.
    riemann::ComposeVelocity(un_star, ut1_star, ut2_star, axis, w_star.u, w_star.v, w_star.w);

    return EulerFlux(w_star, gamma, axis);
}
