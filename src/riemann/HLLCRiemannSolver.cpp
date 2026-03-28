#include "riemann/HLLCRiemannSolver.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

void HLLCRiemannSolver::SplitVelocity(const PrimitiveCell& state,
                                      const FaceNormal& normal,
                                      double& normal_velocity,
                                      double& tangential_x,
                                      double& tangential_y,
                                      double& tangential_z) const {
    normal_velocity = state.u * normal.x + state.v * normal.y + state.w * normal.z;

    tangential_x = state.u - normal_velocity * normal.x;
    tangential_y = state.v - normal_velocity * normal.y;
    tangential_z = state.w - normal_velocity * normal.z;
}

void HLLCRiemannSolver::BuildStarMomentum(const double rho_star,
                                          const double star_normal_velocity,
                                          const double tangential_x,
                                          const double tangential_y,
                                          const double tangential_z,
                                          const FaceNormal& normal,
                                          double& rhoU_star,
                                          double& rhoV_star,
                                          double& rhoW_star) const {
    const double u_star = tangential_x + star_normal_velocity * normal.x;
    const double v_star = tangential_y + star_normal_velocity * normal.y;
    const double w_star = tangential_z + star_normal_velocity * normal.z;

    rhoU_star = rho_star * u_star;
    rhoV_star = rho_star * v_star;
    rhoW_star = rho_star * w_star;
}

ConservativeCell HLLCRiemannSolver::ComputeFlux(const PrimitiveCell& left,
                                                const PrimitiveCell& right,
                                                const double gamma,
                                                const FaceNormal& normal) const {
    if (!IsUnitNormal(normal)) {
        throw std::runtime_error(
            "HLLCRiemannSolver::ComputeFlux: face normal must be unit-length"
        );
    }

    double un_l = 0.0;
    double utx_l = 0.0;
    double uty_l = 0.0;
    double utz_l = 0.0;
    SplitVelocity(left, normal, un_l, utx_l, uty_l, utz_l);

    double un_r = 0.0;
    double utx_r = 0.0;
    double uty_r = 0.0;
    double utz_r = 0.0;
    SplitVelocity(right, normal, un_r, utx_r, uty_r, utz_r);

    const double a_l = SoundSpeed(left, gamma);
    const double a_r = SoundSpeed(right, gamma);

    const double sL = std::min(un_l - a_l, un_r - a_r);
    const double sR = std::max(un_l + a_l, un_r + a_r);

    const ConservativeCell FL = PhysicalFlux(left, gamma, normal);
    const ConservativeCell FR = PhysicalFlux(right, gamma, normal);

    if (sL >= 0.0) {
        return FL;
    }
    if (sR <= 0.0) {
        return FR;
    }

    const ConservativeCell UL = ConservativeFromPrimitive(left, gamma);
    const ConservativeCell UR = ConservativeFromPrimitive(right, gamma);

    const double rhoL = UL.rho;
    const double rhoR = UR.rho;

    const double pL = left.P;
    const double pR = right.P;

    const double denom =
        rhoL * (sL - un_l) -
        rhoR * (sR - un_r);

    const double sM =
        (std::abs(denom) > 1e-14)
            ? (pR - pL + rhoL * un_l * (sL - un_l) - rhoR * un_r * (sR - un_r)) / denom
            : 0.0;

    const double denom_left = sL - sM;
    const double denom_right = sR - sM;

    const double rho_star_l =
        (std::abs(denom_left) > 1e-14)
            ? rhoL * (sL - un_l) / denom_left
            : rhoL;

    const double rho_star_r =
        (std::abs(denom_right) > 1e-14)
            ? rhoR * (sR - un_r) / denom_right
            : rhoR;

    const double e_tot_l = UL.E / rhoL;
    const double e_tot_r = UR.E / rhoR;

    const double e_star_l =
        rho_star_l * (e_tot_l + (sM - un_l) * (sM + pL / (rhoL * (sL - un_l))));

    const double e_star_r =
        rho_star_r * (e_tot_r + (sM - un_r) * (sM + pR / (rhoR * (sR - un_r))));

    double rhoU_star_l = 0.0;
    double rhoV_star_l = 0.0;
    double rhoW_star_l = 0.0;
    BuildStarMomentum(rho_star_l,
                      sM,
                      utx_l,
                      uty_l,
                      utz_l,
                      normal,
                      rhoU_star_l,
                      rhoV_star_l,
                      rhoW_star_l);

    double rhoU_star_r = 0.0;
    double rhoV_star_r = 0.0;
    double rhoW_star_r = 0.0;
    BuildStarMomentum(rho_star_r,
                      sM,
                      utx_r,
                      uty_r,
                      utz_r,
                      normal,
                      rhoU_star_r,
                      rhoV_star_r,
                      rhoW_star_r);

    if (sM >= 0.0) {
        ConservativeCell flux = FL;
        flux.rho += sL * (rho_star_l - UL.rho);
        flux.rhoU += sL * (rhoU_star_l - UL.rhoU);
        flux.rhoV += sL * (rhoV_star_l - UL.rhoV);
        flux.rhoW += sL * (rhoW_star_l - UL.rhoW);
        flux.E += sL * (e_star_l - UL.E);
        return flux;
    }

    ConservativeCell flux = FR;
    flux.rho += sR * (rho_star_r - UR.rho);
    flux.rhoU += sR * (rhoU_star_r - UR.rhoU);
    flux.rhoV += sR * (rhoV_star_r - UR.rhoV);
    flux.rhoW += sR * (rhoW_star_r - UR.rhoW);
    flux.E += sR * (e_star_r - UR.E);
    return flux;
}
