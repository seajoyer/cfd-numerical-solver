#include "riemann/HLLRiemannSolver.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

ConservativeCell HLLRiemannSolver::ComputeFlux(const PrimitiveCell& left,
                                               const PrimitiveCell& right,
                                               const double gamma,
                                               const FaceNormal& normal) const {
    if (!IsUnitNormal(normal)) {
        throw std::runtime_error(
            "HLLRiemannSolver::ComputeFlux: face normal must be unit-length"
        );
    }

    const double un_left = NormalVelocity(left, normal);
    const double un_right = NormalVelocity(right, normal);

    const double a_left = SoundSpeed(left, gamma);
    const double a_right = SoundSpeed(right, gamma);

    const double sL = std::min(un_left - a_left, un_right - a_right);
    const double sR = std::max(un_left + a_left, un_right + a_right);

    const ConservativeCell flux_left = PhysicalFlux(left, gamma, normal);
    const ConservativeCell flux_right = PhysicalFlux(right, gamma, normal);

    const ConservativeCell U_left = ConservativeFromPrimitive(left, gamma);
    const ConservativeCell U_right = ConservativeFromPrimitive(right, gamma);

    if (sL >= 0.0) {
        return flux_left;
    }

    if (sR <= 0.0) {
        return flux_right;
    }

    const double denom = sR - sL;
    if (std::abs(denom) <= 1e-14) {
        throw std::runtime_error(
            "HLLRiemannSolver::ComputeFlux: degenerate wave-speed interval"
        );
    }

    const double inv = 1.0 / denom;

    ConservativeCell flux;
    flux.rho = (sR * flux_left.rho - sL * flux_right.rho + sL * sR * (U_right.rho - U_left.rho)) * inv;
    flux.rhoU = (sR * flux_left.rhoU - sL * flux_right.rhoU + sL * sR * (U_right.rhoU - U_left.rhoU)) * inv;
    flux.rhoV = (sR * flux_left.rhoV - sL * flux_right.rhoV + sL * sR * (U_right.rhoV - U_left.rhoV)) * inv;
    flux.rhoW = (sR * flux_left.rhoW - sL * flux_right.rhoW + sL * sR * (U_right.rhoW - U_left.rhoW)) * inv;
    flux.E = (sR * flux_left.E - sL * flux_right.E + sL * sR * (U_right.E - U_left.E)) * inv;

    return flux;
}
