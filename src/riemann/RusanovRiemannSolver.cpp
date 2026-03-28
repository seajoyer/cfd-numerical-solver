#include "riemann/RusanovRiemannSolver.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

ConservativeCell RusanovRiemannSolver::ComputeFlux(const PrimitiveCell& left,
                                                   const PrimitiveCell& right,
                                                   const double gamma,
                                                   const FaceNormal& normal) const {
    if (!IsUnitNormal(normal)) {
        throw std::runtime_error(
            "RusanovRiemannSolver::ComputeFlux: face normal must be unit-length"
        );
    }

    const double un_left = NormalVelocity(left, normal);
    const double un_right = NormalVelocity(right, normal);

    const double a_left = SoundSpeed(left, gamma);
    const double a_right = SoundSpeed(right, gamma);

    const double s_left = std::abs(un_left) + a_left;
    const double s_right = std::abs(un_right) + a_right;
    const double a_max = std::max(s_left, s_right);

    const ConservativeCell flux_left = PhysicalFlux(left, gamma, normal);
    const ConservativeCell flux_right = PhysicalFlux(right, gamma, normal);

    const ConservativeCell U_left = ConservativeFromPrimitive(left, gamma);
    const ConservativeCell U_right = ConservativeFromPrimitive(right, gamma);

    ConservativeCell flux;
    flux.rho = 0.5 * (flux_left.rho + flux_right.rho - a_max * (U_right.rho - U_left.rho));
    flux.rhoU = 0.5 * (flux_left.rhoU + flux_right.rhoU - a_max * (U_right.rhoU - U_left.rhoU));
    flux.rhoV = 0.5 * (flux_left.rhoV + flux_right.rhoV - a_max * (U_right.rhoV - U_left.rhoV));
    flux.rhoW = 0.5 * (flux_left.rhoW + flux_right.rhoW - a_max * (U_right.rhoW - U_left.rhoW));
    flux.E = 0.5 * (flux_left.E + flux_right.E - a_max * (U_right.E - U_left.E));

    return flux;
}
