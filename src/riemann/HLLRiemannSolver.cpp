#include "riemann/HLLRiemannSolver.hpp"

#include <algorithm>
#include <cmath>


auto HLLRiemannSolver::ComputeFlux(const PrimitiveCell& left,
                                   const PrimitiveCell& right,
                                   const double gamma,
                                   const Axis axis) const -> FluxCell {
    const double un_l = NormalVelocity(left, axis);
    const double un_r = NormalVelocity(right, axis);

    const double a_l = SoundSpeed(left, gamma);
    const double a_r = SoundSpeed(right, gamma);

    // Davis (or Einfeldt-like) wave speed estimates
    const double sL = std::min(un_l - a_l, un_r - a_r);
    const double sR = std::max(un_l + a_l, un_r + a_r);

    const FluxCell FL = EulerFlux(left, gamma, axis);
    const FluxCell FR = EulerFlux(right, gamma, axis);

    double rhoL, rhoUL, rhoVL, rhoWL, EL;
    double rhoR, rhoUR, rhoVR, rhoWR, ER;
    PrimitiveToConservative(left, gamma, rhoL, rhoUL, rhoVL, rhoWL, EL);
    PrimitiveToConservative(right, gamma, rhoR, rhoUR, rhoVR, rhoWR, ER);

    if (sL >= 0.0) {
        return FL;
    }
    if (sR <= 0.0) {
        return FR;
    }

    // HLL flux
    FluxCell F;
    const double inv = 1.0 / (sR - sL);

    F.mass = (sR * FL.mass - sL * FR.mass + sL * sR * (rhoR - rhoL)) * inv;
    F.mom_x = (sR * FL.mom_x - sL * FR.mom_x + sL * sR * (rhoUR - rhoUL)) * inv;
    F.mom_y = (sR * FL.mom_y - sL * FR.mom_y + sL * sR * (rhoVR - rhoVL)) * inv;
    F.mom_z = (sR * FL.mom_z - sL * FR.mom_z + sL * sR * (rhoWR - rhoWL)) * inv;
    F.energy = (sR * FL.energy - sL * FR.energy + sL * sR * (ER - EL)) * inv;

    return F;
}
