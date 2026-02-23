#include "riemann/RusanovRiemannSolver.hpp"

#include <algorithm>
#include <cmath>

auto RusanovRiemannSolver::ComputeFlux(const PrimitiveCell& left,
                                       const PrimitiveCell& right,
                                       const double gamma,
                                       const Axis axis) const -> FluxCell {
    const double un_l = NormalVelocity(left, axis);
    const double un_r = NormalVelocity(right, axis);

    const double a_l = SoundSpeed(left, gamma);
    const double a_r = SoundSpeed(right, gamma);

    const double s_l = std::abs(un_l) + a_l;
    const double s_r = std::abs(un_r) + a_r;
    const double a_max = std::max(s_l, s_r);

    const FluxCell FL = EulerFlux(left, gamma, axis);
    const FluxCell FR = EulerFlux(right, gamma, axis);

    double rhoL, rhoUL, rhoVL, rhoWL, EL;
    double rhoR, rhoUR, rhoVR, rhoWR, ER;
    PrimitiveToConservative(left, gamma, rhoL, rhoUL, rhoVL, rhoWL, EL);
    PrimitiveToConservative(right, gamma, rhoR, rhoUR, rhoVR, rhoWR, ER);

    FluxCell F;
    F.mass = 0.5 * (FL.mass + FR.mass - a_max * (rhoR - rhoL));
    F.mom_x = 0.5 * (FL.mom_x + FR.mom_x - a_max * (rhoUR - rhoUL));
    F.mom_y = 0.5 * (FL.mom_y + FR.mom_y - a_max * (rhoVR - rhoVL));
    F.mom_z = 0.5 * (FL.mom_z + FR.mom_z - a_max * (rhoWR - rhoWL));
    F.energy = 0.5 * (FL.energy + FR.energy - a_max * (ER - EL));
    return F;
}
