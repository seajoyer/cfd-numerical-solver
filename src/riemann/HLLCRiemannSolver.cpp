#include "riemann/HLLCRiemannSolver.hpp"

#include <algorithm>
#include <cmath>

auto HLLCRiemannSolver::ComputeFlux(const PrimitiveCell& left,
                                    const PrimitiveCell& right,
                                    const double gamma,
                                    const Axis axis) const -> FluxCell {
    // Split velocities into normal/tangential for the interface normal axis.
    double un_l, ut1_l, ut2_l;
    double un_r, ut1_r, ut2_r;
    riemann::SplitVelocity(left, axis, un_l, ut1_l, ut2_l);
    riemann::SplitVelocity(right, axis, un_r, ut1_r, ut2_r);

    const double a_l = SoundSpeed(left, gamma);
    const double a_r = SoundSpeed(right, gamma);

    // Davis wave speed estimates (simple and robust)
    const double sL = std::min(un_l - a_l, un_r - a_r);
    const double sR = std::max(un_l + a_l, un_r + a_r);

    const FluxCell FL = EulerFlux(left, gamma, axis);
    const FluxCell FR = EulerFlux(right, gamma, axis);

    // Upwind selection if star region not needed.
    if (sL >= 0.0) return FL;
    if (sR <= 0.0) return FR;

    // Conservative states (global ordering rho, rhoU, rhoV, rhoW, E)
    double rhoL, rhoUL, rhoVL, rhoWL, EL;
    double rhoR, rhoUR, rhoVR, rhoWR, ER;
    PrimitiveToConservative(left, gamma, rhoL, rhoUL, rhoVL, rhoWL, EL);
    PrimitiveToConservative(right, gamma, rhoR, rhoUR, rhoVR, rhoWR, ER);

    // Contact wave speed S*
    const double pL = left.P;
    const double pR = right.P;

    const double denom = rhoL * (sL - un_l) - rhoR * (sR - un_r);
    const double sM =
        (denom != 0.0)
            ? (pR - pL + rhoL * un_l * (sL - un_l) - rhoR * un_r * (sR - un_r)) / denom
            : 0.0;

    // Left star state density
    const double rhoStarL = rhoL * (sL - un_l) / (sL - sM);
    // Right star state density
    const double rhoStarR = rhoR * (sR - un_r) / (sR - sM);

    // Total specific energy
    const double eTotL = EL / rhoL;
    const double eTotR = ER / rhoR;

    // Star-region total energy density (Toro-style)
    const double EStarL = rhoStarL * (eTotL + (sM - un_l) * (sM + pL / (rhoL * (sL - un_l))));
    const double EStarR = rhoStarR * (eTotR + (sM - un_r) * (sM + pR / (rhoR * (sR - un_r))));

    // Star momenta in global ordering: normal momentum uses sM; tangentials copied from each side.
    double rhoUStarL, rhoVStarL, rhoWStarL;
    double rhoUStarR, rhoVStarR, rhoWStarR;
    riemann::ComposeMomentum(rhoStarL, sM, ut1_l, ut2_l, axis, rhoUStarL, rhoVStarL, rhoWStarL);
    riemann::ComposeMomentum(rhoStarR, sM, ut1_r, ut2_r, axis, rhoUStarR, rhoVStarR, rhoWStarR);

    // Flux selection with star correction
    if (sM >= 0.0) {
        FluxCell F = FL;
        F.mass += sL * (rhoStarL - rhoL);
        F.mom_x += sL * (rhoUStarL - rhoUL);
        F.mom_y += sL * (rhoVStarL - rhoVL);
        F.mom_z += sL * (rhoWStarL - rhoWL);
        F.energy += sL * (EStarL - EL);
        return F;
    }

    FluxCell F = FR;
    F.mass += sR * (rhoStarR - rhoR);
    F.mom_x += sR * (rhoUStarR - rhoUR);
    F.mom_y += sR * (rhoVStarR - rhoVR);
    F.mom_z += sR * (rhoWStarR - rhoWR);
    F.energy += sR * (EStarR - ER);
    return F;
}
