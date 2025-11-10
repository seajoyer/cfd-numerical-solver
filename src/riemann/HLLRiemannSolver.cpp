#include "riemann/HLLRiemannSolver.hpp"
#include <algorithm>
#include <cmath>

Flux HLLRiemannSolver::ComputeFlux(const Primitive &left,
                                   const Primitive &right,
                                   double gamma) const {
    const double rhoL = left.rho;
    const double uL = left.u;
    const double pL = left.P;

    const double rhoR = right.rho;
    const double uR = right.u;
    const double pR = right.P;

    const double aL = std::sqrt(gamma * pL / rhoL);
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

    const double invDen = 1.0 / (SR - SL);

    Flux F;
    F.mass = (SR * FL.mass - SL * FR.mass + SL * SR * (UR.rho - UL.rho)) * invDen;
    F.momentum = (SR * FL.momentum - SL * FR.momentum + SL * SR * (UR.rhoU - UL.rhoU)) * invDen;
    F.energy = (SR * FL.energy - SL * FR.energy + SL * SR * (UR.E - UL.E)) * invDen;

    return F;
}
