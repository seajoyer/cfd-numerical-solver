#include "riemann/HLLCRiemannSolver.hpp"
#include <algorithm>
#include <cmath>

Flux HLLCRiemannSolver::ComputeFlux(const Primitive &left,
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

    if (SR - SL < 1e-14) {
        Flux F;
        F.mass = 0.5 * (FL.mass + FR.mass);
        F.momentum = 0.5 * (FL.momentum + FR.momentum);
        F.energy = 0.5 * (FL.energy + FR.energy);
        return F;
    }

    const double numerator = pR - pL
                             + rhoL * (SL - uL) * uL
                             - rhoR * (SR - uR) * uR;
    const double denominator = rhoL * (SL - uL) - rhoR * (SR - uR);

    const double SM = numerator / denominator;

    const double rhoLStar = rhoL * (SL - uL) / (SL - SM);
    const double rhoRStar = rhoR * (SR - uR) / (SR - SM);

    const double pLStar = pL + rhoL * (SL - uL) * (SM - uL);
    const double pRStar = pR + rhoR * (SR - uR) * (SM - uR);

    const double EL = UL.E;
    const double ER = UR.E;

    const double ELStar = ((SL - uL) * EL - pL * uL + pLStar * SM) / (SL - SM);
    const double ERStar = ((SR - uR) * ER - pR * uR + pRStar * SM) / (SR - SM);

    Conservative ULStar;
    ULStar.rho = rhoLStar;
    ULStar.rhoU = rhoLStar * SM;
    ULStar.E = ELStar;

    Conservative URStar;
    URStar.rho = rhoRStar;
    URStar.rhoU = rhoRStar * SM;
    URStar.E = ERStar;

    if (0.0 <= SL) {
        return FL;
    } else if (SL <= 0.0 && 0.0 <= SM) {
        Flux F;
        F.mass = FL.mass + SL * (ULStar.rho - UL.rho);
        F.momentum = FL.momentum + SL * (ULStar.rhoU - UL.rhoU);
        F.energy = FL.energy + SL * (ULStar.E - UL.E);
        return F;
    } else if (SM <= 0.0 && 0.0 <= SR) {
        Flux F;
        F.mass = FR.mass + SR * (URStar.rho - UR.rho);
        F.momentum = FR.momentum + SR * (URStar.rhoU - UR.rhoU);
        F.energy = FR.energy + SR * (URStar.E - UR.E);
        return F;
    } else {
        return FR;
    }
}
