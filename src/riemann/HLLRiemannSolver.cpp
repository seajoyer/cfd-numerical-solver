#include "riemann/HLLRiemannSolver.hpp"

#include <algorithm>
#include <cmath>

#include "solver/EOS.hpp"

auto HLLRiemannSolver::ComputeFlux(const Primitive& left, const Primitive& right,
                                   double gamma) const -> Flux {
    const double rho_L = left.rho;
    const double u_L = left.u;
    const double p_L = left.P;

    const double rho_R = right.rho;
    const double u_R = right.u;
    const double p_R = right.P;

    const double a_L = std::sqrt(gamma * p_L / rho_L);
    const double a_R = std::sqrt(gamma * p_R / rho_R);

    const Conservative U_L = EOS::PrimToCons(left, gamma);
    const Conservative U_R = EOS::PrimToCons(right, gamma);

    const Flux F_L = EulerFlux(left, gamma);
    const Flux F_R = EulerFlux(right, gamma);

    const double S_L = std::min(u_L - a_L, u_R - a_R);
    const double S_R = std::max(u_L + a_L, u_R + a_R);

    if (S_L >= 0.0) {
        return F_L;
    }

    if (S_R <= 0.0) {
        return F_R;
    }

    const double inv_den = 1.0 / (S_R - S_L);

    Flux F;
    F.mass =
        (S_R * F_L.mass - S_L * F_R.mass + S_L * S_R * (U_R.rho - U_L.rho)) * inv_den;
    F.momentum =
        (S_R * F_L.momentum - S_L * F_R.momentum + S_L * S_R * (U_R.rhoU - U_L.rhoU)) *
        inv_den;
    F.energy =
        (S_R * F_L.energy - S_L * F_R.energy + S_L * S_R * (U_R.E - U_L.E)) * inv_den;

    return F;
}
