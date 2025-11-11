#include "riemann/HLLCRiemannSolver.hpp"

#include <algorithm>
#include <cmath>

#include "solver/EOS.hpp"

auto HLLCRiemannSolver::ComputeFlux(const Primitive& left, const Primitive& right,
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

    if (S_R - S_L < 1e-14) {
        Flux F;
        F.mass = 0.5 * (F_L.mass + F_R.mass);
        F.momentum = 0.5 * (F_L.momentum + F_R.momentum);
        F.energy = 0.5 * (F_L.energy + F_R.energy);
        return F;
    }

    const double numerator =
        p_R - p_L + rho_L * (S_L - u_L) * u_L - rho_R * (S_R - u_R) * u_R;
    const double denominator = rho_L * (S_L - u_L) - rho_R * (S_R - u_R);

    const double SM = numerator / denominator;

    const double rho_star_L = rho_L * (S_L - u_L) / (S_L - SM);
    const double rho_star_R = rho_R * (S_R - u_R) / (S_R - SM);

    const double p_star_L = p_L + rho_L * (S_L - u_L) * (SM - u_L);
    const double p_star_R = p_R + rho_R * (S_R - u_R) * (SM - u_R);

    const double E_L = U_L.E;
    const double E_R = U_R.E;

    const double E_star_L = ((S_L - u_L) * E_L - p_L * u_L + p_star_L * SM) / (S_L - SM);
    const double E_star_R = ((S_R - u_R) * E_R - p_R * u_R + p_star_R * SM) / (S_R - SM);

    Conservative U_star_L;
    U_star_L.rho = rho_star_L;
    U_star_L.rhoU = rho_star_L * SM;
    U_star_L.E = E_star_L;

    Conservative U_star_R;
    U_star_R.rho = rho_star_R;
    U_star_R.rhoU = rho_star_R * SM;
    U_star_R.E = E_star_R;

    if (0.0 <= S_L) {
        return F_L;
    } else if (S_L <= 0.0 && 0.0 <= SM) {
        Flux F;
        F.mass = F_L.mass + S_L * (U_star_L.rho - U_L.rho);
        F.momentum = F_L.momentum + S_L * (U_star_L.rhoU - U_L.rhoU);
        F.energy = F_L.energy + S_L * (U_star_L.E - U_L.E);
        return F;
    } else if (SM <= 0.0 && 0.0 <= S_R) {
        Flux F;
        F.mass = F_R.mass + S_R * (U_star_R.rho - U_R.rho);
        F.momentum = F_R.momentum + S_R * (U_star_R.rhoU - U_R.rhoU);
        F.energy = F_R.energy + S_R * (U_star_R.E - U_R.E);
        return F;
    } else {
        return F_R;
    }
}
