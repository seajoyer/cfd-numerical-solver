#include "riemann/HLLCRiemannSolver.hpp"

#include <algorithm>
#include <cmath>

#include "solver/EOS.hpp"

auto HLLCRiemannSolver::ComputeFlux(const Primitive& left, const Primitive& right,
                                    double gamma) const -> Flux {
    const double rho_l = left.rho;
    const double u_l = left.u;
    const double p_l = left.P;

    const double rho_r = right.rho;
    const double u_r = right.u;
    const double p_r = right.P;

    const double a_l = std::sqrt(gamma * p_l / rho_l);
    const double a_r = std::sqrt(gamma * p_r / rho_r);

    const Conservative ul = EOS::PrimToCons(left, gamma);
    const Conservative ur = EOS::PrimToCons(right, gamma);

    const Flux fl = EulerFlux(left, gamma);
    const Flux fr = EulerFlux(right, gamma);

    const double sl = std::min(u_l - a_l, u_r - a_r);
    const double sr = std::max(u_l + a_l, u_r + a_r);

    if (sl >= 0.0) {
        return fl;
    }

    if (sr <= 0.0) {
        return fr;
    }

    if (sr - sl < 1e-14) {
        return 0.5 * (fl + fr);
    }

    const double numerator =
        p_r - p_l + rho_l * (sl - u_l) * u_l - rho_r * (sr - u_r) * u_r;
    const double denominator = rho_l * (sl - u_l) - rho_r * (sr - u_r);

    const double sm = numerator / denominator;

    const double rho_l_star = rho_l * (sl - u_l) / (sl - sm);
    const double rho_r_star = rho_r * (sr - u_r) / (sr - sm);

    const double p_l_star = p_l + rho_l * (sl - u_l) * (sm - u_l);
    const double p_r_star = p_r + rho_r * (sr - u_r) * (sm - u_r);

    const double el = ul.E;
    const double er = ur.E;

    const double el_star = ((sl - u_l) * el - p_l * u_l + p_l_star * sm) / (sl - sm);
    const double er_star = ((sr - u_r) * er - p_r * u_r + p_r_star * sm) / (sr - sm);

    Conservative ul_star;
    ul_star.rho = rho_l_star;
    ul_star.rhoU = rho_l_star * sm;
    ul_star.E = el_star;

    Conservative ur_star;
    ur_star.rho = rho_r_star;
    ur_star.rhoU = rho_r_star * sm;
    ur_star.E = er_star;

    if (0.0 <= sl) {
        return fl;
    }
    if (sl <= 0.0 && 0.0 <= sm) {
        return fl + sl * (ul_star - ul);
    }
    if (sm <= 0.0 && 0.0 <= sr) {
        return fr + sr * (ur_star - ur);
    }
    return fr;
}
