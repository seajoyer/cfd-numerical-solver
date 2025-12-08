#include "riemann/RusanovRiemannSolver.hpp"

#include <algorithm>
#include <cmath>

#include "solver/EOS.hpp"

auto RusanovRiemannSolver::ComputeFlux(const Primitive& left,
                                       const Primitive& right,
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

    const double s_l = std::abs(u_l) + a_l;
    const double s_r = std::abs(u_r) + a_r;
    const double a_max = std::max(s_l, s_r);

    const Flux f = 0.5 * (fl + fr + a_max * (ul - ur));

    return f;
}