#include "riemann/HLLRiemannSolver.hpp"
#include <algorithm>
#include <cmath>

auto HLLRiemannSolver::ComputeFlux(const Primitive& left,
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

    const double sl = std::min(u_l - a_l, u_r - a_r);
    const double sr = std::max(u_l + a_l, u_r + a_r);

    if (sl >= 0.0) {
        return fl;
    }

    if (sr <= 0.0) {
        return fr;
    }

    const double inv_den = 1.0 / (sr - sl);

    Flux f = (sr * fl - sl * fr + sl * sr * (ur - ul)) * inv_den;

    return f;
}