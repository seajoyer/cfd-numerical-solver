#include "solver/EOS.hpp"

#include <cmath>

auto EOS::PrimToCons(const Primitive& w, double gamma) -> Conservative {
    Conservative uc;
    uc.rho = w.rho;
    uc.rhoU = w.rho * w.u;
    uc.rhoV = w.rho * w.v;
    uc.E = w.P / (gamma - 1.0) + 0.5 * w.rho * (w.u * w.u + w.v * w.v);
    return uc;
}

auto EOS::ConsToPrim(const Conservative& uc, double gamma) -> Primitive {
    Primitive w;
    w.rho = uc.rho;
    if (uc.rho <= 0.0) {
        return {0.0, 0.0, 0.0, 0.0};
    }
    const double inv_rho = 1.0 / uc.rho;
    w.u = uc.rhoU * inv_rho;
    w.v = uc.rhoV * inv_rho;
    w.P = Pressure(uc, gamma);
    return w;
}

auto EOS::Pressure(const Conservative& uc, double gamma) -> double {
    const double rho = uc.rho;
    if (rho <= 0.0) return 0.0;
    const double u_vel = uc.rhoU / rho;
    const double v_vel = uc.rhoV / rho;
    const double kinetic = 0.5 * rho * (u_vel * u_vel + v_vel * v_vel);
    const double internal = uc.E - kinetic;
    return (gamma - 1.0) * internal;
}

auto EOS::SoundSpeed(const Primitive& w, double gamma) -> double {
    if (w.rho <= 0.0 || w.P <= 0.0) return 0.0;
    return std::sqrt(gamma * w.P / w.rho);
}
