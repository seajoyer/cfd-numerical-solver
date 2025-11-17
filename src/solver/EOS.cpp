#include "solver/EOS.hpp"

#include <cmath>

auto EOS::PrimToCons(const Primitive& w, double gamma) -> Conservative {
    Conservative u;
    u.rho = w.rho;
    u.rhoU = w.rho * w.u;
    u.E = w.P / (gamma - 1.0) + 0.5 * w.rho * w.u * w.u;
    return u;
}

auto EOS::ConsToPrim(const Conservative& u, double gamma) -> Primitive {
    Primitive w;
    w.rho = u.rho;
    const double inv_rho = 1.0 / u.rho;
    w.u = u.rhoU * inv_rho;
    w.P = Pressure(u, gamma);
    return w;
}

auto EOS::Pressure(const Conservative& u, double gamma) -> double {
    const double rho = u.rho;
    const double u_vel = u.rhoU / rho;
    const double kinetic = 0.5 * rho * u_vel * u_vel;
    const double internal = u.E - kinetic;
    return (gamma - 1.0) * internal;
}

auto EOS::SoundSpeed(const Primitive& w, double gamma) -> double {
    return std::sqrt(gamma * w.P / w.rho);
}
