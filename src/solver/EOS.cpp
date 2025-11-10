#include "solver/EOS.hpp"
#include <cmath>

Conservative EOS::PrimToCons(const Primitive &w, double gamma) {
    Conservative u;
    u.rho = w.rho;
    u.rhoU = w.rho * w.u;
    u.E = w.P / (gamma - 1.0) + 0.5 * w.rho * w.u * w.u;
    return u;
}

Primitive EOS::ConsToPrim(const Conservative &u, double gamma) {
    Primitive w;
    w.rho = u.rho;
    const double inv_rho = 1.0 / u.rho;
    w.u = u.rhoU * inv_rho;
    w.P = Pressure(u, gamma);
    return w;
}

double EOS::Pressure(const Conservative &u, double gamma) {
    const double rho = u.rho;
    const double u_vel = u.rhoU / rho;
    const double kinetic = 0.5 * rho * u_vel * u_vel;
    const double internal = u.E - kinetic;
    return (gamma - 1.0) * internal;
}

double EOS::SoundSpeed(const Primitive &w, double gamma) {
    return std::sqrt(gamma * w.P / w.rho);
}
