#include "data/Variables.hpp"

#include <cmath>

// ---------- Primitive ----------

Primitive::Primitive() : rho(0.0), u(0.0), P(0.0) {}

Primitive::Primitive(double rho_, double u_, double P_) : rho(rho_), u(u_), P(P_) {}

// ---------- Conservative ----------

Conservative::Conservative() : rho(0.0), rhoU(0.0), E(0.0) {}

Conservative::Conservative(double rho_, double rhoU_, double E_)
    : rho(rho_), rhoU(rhoU_), E(E_) {}

// ---------- Flux ----------

Flux::Flux() : mass(0.0), momentum(0.0), energy(0.0) {}

Flux::Flux(double m, double mu, double e) : mass(m), momentum(mu), energy(e) {}

auto Flux::Diff(const Flux& Fplus, const Flux& Fminus) -> Flux { return Fplus - Fminus; }

// ---------- Euler flux ----------

auto EulerFlux(const Primitive& state, double gamma) -> Flux {
    const double rho = state.rho;
    const double u = state.u;
    const double P = state.P;

    const double mass_flux = rho * u;
    const double momentum_flux = rho * u * u + P;

    const double E = P / (gamma - 1.0) + 0.5 * rho * u * u;
    const double energy_flux = u * (E + P);

    Flux flux;
    flux.mass = mass_flux;
    flux.momentum = momentum_flux;
    flux.energy = energy_flux;

    return flux;
}

// ---------- Primitive <-> Conservative conversions ----------

auto ToConservative(const Primitive& w, double gamma) -> Conservative {
    const double rho = w.rho;
    const double u = w.u;
    const double P = w.P;

    const double rhoU = rho * u;
    const double E = P / (gamma - 1.0) + 0.5 * rho * u * u;

    return {rho, rhoU, E};
}

auto ToPrimitive(const Conservative& U, double gamma) -> Primitive {
    const double rho = U.rho;
    const double rhoU = U.rhoU;
    const double E = U.E;

    if (rho <= 0.0) {
        return {0.0, 0.0, 0.0};
    }

    const double u = rhoU / rho;
    const double kinetic = 0.5 * rho * u * u;
    const double P = (gamma - 1.0) * (E - kinetic);

    return {rho, u, P};
}

// ---------- Small algebra helpers for Conservative ----------

auto operator+=(Conservative& a, const Conservative& b) -> Conservative& {
    a.rho += b.rho;
    a.rhoU += b.rhoU;
    a.E += b.E;
    return a;
}

auto operator-=(Conservative& a, const Conservative& b) -> Conservative& {
    a.rho -= b.rho;
    a.rhoU -= b.rhoU;
    a.E -= b.E;
    return a;
}

auto operator+(Conservative a, const Conservative& b) -> Conservative {
    a += b;
    return a;
}

auto operator-(Conservative a, const Conservative& b) -> Conservative {
    a -= b;
    return a;
}

auto operator*=(Conservative& a, double s) -> Conservative& {
    a.rho *= s;
    a.rhoU *= s;
    a.E *= s;
    return a;
}

auto operator*(Conservative a, double s) -> Conservative {
    a *= s;
    return a;
}

auto operator*(double s, Conservative a) -> Conservative {
    a *= s;
    return a;
}

auto operator-=(Conservative& u, const Flux& f) -> Conservative& {
    u.rho -= f.mass;
    u.rhoU -= f.momentum;
    u.E -= f.energy;
    return u;
}

auto operator-(Conservative u, const Flux& f) -> Conservative {
    u -= f;
    return u;
}

// ---------- Small algebra helpers for Flux ----------

auto operator+=(Flux& a, const Flux& b) -> Flux& {
    a.mass += b.mass;
    a.momentum += b.momentum;
    a.energy += b.energy;
    return a;
}

auto operator-=(Flux& a, const Flux& b) -> Flux& {
    a.mass -= b.mass;
    a.momentum -= b.momentum;
    a.energy -= b.energy;
    return a;
}

auto operator+(Flux a, const Flux& b) -> Flux {
    a += b;
    return a;
}

auto operator-(Flux a, const Flux& b) -> Flux {
    a -= b;
    return a;
}

auto operator*=(Flux& a, double s) -> Flux& {
    a.mass *= s;
    a.momentum *= s;
    a.energy *= s;
    return a;
}

auto operator*(Flux a, double s) -> Flux {
    a *= s;
    return a;
}

auto operator*(double s, Flux a) -> Flux {
    a *= s;
    return a;
}

auto operator+=(Flux& f, const Conservative& u) -> Flux& {
    f.mass += u.rho;
    f.momentum += u.rhoU;
    f.energy += u.E;
    return f;
}

auto operator+(Flux f, const Conservative& u) -> Flux {
    f += u;
    return f;
}

// ---------- Optional algebra for Primitive ----------

auto operator+=(Primitive& a, const Primitive& b) -> Primitive& {
    a.rho += b.rho;
    a.u += b.u;
    a.P += b.P;
    return a;
}

auto operator-=(Primitive& a, const Primitive& b) -> Primitive& {
    a.rho -= b.rho;
    a.u -= b.u;
    a.P -= b.P;
    return a;
}

auto operator+(Primitive a, const Primitive& b) -> Primitive {
    a += b;
    return a;
}

auto operator-(Primitive a, const Primitive& b) -> Primitive {
    a -= b;
    return a;
}

auto operator*=(Primitive& a, double s) -> Primitive& {
    a.rho *= s;
    a.u *= s;
    a.P *= s;
    return a;
}

auto operator*(Primitive a, double s) -> Primitive {
    a *= s;
    return a;
}

auto operator*(double s, Primitive a) -> Primitive {
    a *= s;
    return a;
}
