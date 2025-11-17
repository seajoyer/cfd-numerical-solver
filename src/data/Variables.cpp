#include "data/Variables.hpp"

#include <cmath>
#include <ostream>

// ---------- Primitive ----------

Primitive::Primitive()
    : rho(0.0), u(0.0), P(0.0) {
}

Primitive::Primitive(double rho_, double u_, double P_)
    : rho(rho_), u(u_), P(P_) {
}

// ---------- Conservative ----------

Conservative::Conservative()
    : rho(0.0), rhoU(0.0), E(0.0) {
}

Conservative::Conservative(double rho_, double rhoU_, double E_)
    : rho(rho_), rhoU(rhoU_), E(E_) {
}

// ---------- Flux ----------

Flux::Flux()
    : mass(0.0), momentum(0.0), energy(0.0) {
}

Flux::Flux(double m, double mu, double e)
    : mass(m), momentum(mu), energy(e) {
}

Flux Flux::Diff(const Flux& Fplus, const Flux& Fminus) {
    return Fplus - Fminus;
}

// ---------- Euler flux ----------

Flux EulerFlux(const Primitive& state, double gamma) {
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

Conservative ToConservative(const Primitive& w, double gamma) {
    const double rho = w.rho;
    const double u = w.u;
    const double P = w.P;

    const double rhoU = rho * u;
    const double E = P / (gamma - 1.0) + 0.5 * rho * u * u;

    return Conservative(rho, rhoU, E);
}

Primitive ToPrimitive(const Conservative& U, double gamma) {
    const double rho = U.rho;
    const double rhoU = U.rhoU;
    const double E = U.E;

    if (rho <= 0.0) {
        return Primitive(0.0, 0.0, 0.0);
    }

    const double u = rhoU / rho;
    const double kinetic = 0.5 * rho * u * u;
    const double P = (gamma - 1.0) * (E - kinetic);

    return Primitive(rho, u, P);
}

// ---------- Small algebra helpers for Conservative ----------

Conservative& operator+=(Conservative& a, const Conservative& b) {
    a.rho += b.rho;
    a.rhoU += b.rhoU;
    a.E += b.E;
    return a;
}

Conservative& operator-=(Conservative& a, const Conservative& b) {
    a.rho -= b.rho;
    a.rhoU -= b.rhoU;
    a.E -= b.E;
    return a;
}

Conservative operator+(Conservative a, const Conservative& b) {
    a += b;
    return a;
}

Conservative operator-(Conservative a, const Conservative& b) {
    a -= b;
    return a;
}

Conservative& operator*=(Conservative& a, double s) {
    a.rho *= s;
    a.rhoU *= s;
    a.E *= s;
    return a;
}

Conservative operator*(Conservative a, double s) {
    a *= s;
    return a;
}

Conservative operator*(double s, Conservative a) {
    a *= s;
    return a;
}

Conservative& operator-=(Conservative& u, const Flux& f) {
    u.rho -= f.mass;
    u.rhoU -= f.momentum;
    u.E -= f.energy;
    return u;
}

Conservative operator-(Conservative u, const Flux& f) {
    u -= f;
    return u;
}

// ---------- Small algebra helpers for Flux ----------

Flux& operator+=(Flux& a, const Flux& b) {
    a.mass += b.mass;
    a.momentum += b.momentum;
    a.energy += b.energy;
    return a;
}

Flux& operator-=(Flux& a, const Flux& b) {
    a.mass -= b.mass;
    a.momentum -= b.momentum;
    a.energy -= b.energy;
    return a;
}

Flux operator+(Flux a, const Flux& b) {
    a += b;
    return a;
}

Flux operator-(Flux a, const Flux& b) {
    a -= b;
    return a;
}

Flux& operator*=(Flux& a, double s) {
    a.mass *= s;
    a.momentum *= s;
    a.energy *= s;
    return a;
}

Flux operator*(Flux a, double s) {
    a *= s;
    return a;
}

Flux operator*(double s, Flux a) {
    a *= s;
    return a;
}

Flux& operator+=(Flux& f, const Conservative& u) {
    f.mass += u.rho;
    f.momentum += u.rhoU;
    f.energy += u.E;
    return f;
}

Flux operator+(Flux f, const Conservative& u) {
    f += u;
    return f;
}

// ---------- Optional algebra for Primitive ----------

Primitive& operator+=(Primitive& a, const Primitive& b) {
    a.rho += b.rho;
    a.u += b.u;
    a.P += b.P;
    return a;
}

Primitive& operator-=(Primitive& a, const Primitive& b) {
    a.rho -= b.rho;
    a.u -= b.u;
    a.P -= b.P;
    return a;
}

Primitive operator+(Primitive a, const Primitive& b) {
    a += b;
    return a;
}

Primitive operator-(Primitive a, const Primitive& b) {
    a -= b;
    return a;
}

Primitive& operator*=(Primitive& a, double s) {
    a.rho *= s;
    a.u *= s;
    a.P *= s;
    return a;
}

Primitive operator*(Primitive a, double s) {
    a *= s;
    return a;
}

Primitive operator*(double s, Primitive a) {
    a *= s;
    return a;
}