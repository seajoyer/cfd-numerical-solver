#include "data/Variables.hpp"

// ---------- Primitive ----------

Primitive::Primitive() : rho(0.0), u(0.0), v(0.0), P(0.0) {}

Primitive::Primitive(double rho_, double u_, double P_)
    : rho(rho_), u(u_), v(0.0), P(P_) {}

Primitive::Primitive(double rho_, double u_, double v_, double P_)
    : rho(rho_), u(u_), v(v_), P(P_) {}

// ---------- Conservative ----------

Conservative::Conservative() : rho(0.0), rhoU(0.0), rhoV(0.0), E(0.0) {}

Conservative::Conservative(double rho_, double rhoU_, double E_)
    : rho(rho_), rhoU(rhoU_), rhoV(0.0), E(E_) {}

Conservative::Conservative(double rho_, double rhoU_, double rhoV_, double E_)
    : rho(rho_), rhoU(rhoU_), rhoV(rhoV_), E(E_) {}

auto Conservative::operator=(const Flux& f) -> Conservative& {
    this->rho = f.mass;
    this->rhoU = f.momentum_x;
    this->rhoV = f.momentum_y;
    this->E = f.energy;
    return *this;
}

// ---------- Flux ----------

Flux::Flux() : mass(0.0), momentum_x(0.0), momentum_y(0.0), energy(0.0),
               momentum(momentum_x) {}

Flux::Flux(double m, double mu, double e)
    : mass(m), momentum_x(mu), momentum_y(0.0), energy(e),
      momentum(momentum_x) {}

Flux::Flux(double m, double mu_x, double mu_y, double e)
    : mass(m), momentum_x(mu_x), momentum_y(mu_y), energy(e),
      momentum(momentum_x) {}

Flux::Flux(const Flux& other)
    : mass(other.mass), momentum_x(other.momentum_x),
      momentum_y(other.momentum_y), energy(other.energy),
      momentum(momentum_x) {}

auto Flux::operator=(const Flux& other) -> Flux& {
    if (this != &other) {
        mass = other.mass;
        momentum_x = other.momentum_x;
        momentum_y = other.momentum_y;
        energy = other.energy;
    }
    return *this;
}

auto Flux::Diff(const Flux& Fplus, const Flux& Fminus) -> Flux {
    return Fplus - Fminus;
}

// ---------- Euler flux (x-direction) ----------

auto EulerFlux(const Primitive& state, double gamma) -> Flux {
    const double rho = state.rho;
    const double u = state.u;
    const double v = state.v;
    const double P = state.P;

    const double mass_flux = rho * u;
    const double momentum_x_flux = rho * u * u + P;
    const double momentum_y_flux = rho * u * v;

    const double E = P / (gamma - 1.0) + 0.5 * rho * (u * u + v * v);
    const double energy_flux = u * (E + P);

    return {mass_flux, momentum_x_flux, momentum_y_flux, energy_flux};
}

// ---------- Euler flux (y-direction) ----------

auto EulerFluxY(const Primitive& state, double gamma) -> Flux {
    const double rho = state.rho;
    const double u = state.u;
    const double v = state.v;
    const double P = state.P;

    const double mass_flux = rho * v;
    const double momentum_x_flux = rho * u * v;
    const double momentum_y_flux = rho * v * v + P;

    const double E = P / (gamma - 1.0) + 0.5 * rho * (u * u + v * v);
    const double energy_flux = v * (E + P);

    return {mass_flux, momentum_x_flux, momentum_y_flux, energy_flux};
}

// ---------- Primitive <-> Conservative conversions ----------

auto ToConservative(const Primitive& w, double gamma) -> Conservative {
    const double rho = w.rho;
    const double u = w.u;
    const double v = w.v;
    const double P = w.P;

    const double rhoU = rho * u;
    const double rhoV = rho * v;
    const double E = P / (gamma - 1.0) + 0.5 * rho * (u * u + v * v);

    return {rho, rhoU, rhoV, E};
}

auto ToPrimitive(const Conservative& U, double gamma) -> Primitive {
    const double rho = U.rho;
    const double rhoU = U.rhoU;
    const double rhoV = U.rhoV;
    const double E = U.E;

    if (rho <= 0.0) {
        return {0.0, 0.0, 0.0, 0.0};
    }

    const double u = rhoU / rho;
    const double v = rhoV / rho;
    const double kinetic = 0.5 * rho * (u * u + v * v);
    const double P = (gamma - 1.0) * (E - kinetic);

    return {rho, u, v, P};
}

// ---------- Small algebra helpers for Conservative ----------

auto operator+=(Conservative& a, const Conservative& b) -> Conservative& {
    a.rho += b.rho;
    a.rhoU += b.rhoU;
    a.rhoV += b.rhoV;
    a.E += b.E;
    return a;
}

auto operator-=(Conservative& a, const Conservative& b) -> Conservative& {
    a.rho -= b.rho;
    a.rhoU -= b.rhoU;
    a.rhoV -= b.rhoV;
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
    a.rhoV *= s;
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
    u.rhoU -= f.momentum_x;
    u.rhoV -= f.momentum_y;
    u.E -= f.energy;
    return u;
}

auto operator-(Conservative u, const Flux& f) -> Conservative {
    u -= f;
    return u;
}

auto operator+=(Conservative& u, const Flux& f) -> Conservative& {
    u.rho += f.mass;
    u.rhoU += f.momentum_x;
    u.rhoV += f.momentum_y;
    u.E += f.energy;
    return u;
}

auto operator+(Conservative u, const Flux& f) -> Conservative {
    u += f;
    return u;
}

// ---------- Small algebra helpers for Flux ----------

auto operator+=(Flux& a, const Flux& b) -> Flux& {
    a.mass += b.mass;
    a.momentum_x += b.momentum_x;
    a.momentum_y += b.momentum_y;
    a.energy += b.energy;
    return a;
}

auto operator-=(Flux& a, const Flux& b) -> Flux& {
    a.mass -= b.mass;
    a.momentum_x -= b.momentum_x;
    a.momentum_y -= b.momentum_y;
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
    a.momentum_x *= s;
    a.momentum_y *= s;
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
    f.momentum_x += u.rhoU;
    f.momentum_y += u.rhoV;
    f.energy += u.E;
    return f;
}

auto operator+(Flux f, const Conservative& u) -> Flux {
    f += u;
    return f;
}

auto operator-(Flux f, const Conservative& u) -> Flux {
    f.mass -= u.rho;
    f.momentum_x -= u.rhoU;
    f.momentum_y -= u.rhoV;
    f.energy -= u.E;
    return f;
}

// ---------- Optional algebra for Primitive ----------

auto operator+=(Primitive& a, const Primitive& b) -> Primitive& {
    a.rho += b.rho;
    a.u += b.u;
    a.v += b.v;
    a.P += b.P;
    return a;
}

auto operator-=(Primitive& a, const Primitive& b) -> Primitive& {
    a.rho -= b.rho;
    a.u -= b.u;
    a.v -= b.v;
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
    a.v *= s;
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
