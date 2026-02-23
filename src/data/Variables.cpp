#include "data/Variables.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

PrimitiveCell PrimitiveFromConservative(
    const xt::xtensor<double, 4>& U,
    const int i, const int j, const int k,
    const double gamma,
    const double rho_floor,
    const double p_floor
) {
    PrimitiveCell w;

    const double rho_in  = U(var::rho,  i, j, k);
    const double rhoU_in = U(var::rhoU, i, j, k);
    const double rhoV_in = U(var::rhoV, i, j, k);
    const double rhoW_in = U(var::rhoW, i, j, k);
    const double E_in    = U(var::E,    i, j, k);

    const double rho = (rho_in > rho_floor) ? rho_in : rho_floor;
    const double inv_rho = 1.0 / rho;

    const double u = rhoU_in * inv_rho;
    const double v = rhoV_in * inv_rho;
    const double wv = rhoW_in * inv_rho;

    const double kinetic = 0.5 * rho * (u*u + v*v + wv*wv);
    const double eint = E_in - kinetic;
    const double P_raw = (gamma - 1.0) * eint;
    const double P = (P_raw > p_floor) ? P_raw : p_floor;

    w.rho = rho;
    w.u = u;
    w.v = v;
    w.w = wv;
    w.P = P;
    return w;
}

double SoundSpeed(const PrimitiveCell& w, const double gamma, const double p_floor) {
    const double P = (w.P > p_floor) ? w.P : p_floor;
    const double rho = (w.rho > 0.0) ? w.rho : 1e-14;
    const double a2 = gamma * P / rho;
    return std::sqrt(std::max(a2, 0.0));
}

static double TotalEnergyFromPrimitive(const PrimitiveCell& w, const double gamma) {
    const double kinetic = 0.5 * w.rho * (w.u*w.u + w.v*w.v + w.w*w.w);
    return w.P / (gamma - 1.0) + kinetic;
}

FluxCell EulerFlux(const PrimitiveCell& w, const double gamma, const Axis axis) {
    FluxCell f;
    const double E = TotalEnergyFromPrimitive(w, gamma);

    if (axis == Axis::X) {
        f.mass   = w.rho * w.u;
        f.mom_x  = w.rho * w.u * w.u + w.P;
        f.mom_y  = w.rho * w.u * w.v;
        f.mom_z  = w.rho * w.u * w.w;
        f.energy = w.u * (E + w.P);
        return f;
    }

    if (axis == Axis::Y) {
        f.mass   = w.rho * w.v;
        f.mom_x  = w.rho * w.u * w.v;
        f.mom_y  = w.rho * w.v * w.v + w.P;
        f.mom_z  = w.rho * w.v * w.w;
        f.energy = w.v * (E + w.P);
        return f;
    }

    // Axis::Z
    f.mass   = w.rho * w.w;
    f.mom_x  = w.rho * w.u * w.w;
    f.mom_y  = w.rho * w.v * w.w;
    f.mom_z  = w.rho * w.w * w.w + w.P;
    f.energy = w.w * (E + w.P);
    return f;
}

void ConvertUtoW(
    const xt::xtensor<double, 4>& U,
    xt::xtensor<double, 4>& W,
    const double gamma,
    const int i0, const int i1,
    const int j0, const int j1,
    const int k0, const int k1,
    const double rho_floor,
    const double p_floor
) {
    // Basic shape sanity: only check var dimension to avoid overhead.
    if (U.shape().size() != 4 || W.shape().size() != 4) {
        throw std::invalid_argument("ConvertUtoW expects 4D tensors");
    }
    if (U.shape()[0] != var::nvar || W.shape()[0] != var::nvar) {
        throw std::invalid_argument("ConvertUtoW expects first dimension size = 5");
    }

    for (int k = k0; k < k1; ++k) {
        for (int j = j0; j < j1; ++j) {
            for (int i = i0; i < i1; ++i) {
                const PrimitiveCell w = PrimitiveFromConservative(U, i, j, k, gamma, rho_floor, p_floor);
                W(var::u_rho, i, j, k) = w.rho;
                W(var::u_u,   i, j, k) = w.u;
                W(var::u_v,   i, j, k) = w.v;
                W(var::u_w,   i, j, k) = w.w;
                W(var::u_P,   i, j, k) = w.P;
            }
        }
    }
}

double NormalVelocity(const PrimitiveCell& w, const Axis axis) {
    if (axis == Axis::X) return w.u;
    if (axis == Axis::Y) return w.v;
    return w.w;
}

void PrimitiveToConservative(const PrimitiveCell& w,
                             const double gamma,
                             double& rho,
                             double& rhoU,
                             double& rhoV,
                             double& rhoW,
                             double& E) {
    rho  = w.rho;
    rhoU = w.rho * w.u;
    rhoV = w.rho * w.v;
    rhoW = w.rho * w.w;

    const double kinetic = 0.5 * w.rho * (w.u*w.u + w.v*w.v + w.w*w.w);
    E = w.P / (gamma - 1.0) + kinetic;
}

ConservativeCell& ConservativeCell::operator+=(const ConservativeCell& o) {
    rho  += o.rho;
    rhoU += o.rhoU;
    rhoV += o.rhoV;
    rhoW += o.rhoW;
    E    += o.E;
    return *this;
}

ConservativeCell& ConservativeCell::operator-=(const ConservativeCell& o) {
    rho  -= o.rho;
    rhoU -= o.rhoU;
    rhoV -= o.rhoV;
    rhoW -= o.rhoW;
    E    -= o.E;
    return *this;
}

ConservativeCell operator+(ConservativeCell a, const ConservativeCell& b) {
    a += b;
    return a;
}

ConservativeCell operator-(ConservativeCell a, const ConservativeCell& b) {
    a -= b;
    return a;
}

ConservativeCell operator*(const double s, ConservativeCell a) {
    a.rho  *= s;
    a.rhoU *= s;
    a.rhoV *= s;
    a.rhoW *= s;
    a.E    *= s;
    return a;
}

ConservativeCell operator*(ConservativeCell a, const double s) {
    return s * a;
}

ConservativeCell ConservativeFromPrimitive(const PrimitiveCell& w, const double gamma) {
    ConservativeCell U;
    U.rho  = w.rho;
    U.rhoU = w.rho * w.u;
    U.rhoV = w.rho * w.v;
    U.rhoW = w.rho * w.w;

    const double kinetic = 0.5 * w.rho * (w.u*w.u + w.v*w.v + w.w*w.w);
    U.E = w.P / (gamma - 1.0) + kinetic;
    return U;
}

PrimitiveCell PrimitiveFromConservativeCell(const ConservativeCell& U,
                                            const double gamma,
                                            const double rho_floor,
                                            const double p_floor) {
    PrimitiveCell w;

    const double rho = (U.rho > rho_floor) ? U.rho : rho_floor;
    const double inv_rho = 1.0 / rho;

    const double u = U.rhoU * inv_rho;
    const double v = U.rhoV * inv_rho;
    const double ww = U.rhoW * inv_rho;

    const double kinetic = 0.5 * rho * (u*u + v*v + ww*ww);
    const double eint = U.E - kinetic;
    const double P_raw = (gamma - 1.0) * eint;
    const double P = (P_raw > p_floor) ? P_raw : p_floor;

    w.rho = rho;
    w.u = u;
    w.v = v;
    w.w = ww;
    w.P = P;
    return w;
}