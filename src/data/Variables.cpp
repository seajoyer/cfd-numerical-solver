#include "data/Variables.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

ConservativeCell& ConservativeCell::operator+=(const ConservativeCell& other) {
    rho += other.rho;
    rhoU += other.rhoU;
    rhoV += other.rhoV;
    rhoW += other.rhoW;
    E += other.E;
    return *this;
}

ConservativeCell& ConservativeCell::operator-=(const ConservativeCell& other) {
    rho -= other.rho;
    rhoU -= other.rhoU;
    rhoV -= other.rhoV;
    rhoW -= other.rhoW;
    E -= other.E;
    return *this;
}

ConservativeCell operator+(ConservativeCell lhs, const ConservativeCell& rhs) {
    lhs += rhs;
    return lhs;
}

ConservativeCell operator-(ConservativeCell lhs, const ConservativeCell& rhs) {
    lhs -= rhs;
    return lhs;
}

ConservativeCell operator*(const double scalar, ConservativeCell value) {
    value.rho *= scalar;
    value.rhoU *= scalar;
    value.rhoV *= scalar;
    value.rhoW *= scalar;
    value.E *= scalar;
    return value;
}

ConservativeCell operator*(ConservativeCell value, const double scalar) {
    return scalar * value;
}

ConservativeCell ConservativeFromPrimitive(const PrimitiveCell& w, const double gamma) {
    ConservativeCell U;

    U.rho = w.rho;
    U.rhoU = w.rho * w.u;
    U.rhoV = w.rho * w.v;
    U.rhoW = w.rho * w.w;
    U.E = TotalEnergyDensity(w, gamma);

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

    const double kinetic = 0.5 * rho * (u * u + v * v + ww * ww);
    const double internal_energy_density = U.E - kinetic;
    const double P_raw = (gamma - 1.0) * internal_energy_density;
    const double P = (P_raw > p_floor) ? P_raw : p_floor;

    w.rho = rho;
    w.u = u;
    w.v = v;
    w.w = ww;
    w.P = P;

    return w;
}

void PrimitiveToConservative(const PrimitiveCell& w,
                             const double gamma,
                             double& rho,
                             double& rhoU,
                             double& rhoV,
                             double& rhoW,
                             double& E) {
    rho = w.rho;
    rhoU = w.rho * w.u;
    rhoV = w.rho * w.v;
    rhoW = w.rho * w.w;
    E = TotalEnergyDensity(w, gamma);
}

double SoundSpeed(const PrimitiveCell& w,
                  const double gamma,
                  const double p_floor,
                  const double rho_floor) {
    const double P = (w.P > p_floor) ? w.P : p_floor;
    const double rho = (w.rho > rho_floor) ? w.rho : rho_floor;
    const double a2 = gamma * P / rho;
    return std::sqrt(std::max(a2, 0.0));
}

double KineticEnergyDensity(const PrimitiveCell& w) {
    return 0.5 * w.rho * (w.u * w.u + w.v * w.v + w.w * w.w);
}

double TotalEnergyDensity(const PrimitiveCell& w, const double gamma) {
    return w.P / (gamma - 1.0) + KineticEnergyDensity(w);
}

bool IsUnitNormal(const FaceNormal& normal, const double tolerance) {
    const double norm =
        std::sqrt(normal.x * normal.x +
                  normal.y * normal.y +
                  normal.z * normal.z);

    return std::abs(norm - 1.0) <= tolerance;
}

double NormalVelocity(const PrimitiveCell& w, const FaceNormal& normal) {
    return w.u * normal.x + w.v * normal.y + w.w * normal.z;
}

ConservativeCell PhysicalFlux(const PrimitiveCell& w,
                              const double gamma,
                              const FaceNormal& normal) {
    const double vn = NormalVelocity(w, normal);
    const double E = TotalEnergyDensity(w, gamma);

    ConservativeCell flux;
    flux.rho = w.rho * vn;
    flux.rhoU = w.rho * w.u * vn + w.P * normal.x;
    flux.rhoV = w.rho * w.v * vn + w.P * normal.y;
    flux.rhoW = w.rho * w.w * vn + w.P * normal.z;
    flux.E = (E + w.P) * vn;

    return flux;
}

FluxCell EulerFlux(const PrimitiveCell& w,
                   const double gamma,
                   const FaceNormal& normal) {
    const ConservativeCell flux = PhysicalFlux(w, gamma, normal);

    FluxCell result;
    result.mass = flux.rho;
    result.mom_x = flux.rhoU;
    result.mom_y = flux.rhoV;
    result.mom_z = flux.rhoW;
    result.energy = flux.E;

    return result;
}
