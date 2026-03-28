#ifndef VARIABLES_HPP
#define VARIABLES_HPP

#include <cstddef>

/**
 * @file Variables.hpp
 * @brief Minimal variable utilities for Euler equations (ideal gas).
 *
 * Conventions:
 *  - Conservative state: (rho, rhoU, rhoV, rhoW, E)
 *  - Primitive state:    (rho, u, v, w, P)
 *
 * This module provides:
 *  - Cell-level conversion between conservative and primitive variables.
 *  - Ideal-gas sound speed.
 *  - Euler flux projected onto an arbitrary face normal.
 */

namespace var {
    // Conservative / primitive variable ordering
    static constexpr std::size_t rho = 0;
    static constexpr std::size_t u1 = 1;
    static constexpr std::size_t u2 = 2;
    static constexpr std::size_t u3 = 3;
    static constexpr std::size_t p_or_E = 4;

    static constexpr std::size_t nvar = 5;
} // namespace var

/**
 * @brief Primitive state in one cell.
 */
struct PrimitiveCell final {
    double rho = 0.0;
    double u = 0.0;
    double v = 0.0;
    double w = 0.0;
    double P = 0.0;
};

/**
 * @brief Conservative state in one cell.
 * @details Order matches U = (rho, rhoU, rhoV, rhoW, E).
 */
struct ConservativeCell final {
    double rho = 0.0;
    double rhoU = 0.0;
    double rhoV = 0.0;
    double rhoW = 0.0;
    double E = 0.0;

    ConservativeCell& operator+=(const ConservativeCell& other);
    ConservativeCell& operator-=(const ConservativeCell& other);
};

[[nodiscard]] ConservativeCell operator+(ConservativeCell lhs, const ConservativeCell& rhs);
[[nodiscard]] ConservativeCell operator-(ConservativeCell lhs, const ConservativeCell& rhs);
[[nodiscard]] ConservativeCell operator*(double scalar, ConservativeCell value);
[[nodiscard]] ConservativeCell operator*(ConservativeCell value, double scalar);

/**
 * @brief Euler flux vector (5 components) in Cartesian form.
 * @details Order matches conservative variables: (mass, mom_x, mom_y, mom_z, energy).
 */
struct FluxCell final {
    double mass = 0.0;
    double mom_x = 0.0;
    double mom_y = 0.0;
    double mom_z = 0.0;
    double energy = 0.0;
};

/**
 * @brief One unit face normal.
 */
struct FaceNormal final {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

/**
 * @brief Convert primitive state to conservative state (ideal gas).
 */
[[nodiscard]] ConservativeCell ConservativeFromPrimitive(const PrimitiveCell& w, double gamma);

/**
 * @brief Convert conservative state to primitive state (ideal gas).
 */
[[nodiscard]] PrimitiveCell PrimitiveFromConservativeCell(const ConservativeCell& U,
                                                          double gamma,
                                                          double rho_floor = 1e-14,
                                                          double p_floor = 1e-14);

/**
 * @brief Convert primitive state to conservative components.
 */
void PrimitiveToConservative(const PrimitiveCell& w,
                             double gamma,
                             double& rho,
                             double& rhoU,
                             double& rhoV,
                             double& rhoW,
                             double& E);

/**
 * @brief Speed of sound for ideal gas from primitive state.
 */
[[nodiscard]] double SoundSpeed(const PrimitiveCell& w,
                                double gamma,
                                double p_floor = 1e-14,
                                double rho_floor = 1e-14);

/**
 * @brief Total specific kinetic energy density.
 * @details Returns 0.5 * rho * (u^2 + v^2 + w^2).
 */
[[nodiscard]] double KineticEnergyDensity(const PrimitiveCell& w);

/**
 * @brief Total energy density from primitive state.
 */
[[nodiscard]] double TotalEnergyDensity(const PrimitiveCell& w, double gamma);

/**
 * @brief Check whether face normal is unit-length within tolerance.
 */
[[nodiscard]] bool IsUnitNormal(const FaceNormal& normal, double tolerance = 1e-10);

/**
 * @brief Return velocity component projected on the given face normal.
 */
[[nodiscard]] double NormalVelocity(const PrimitiveCell& w, const FaceNormal& normal);

/**
 * @brief Cartesian Euler flux tensor contracted with the given face normal.
 * @details Returns physical flux F(U) · n in conservative-variable ordering.
 */
[[nodiscard]] ConservativeCell PhysicalFlux(const PrimitiveCell& w,
                                            double gamma,
                                            const FaceNormal& normal);

/**
 * @brief Same as PhysicalFlux, returned in FluxCell type.
 */
[[nodiscard]] FluxCell EulerFlux(const PrimitiveCell& w,
                                 double gamma,
                                 const FaceNormal& normal);

#endif  // VARIABLES_HPP
