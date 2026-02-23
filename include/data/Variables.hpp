#ifndef VARIABLES_HPP
#define VARIABLES_HPP

#include <cstddef>
#include <cstdint>
#include <xtensor.hpp>

/**
 * @file Variables.hpp
 * @brief Minimal variable utilities for Euler equations (ideal gas), optimized for clarity.
 *
 * Conventions:
 *  - Conservative field: U(var, i, j, k), var = (rho, rhoU, rhoV, rhoW, E)
 *  - Primitive field:    W(var, i, j, k), var = (rho, u, v, w, P)
 *
 * This module provides:
 *  - Cell-level conversion U -> W (ideal gas).
 *  - Cell-level Euler flux from primitive variables.
 *  - Field conversion U -> W over a region (no allocations inside).
 */

enum class Axis : std::uint8_t { X = 0, Y = 1, Z = 2 };

/**
 * @brief Integer stride (di,dj,dk) for stepping along a given axis.
 *
 * Axis::X -> (1,0,0), Axis::Y -> (0,1,0), Axis::Z -> (0,0,1)
 */
struct AxisStride final {
    int di = 0;
    int dj = 0;
    int dk = 0;

    [[nodiscard]] static constexpr AxisStride FromAxis(Axis axis) {
        return axis == Axis::X
                   ? AxisStride{.di = 1, .dj = 0, .dk = 0}
                   : axis == Axis::Y
                   ? AxisStride{.di = 0, .dj = 1, .dk = 0}
                   : AxisStride{.di = 0, .dj = 0, .dk = 1};
    }
};

namespace var {
    // Conservative indices in U(var, i, j, k)
    static constexpr std::size_t rho = 0;
    static constexpr std::size_t rhoU = 1;
    static constexpr std::size_t rhoV = 2;
    static constexpr std::size_t rhoW = 3;
    static constexpr std::size_t E = 4;

    // Primitive indices in W(var, i, j, k)
    static constexpr std::size_t u_rho = 0;
    static constexpr std::size_t u_u = 1;
    static constexpr std::size_t u_v = 2;
    static constexpr std::size_t u_w = 3;
    static constexpr std::size_t u_P = 4;

    // Number of variables
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
 * @brief Euler flux vector (5 components) at a face/cell.
 * Order matches conservative variables: (mass, mom_x, mom_y, mom_z, energy).
 */
struct FluxCell final {
    double mass = 0.0;
    double mom_x = 0.0;
    double mom_y = 0.0;
    double mom_z = 0.0;
    double energy = 0.0;
};

/**
 * @brief Convert conservative state U at (i,j,k) to primitive state W (ideal gas).
 * @param U Conservative field U(var, i, j, k).
 * @param i,j,k Cell indices.
 * @param gamma Ratio of specific heats.
 * @param rho_floor Minimum density used if rho <= 0 (safety).
 * @param p_floor Minimum pressure used if computed P <= 0 (safety).
 */
[[nodiscard]] PrimitiveCell PrimitiveFromConservative(
    const xt::xtensor<double, 4>& U,
    int i, int j, int k,
    double gamma,
    double rho_floor = 1e-14,
    double p_floor = 1e-14
);

/**
 * @brief Speed of sound for ideal gas from primitive state.
 */
[[nodiscard]] double SoundSpeed(const PrimitiveCell& w, double gamma, double p_floor = 1e-14);

/**
 * @brief Euler flux from primitive state along a given axis (ideal gas).
 */
[[nodiscard]] FluxCell EulerFlux(const PrimitiveCell& w, double gamma, Axis axis);

/**
 * @brief Convert U -> W over a rectangular region [i0,i1) x [j0,j1) x [k0,k1).
 * @details No allocations inside; W must have shape (5, sx, sy, sz).
 */
void ConvertUtoW(
    const xt::xtensor<double, 4>& U,
    xt::xtensor<double, 4>& W,
    double gamma,
    int i0, int i1,
    int j0, int j1,
    int k0, int k1,
    double rho_floor = 1e-14,
    double p_floor = 1e-14
);

/**
 * @brief Normal velocity component for a primitive state along the given axis.
 */
[[nodiscard]] double NormalVelocity(const PrimitiveCell& w, Axis axis);

/**
 * @brief Convert primitive to conservative components (rho, rhoU, rhoV, rhoW, E).
 */
void PrimitiveToConservative(const PrimitiveCell& w,
                             double gamma,
                             double& rho,
                             double& rhoU,
                             double& rhoV,
                             double& rhoW,
                             double& E);

/**
 * @brief Conservative state in one cell.
 * Order matches U = (rho, rhoU, rhoV, rhoW, E).
 */
struct ConservativeCell final {
    double rho  = 0.0;
    double rhoU = 0.0;
    double rhoV = 0.0;
    double rhoW = 0.0;
    double E    = 0.0;

    // Lightweight arithmetic helpers (no allocations)
    ConservativeCell& operator+=(const ConservativeCell& other);
    ConservativeCell& operator-=(const ConservativeCell& other);
};

[[nodiscard]] ConservativeCell operator+(ConservativeCell a, const ConservativeCell& b);
[[nodiscard]] ConservativeCell operator-(ConservativeCell a, const ConservativeCell& b);
[[nodiscard]] ConservativeCell operator*(double s, ConservativeCell a);
[[nodiscard]] ConservativeCell operator*(ConservativeCell a, double s);

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

#endif  // VARIABLES_HPP
