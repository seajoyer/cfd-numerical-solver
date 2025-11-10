#ifndef VARIABLES_HPP
#define VARIABLES_HPP

#include <cstddef>

/**
 * @struct Primitive
 * @brief Primitive (physical) variables for 1D Euler equations.
 *
 * This structure stores the standard primitive state used in
 * reconstruction and Riemann solvers:
 * - rho: density
 * - u:   velocity (x-component)
 * - P:   thermodynamic pressure
 */
struct Primitive {
    double rho;  ///< Density
    double u;    ///< Velocity in x-direction
    double P;    ///< Thermodynamic pressure
};

/**
 * @struct Conservative
 * @brief Conservative variables for 1D Euler equations.
 *
 * This structure represents the conserved quantities:
 * - rho:  mass density
 * - rhoU: momentum density in x-direction
 * - E:    total energy density
 */
struct Conservative {
    double rho;   ///< Mass density
    double rhoU;  ///< Momentum density in x-direction
    double E;     ///< Total energy density
};

/**
 * @struct Flux
 * @brief Numerical flux vector for 1D Euler equations.
 *
 * Represents the fluxes of the conservative variables across
 * a cell interface:
 * - mass:    mass flux
 * - momentum: momentum flux
 * - energy:   energy flux
 */
struct Flux {
    double mass;     ///< Mass flux
    double momentum; ///< Momentum flux
    double energy;   ///< Energy flux
};

/**
 * @brief Computes the physical Euler flux for given primitive state.
 *
 * This function returns the ideal-gas Euler flux vector corresponding
 * to the provided primitive variables. It is typically used inside
 * Riemann solvers and for consistency checks.
 *
 * @param state Primitive state (rho, u, P) at the interface.
 * @param gamma Ratio of specific heats.
 * @return Flux vector (mass, momentum, energy).
 */
Flux EulerFlux(const Primitive &state, double gamma);

#endif  // VARIABLES_HPP
