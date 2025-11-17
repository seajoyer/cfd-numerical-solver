#ifndef VARIABLES_HPP
#define VARIABLES_HPP

#include <cstddef>

/**
 * @class Primitive
 * @brief Primitive (physical) state for the 1D Euler equations.
 *
 * This structure stores the standard primitive variables used in
 * reconstruction and Riemann solvers:
 *
 *  - rho: density
 *  - u:   velocity in x-direction
 *  - P:   thermodynamic pressure
 */
struct Primitive {
    double rho; ///< Density
    double u; ///< Velocity in x-direction
    double P; ///< Thermodynamic pressure

    /**
     * @brief Constructs a zero-initialized primitive state.
     */
    Primitive();

    /**
     * @brief Constructs a primitive state from components.
     *
     * @param rho_ Density.
     * @param u_   Velocity.
     * @param P_   Thermodynamic pressure.
     */
    Primitive(double rho_, double u_, double P_);
};

/**
 * @class Conservative
 * @brief Conservative state for the 1D Euler equations.
 *
 * This structure represents the conserved quantities:
 *
 *  - rho:  mass density
 *  - rhoU: momentum density in x-direction
 *  - E:    total energy density
 */
struct Conservative {
    double rho; ///< Mass density
    double rhoU; ///< Momentum density in x-direction
    double E; ///< Total energy density

    /**
     * @brief Constructs a zero-initialized conservative state.
     */
    Conservative();

    /**
     * @brief Constructs a conservative state from components.
     *
     * @param rho_  Mass density.
     * @param rhoU_ Momentum density.
     * @param E_    Total energy density.
     */
    Conservative(double rho_, double rhoU_, double E_);
};

/**
 * @class Flux
 * @brief Numerical flux for the 1D Euler equations.
 *
 * This structure stores the fluxes of the conservative variables
 * across a cell interface:
 *
 *  - mass:     mass flux
 *  - momentum: momentum flux
 *  - energy:   energy flux
 */
struct Flux {
    double mass; ///< Mass flux
    double momentum; ///< Momentum flux
    double energy; ///< Energy flux

    /**
     * @brief Constructs a zero-initialized flux.
     */
    Flux();

    /**
     * @brief Constructs a flux from components.
     *
     * @param m  Mass flux.
     * @param mu Momentum flux.
     * @param e  Energy flux.
     */
    Flux(double m, double mu, double e);


    /**
     * @brief Component-wise difference of two fluxes.
     *
     * Computes Fplus - Fminus for each component.
     *
     * @param Fplus Upwind / right flux.
     * @param Fminus Downwind / left flux.
     * @return Component-wise difference Fplus - Fminus.
     */
    static Flux Diff(const Flux& Fplus, const Flux& Fminus);
};

/**
 * @brief Computes the physical Euler flux for a primitive state.
 *
 * Given a primitive state (rho, u, P), this function computes the
 * corresponding Euler flux for an ideal gas:
 *
 *  - mass flux     = rho * u
 *  - momentum flux = rho * u^2 + P
 *  - energy flux   = u * (E + P)
 *
 * where
 *
 *  E = P / (gamma - 1) + 0.5 * rho * u^2
 *
 * @param state Primitive state at the interface.
 * @param gamma Ratio of specific heats.
 * @return Flux vector (mass, momentum, energy).
 */
Flux EulerFlux(const Primitive& state, double gamma);

/**
 * @brief Converts primitive variables to conservative variables.
 *
 * Uses the ideal-gas relation
 *
 *  E = P / (gamma - 1) + 0.5 * rho * u^2
 *
 * @param w     Primitive state (rho, u, P).
 * @param gamma Ratio of specific heats.
 * @return Conservative state (rho, rhoU, E).
 */
Conservative ToConservative(const Primitive& w, double gamma);

/**
 * @brief Converts conservative variables to primitive variables.
 *
 * Uses the ideal-gas relations
 *
 *  u = rhoU / rho
 *  P = (gamma - 1) * (E - 0.5 * rho * u^2)
 *
 * If rho <= 0, the function returns a zero state (rho = u = P = 0),
 * which should normally be prevented by positivity limiters.
 *
 * @param U     Conservative state (rho, rhoU, E).
 * @param gamma Ratio of specific heats.
 * @return Primitive state (rho, u, P).
 */
Primitive ToPrimitive(const Conservative& U, double gamma);

/// --- Small algebra helpers for Conservative ---

Conservative& operator+=(Conservative& a, const Conservative& b);
Conservative& operator-=(Conservative& a, const Conservative& b);
Conservative operator+(Conservative a, const Conservative& b);
Conservative operator-(Conservative a, const Conservative& b);
Conservative& operator*=(Conservative& a, double s);
Conservative operator*(Conservative a, double s);
Conservative operator*(double s, Conservative a);
Conservative& operator-=(Conservative& u, const Flux& f);
Conservative operator-(Conservative u, const Flux& f);

/// --- Small algebra helpers for Flux ---

Flux& operator+=(Flux& a, const Flux& b);
Flux& operator-=(Flux& a, const Flux& b);
Flux operator+(Flux a, const Flux& b);
Flux operator-(Flux a, const Flux& b);
Flux& operator*=(Flux& a, double s);
Flux operator*(Flux a, double s);
Flux operator*(double s, Flux a);
Flux& operator+=(Flux& f, const Conservative& u);
Flux operator+(Flux f, const Conservative& u);

/// --- Small algebra helpers for Primitive (useful in reconstruction) ---

Primitive& operator+=(Primitive& a, const Primitive& b);
Primitive& operator-=(Primitive& a, const Primitive& b);
Primitive operator+(Primitive a, const Primitive& b);
Primitive operator-(Primitive a, const Primitive& b);
Primitive& operator*=(Primitive& a, double s);
Primitive operator*(Primitive a, double s);
Primitive operator*(double s, Primitive a);

#endif  // VARIABLES_HPP