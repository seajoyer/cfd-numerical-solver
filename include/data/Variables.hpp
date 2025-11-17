#ifndef VARIABLES_HPP
#define VARIABLES_HPP

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
    double rho;  ///< Density
    double u;    ///< Velocity in x-direction
    double P;    ///< Thermodynamic pressure

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
    double rho;   ///< Mass density
    double rhoU;  ///< Momentum density in x-direction
    double E;     ///< Total energy density

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
    double mass;      ///< Mass flux
    double momentum;  ///< Momentum flux
    double energy;    ///< Energy flux

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
    static auto Diff(const Flux& Fplus, const Flux& Fminus) -> Flux;
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
auto EulerFlux(const Primitive& state, double gamma) -> Flux;

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
auto ToConservative(const Primitive& w, double gamma) -> Conservative;

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
auto ToPrimitive(const Conservative& U, double gamma) -> Primitive;

/// --- Small algebra helpers for Conservative ---

auto operator+=(Conservative& a, const Conservative& b) -> Conservative&;
auto operator-=(Conservative& a, const Conservative& b) -> Conservative&;
auto operator+(Conservative a, const Conservative& b) -> Conservative;
auto operator-(Conservative a, const Conservative& b) -> Conservative;
auto operator*=(Conservative& a, double s) -> Conservative&;
auto operator*(Conservative a, double s) -> Conservative;
auto operator*(double s, Conservative a) -> Conservative;
auto operator-=(Conservative& u, const Flux& f) -> Conservative&;
auto operator-(Conservative u, const Flux& f) -> Conservative;

/// --- Small algebra helpers for Flux ---

auto operator+=(Flux& a, const Flux& b) -> Flux&;
auto operator-=(Flux& a, const Flux& b) -> Flux&;
auto operator+(Flux a, const Flux& b) -> Flux;
auto operator-(Flux a, const Flux& b) -> Flux;
auto operator*=(Flux& a, double s) -> Flux&;
auto operator*(Flux a, double s) -> Flux;
auto operator*(double s, Flux a) -> Flux;
auto operator+=(Flux& f, const Conservative& u) -> Flux&;
auto operator+(Flux f, const Conservative& u) -> Flux;

/// --- Small algebra helpers for Primitive (useful in reconstruction) ---

auto operator+=(Primitive& a, const Primitive& b) -> Primitive&;
auto operator-=(Primitive& a, const Primitive& b) -> Primitive&;
auto operator+(Primitive a, const Primitive& b) -> Primitive;
auto operator-(Primitive a, const Primitive& b) -> Primitive;
auto operator*=(Primitive& a, double s) -> Primitive&;
auto operator*(Primitive a, double s) -> Primitive;
auto operator*(double s, Primitive a) -> Primitive;

#endif  // VARIABLES_HPP
