#ifndef VARIABLES_HPP
#define VARIABLES_HPP

/**
 * @class Primitive
 * @brief Primitive (physical) state for the 1D/2D Euler equations.
 *
 *  - rho: density
 *  - u:   velocity in x-direction
 *  - v:   velocity in y-direction (used only in 2D; 0 in 1D)
 *  - P:   thermodynamic pressure
 */
struct Primitive {
    double rho;  ///< Density
    double u;    ///< Velocity in x-direction
    double v;    ///< Velocity in y-direction (0 in 1D)
    double P;    ///< Thermodynamic pressure

    Primitive();
    Primitive(double rho_, double u_, double P_);
    Primitive(double rho_, double u_, double v_, double P_);
};

/**
 * @class Flux
 * @brief Numerical flux for the 1D/2D Euler equations.
 *
 *  - mass:       mass flux
 *  - momentum_x: x-momentum flux
 *  - momentum_y: y-momentum flux (transverse; 0 in 1D)
 *  - energy:     energy flux
 */
struct Flux {
    double mass;        ///< Mass flux
    double momentum_x;  ///< x-Momentum flux (alias: momentum)
    double momentum_y;  ///< y-Momentum flux (transverse advection in x-sweep)
    double energy;      ///< Energy flux

    // Legacy alias for backward compatibility with 1D code
    double& momentum;   ///< Reference to momentum_x for 1D compatibility

    Flux();
    Flux(double m, double mu, double e);
    Flux(double m, double mu_x, double mu_y, double e);

    Flux(const Flux& other);
    auto operator=(const Flux& other) -> Flux&;

    static auto Diff(const Flux& Fplus, const Flux& Fminus) -> Flux;
};

/**
 * @class Conservative
 * @brief Conservative state for the 1D/2D Euler equations.
 *
 *  - rho:  mass density
 *  - rhoU: x-momentum density
 *  - rhoV: y-momentum density (0 in 1D)
 *  - E:    total energy density
 */
struct Conservative {
    double rho;   ///< Mass density
    double rhoU;  ///< x-Momentum density
    double rhoV;  ///< y-Momentum density (0 in 1D)
    double E;     ///< Total energy density

    Conservative();
    Conservative(double rho_, double rhoU_, double E_);
    Conservative(double rho_, double rhoU_, double rhoV_, double E_);
    auto operator=(const Flux& f) -> Conservative&;
};

/**
 * @brief Computes the physical Euler flux in the x-direction for a primitive state.
 *
 * For the 2D Euler equations, the x-direction flux is:
 *  - mass flux       = rho * u
 *  - x-momentum flux = rho * u^2 + P
 *  - y-momentum flux = rho * u * v
 *  - energy flux     = u * (E + P)
 *
 * where E = P / (gamma - 1) + 0.5 * rho * (u^2 + v^2)
 */
auto EulerFlux(const Primitive& state, double gamma) -> Flux;

/**
 * @brief Computes the physical Euler flux in the y-direction for a primitive state.
 *
 * For the 2D Euler equations, the y-direction flux is:
 *  - mass flux       = rho * v
 *  - x-momentum flux = rho * u * v
 *  - y-momentum flux = rho * v^2 + P
 *  - energy flux     = v * (E + P)
 */
auto EulerFluxY(const Primitive& state, double gamma) -> Flux;

auto ToConservative(const Primitive& w, double gamma) -> Conservative;
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
auto operator+=(Conservative& u, const Flux& f) -> Conservative&;
auto operator+(Conservative u, const Flux& f) -> Conservative;

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
auto operator-(Flux f, const Conservative& u) -> Flux;

/// --- Small algebra helpers for Primitive ---

auto operator+=(Primitive& a, const Primitive& b) -> Primitive&;
auto operator-=(Primitive& a, const Primitive& b) -> Primitive&;
auto operator+(Primitive a, const Primitive& b) -> Primitive;
auto operator-(Primitive a, const Primitive& b) -> Primitive;
auto operator*=(Primitive& a, double s) -> Primitive&;
auto operator*(Primitive a, double s) -> Primitive;
auto operator*(double s, Primitive a) -> Primitive;

#endif  // VARIABLES_HPP
