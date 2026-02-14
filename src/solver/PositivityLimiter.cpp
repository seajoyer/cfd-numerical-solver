#include "solver/PositivityLimiter.hpp"

#include <algorithm>
#include <cmath>

void PositivityLimiter::Apply(Conservative& U, double gamma,
                              double rho_min, double p_min) {
    // Enforce minimum density
    if (U.rho < rho_min) {
        U.rho = rho_min;
    }

    // Recompute velocity from (potentially clamped) density
    const double u = U.rhoU / U.rho;
    const double v = U.rhoV / U.rho;

    // Kinetic energy (includes both u and v components for 2D)
    const double kinetic = 0.5 * U.rho * (u * u + v * v);

    // Internal energy
    const double e_int = U.E - kinetic;

    // Minimum internal energy from minimum pressure
    const double e_int_min = p_min / (gamma - 1.0);

    // Enforce minimum pressure
    if (e_int < e_int_min) {
        U.E = kinetic + e_int_min;
    }
}
