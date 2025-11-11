#include "solver/Variables.hpp"

auto EulerFlux(const Primitive &state, double gamma) -> Flux {
    const double rho = state.rho;
    const double u = state.u;
    const double P = state.P;

    const double mass_flux = rho * u;
    const double momentum_flux = rho * u * u + P;

    const double E = P / (gamma - 1.0) + 0.5 * rho * u * u;
    const double energy_flux = u * (E + P);

    Flux flux;
    flux.mass = mass_flux;
    flux.momentum = momentum_flux;
    flux.energy = energy_flux;
    return flux;
}
