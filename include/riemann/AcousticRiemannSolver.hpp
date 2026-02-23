#ifndef ACOUSTICRIEMANNSOLVER_HPP
#define ACOUSTICRIEMANNSOLVER_HPP

#include "riemann/RiemannSolver.hpp"

/**
 * @class AcousticRiemannSolver
 * @brief Linearized (acoustic) Riemann solver for Euler equations (axis-aligned).
 *
 * Uses acoustic impedance matching to compute interface (u*, p*) for the normal velocity,
 * recovers rho* from p' = c^2 rho', and returns Euler flux at W*.
 *
 * Intended for weak perturbations / educational use.
 */
class AcousticRiemannSolver final : public RiemannSolver {
public:
    AcousticRiemannSolver();

    [[nodiscard]] auto ComputeFlux(const PrimitiveCell& left,
                                   const PrimitiveCell& right,
                                   double gamma,
                                   Axis axis) const -> FluxCell override;

private:
    double rho_min_;
    double p_min_;
};

#endif  // ACOUSTICRIEMANNSOLVER_HPP