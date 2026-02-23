#ifndef ROERIEMANNSOLVER_HPP
#define ROERIEMANNSOLVER_HPP

#include "riemann/RiemannSolver.hpp"

/**
 * @class RoeRiemannSolver
 * @brief Roe approximate Riemann solver for Euler equations (ideal gas), axis-aligned.
 *
 * Uses Roe-averaged states to linearize the Jacobian along the interface normal axis.
 * Includes a Harten–Hyman type entropy fix for the acoustic eigenvalues.
 */
class RoeRiemannSolver final : public RiemannSolver {
public:
    RoeRiemannSolver() = default;

    [[nodiscard]] auto ComputeFlux(const PrimitiveCell& left,
                                   const PrimitiveCell& right,
                                   double gamma,
                                   Axis axis) const -> FluxCell override;
};

#endif  // ROERIEMANNSOLVER_HPP