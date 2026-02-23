#ifndef RUSANOVRIEMANNSOLVER_HPP
#define RUSANOVRIEMANNSOLVER_HPP

#include "riemann/RiemannSolver.hpp"

/**
 * @class RusanovRiemannSolver
 * @brief Local Lax–Friedrichs (Rusanov) approximate Riemann solver.
 */
class RusanovRiemannSolver final : public RiemannSolver {
public:
    RusanovRiemannSolver() = default;

    [[nodiscard]] auto ComputeFlux(const PrimitiveCell& left,
                                   const PrimitiveCell& right,
                                   double gamma,
                                   Axis axis) const -> FluxCell override;
};

#endif  // RUSANOVRIEMANNSOLVER_HPP