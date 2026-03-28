#ifndef RUSANOVRIEMANNSOLVER_HPP
#define RUSANOVRIEMANNSOLVER_HPP

#include "riemann/RiemannSolver.hpp"

/**
 * @class RusanovRiemannSolver
 * @brief Local Lax-Friedrichs (Rusanov) approximate Riemann solver.
 *
 * Computes conservative numerical flux through one face with an arbitrary
 * unit normal.
 */
class RusanovRiemannSolver final : public RiemannSolver {
public:
    RusanovRiemannSolver() = default;

    [[nodiscard]] ConservativeCell ComputeFlux(const PrimitiveCell& left,
                                               const PrimitiveCell& right,
                                               double gamma,
                                               const FaceNormal& normal) const override;
};

#endif  // RUSANOVRIEMANNSOLVER_HPP
