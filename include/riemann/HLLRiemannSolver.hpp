#ifndef HLLRIEMANNSOLVER_HPP
#define HLLRIEMANNSOLVER_HPP

#include "riemann/RiemannSolver.hpp"

/**
 * @class HLLRiemannSolver
 * @brief HLL approximate Riemann solver (two-wave model).
 *
 * Computes conservative numerical flux through one face with an arbitrary
 * unit normal.
 */
class HLLRiemannSolver final : public RiemannSolver {
public:
    HLLRiemannSolver() = default;

    [[nodiscard]] ConservativeCell ComputeFlux(const PrimitiveCell& left,
                                               const PrimitiveCell& right,
                                               double gamma,
                                               const FaceNormal& normal) const override;
};

#endif  // HLLRIEMANNSOLVER_HPP
