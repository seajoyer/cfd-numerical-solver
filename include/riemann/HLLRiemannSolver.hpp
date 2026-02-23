#ifndef HLLRIEMANNSOLVER_HPP
#define HLLRIEMANNSOLVER_HPP

#include "riemann/RiemannSolver.hpp"

/**
 * @class HLLRiemannSolver
 * @brief HLL approximate Riemann solver (two-wave model).
 */
class HLLRiemannSolver final : public RiemannSolver {
public:
    HLLRiemannSolver() = default;

    [[nodiscard]] auto ComputeFlux(const PrimitiveCell& left,
                                   const PrimitiveCell& right,
                                   double gamma,
                                   Axis axis) const -> FluxCell override;
};

#endif  // HLLRIEMANNSOLVER_HPP