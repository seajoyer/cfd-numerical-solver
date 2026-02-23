#ifndef HLLCRIEMANNSOLVER_HPP
#define HLLCRIEMANNSOLVER_HPP

#include "riemann/RiemannSolver.hpp"

/**
 * @class HLLCRiemannSolver
 * @brief HLLC approximate Riemann solver (restores contact wave).
 *
 * Computes numerical flux for Euler equations (ideal gas) at an interface aligned
 * with a coordinate axis.
 */
class HLLCRiemannSolver final : public RiemannSolver {
public:
    HLLCRiemannSolver() = default;

    [[nodiscard]] auto ComputeFlux(const PrimitiveCell& left,
                                   const PrimitiveCell& right,
                                   double gamma,
                                   Axis axis) const -> FluxCell override;
};

#endif  // HLLCRIEMANNSOLVER_HPP