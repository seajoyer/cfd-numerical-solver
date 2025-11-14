#ifndef HLLRIEMANNSOLVER_HPP
#define HLLRIEMANNSOLVER_HPP

#include "RiemannSolver.hpp"

/**
 * @class HLLRiemannSolver
 * @brief HLL (Harten–Lax–van Leer) approximate Riemann solver.
 *
 * Uses a two-wave approximation with estimates of the minimal and maximal
 * signal velocities. It is robust and diffusive and does not resolve
 * contact discontinuities exactly.
 */
class HLLRiemannSolver : public RiemannSolver {
public:
    /**
     * @brief Default constructor.
     */
    HLLRiemannSolver() = default;

    /**
     * @brief Computes the HLL numerical flux.
     *
     * @param left Left primitive state.
     * @param right Right primitive state.
     * @param gamma Ratio of specific heats.
     * @return HLL flux at the interface.
     */
    [[nodiscard]] auto ComputeFlux(const Primitive& left, const Primitive& right,
                                   double gamma) const -> Flux override;
};

#endif  // HLLRIEMANNSOLVER_HPP
