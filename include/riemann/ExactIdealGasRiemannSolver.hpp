#ifndef EXACTIDEALGASRIEMANNSOLVER_HPP
#define EXACTIDEALGASRIEMANNSOLVER_HPP

#include "RiemannSolver.hpp"
#include "solver/EOS.hpp"

/**
 * @class ExactIdealGasRiemannSolver
 * @brief Exact Riemann solver for the 1D ideal-gas Euler equations.
 *
 * Computes the exact solution of the Riemann problem for a calorically
 * perfect ideal gas and evaluates the corresponding interface flux.
 *
 * @note Intended mainly for validation and reference; more expensive
 *       than approximate solvers such as HLL/HLLC.
 */
class ExactIdealGasRiemannSolver : public RiemannSolver {
public:
    /**
     * @brief Default constructor.
     */
    ExactIdealGasRiemannSolver() {}

    /**
     * @brief Computes the exact Riemann flux for ideal gas.
     *
     * @param left Left primitive state.
     * @param right Right primitive state.
     * @param gamma Ratio of specific heats.
     * @return Exact flux at the interface.
     */
    Flux ComputeFlux(const Primitive &left,
                     const Primitive &right,
                     double gamma) const override;
};

#endif  // EXACTIDEALGASRIEMANNSOLVER_HPP
