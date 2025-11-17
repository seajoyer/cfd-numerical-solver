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
     * @brief Constructs the solver with xi = 0 and Q = 2.
     */
    ExactIdealGasRiemannSolver();

    /**
     * @brief Constructs the solver with specific xi and Q_user.
     *
     * @param xi Similarity coordinate (x - x0) / t.
     * @param Q_user boundary for P_max / P_min.
     */
    ExactIdealGasRiemannSolver(const double xi,
                               const double Q_user) : xi_(xi), Q_user_(Q_user) {}

    /**
     * @brief Sets the similarity coordinate xi used for sampling.
     *
     * For the Godunov method, xi = 0 corresponds to the cell interface.
     *
     * @param xi Similarity coordinate (x - x0) / t.
     */
    void SetXi(double xi);

    /**
     * @brief Sets the Q
     * @param Q .
     */
    void SetQ(double Q);

    /**
     * @brief Samples the exact Riemann solution at given xi.
     *
     * Uses the same internal exact solution as ComputeFlux.
     *
     * @param left Left primitive state.
     * @param right Right primitive state.
     * @param gamma Ratio of specific heats.
     * @param xi Similarity coordinate (x - x0) / t.
     * @return Primitive state of the exact solution at xi.
     */
    [[nodiscard]] auto Sample(const Primitive& left,
                              const Primitive& right,
                              double gamma,
                              double xi) const -> Primitive;

    /**
     * @brief Computes the exact Riemann flux for ideal gas.
     *
     * @param left Left primitive state.
     * @param right Right primitive state.
     * @param gamma Ratio of specific heats.
     * @return Exact flux at the interface.
     */
    [[nodiscard]] auto ComputeFlux(const Primitive& left,
                                   const Primitive& right,
                                   double gamma) const -> Flux override;

private:
    double xi_;
    double Q_user_{2.};
};

#endif  // EXACTIDEALGASRIEMANNSOLVER_HPP