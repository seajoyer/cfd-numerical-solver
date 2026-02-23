#ifndef EXACTIDEALGASRIEMANNSOLVER_HPP
#define EXACTIDEALGASRIEMANNSOLVER_HPP

#include "riemann/RiemannSolver.hpp"

/**
 * @class ExactIdealGasRiemannSolver
 * @brief Exact Riemann solver for the 1D ideal-gas Euler equations, embedded in axis-aligned 3D flux.
 *
 * Solves the exact 1D Riemann problem along the interface normal (Axis) for (rho, u_n, p),
 * then reconstructs a 3D primitive state by taking tangential velocities from the appropriate side
 * of the contact discontinuity.
 *
 * Intended mainly for validation; more expensive than approximate solvers.
 */
class ExactIdealGasRiemannSolver final : public RiemannSolver {
public:
    /** @brief Constructs the solver with xi = 0 and Q = 2. */
    ExactIdealGasRiemannSolver();

    /** @brief Constructs the solver with specific xi and Q_user (PVRS switching threshold). */
    ExactIdealGasRiemannSolver(double xi, double Q_user) : xi_(xi), Q_user_(Q_user) {}

    /** @brief Sets the similarity coordinate xi used for sampling (interface: xi = 0). */
    void SetXi(double xi);

    /** @brief Sets the pressure ratio threshold Q used in initial guess selection. */
    void SetQ(double Q);

    /**
     * @brief Sample the exact solution at given xi and interface axis.
     *
     * @param left  Left primitive state (rho,u,v,w,P).
     * @param right Right primitive state (rho,u,v,w,P).
     * @param gamma Ratio of specific heats.
     * @param xi    Similarity coordinate (x-x0)/t in the 1D normal direction.
     * @param axis  Interface normal axis.
     * @return PrimitiveCell sample (rho,u,v,w,P).
     */
    [[nodiscard]] PrimitiveCell Sample(const PrimitiveCell& left,
                                       const PrimitiveCell& right,
                                       double gamma,
                                       double xi,
                                       Axis axis) const;

    /** @brief Compute the exact numerical flux at the interface (xi = xi_). */
    [[nodiscard]] auto ComputeFlux(const PrimitiveCell& left,
                                   const PrimitiveCell& right,
                                   double gamma,
                                   Axis axis) const -> FluxCell override;

private:
    double xi_ = 0.0;
    double Q_user_ = 2.0;
};

#endif  // EXACTIDEALGASRIEMANNSOLVER_HPP