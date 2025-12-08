#ifndef OSHERRIEMANNSOLVER_HPP
#define OSHERRIEMANNSOLVER_HPP

#include "RiemannSolver.hpp"

/**
 * @class OsherRiemannSolver
 * @brief Osher–Solomon–type approximate Riemann solver for 1D Euler equations.
 *
 * This solver computes the numerical flux as
 * \f[
 *   F^{OS}(U_L, U_R)
 *     = \frac{1}{2}\big(F(U_L) + F(U_R)\big)
 *       - \frac{1}{2} \int_0^1 |A(U(\theta))|\,(U_R - U_L)\,d\theta,
 * \f]
 * where \f$A = \partial F/\partial U\f$ is the flux Jacobian and
 * the path in state space is chosen as a straight segment
 * \f$U(\theta) = (1-\theta)U_L + \theta U_R\f$.
 *
 * The path integral of the dissipation term is approximated by a
 * three–point Gauss–Legendre quadrature.  At each quadrature point
 * the Jacobian is diagonalized via its local eigenstructure and the
 * action of \f$|A|\f$ on \f$\Delta U = U_R-U_L\f$ is computed as
 * \f[
 *   |A|\,\Delta U = \sum_{k=1}^3 |\lambda_k|\,\alpha_k\,r_k,
 * \f]
 * with \f$\lambda_k\f$ and \f$r_k\f$ being eigenvalues and right
 * eigenvectors of the Euler system, and coefficients \f$\alpha_k\f$
 * obtained from solving
 * \f$R \alpha = \Delta U\f$ for \f$\alpha\f$.
 *
 * If the local eigenmatrix becomes (almost) singular or the primitive
 * variables are non-physical, the solver falls back to a simple
 * Rusanov (local Lax–Friedrichs) flux to maintain robustness.
 */
class OsherRiemannSolver : public RiemannSolver {
public:
    /**
     * @brief Default constructor.
     */
    OsherRiemannSolver() = default;

    /**
     * @brief Computes the Osher–Solomon numerical flux.
     *
     * @param left  Left primitive state.
     * @param right Right primitive state.
     * @param gamma Ratio of specific heats.
     * @return Numerical flux at the interface.
     */
    [[nodiscard]] auto ComputeFlux(const Primitive& left,
                                   const Primitive& right,
                                   double gamma) const -> Flux override;

private:
    /**
     * @brief Computes the Euler Jacobian action |A(U)|·dU via eigen-decomposition.
     *
     * Given a local primitive state @p w and a conservative jump @p dU,
     * this routine constructs the right eigenvectors of the 1D Euler system,
     * solves for wave strengths alpha from R * alpha = dU, and then forms
     *
     *   |A| dU = sum_k |lambda_k| alpha_k r_k.
     *
     * @param w      Local primitive state.
     * @param gamma  Ratio of specific heats.
     * @param dU     Conservative jump U_R - U_L.
     * @param[out] result  Vector |A| dU.
     * @return false if the eigenmatrix is nearly singular or state is invalid
     *         (in which case the caller should fall back to a more diffusive flux).
     */
    static auto ApplyAbsJacobian(const Primitive& w,
                                 double gamma,
                                 const Conservative& dU,
                                 Conservative& result) -> bool;

    /**
     * @brief Simple Rusanov (local Lax–Friedrichs) fallback flux.
     *
     * Used when the Osher eigen-decomposition becomes ill-conditioned
     * (e.g. near vacuum).  It uses
     *
     *   F = 0.5 (F_L + F_R) - 0.5 * a_max * (U_R - U_L),
     *
     * where a_max is the maximal signal speed estimated from both sides.
     *
     * @param left  Left primitive state.
     * @param right Right primitive state.
     * @param gamma Ratio of specific heats.
     * @return Rusanov flux.
     */
    [[nodiscard]] static auto RusanovFlux(const Primitive& left,
                                          const Primitive& right,
                                          double gamma) -> Flux;
};

#endif  // OSHERRIEMANNSOLVER_HPP
