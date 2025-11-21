#ifndef ACOUSTICRIEMANNSOLVER_HPP
#define ACOUSTICRIEMANNSOLVER_HPP

#include <cmath>

#include "data/Variables.hpp"
#include "riemann/RiemannSolver.hpp"

/**
 * @class AcousticRiemannSolver
 * @brief Linearized (acoustic) Riemann solver for the Euler equations.
 *
 * This solver implements a simple acoustic approximation of the
 * Euler Riemann problem by linearizing the equations in terms of
 * velocity @f$u@f$ and pressure @f$p@f$ near the interface.
 *
 * For given left/right states
 *   - @f$(\rho_L, u_L, p_L)@f$,
 *   - @f$(\rho_R, u_R, p_R)@f$,
 *
 * we introduce sound speeds and acoustic impedances
 *
 *   @f[
 *     c_L = \sqrt{\gamma p_L / \rho_L}, \qquad
 *     c_R = \sqrt{\gamma p_R / \rho_R},
 *   @f]
 *
 *   @f[
 *     Z_L = \rho_L c_L, \qquad Z_R = \rho_R c_R .
 *   @f]
 *
 * The linearized (acoustic) Riemann solver for the system
 * @f[
 *   p_t + Z_0 c_0 u_x = 0, \qquad
 *   u_t + \frac{1}{Z_0} p_x = 0
 * @f]
 * yields values at the interface point @f$\xi = x/t = 0@f$
 *
 *   @f[
 *     u^* = \frac{Z_L u_L + Z_R u_R + (p_L - p_R)}{Z_L + Z_R},
 *   @f]
 *
 *   @f[
 *     p^* = \frac{Z_L p_R + Z_R p_L + Z_L Z_R (u_L - u_R)}
 *                {Z_L + Z_R}.
 *   @f]
 *
 * Next, the density on each side is recovered from the linear
 * relationship @f$p' = c^2 \rho'@f$:
 *
 *   @f[
 *     \rho_L^* = \rho_L + \frac{p^* - p_L}{c_L^2}, \qquad
 *     \rho_R^* = \rho_R + \frac{p^* - p_R}{c_R^2},
 *   @f]
 *
 * and the interface density is taken as the average
 *
 *   @f[
 *     \rho^* = \frac12(\rho_L^* + \rho_R^*).
 *   @f]
 *
 * Finally, the numerical flux is computed as the usual Euler
 * flux at the point @f$(\rho^*, u^*, p^*)@f$:
 *
 *   @f[
 *     F^* = F_{\text{Euler}}(\rho^*, u^*, p^*).
 *   @f]
 *
 * This approach is only correct for "weak" perturbations and serves
 * as an educational example; for strong discontinuities (Sod) it
 * provides a rough approximation.
 */
class AcousticRiemannSolver : public RiemannSolver {
   public:
    /**
     * @brief Default constructor.
     *
     * Initializes internal density / pressure floors used to safeguard
     * the computation of sound speeds and impedances.
     */
    AcousticRiemannSolver();

    /**
     * @brief Computes the numerical flux using an acoustic approximation.
     *
     * Algorithm:
     *  1. Compute @f$c_L, c_R@f$ and impedances @f$Z_L, Z_R@f$.
     *  2. Find interface values @f$u^*, p^*@f$ using the formulas above.
     *  3. Recover @f$\rho_L^*, \rho_R^*@f$ from the linear EOS
     *     @f$p' = c^2\rho'@f$ and take
     *        @f$\rho^* = 0.5(\rho_L^* + \rho_R^*)@f$.
     *  4. Assemble the primitive state @f$W^* = (\rho^*, u^*, p^*)@f$.
     *  5. Return @f$F^* = \text{EulerFlux}(W^*, \gamma)@f$.
     *
     * If the averaged state is degenerate (non-positive density
     * or pressure, zero impedances, etc.), the method falls back to
     * a simple central flux
     *
     *   @f[
     *     F = \frac12 \bigl( F(W_L) + F(W_R) \bigr).
     *   @f]
     *
     * @param left  Left primitive state at the interface.
     * @param right Right primitive state at the interface.
     * @param gamma Ratio of specific heats.
     * @return Numerical flux vector through the interface.
     */
    [[nodiscard]] auto ComputeFlux(const Primitive& left, const Primitive& right,
                                   double gamma) const -> Flux override;

   private:
    double rho_min_;  ///< Density floor used when computing sound speed.
    double p_min_;    ///< Pressure floor used when computing sound speed.
};

#endif  // ACOUSTICRIEMANNSOLVER_HPP
