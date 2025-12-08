#ifndef ROERIEMANNSOLVER_HPP
#define ROERIEMANNSOLVER_HPP

#include "RiemannSolver.hpp"

/**
 * @class RoeRiemannSolver
 * @brief Roe approximate Riemann solver for the 1D Euler equations.
 *
 * Uses Roe-averaged states to linearize the flux Jacobian and constructs
 * the numerical flux as
 *
 * \f[
 *   F_{\text{roe}} =
 *   \frac{1}{2}\left(F_L + F_R\right)
 *   - \frac{1}{2} \sum_{k=1}^{3} |\lambda_k| \alpha_k \mathbf{r}_k,
 * \f]
 *
 * where \f$\lambda_k\f$ and \f$\mathbf{r}_k\f$ are eigenvalues and
 * right eigenvectors of the Roe matrix, and \f$\alpha_k\f$ are the
 * wave strengths for the jump in conservative variables.
 *
 * An entropy fix of Harten–Hyman type is applied to the acoustic
 * eigenvalues to avoid non-physical expansion shocks.
 */
class RoeRiemannSolver : public RiemannSolver {
public:
    /**
     * @brief Default constructor.
     */
    RoeRiemannSolver() = default;

    /**
     * @brief Computes the Roe numerical flux with entropy fix.
     *
     * @param left  Left primitive state.
     * @param right Right primitive state.
     * @param gamma Ratio of specific heats.
     * @return Roe flux at the interface.
     */
    [[nodiscard]] auto ComputeFlux(const Primitive& left,
                                   const Primitive& right,
                                   double gamma) const -> Flux override;
};

#endif  // ROERIEMANNSOLVER_HPP