#ifndef RUSANOVRIEMANNSOLVER_HPP
#define RUSANOVRIEMANNSOLVER_HPP

#include "RiemannSolver.hpp"

/**
 * @class RusanovRiemannSolver
 * @brief Local Lax–Friedrichs (Rusanov) approximate Riemann solver.
 *
 * This solver uses a single dissipation term based on the maximal
 * local signal speed. The numerical flux is given by
 *
 * \f[
 *   F_{\text{rus}} =
 *   \frac{1}{2}\left(F_L + F_R\right)
 *   - \frac{1}{2}\,a_{\max}\left(U_R - U_L\right),
 * \f]
 *
 * where \f$a_{\max} = \max(|u_L| + c_L,\; |u_R| + c_R)\f$ and
 * \f$c = \sqrt{\gamma P / \rho}\f is the sound speed.
 *
 * It is very robust but relatively diffusive and does not resolve
 * contact discontinuities sharply.
 */
class RusanovRiemannSolver : public RiemannSolver {
public:
    /**
     * @brief Default constructor.
     */
    RusanovRiemannSolver() = default;

    /**
     * @brief Computes the Rusanov (local Lax–Friedrichs) numerical flux.
     *
     * @param left  Left primitive state.
     * @param right Right primitive state.
     * @param gamma Ratio of specific heats.
     * @return Rusanov flux at the interface.
     */
    [[nodiscard]] auto ComputeFlux(const Primitive& left,
                                   const Primitive& right,
                                   double gamma) const -> Flux override;
};

#endif  // RUSANOVRIEMANNSOLVER_HPP