#ifndef ACOUSTICRIEMANNSOLVER_HPP
#define ACOUSTICRIEMANNSOLVER_HPP

#include "riemann/RiemannSolver.hpp"
#include "data/Variables.hpp"

#include <cmath>

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
 * Линеаризованный (акустический) Riemann-решатель для системы
 * @f[
 *   p_t + Z_0 c_0 u_x = 0, \qquad
 *   u_t + \frac{1}{Z_0} p_x = 0
 * @f]
 * даёт значения в точке интерфейса @f$\xi = x/t = 0@f$
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
 * Далее плотность на каждой стороне восстанавливается из линейной
 * связи @f$p' = c^2 \rho'@f$:
 *
 *   @f[
 *     \rho_L^* = \rho_L + \frac{p^* - p_L}{c_L^2}, \qquad
 *     \rho_R^* = \rho_R + \frac{p^* - p_R}{c_R^2},
 *   @f]
 *
 * и интерфейсная плотность берётся как среднее
 *
 *   @f[
 *     \rho^* = \tfrac12(\rho_L^* + \rho_R^*).
 *   @f]
 *
 * Наконец, численный поток вычисляется как обычный эйлеровский
 * поток в точке @f$(\rho^*, u^*, p^*)@f$:
 *
 *   @f[
 *     F^* = F_{\text{Euler}}(\rho^*, u^*, p^*).
 *   @f]
 *
 * Этот подход корректен только для «слабых» возмущений и служит
 * учебным примером; для сильных разрывов (Sod) даёт грубую
 * аппроксимацию.
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
     * Алгоритм:
     *  1. Вычислить @f$c_L, c_R@f$ и импедансы @f$Z_L, Z_R@f$.
     *  2. Найти интерфейсные значения @f$u^*, p^*@f$ по формулам выше.
     *  3. Восстановить @f$\rho_L^*, \rho_R^*@f$ из линейной EOS
     *     @f$p' = c^2\rho'@f$ и взять
     *        @f$\rho^* = 0.5(\rho_L^* + \rho_R^*)@f$.
     *  4. Собрать примитивное состояние @f$W^* = (\rho^*, u^*, p^*)@f$.
     *  5. Вернуть @f$F^* = \text{EulerFlux}(W^*, \gamma)@f$.
     *
     * Если усреднённое состояние вырождено (неположительные плотность
     * или давление, нулевые импедансы и т.п.), метод переходит на
     * простой центральный поток
     *
     *   @f[
     *     F = \tfrac12 \bigl( F(W_L) + F(W_R) \bigr).
     *   @f]
     *
     * @param left  Left primitive state at the interface.
     * @param right Right primitive state at the interface.
     * @param gamma Ratio of specific heats.
     * @return Numerical flux vector through the interface.
     */
    [[nodiscard]] auto ComputeFlux(const Primitive& left,
                                   const Primitive& right,
                                   double gamma) const -> Flux override;

private:
    double rho_min_; ///< Density floor used when computing sound speed.
    double p_min_;   ///< Pressure floor used when computing sound speed.
};

#endif  // ACOUSTICRIEMANNSOLVER_HPP
