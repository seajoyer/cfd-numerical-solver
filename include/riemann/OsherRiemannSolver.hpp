#ifndef OSHERRIEMANNSOLVER_HPP
#define OSHERRIEMANNSOLVER_HPP

#include "riemann/RiemannSolver.hpp"

/**
 * @class OsherRiemannSolver
 * @brief Osher–Solomon approximate Riemann solver (axis-aligned Euler, ideal gas).
 *
 * Uses straight-line path in primitive space and 3-point Gauss–Legendre quadrature:
 *   F = 0.5(F_L + F_R) - 0.5 * \int_0^1 |A(W(theta))| (U_R - U_L) d theta.
 *
 * The action |A|·dU is computed "as before":
 *  - 1D subsystem (rho, rho*u_n, E): eigen-decomposition with 3x3 solve (R alpha = dU).
 *  - Tangential momenta: treated as passive advection with eigenvalue u_n:
 *      (|A| dU)_{rho*u_t} = |u_n| * d(rho*u_t).
 *
 * Falls back to Rusanov (LLF) if the local state becomes non-physical or the 3x3 solve fails.
 */
class OsherRiemannSolver final : public RiemannSolver {
public:
    OsherRiemannSolver() = default;

    [[nodiscard]] auto ComputeFlux(const PrimitiveCell& left,
                                   const PrimitiveCell& right,
                                   double gamma,
                                   Axis axis) const -> FluxCell override;

private:
    [[nodiscard]] static bool IsPhysical(const PrimitiveCell& w);

    [[nodiscard]] static PrimitiveCell ToNormalFrame(const PrimitiveCell& w, Axis axis);

    [[nodiscard]] static ConservativeCell Diff(const ConservativeCell& a, const ConservativeCell& b);

    [[nodiscard]] static bool Solve3x3(double r1_0, double r1_1, double r1_2,
                                       double r2_0, double r2_1, double r2_2,
                                       double r3_0, double r3_1, double r3_2,
                                       double b0, double b1, double b2,
                                       double& x1, double& x2, double& x3);

    // Apply |A(W)| * dU in 1D (normal frame) for components (rho, rho*un, E),
    // tangential momenta as passive advection with |un|.
    [[nodiscard]] static bool ApplyAbsJacobian(const PrimitiveCell& w,
                                               double gamma,
                                               const ConservativeCell& dU,
                                               ConservativeCell& out);

    [[nodiscard]] static FluxCell RusanovFluxRot(const PrimitiveCell& WL,
                                                 const PrimitiveCell& WR,
                                                 double gamma);

    static void ComposeMomentumFluxToGlobal(const FluxCell& F_rot, Axis axis, FluxCell& F);
};

#endif  // OSHERRIEMANNSOLVER_HPP
