#ifndef ANALYTICALSOLVER_HPP
#define ANALYTICALSOLVER_HPP

#include <memory>

#include "../data/Variables.hpp"
#include "Solver.hpp"
#include "config/Settings.hpp"
#include "data/DataLayer.hpp"
#include "riemann/ExactIdealGasRiemannSolver.hpp"

/**
 * @class AnalyticalSolver
 * @brief Exact 1D Riemann problem solver for verification and test problems.
 *
 * This solver does not perform a numerical Godunov update. Instead, it:
 *  - assumes the initial condition is a single Riemann problem with a
 *    discontinuity at x = x0 (from Settings),
 *  - infers left/right primitive states from the initial DataLayer on the
 *    first call to Step(),
 *  - uses ExactIdealGasRiemannSolver to evaluate the self-similar solution
 *    at xi = (x - x0) / t for each core cell.
 *
 * Design:
 *  - Constructor takes only Settings (uniform with other solvers).
 *  - Initial left/right states are detected once from DataLayer on first Step.
 *  - Subsequent calls reuse these stored states.
 *  - Boundary conditions are ignored; the exact solution is defined globally.
 *
 * Time handling:
 *  - If no external dt is set (dt_external_ <= 0), Step() jumps directly
 *    to t_end in a single call.
 *  - If dt_external_ > 0, each Step() advances by dt (clipped to t_end),
 *    and t_cur is incremented accordingly.
 */
class AnalyticalSolver : public Solver {
   public:
    /**
     * @brief Constructs an analytical Riemann solver from global settings.
     *
     * Uses:
     *  - settings.gamma for EOS,
     *  - settings.x0 for discontinuity position,
     *  - settings.t_end for final time.
     *
     * Left/right primitive states will be inferred from DataLayer on the
     * first Step() call.
     *
     * @param settings Global simulation settings (copied internally).
     */
    explicit AnalyticalSolver(const Settings& settings);

    /**
     * @brief Advances (or jumps) the analytical solution in time.
     *
     * On the first call:
     *  - infers left/right states from the current DataLayer using x0,
     *    assuming a piecewise-constant Riemann initial condition.
     *
     * On every call:
     *  - chooses dt (external or jump to t_end),
     *  - evaluates the exact Riemann solution at t = t_cur + dt,
     *  - overwrites core cells in DataLayer with the exact solution,
     *  - updates derived quantities (momentum, energies, etc.),
     *  - increments t_cur by dt.
     *
     * @param layer Data layer (core cells are overwritten).
     * @param t_cur Current simulation time (incremented by dt on return).
     * @return Actual time step taken (dt), or 0.0 if nothing is done.
     */
    auto Step(DataLayer& layer, double& t_cur) -> double override;

    /**
     * @brief Sets an external time step.
     *
     * If dt <= 0, the solver will use a single step to t_end.
     *
     * @param dt External time step size.
     */
    void SetDt(double dt);

    /**
     * @brief Ignored for the analytical solver.
     */
    void SetCfl(double cfl) override;

    /**
     * @brief Ignored for the analytical solver.
     *
     * Boundary conditions are not applied; the solution is defined analytically.
     */
    void AddBoundary(int axis, std::shared_ptr<BoundaryCondition> left_bc,
                     std::shared_ptr<BoundaryCondition> right_bc) override;

   private:
    Settings settings_;
    double gamma_;
    double x0_;

    Primitive left_;
    Primitive right_;
    bool has_initial_states_;

    ExactIdealGasRiemannSolver exact_;
    double dt_external_;

    /**
     * @brief Infer left/right initial states from the current DataLayer.
     *
     * Strategy:
     *  - Find the first core index i such that xc(i) >= x0.
     *  - Take left state from i-1 (clamped to coreStart),
     *    right state from i   (clamped to coreEnd-1).
     *  - Read (rho, u, P) from DataLayer at those indices.
     *
     * Assumes the initial condition is a Riemann problem with piecewise
     * constant states separated at x0.
     */
    void InitializeStatesFromLayer(const DataLayer& layer);
};

#endif  // ANALYTICALSOLVER_HPP
