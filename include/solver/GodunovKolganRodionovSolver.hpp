#ifndef GODUNOVKOLGANRODIONOVSOLVER_HPP
#define GODUNOVKOLGANRODIONOVSOLVER_HPP

#include <cctype>
#include <memory>

#include "Solver.hpp"
#include "bc/BoundaryManager.hpp"
#include "config/Settings.hpp"
#include "data/DataLayer.hpp"
#include "data/Variables.hpp"
#include "reconstruction/Reconstruction.hpp"
#include "riemann/RiemannSolver.hpp"

/**
 * @class GodunovKolganRodionovSolver
 * @brief Second-order Godunov-type solver (Kolgan–Rodionov / MUSCL–Hancock).
 *
 * This solver extends the first-order Godunov method to second-order accuracy
 * in smooth regions using:
 *  - piecewise-linear reconstruction of primitive variables with a
 *    TVD minmod-type limiter (P1Reconstruction),
 *  - a MUSCL–Hancock predictor step to advance interface states to half time,
 *  - a configurable Riemann solver (HLL / HLLC / Exact) evaluated at the
 *    predicted half-time interface states.
 *
 * The solver operates directly on the xtensor-based DataLayer storage and
 * uses only small stack-allocated helper objects (Primitive, Conservative,
 * Flux) for per-cell and per-interface computations. A single xt::xarray of
 * Conservative states is used to buffer the updated core states before
 * writing them back to DataLayer.
 */
class GodunovKolganRodionovSolver : public Solver {
   public:
    /**
     * @brief Constructs the solver using global simulation settings.
     *
     * The constructor:
     *  - stores the Settings object locally,
     *  - initializes the BoundaryManager,
     *  - selects and constructs the reconstruction and Riemann solver,
     *  - sets CFL and positivity thresholds.
     *
     * @param settings Global simulation settings (copied internally).
     */
    explicit GodunovKolganRodionovSolver(const Settings& settings);

    /**
     * @brief Performs one second-order Godunov–Kolgan–Rodionov timestep.
     *
     * The method:
     *  - applies boundary conditions,
     *  - computes a stable dt from the CFL condition (and clamps it so that
     *    the final time t_end is not exceeded),
     *  - for each physical cell:
     *      * uses P1Reconstruction to obtain left/right primitive states
     *        inside the cell,
     *      * performs a MUSCL–Hancock predictor to half time,
     *      * builds interface states at j-1/2 and j+1/2 from predicted
     *        cell-edge values,
     *      * computes Riemann fluxes and updates the conservative state,
     *      * applies the positivity limiter,
     *  - writes updated states back into DataLayer core cells,
     *  - updates the simulation time t_cur.
     *
     * If dt <= 0 or the physical domain has fewer than two cells, no update
     * is performed and 0.0 is returned.
     *
     * @param layer Data layer containing the current state (modified in-place).
     * @param t_cur Current simulation time; incremented by dt on return.
     * @return The actual timestep taken (dt), or 0.0 if no step is performed.
     */
    auto Step(DataLayer& layer, double& t_cur) -> double override;

    /**
     * @brief Sets the CFL number used for timestep computation.
     *
     * @param cfl CFL number (0 < cfl <= 1).
     */
    void SetCfl(double cfl) override { cfl_ = cfl; }

    /**
     * @brief Assigns boundary conditions for a specific axis.
     *
     * For 1D problems axis = 0 is expected.
     *
     * @param axis     Spatial axis index.
     * @param left_bc  Boundary condition for the left / lower boundary.
     * @param right_bc Boundary condition for the right / upper boundary.
     */
    void AddBoundary(int axis, std::shared_ptr<BoundaryCondition> left_bc,
                     std::shared_ptr<BoundaryCondition> right_bc) override;

   private:
    /** @brief Local copy of simulation settings. */
    Settings settings_;

    /** @brief Manager for boundary conditions. */
    BoundaryManager boundary_manager_;

    /** @brief Reconstruction strategy (P1 piecewise-linear). */
    std::shared_ptr<Reconstruction> reconstruction_;

    /** @brief Selected Riemann solver implementation. */
    std::shared_ptr<RiemannSolver> riemann_solver_;

    /** @brief Minimal allowed density for positivity limiting. */
    double rho_min_;

    /** @brief Minimal allowed pressure for positivity limiting. */
    double p_min_;

    /**
     * @brief Initializes the reconstruction scheme.
     *
     * For this solver, a piecewise-linear reconstruction (P1) is required
     * to achieve second-order accuracy. If settings_.reconstruction
     * contains "p1" (case-insensitive), P1Reconstruction is used; otherwise
     * the solver still defaults to P1Reconstruction.
     */
    void InitializeReconstruction();

    /**
     * @brief Initializes the Riemann solver according to settings_.riemann_solver.
     *
     * Case-insensitive selection:
     *  - contains "hllc"  → HLLCRiemannSolver
     *  - contains "hll"   → HLLRiemannSolver
     *  - contains "exact" → ExactIdealGasRiemannSolver
     *  - otherwise        → HLLRiemannSolver (default)
     */
    void InitializeRiemannSolver();

    /**
     * @brief Computes the spatial step size for a uniform 1D grid.
     *
     * Preference order:
     *  1. If settings_.L_x > 0 and settings_.N > 0:
     *       dx = L_x / N
     *  2. Otherwise, attempt to deduce dx from DataLayer cell centers xc.
     *  3. Fall back to dx = 1.0 if none of the above are available.
     *
     * @param layer Data layer with mesh coordinates.
     * @return Cell size dx.
     */
    [[nodiscard]] auto ComputeDx(const DataLayer& layer) const -> double;

    /**
     * @brief Writes a single conservative state back into the DataLayer.
     *
     * Converts (rho, rhoU, E) into primitive and auxiliary quantities and
     * stores them in:
     *  - rho, u, P
     *  - p (momentum)
     *  - V (specific volume)
     *  - U (specific internal energy)
     *  - e (total energy density)
     *  - m (cell mass = rho * dx)
     *
     * Pressure is computed via EOS::Pressure().
     *
     * @param uc    Conservative state to store.
     * @param i     Cell index in the 1D layout (including ghosts).
     * @param dx    Cell size (used for mass computation).
     * @param layer Data layer to modify.
     */
    void StoreConservativeCell(const Conservative& uc, int i, double dx,
                               DataLayer& layer) const;

    /**
     * @brief Computes MUSCL–Hancock predicted states in a given cell using P1
     * reconstruction.
     *
     * For cell index i, this function:
     *  - calls P1-based reconstruction at interfaces i-1/2 and i+1/2
     *    through Reconstruction::ComputeInterfaceStates,
     *  - extracts the left and right primitive states inside the cell from
     *    the corresponding interface states,
     *  - converts them to conservative form,
     *  - computes Euler fluxes for these left/right states,
     *  - performs a half-step MUSCL–Hancock predictor:
     *
     *       U_i^{L,*} = U_i^L - (dt/(2 dx)) [F(U_i^R) - F(U_i^L)]
     *       U_i^{R,*} = U_i^R - (dt/(2 dx)) [F(U_i^R) - F(U_i^L)].
     *
     * Interface indices are clamped to the valid range to safely work
     * near the domain boundaries using ghost cells.
     *
     * @param layer        Data layer with primitive fields.
     * @param i            Cell index.
     * @param halfDtOverDx dt / (2 * dx) factor.
     * @param U_L_star_out Output: predicted left state U_i^{L,*}.
     * @param U_R_star_out Output: predicted right state U_i^{R,*}.
     */
    void ComputePredictedStatesAtCell(const DataLayer& layer, int i,
                                      double half_dt_over_dx, Conservative& U_L_star_out,
                                      Conservative& U_R_star_out) const;
};

#endif  // GODUNOVKOLGANRODIONOVSOLVER_HPP
