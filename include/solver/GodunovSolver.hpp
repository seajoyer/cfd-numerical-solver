#ifndef GODUNOVSOLVER_HPP
#define GODUNOVSOLVER_HPP

#include <memory>
#include <algorithm>
#include <cctype>
#include <string>

#include "Solver.hpp"
#include "config/Settings.hpp"
#include "data/DataLayer.hpp"
#include "bc/BoundaryManager.hpp"

#include "reconstruction/Reconstruction.hpp"

#include "solver/EOS.hpp"
#include "data/Variables.hpp"

#include "riemann/RiemannSolver.hpp"

/**
 * @class GodunovSolver
 * @brief Explicit finite-volume Godunov solver for the 1D Euler equations.
 *
 * This class implements a first-order finite-volume Godunov scheme for
 * solving the 1D Euler equations of gas dynamics. The solver operates
 * directly on the xtensor-based storage provided by DataLayer and uses
 * small stack-allocated helper structures (Primitive, Conservative, Flux)
 * for all per-cell and per-interface computations.
 *
 * ## Numerical algorithm per timestep:
 *
 * 1. Apply boundary conditions via BoundaryManager.
 * 2. Compute a stable time step `dt` using the CFL condition through
 *    TimeStepCalculator.
 * 3. Optionally clamp `dt` so that the simulation does not exceed `t_end`.
 * 4. For every physical cell `j`:
 *      - Convert its primitive variables into conservative form.
 *      - Reconstruct left/right states at interfaces `j-1/2` and `j+1/2`
 *        via the configured reconstruction scheme (P0 or P1).
 *      - Compute interface fluxes with the selected Riemann solver
 *        (HLL, HLLC, or Exact).
 *      - Update the conservative state using the Godunov flux difference:
 *
 *            U^{n+1}_j = U^n_j - (dt/dx) * (F_{j+1/2} - F_{j-1/2})
 *
 *      - Apply the positivity limiter (density & pressure floors).
 *      - Write back the updated conservative and primitive variables
 *        into DataLayer.
 *
 * ## Design notes:
 * - No std::vector is used for entire fields. The solver works
 *   interface-by-interface and cell-by-cell, keeping the memory footprint
 *   minimal and eliminating unnecessary copies.
 * - Reconstruction and Riemann solver types are selected based on the
 *   text configuration stored in Settings.
 * - The solver is written with future extension in mind:
 *   - 2D/3D grids,
 *   - GPU/OMP/MPI parallelization,
 *   - higher-order schemes,
 *   - characteristic-based reconstruction.
 */
class GodunovSolver : public Solver {
public:
    /**
     * @brief Constructs a Godunov solver using the global simulation settings.
     *
     * The constructor:
     *  - stores the given Settings instance,
     *  - initializes the BoundaryManager,
     *  - selects the reconstruction scheme (P0/P1),
     *  - selects the Riemann solver (HLL/HLLC/Exact),
     *  - sets default positivity-limiting thresholds.
     *
     * @param settings Global simulation settings (copied internally).
     */
    explicit GodunovSolver(const Settings& settings);

    /**
     * @brief Advances the solution by a single explicit Godunov timestep.
     *
     * Performs a complete update of the conservative variables stored in
     * the provided DataLayer. The timestep `dt` is selected using the CFL
     * condition, and the simulation time `t_cur` is incremented accordingly.
     *
     * If the physical domain contains fewer than two cells or if the computed
     * timestep is non-positive, no update is performed and `0.0` is returned.
     *
     * @param layer DataLayer containing the current state (modified in-place).
     * @param t_cur Current simulation time (incremented by `dt` on return).
     * @return The actual timestep taken (`dt`), or `0.0` if no step is performed.
     */
    double Step(DataLayer& layer, double& t_cur) override;

    /**
     * @brief Overrides the CFL number used for timestep computation.
     *
     * Allows dynamically changing the CFL value after the solver is created.
     *
     * @param cfl CFL number (0 < cfl <= 1).
     */
    void SetCfl(double cfl) override;

    /**
     * @brief Assigns boundary conditions for a specific axis.
     *
     * For 1D simulations, `axis = 0` is expected.
     *
     * @param axis Spatial axis index.
     * @param left_bc  Boundary condition at the lower/left boundary.
     * @param right_bc Boundary condition at the upper/right boundary.
     */
    void AddBoundary(int axis,
                     std::shared_ptr<BoundaryCondition> left_bc,
                     std::shared_ptr<BoundaryCondition> right_bc) override;

private:
    /** @brief Local copy of global simulation settings. */
    Settings settings_;

    /** @brief Manager responsible for applying all boundary conditions. */
    BoundaryManager boundaryManager_;

    /** @brief Selected reconstruction strategy (P0 or P1). */
    std::shared_ptr<Reconstruction> reconstruction_;

    /** @brief Selected Riemann solver (HLL, HLLC, or Exact). */
    std::shared_ptr<RiemannSolver> riemannSolver_;

    /** @brief Minimum allowed density (for positivity limiting). */
    double rhoMin_;

    /** @brief Minimum allowed pressure (for positivity limiting). */
    double pMin_;

    /**
     * @brief Selects and constructs the reconstruction scheme.
     *
     * Based on the string in `settings_.reconstruction` (case-insensitive):
     *  - contains "p1" → P1Reconstruction
     *  - contains "p0" → P0Reconstruction
     *  - otherwise     → defaults to P0Reconstruction
     */
    void InitializeReconstruction();

    /**
     * @brief Selects and constructs the Riemann solver.
     *
     * Based on the string in `settings_.riemann_solver` (case-insensitive):
     *  - contains "hllc" → HLLCRiemannSolver
     *  - contains "hll"  → HLLRiemannSolver
     *  - contains "exact" → ExactIdealGasRiemannSolver
     *  - otherwise       → defaults to HLLRiemannSolver
     */
    void InitializeRiemannSolver();

    /**
     * @brief Computes the spatial step size `dx` for a uniform 1D grid.
     *
     * Preference order:
     *  1. If Settings specify both domain length `L_x` and number of
     *     physical cells `N`, use `dx = L_x / N`.
     *  2. Otherwise, attempt to deduce `dx` from the DataLayer cell centers
     *     using `xc[i+1] - xc[i]`.
     *  3. Fallback to `dx = 1.0` if both methods fail.
     *
     * @param layer DataLayer storing mesh coordinates.
     * @return Cell size `dx`.
     */
    double ComputeDx(const DataLayer& layer) const;

    /**
     * @brief Writes a conservative state back into the DataLayer.
     *
     * Converts a conservative state `(rho, rhoU, E)` to primitive and
     * auxiliary quantities and stores them in the DataLayer arrays:
     *
     *  - rho, u, P (primitive)
     *  - p (momentum)
     *  - V (specific volume)
     *  - U (specific internal energy)
     *  - e (total energy)
     *  - m (cell mass `rho * dx`)
     *
     * Pressure is computed via `EOS::Pressure()`.
     *
     * @param uc Updated conservative state.
     * @param i  Cell index in the linear layout (including ghosts).
     * @param dx Spatial step size.
     * @param layer DataLayer to modify.
     */
    void StoreConservativeCell(const Conservative& uc,
                               int i,
                               double dx,
                               DataLayer& layer) const;
};

#endif  // GODUNOVSOLVER_HPP
