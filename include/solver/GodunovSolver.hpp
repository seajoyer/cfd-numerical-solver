#ifndef GODUNOVSOLVER_HPP
#define GODUNOVSOLVER_HPP

#include <memory>
#include <vector>

#include "Solver.hpp"
#include "config/Settings.hpp"
#include "data/DataLayer.hpp"
#include "bc/BoundaryManager.hpp"

#include "reconstruction/Reconstruction.hpp"
#include "reconstruction/P0Reconstruction.hpp"
#include "TimeStepCalculator.hpp"
#include "PositivityLimiter.hpp"

// Подразумевается, что Riemann-солверы лежат в /include/reimann
#include "riemann/RiemannSolver.hpp"
#include "riemann/HLLRiemannSolver.hpp"
#include "riemann/HLLCRiemannSolver.hpp"
#include "riemann/ExactIdealGasRiemannSolver.hpp"

/**
 * @class GodunovSolver
 * @brief First-order finite-volume Godunov solver for 1D Euler equations.
 *
 * This solver implements a classical Godunov scheme with:
 * - piecewise-constant (P0) reconstruction,
 * - configurable Riemann solver (HLL / HLLC / Exact),
 * - CFL-based time step selection,
 * - optional positivity limiting for density and pressure.
 *
 * The solver operates on a provided DataLayer instance and uses
 * BoundaryManager to update ghost cells before each step.
 *
 * Responsibilities per Step():
 *  1. Apply all configured boundary conditions.
 *  2. Convert conservative variables to primitive variables.
 *  3. Compute the stable time step from CFL condition.
 *  4. Reconstruct left/right states at interfaces.
 *  5. Evaluate numerical fluxes via the chosen Riemann solver.
 *  6. Update conservative variables (Godunov FV update).
 *  7. Enforce positivity constraints on updated states.
 *
 * Output is handled externally (e.g. via StepWriter in Simulation).
 */
class GodunovSolver : public Solver {
public:
    /**
     * @brief Constructs a Godunov solver using global simulation settings.
     *
     * The constructor initializes internal parameters (CFL, gamma, domain size)
     * and allocates helper objects (BoundaryManager, reconstruction, Riemann solver).
     * Selection of a specific Riemann solver implementation is based on settings
     * (e.g. configuration key, to be interpreted in the implementation).
     *
     * @param settings Global simulation settings (copied internally).
     */
    explicit GodunovSolver(const Settings &settings);

    /**
     * @brief Advances the solution by one explicit Godunov time step.
     *
     * Performs a full update of conservative variables stored in the given
     * DataLayer, starting from the current time t_cur. The actual time step
     * size is computed from the CFL condition and returned.
     *
     * @param layer Data layer containing current state (modified in-place).
     * @param t_cur Current simulation time (incremented by dt on return).
     * @return Actual time step taken (dt).
     */
    double Step(DataLayer &layer, double &t_cur) override;

    /**
     * @brief Sets the CFL number used for time step computation.
     *
     * Overrides the value taken from Settings, if needed.
     *
     * @param cfl CFL number (0 < cfl <= 1).
     */
    void SetCfl(double cfl) override;

    /**
     * @brief Assigns boundary conditions for a specific axis.
     *
     * This method forwards boundary configuration to the internal
     * BoundaryManager. For 1D simulations, axis = 0 is expected.
     *
     * @param axis Spatial axis index (0=x, 1=y, 2=z).
     * @param left_bc Boundary condition for the left/lower boundary.
     * @param right_bc Boundary condition for the right/upper boundary.
     */
    void AddBoundary(int axis,
                     std::shared_ptr<BoundaryCondition> left_bc,
                     std::shared_ptr<BoundaryCondition> right_bc) override;

private:
    /**
     * @brief Global simulation settings (copied from input).
     */
    Settings settings_;

    /**
     * @brief Manager for all boundary conditions.
     *
     * Used to update ghost cells before each Godunov step.
     */
    BoundaryManager boundaryManager_;

    /**
     * @brief Reconstruction scheme (currently P0 for first-order Godunov).
     */
    std::shared_ptr<Reconstruction> reconstruction_;

    /**
     * @brief Selected Riemann solver implementation.
     */
    std::shared_ptr<RiemannSolver> riemannSolver_;

    /**
     * @brief Minimal allowed density for positivity limiting.
     */
    double rhoMin_;

    /**
     * @brief Minimal allowed pressure for positivity limiting.
     */
    double pMin_;

    /**
     * @brief Initializes the reconstruction strategy.
     *
     * Currently selects piecewise-constant reconstruction.
     */
    void InitializeReconstruction();

    /**
     * @brief Initializes the Riemann solver according to settings.
     *
     * Chooses between HLL, HLLC, or exact ideal-gas solver depending
     * on configuration (details implemented in source file).
     */
    void InitializeRiemannSolver();

    /**
     * @brief Computes spatial step size for 1D uniform grid.
     *
     * Uses either configuration parameters (domain length and N)
     * or coordinates stored in the DataLayer in the implementation.
     *
     * @param layer Data layer with grid description.
     * @return Cell size dx.
     */
    double ComputeDx(const DataLayer &layer) const;

    /**
     * @brief Converts conservative variables in DataLayer to primitives.
     *
     * Extracts 1D conservative fields from DataLayer, converts them into
     * a contiguous array of Primitive states for the physical cells.
     *
     * @param layer Data layer with current conservative variables.
     * @param primitives Output vector of primitive states.
     */
    void BuildPrimitiveArray(const DataLayer &layer,
                             std::vector<Primitive> &primitives) const;

    /**
     * @brief Writes updated conservative variables back into DataLayer.
     *
     * Copies the updated conservative states of physical cells into
     * the corresponding arrays (U, V, etc.) stored in the DataLayer.
     *
     * @param updatedConservative Updated conservative states.
     * @param layer Data layer to be modified.
     */
    void StoreConservativeArray(const std::vector<Conservative> &updatedConservative,
                                DataLayer &layer) const;
};

#endif  // GODUNOVSOLVER_HPP
