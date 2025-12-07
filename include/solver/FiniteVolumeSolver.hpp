#ifndef FINITEVOLUMESOLVER_HPP
#define FINITEVOLUMESOLVER_HPP

#include <memory>

#include "Solver.hpp"
#include "bc/BoundaryManager.hpp"
#include "config/Settings.hpp"
#include "data/DataLayer.hpp"
#include "data/Variables.hpp"
#include "spatial/SpatialOperator.hpp"
#include "time/TimeIntegrator.hpp"

/**
 * @class FiniteVolumeSolver
 * @brief Generic finite-volume solver using pluggable spatial and time integrators.
 *
 * This class implements the high-level explicit FV time-stepping pipeline:
 *
 *  1. Apply boundary conditions through BoundaryManager.
 *  2. Compute spatial step size dx for a uniform 1D grid.
 *  3. Compute a CFL-limited time step dt via TimeStepCalculator.
 *  4. Optionally clamp dt to avoid overshooting t_end.
 *  5. Delegate the actual update of the conservative variables to a
 *     TimeIntegrator, which in turn repeatedly invokes a SpatialOperator
 *     to compute the semi-discrete RHS dU/dt.
 *
 * The FiniteVolumeSolver does not prescribe the numerical scheme:
 * - First-order Godunov corresponds to:
 *       GodunovSpatialOperator + ForwardEulerTimeIntegrator
 * - Second-order MUSCL–Hancock corresponds to:
 *       MusclHancockSpatialOperator + ForwardEulerTimeIntegrator
 * - MacCormack can be implemented as:
 *       any suitable SpatialOperator + MacCormackTimeIntegrator
 *
 * AnalyticalSolver remains a separate Solver and is not forced
 * into this architecture.
 */
class FiniteVolumeSolver : public Solver {
public:
    /**
     * @brief Constructs a finite-volume solver with given settings and components.
     *
     * The constructor:
     *  - stores a local copy of Settings,
     *  - initializes BoundaryManager with the configured dimension,
     *  - stores shared pointers to a SpatialOperator and a TimeIntegrator,
     *  - initializes positivity-limiting thresholds.
     *
     * @param settings         Global simulation settings (copied internally).
     * @param spatial_operator Spatial discretization operator (Godunov, MUSCL, ...).
     */
    FiniteVolumeSolver(const Settings& settings,
                       std::shared_ptr<SpatialOperator> spatial_operator);

    /**
     * @brief Advances the solution by one explicit finite-volume time step.
     *
     * The method:
     *  - applies boundary conditions for all axes,
     *  - computes dx,
     *  - computes a stable dt from CFL (and optionally clamps to t_end),
     *  - invokes the TimeIntegrator to update the conservative variables
     *    stored in DataLayer in-place,
     *  - updates the simulation time t_cur.
     *
     * If the physical domain has fewer than two cells or if dt <= 0,
     * no update is performed and 0.0 is returned.
     *
     * @param layer DataLayer containing the current state (modified in-place).
     * @param t_cur Current simulation time (incremented by dt on return).
     * @return Actual time step taken dt, or 0.0 if no step is performed.
     */
    auto Step(DataLayer& layer, double& t_cur) -> double override;

    /**
     * @brief Sets the CFL number used for dt computation.
     *
     * @param cfl CFL number (0 < cfl <= 1).
     */
    void SetCfl(double cfl) override;

    /**
     * @brief Assigns boundary conditions for a specific axis.
     *
     * For 1D problems axis = 0 is expected.
     *
     * @param axis     Spatial axis index.
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
    BoundaryManager boundary_manager_;

    /** @brief Spatial semi-discrete operator dU/dt = L(U). */
    std::shared_ptr<SpatialOperator> spatial_operator_;

    /** @brief Time integration scheme operating on the semi-discrete system. */
    std::shared_ptr<TimeIntegrator> time_integrator_;

    /** @brief Minimal allowed density for positivity limiting. */
    double rho_min_;

    /** @brief Minimal allowed pressure for positivity limiting. */
    double p_min_;

    /**
     * @brief Computes spatial step size dx for a uniform 1D grid.
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
     * @brief Selects and constructs the time integrator.
     *
     * Based on the string in `settings_.time_integrator` (case-insensitive):
     *  - contains "euler"       → ForwardEulerTimeIntegrator
     *  - contains "ssprk2"      → SSPRK2TimeIntegrator
     *  - contains "ssprk3"      → SSPRK3TimeIntegrator
     *  - contains "maccormack"  → MacCormackTimeIntegrator (works only if solver set to be MacCormack)
     *  - otherwise              → defaults to ForwardEulerTimeIntegrator
     */
    void InitializeTimeIntegrator();
};

#endif  // FINITEVOLUMESOLVER_HPP