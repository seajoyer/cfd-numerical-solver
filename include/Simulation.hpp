#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include <cstddef>
#include <memory>

#include "config/InitialConditions.hpp"
#include "config/Settings.hpp"
#include "data/DataLayer.hpp"
#include "output/StepWriter.hpp"
#include "solver/AnalyticalSolver.hpp"
#include "solver/Solver.hpp"

/**
 * @file Simulation.hpp
 * @brief Main simulation orchestrator and workflow manager
 */

/**
 * @class Simulation
 * @brief Orchestrates CFD simulation workflow from initialization to completion
 * 
 * The Simulation class is the top-level controller that coordinates all components
 * of the finite-volume CFD solver:
 * 
 * ## Responsibilities
 * - Initializes computational grid and data structures
 * - Creates and configures solver based on settings
 * - Applies initial and boundary conditions
 * - Manages time integration loop
 * - Controls output writing at specified intervals
 * - Optionally runs analytical solutions for verification
 * 
 * ## Workflow
 * 1. **Construction**: Store settings and initial conditions
 * 2. **Initialization** (via Initialize()):
 *    - Allocate DataLayer with ghost cells
 *    - Create solver (Godunov, MUSCL-Hancock, etc.)
 *    - Create boundary conditions
 *    - Create output writer (VTK)
 *    - Apply initial conditions
 *    - Optionally initialize analytical solver
 * 3. **Time loop** (via Run()):
 *    - Check stopping criteria (t_end, step_end)
 *    - Advance solution one time step
 *    - Update analytical solution (if enabled)
 *    - Write output files at specified intervals
 *    - Log progress to console
 * 4. **Completion**: Final output and statistics
 * 
 * ## Usage Example
 * ```cpp
 * Settings settings = \/* loaded from config *\/;
 * InitialConditions ic = \/* loaded from config *\/;
 * 
 * Simulation sim(settings, ic);
 * sim.Run();  // Complete simulation
 * 
 * // Access results
 * const DataLayer& data = sim.GetDataLayer();
 * double final_time = sim.GetCurrentTime();
 * ```
 * 
 * ## Output Organization
 * Numerical solutions:
 * ```
 * output_dir/solver__R_reconstruction__N_gridsize__CFL_value/
 *   ├── step_0000.vtk
 *   ├── step_0050.vtk
 *   └── ...
 * ```
 * 
 * Analytical solutions (if enabled):
 * ```
 * output_dir/analytical/
 *   ├── step_0000.vtk
 *   ├── step_0050.vtk
 *   └── ...
 * ```
 * 
 * @note Simulation owns all dynamically allocated components
 * @note Not copyable or movable (contains unique_ptr members)
 * @see Settings, Solver, DataLayer, StepWriter
 */
class Simulation {
public:
    // ==================== Construction ====================
    
    /**
     * @brief Constructs simulation with specified configuration
     * 
     * Stores settings and initial conditions but does not allocate
     * memory or initialize components. Call Run() to execute.
     * 
     * @param settings Configuration parameters (copied internally)
     * @param initial_conditions Initial state definition (copied internally)
     * @note Actual initialization deferred until Run() is called
     */
    explicit Simulation(Settings settings, const InitialConditions& initial_conditions);

    // ==================== Main Interface ====================
    
    /**
     * @brief Executes complete simulation workflow
     * 
     * Performs the following sequence:
     * 1. Initialize() - Set up all components
     * 2. WriteInitialState() - Output t=0 state
     * 3. Time loop until stopping criterion met:
     *    - Solver::Step() advances solution
     *    - Update analytical solution if enabled
     *    - WriteStepState() at output intervals
     *    - PrintLog() at logging intervals
     * 4. Print completion statistics
     * 
     * Stopping criteria (first satisfied):
     * - t_cur >= settings.t_end (if t_end > 0)
     * - step_cur >= settings.step_end (if step_end > 0)
     * 
     * @note This is the primary entry point for running simulations
     */
    void Run();

    // ==================== State Access ====================
    
    /**
     * @brief Accesses computational data layer
     * 
     * Provides read/write access to all simulation fields (ρ, u, P, etc.)
     * including ghost cells.
     * 
     * @return Reference to DataLayer
     * @throw std::runtime_error if called before initialization
     */
    auto GetDataLayer() -> DataLayer&;

    /**
     * @brief Returns current time step index
     * @return Step number (0 = initial state)
     */
    [[nodiscard]] auto GetCurrentStep() const -> std::size_t;

    /**
     * @brief Returns current physical simulation time
     * @return Time in seconds (or dimensionless)
     */
    [[nodiscard]] auto GetCurrentTime() const -> double;

private:
    // ==================== Initialization ====================
    
    /**
     * @brief Initializes all simulation components
     * 
     * Called once at start of Run(). Performs:
     * - DataLayer allocation (N + 2·padding cells)
     * - Solver creation via factory
     * - Boundary condition setup
     * - StepWriter creation (VTK)
     * - Initial condition application
     * - Analytical solver initialization (if enabled)
     * 
     * Prints configuration summary to stdout.
     * 
     * @note Automatically called by Run()
     */
    void Initialize();

    /**
     * @brief Creates solver instance based on settings
     * 
     * Uses SolverFactory to instantiate appropriate solver type:
     * - "godunov" → GodunovSolver
     * - "godunov-kolgan" → GodunovSolver with P1
     * - "godunov-kolgan-rodionov" → GodunovKolganRodionovSolver
     * - "analytical" → AnalyticalSolver
     * 
     * @return Unique pointer to created solver
     * @throw std::runtime_error if solver type unknown
     * @see SolverFactory
     */
    auto CreateSolver() -> std::unique_ptr<Solver>;

    /**
     * @brief Creates boundary condition instance
     * 
     * Uses BoundaryFactory to instantiate based on type string.
     * For inlet/free_stream types, external state (rho_inf, u_inf, p_inf)
     * is applied from initial conditions.
     * 
     * @param boundary_type Type identifier ("inlet", "outlet", etc.)
     * @param rho_inf External density
     * @param u_inf External velocity
     * @param p_inf External pressure
     * @return Shared pointer to boundary condition
     * @throw std::runtime_error if type unknown
     * @see BoundaryFactory, BoundaryCondition
     */
    auto CreateBoundaryCondition(const std::string& boundary_type,
                                 double rho_inf, double u_inf, double p_inf)
        -> std::shared_ptr<BoundaryCondition>;

    /**
     * @brief Creates output writer instance
     * 
     * Uses WriterFactory to instantiate appropriate writer.
     * Currently only "vtk" format supported.
     * 
     * @param output_format Format identifier ("vtk")
     * @param output_dir Directory for output files
     * @return Unique pointer to writer
     * @throw std::runtime_error if format unsupported
     * @see WriterFactory, VTKWriter
     */
    auto CreateWriter(const std::string& output_format,
                      const std::string& output_dir)
        -> std::unique_ptr<StepWriter>;

    /**
     * @brief Applies initial conditions to data layer
     * 
     * Sets up Riemann problem initial state:
     * - Left state (x < x0): (ρ_L, u_L, P_L)
     * - Right state (x ≥ x0): (ρ_R, u_R, P_R)
     * 
     * Also computes derived quantities:
     * - Cell mass m = ρ·Δx
     * - Momentum p = ρ·u
     * - Energies (internal, kinetic, total)
     * - Specific volume V = 1/ρ
     * 
     * @param layer Data layer to initialize
     * @note Only initializes core cells; ghosts filled by first BC application
     */
    void ApplyInitialConditions(DataLayer& layer);

    // ==================== Output Control ====================
    
    /**
     * @brief Checks if current step should produce output
     * 
     * Returns true if either condition met:
     * - Step-based: step_cur % output_every_steps == 0
     * - Time-based: crossed output_every_time boundary
     * - Always true at t ≥ t_end (final state)
     * 
     * @return true if output should be written
     */
    [[nodiscard]] auto ShouldWrite() const -> bool;

    /**
     * @brief Checks if current step should print log
     * 
     * Similar logic to ShouldWrite() but for console logging.
     * 
     * @return true if progress should be logged
     */
    [[nodiscard]] auto ShouldLog() const -> bool;

    /**
     * @brief Checks if time loop should continue
     * 
     * Returns false if any stopping criterion met:
     * - t_cur >= t_end (if t_end > 0)
     * - step_cur >= step_end (if step_end > 0)
     * 
     * @return true if next step should be computed
     */
    [[nodiscard]] auto ShouldRun() const -> bool;

    /**
     * @brief Writes initial state (t=0) to disk
     * 
     * Outputs both numerical and analytical solutions (if enabled)
     * as step 0 files.
     */
    void WriteInitialState() const;

    /**
     * @brief Writes current state to disk if output interval reached
     * 
     * @param t_cur Current simulation time
     * @param step_cur Current step number
     */
    void WriteStepState(double t_cur, std::size_t step_cur) const;

    /**
     * @brief Prints progress log if logging interval reached
     * 
     * Format: "[PROGRESS]: Step N, X% processed, time: t of t_end"
     * Uses carriage return for in-place console updates.
     */
    void PrintLog() const;

    // ==================== Member Data ====================
    
    /** @brief Numerical solver configuration */
    Settings settings_;
    
    /** @brief Analytical solver configuration (modified copy of settings_) */
    Settings analytical_settings_;
    
    /** @brief Initial condition definition */
    InitialConditions initial_conditions_;
    
    /** @brief Numerical solver instance */
    std::unique_ptr<Solver> solver_;
    
    /** @brief Analytical solver instance (if enabled) */
    std::unique_ptr<AnalyticalSolver> analytical_solver_;
    
    /** @brief Output writer for numerical solution */
    std::unique_ptr<StepWriter> writer_;
    
    /** @brief Output writer for analytical solution */
    std::unique_ptr<StepWriter> analytical_writer_;
    
    /** @brief Numerical solution data */
    std::unique_ptr<DataLayer> layer_;
    
    /** @brief Analytical solution data (if enabled) */
    std::unique_ptr<DataLayer> analytical_layer_;

    /** @brief Current simulation time */
    double t_cur_ = 0.0;
    
    /** @brief Current time step index */
    std::size_t step_cur_ = 0;
    
    /** @brief Last computed time step size */
    double dt_ = 1.0;
};

#endif  // SIMULATION_HPP
