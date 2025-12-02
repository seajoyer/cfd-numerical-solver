#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include <cstddef>
#include <memory>
#include <vector>

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
 * ## Output Organization
 * ```
 * output_dir/run_<timestamp>/
 *   └── <case_name>/
 *       ├── png/                              (if PNG format enabled)
 *       │   └── step_NNNN.png
 *       └── vtk/                              (if VTK format enabled)
 *           ├── analytical/                   (if analytical enabled)
 *           │   └── step_NNNN.vtk
 *           └── solver__R_recon__N_size__CFL_value/
 *               └── solver__...__step_NNNN.vtk
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
     * @param settings Configuration parameters (copied internally)
     * @param initial_conditions Initial state definition (copied internally)
     */
    explicit Simulation(Settings settings, const InitialConditions& initial_conditions);

    // ==================== Main Interface ====================
    
    /**
     * @brief Executes complete simulation workflow
     */
    void Run();

    // ==================== State Access ====================
    
    /**
     * @brief Accesses computational data layer
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
     */
    void Initialize();

    /**
     * @brief Creates solver instance based on settings
     * @return Unique pointer to created solver
     */
    auto CreateSolver() -> std::unique_ptr<Solver>;

    /**
     * @brief Creates boundary condition instance
     * @return Shared pointer to boundary condition
     */
    auto CreateBoundaryCondition(const std::string& boundary_type,
                                 double rho_inf, double u_inf, double p_inf)
        -> std::shared_ptr<BoundaryCondition>;

    /**
     * @brief Creates output writers for all enabled formats
     * 
     * Creates separate writers for each output format (vtk, png).
     * VTK format creates separate writers for numerical and analytical.
     * PNG format creates a single writer that handles both.
     */
    void CreateWriters();

    /**
     * @brief Applies initial conditions to data layer
     * @param layer Data layer to initialize
     */
    void ApplyInitialConditions(DataLayer& layer);

    // ==================== Output Control ====================
    
    /**
     * @brief Checks if current step should produce output
     * @return true if output should be written
     */
    [[nodiscard]] auto ShouldWrite() const -> bool;

    /**
     * @brief Checks if current step should print log
     * @return true if progress should be logged
     */
    [[nodiscard]] auto ShouldLog() const -> bool;

    /**
     * @brief Checks if time loop should continue
     * @return true if next step should be computed
     */
    [[nodiscard]] auto ShouldRun() const -> bool;

    /**
     * @brief Writes initial state (t=0) to disk
     */
    void WriteInitialState() const;

    /**
     * @brief Writes current state to disk if output interval reached
     * @param t_cur Current simulation time
     * @param step_cur Current step number
     */
    void WriteStepState(double t_cur, std::size_t step_cur) const;

    /**
     * @brief Prints progress log if logging interval reached
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
    
    /** @brief VTK writer for numerical solution */
    std::unique_ptr<StepWriter> vtk_writer_;
    
    /** @brief VTK writer for analytical solution */
    std::unique_ptr<StepWriter> vtk_analytical_writer_;
    
    /** @brief PNG writer (handles both numerical and analytical) */
    std::unique_ptr<StepWriter> png_writer_;
    
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
    
    /** @brief Base output directory for this case */
    std::string case_output_dir_;
};

#endif  // SIMULATION_HPP
