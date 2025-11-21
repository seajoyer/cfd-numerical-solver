#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include <cstddef>
#include <memory>

#include "config/InitialConditions.hpp"
#include "config/Settings.hpp"
#include "data/DataLayer.hpp"
#include "output/StepWriter.hpp"
#include "solver/Solver.hpp"
#include "solver/AnalyticalSolver.hpp"

/**
 * @class Simulation
 * @brief Central orchestrator that coordinates all major components of the solver.
 *
 * The Simulation class encapsulates the entire CFD workflow:
 *  - reading configuration parameters,
 *  - initializing data structures,
 *  - creating solver and boundary conditions,
 *  - controlling time integration,
 *  - periodically writing simulation results to disk.
 *
 * The main goal of this class is to minimize logic inside main()
 * and provide a single entry point for running computations.
 *
 * @note Simulation owns all dynamically allocated components (Solver, StepWriter,
 * DataLayer) and is responsible for their initialization and lifetime management.
 */
class Simulation {
   public:
    /**
     * @brief Constructs the Simulation object with the specified configuration.
     *
     * Initializes internal parameters and stores the configuration structure,
     * but does not allocate memory or initialize physical data yet.
     *
     * @param settings Configuration structure containing numerical parameters and file
     * paths.
     * @param initial_conditions Initial conditions for the simulation (optional).
     */
    explicit Simulation(Settings settings, const InitialConditions& initial_conditions);

    /**
     * @brief Runs the full simulation loop.
     *
     * Performs all necessary setup steps (grid initialization, solver construction,
     * etc.), then advances the physical fields in time until the end time is reached.
     *
     * @note This function represents the highest-level execution entry point.
     */
    void Run();

    /**
     * @brief Accesses the underlying data layer for read/write operations.
     *
     * @return Reference to the DataLayer object.
     */
    auto GetDataLayer() -> DataLayer&;

    /**
     * @brief Gets the current time step index.
     *
     * @return The current time step index (std::size_t).
     */
    [[nodiscard]] auto GetCurrentStep() const -> std::size_t;

    /**
     * @brief Gets the current physical time of the simulation.
     *
     * @return The current physical time (double).
     */
    [[nodiscard]] auto GetCurrentTime() const -> double;

   private:
    /**
     * @brief Initializes all simulation components before the first time step.
     *
     * Allocates the DataLayer, creates the Solver, assigns boundary conditions,
     * and prepares the StepWriter for output.
     *
     * @note Called automatically inside Run(), typically not used externally.
     */
    void Initialize();

    /**
     * @brief Creates the appropriate solver based on the settings.
     *
     * @return Unique pointer to the created solver.
     * @throws std::runtime_error if solver type is not recognized.
     */
    auto CreateSolver() -> std::unique_ptr<Solver>;

    /**
     * @brief Creates a boundary condition based on the boundary type string.
     *
     * @param boundary_type String identifier for the boundary type.
     * @param rho_inf External/inlet density
     * @param u_inf External/inlet velocity
     * @param p_inf External/inlet pressure
     * @return Shared pointer to the created boundary condition.
     * @throws std::runtime_error if boundary type is not recognized.
     */
    auto CreateBoundaryCondition(const std::string& boundary_type, double rho_inf,
                                 double u_inf, double p_inf)
        -> std::shared_ptr<BoundaryCondition>;

    /**
     * @brief Creates an output writer based on the configured format and directory.
     *
     * Delegates to WriterFactory for instantiation.
     *
     * @param output_format String identifier for the output format (e.g., "vtk").
     * @param output_dir Directory path for output files.
     * @return Unique pointer to the created StepWriter.
     * @throws std::runtime_error if format is unsupported (propagated from factory).
     */
    auto CreateWriter(const std::string& output_format, const std::string& output_dir)
        -> std::unique_ptr<StepWriter>;

    /**
     * @brief Applies initial conditions to the data layer.
     *
     * Sets up the initial state of the physical fields based on the
     * configured initial conditions (e.g., Sod shock tube problems).
     */
    void ApplyInitialConditions(DataLayer& layer);

    /**
     * @brief Determines whether the current time step should be written to disk.
     *
     * @return true if the step should be written according to configuration.
     */
    [[nodiscard]] auto ShouldWrite() const -> bool;

    /**
     * @brief Determines whether the current time step should be logged.
     *
     * @return true if the step should be logged.
     */
    [[nodiscard]] auto ShouldLog() const -> bool;

    /**
     * @brief Determines whether the simulation shoult calculate the next step.
     *
     * @return true if the next step should be calculated.
     */
    [[nodiscard]] auto ShouldRun() const -> bool;

    /**
     * @brief Writes the initial state of the system before any time advancement.
     *
     * The output is handled by the StepWriter component.
     */
    void WriteInitialState() const;

    /**
     * @brief Writes simulation data for a specific time step.
     *
     * @param t_cur Physical time corresponding to this step.
     * @param step_cur Current time step index.
     *
     * @note The method checks output frequency using ShouldWrite() before writing.
     */
    void WriteStepState(double t_cur, std::size_t step_cur) const;

    /**
     * @brief Prints simulation log for a specific time step.
     *
     * @note The method checks output frequency using ShouldLog() before writing.
     */
    void PrintLog() const;

    Settings settings_;
    Settings analytical_settings_;
    InitialConditions initial_conditions_;
    bool log_progress_;
    std::unique_ptr<Solver> solver_;
    std::unique_ptr<AnalyticalSolver> analytical_solver_;
    std::unique_ptr<StepWriter> writer_;
    std::unique_ptr<StepWriter> analytical_writer_;
    std::unique_ptr<DataLayer> layer_;
    std::unique_ptr<DataLayer> analytical_layer_;

    double t_cur_ = 0.0;
    std::size_t step_cur_ = 0;
    double dt_ = 1.0;
};

#endif  // SIMULATION_HPP
