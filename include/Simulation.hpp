#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include <memory>
#include <string>
#include "config/Config.hpp"
#include "data/DataLayer.hpp"
#include "solver/Solver.hpp"
#include "config/StepWriter.hpp"


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
 * @note Simulation owns all dynamically allocated components (Solver, StepWriter, DataLayer)
 *       and is responsible for their initialization and lifetime management.
 */
class Simulation {
public:
    /**
     * @brief Constructs the Simulation object with the specified configuration.
     *
     * Initializes internal parameters and stores the configuration structure,
     * but does not allocate memory or initialize physical data yet.
     *
     * @param config Configuration structure containing numerical parameters and file paths.
     */
    explicit Simulation(Config config);

    /**
     * @brief Runs the full simulation loop.
     *
     * Performs all necessary setup steps (grid initialization, solver construction, etc.),
     * then advances the physical fields in time until the end time is reached.
     *
     * @note This function represents the highest-level execution entry point.
     */
    void Run();

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
     * @brief Writes the initial state of the system before any time advancement.
     *
     * The output is handled by the StepWriter component.
     */
    void WriteInitialState() const;

    /**
     * @brief Writes simulation data for a specific time step.
     *
     * @param step Current time step index.
     * @param time Physical time corresponding to this step.
     *
     * @note The method checks output frequency using ShouldWrite() before writing.
     */
    void WriteStepState(std::size_t step, double time) const;

    /**
     * @brief Determines whether the current time step should be written to disk.
     *
     * @param step Current time step index.
     * @return true if the step should be written according to configuration.
     */
    [[nodiscard]] bool ShouldWrite(std::size_t step) const;

    Config config;
    std::unique_ptr<Solver> solver;
    std::unique_ptr<StepWriter> writer;
    std::unique_ptr<DataLayer> layer;

    double time = 0.0;
    std::size_t step = 0;
};

#endif // SIMULATION_HPP
