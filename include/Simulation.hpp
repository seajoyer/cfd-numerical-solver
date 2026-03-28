#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include <cstddef>
#include <memory>
#include <string>

#include "config/InitialConditions.hpp"
#include "config/Settings.hpp"

class BoundaryManager;
class DataLayer;
class Mesh;
class Solver;
class StepWriter;
class Workspace;

/**
 * @file Simulation.hpp
 * @brief Single-case simulation orchestrator for generic meshes.
 */

/**
 * @class Simulation
 * @brief Executes one simulation case from initialization to completion.
 *
 * Responsibilities:
 * - validate case settings
 * - build mesh
 * - allocate runtime state
 * - initialize solution fields
 * - initialize boundary handling
 * - initialize solver and writer
 * - run time loop
 */
class Simulation final {
public:
    explicit Simulation(Settings settings, InitialConditions initial_conditions);

    ~Simulation();

    /**
     * @brief Execute full simulation workflow.
     */
    void Run();

    /**
     * @brief Access main data layer.
     */
    [[nodiscard]] DataLayer& GetDataLayer();

    /**
     * @brief Current step index.
     */
    [[nodiscard]] std::size_t GetCurrentStep() const;

    /**
     * @brief Current physical time.
     */
    [[nodiscard]] double GetCurrentTime() const;

private:
    void Initialize();
    void ValidateConfiguration() const;
    void ValidateBoundaryCoverage() const;

    void BuildMesh();
    void AllocateState();
    void InitializeFields();
    void InitializeBoundaryConditions();
    void InitializeSolver();
    void InitializeWriter();

    [[nodiscard]] std::unique_ptr<Mesh> CreateMesh() const;
    [[nodiscard]] std::unique_ptr<Solver> CreateSolver();

    [[nodiscard]] bool IsKnownSolver(const std::string& solver) const;
    [[nodiscard]] bool IsKnownTimeIntegrator(const std::string& time_integrator) const;
    [[nodiscard]] bool IsKnownReconstruction(const std::string& reconstruction) const;
    [[nodiscard]] bool IsKnownRiemannSolver(const std::string& riemann_solver) const;
    [[nodiscard]] bool IsKnownOutputFormat(const std::string& format) const;

    [[nodiscard]] bool ShouldWrite() const;
    [[nodiscard]] bool ShouldLog() const;
    [[nodiscard]] bool ShouldRun() const;

    void WriteInitialState() const;
    void WriteStepState() const;
    void PrintLog() const;
    void FinalizeWriter();

    Settings settings_;
    InitialConditions initial_conditions_;

    std::unique_ptr<Mesh> mesh_;
    std::unique_ptr<DataLayer> layer_;
    std::unique_ptr<Workspace> workspace_;
    std::unique_ptr<Solver> solver_;
    std::unique_ptr<StepWriter> vtk_writer_;

    std::shared_ptr<BoundaryManager> boundary_manager_;

    double t_cur_ = 0.0;
    std::size_t step_cur_ = 0;
    double dt_ = 0.0;
};

#endif  // SIMULATION_HPP
