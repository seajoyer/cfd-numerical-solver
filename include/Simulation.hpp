#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include <cstddef>
#include <memory>
#include <string>

#include "bc/BoundaryManager.hpp"
#include "config/InitialConditions.hpp"
#include "config/Settings.hpp"
#include "data/DataLayer.hpp"
#include "data/Mesh.hpp"
#include "output/StepWriter.hpp"
#include "parallel/DomainDecomposition.hpp"
#include "parallel/MPIContext.hpp"
#include "solver/Solver.hpp"

/**
 * @file Simulation.hpp
 * @brief Single-case simulation orchestrator
 */

/**
 * @class Simulation
 * @brief Executes one simulation case from initialization to completion
 *
 * Responsibilities:
 * - build runtime objects for one case
 * - initialize mesh and data
 * - apply initial and boundary conditions
 * - run time loop
 * - write VTK output
 */
class Simulation {
public:
    explicit Simulation(Settings settings, const InitialConditions& initial_conditions);

    /**
     * @brief Executes full simulation workflow
     */
    void Run();

    /**
     * @brief Accesses main data layer
     * @return Reference to initialized data layer
     */
    auto GetDataLayer() -> DataLayer&;

    /**
     * @brief Returns current step index
     */
    [[nodiscard]] auto GetCurrentStep() const -> std::size_t;

    /**
     * @brief Returns current physical time
     */
    [[nodiscard]] auto GetCurrentTime() const -> double;

private:
    void Initialize();
    void InitializeParallel();
    void InitializeMesh();
    void InitializeCoordinates();
    void InitializeDataLayer();
    void InitializeGeometry();
    void InitializeBoundaryConditions();
    void InitializeSolver();
    void InitializeWriter();

    void ApplyInitialConditions(DataLayer& layer, Mesh& mesh);

    [[nodiscard]] auto DeterminePadding() const -> int;
    void ValidateConfiguration() const;

    [[nodiscard]] auto IsKnownSolver(const std::string& solver) const -> bool;
    [[nodiscard]] auto IsKnownTimeIntegrator(const std::string& time_integrator) const -> bool;
    [[nodiscard]] auto IsKnownReconstruction(const std::string& reconstruction) const -> bool;
    [[nodiscard]] auto IsKnownRiemannSolver(const std::string& riemann_solver) const -> bool;
    [[nodiscard]] auto IsKnownBoundaryCondition(const std::string& bc) const -> bool;
    [[nodiscard]] auto IsKnownOutputFormat(const std::string& format) const -> bool;

    auto CreateSolver() -> std::unique_ptr<Solver>;

    [[nodiscard]] auto ShouldWrite() const -> bool;
    [[nodiscard]] auto ShouldLog() const -> bool;
    [[nodiscard]] auto ShouldRun() const -> bool;

    void WriteInitialState() const;
    void WriteStepState(double t_cur, std::size_t step_cur) const;
    void PrintLog() const;
    void FinalizeWriter();

private:
    Settings settings_;
    InitialConditions initial_conditions_;

    std::unique_ptr<Solver> solver_;
    std::unique_ptr<StepWriter> vtk_writer_;

    std::unique_ptr<DataLayer> layer_;
    std::unique_ptr<Mesh> mesh_;

    std::unique_ptr<MPIContext> mpi_context_;
    std::unique_ptr<DomainDecomposition> decomposition_;

    std::shared_ptr<BoundaryManager> boundary_manager_;

    double t_cur_ = 0.0;
    std::size_t step_cur_ = 0;
    double dt_ = 1.0;

    std::string case_output_dir_;
};

#endif  // SIMULATION_HPP
