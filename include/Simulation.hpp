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

class Simulation {
public:
    explicit Simulation(Settings settings, const InitialConditions& initial_conditions);

    void Run();

    auto GetDataLayer() -> DataLayer&;
    [[nodiscard]] auto GetCurrentStep() const -> std::size_t;
    [[nodiscard]] auto GetCurrentTime() const -> double;

private:
    void Initialize();
    auto CreateSolver() -> std::unique_ptr<Solver>;
    auto CreateBoundaryCondition(const std::string& boundary_type,
                                 double rho_inf, double u_inf, double p_inf)
        -> std::shared_ptr<BoundaryCondition>;
    auto CreateBoundaryCondition2D(const std::string& boundary_type,
                                   double rho_inf, double u_inf,
                                   double v_inf, double p_inf)
        -> std::shared_ptr<BoundaryCondition>;
    void CreateWriters();
    void ApplyInitialConditions(DataLayer& layer);
    void ApplyInitialConditions2D(DataLayer& layer);

    [[nodiscard]] auto ShouldWrite() const -> bool;
    [[nodiscard]] auto ShouldLog() const -> bool;
    [[nodiscard]] auto ShouldRun() const -> bool;
    void WriteInitialState() const;
    void WriteStepState(double t_cur, std::size_t step_cur) const;
    void PrintLog() const;
    void FinalizeWriters();

    Settings settings_;
    Settings analytical_settings_;
    InitialConditions initial_conditions_;
    std::unique_ptr<Solver> solver_;
    std::unique_ptr<AnalyticalSolver> analytical_solver_;
    std::unique_ptr<StepWriter> vtk_writer_;
    std::unique_ptr<StepWriter> vtk_analytical_writer_;
    std::unique_ptr<StepWriter> png_writer_;
    std::unique_ptr<StepWriter> gif_writer_;
    std::unique_ptr<DataLayer> layer_;
    std::unique_ptr<DataLayer> analytical_layer_;

    double t_cur_ = 0.0;
    std::size_t step_cur_ = 0;
    double dt_ = 1.0;
    std::string case_output_dir_;
};

#endif  // SIMULATION_HPP
