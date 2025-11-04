#include "Simulation.hpp"
#include "solver/GodunovSolver.hpp"
#include "bc/BoundaryManager.hpp"
#include "bc/PeriodicBoundary.hpp"
#include "output/VTKWriter.hpp"
#include <utility>
#include <iostream>

Simulation::Simulation(Settings settings)
    : settings_(std::move(settings)) {
}

void Simulation::Initialize() {
    // Initialize data layer
    layer_ = std::make_unique<DataLayer>(settings_.N, settings_.padding, settings_.dim);

    // Initialize solver
    solver_ = std::make_unique<GodunovSolver>(settings_.dim);
    solver_->SetCfl(settings_.CFL);

    // Set up boundary conditions
    auto periodic = std::make_shared<PeriodicBoundary>();
    solver_->AddBoundary(0, periodic, periodic);

    // Initialize output writer
    writer_ = std::make_unique<VTKWriter>(settings_.output_dir);
    
    std::cout << "Simulation initialized:" << '\n';
    std::cout << "  Grid size (N): " << settings_.N << '\n';
    std::cout << "  Dimension: " << settings_.dim << '\n';
    std::cout << "  CFL: " << settings_.CFL << '\n';
    std::cout << "  End time: " << settings_.t_end << '\n';
    std::cout << "  Output directory: " << settings_.output_dir << '\n';
}

auto Simulation::GetDataLayer() -> DataLayer& {
    if (!layer_) {
        throw std::runtime_error("DataLayer is not initialized.");
    }
    return *layer_;
}

auto Simulation::GetCurrentStep() const -> std::size_t {
    return step_;
}

auto Simulation::GetCurrentTime() const -> double {
    return t_cur_;
}

auto Simulation::ShouldWrite(std::size_t step) const -> bool {
    if (settings_.output_every_steps == 0) return false;
    return (step % settings_.output_every_steps) == 0;
}

void Simulation::WriteInitialState() const {
    if (writer_) {
        writer_->Write(*layer_, 0, 0.0);
    }
}

void Simulation::WriteStepState(std::size_t step, double t_cur) const {
    if (writer_ && ShouldWrite(step)) {
        writer_->Write(*layer_, step, t_cur);
    }
}

void Simulation::Run() {
    Initialize();
    WriteInitialState();

    t_cur_ = 0.0;
    step_ = 0;

    std::cout << "\nStarting simulation..." << '\n';

    while (t_cur_ < settings_.t_end) {
        // Advance one time step and get the actual dt used
        double dt = solver_->Step(*layer_, t_cur_);
        
        t_cur_ += dt;
        ++step_;
        
        // Write output if needed
        WriteStepState(step_, t_cur_);
        
        // Progress output
        if (0 == step_ % settings_.output_every_steps) {
            std::cout << "Step " << step_ << ", time = " << t_cur_  << '\n';
        }
    }

    std::cout << "Simulation completed!" << '\n';
    std::cout << "  Final time: " << t_cur_ << '\n';
    std::cout << "  Total steps: " << step_ << '\n';
}
