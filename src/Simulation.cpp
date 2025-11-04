#include "Simulation.hpp"
#include "solver/GodunovSolver.hpp"
#include "bc/BoundaryManager.hpp"
#include "bc/PeriodicBoundary.hpp"
#include "output/VTKWriter.hpp"
#include <utility>

Simulation::Simulation(Settings settings)
    : settings_(std::move(settings)) {
}

void Simulation::Initialize() {
    layer_ = std::make_unique<DataLayer>(settings_.N, settings_.padding, settings_.dim);

    solver_ = std::make_unique<GodunovSolver>();
    solver_->SetCfl(settings_.CFL);

    auto periodic = std::make_shared<PeriodicBoundary>();
    solver_->AddBoundary(0, periodic, periodic);

    writer_ = std::make_unique<VTKWriter>(settings_.output_dir);
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
    return time_;
}

auto Simulation::ShouldWrite(std::size_t s) const -> bool {
    if (settings_.output_every_steps == 0) return false;
    return (s % settings_.output_every_steps) == 0;
}

void Simulation::WriteInitialState() const {
    if (writer_) writer_->Write(*layer_, /*step*/0, /*time*/0.0);
}

void Simulation::WriteStepState(std::size_t s, double t) const {
    if (writer_ && ShouldWrite(s)) {
        writer_->Write(*layer_, s, t);
    }
}

void Simulation::Run() {
    Initialize();
    WriteInitialState();

    time_ = 0.0;
    step_ = 0;

    while (time_ < settings_.t_end) {
        solver_->Step(*layer_, time_, settings_.t_end);
        ++step_;
        WriteStepState(step_, time_);
        time_ += 0.01;
    }
}
