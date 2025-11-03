#include "Simulation.hpp"
#include "solver/GodunovSolver.hpp"
#include "bc/BoundaryManager.hpp"
#include "bc/PeriodicBoundary.hpp"
#include "visualization/VTKWriter.hpp"
#include <utility>

Simulation::Simulation(Config config)
    : config_(std::move(config)) {
}

void Simulation::Initialize() {
    layer_ = std::make_unique<DataLayer>(config_.N, config_.padding, config_.dim);


    solver_ = std::make_unique<GodunovSolver>();
    solver_->SetCfl(config_.CFL);


    auto periodic = std::make_shared<PeriodicBoundary>();
    solver_->AddBoundary(0, periodic, periodic);


    writer_ = std::make_unique<VTKWriter>(config_.output_dir);
}

void Simulation::WriteInitialState() const {
    if (writer_) writer_->Write(*layer_, /*step*/0, /*time*/0.0);
}

auto Simulation::ShouldWrite(std::size_t s) const -> bool {
    if (config_.output_every == 0) return false;
    return (s % config_.output_every) == 0;
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

    while (time_ < config_.t_end) {
        solver_->Step(*layer_, time_, config_.t_end);
        ++step_;
        WriteStepState(step_, time_);
        time_ += 0.01;
    }
}
