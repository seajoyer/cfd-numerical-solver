#include "Simulation.hpp"
#include "solver/GodunovSolver.hpp"
#include "bc/BoundaryManager.hpp"
#include "bc/PeriodicBoundary.hpp"
#include "visualization/VTKWriter.hpp"
#include <utility>

Simulation::Simulation(Config config)
    : config(std::move(config)) {
}

void Simulation::Initialize() {
    layer = std::make_unique<DataLayer>(config.n, config.padding, config.dim);


    solver = std::make_unique<GodunovSolver>();
    solver->SetCfl(config.cfl);


    auto periodic = std::make_shared<PeriodicBoundary>();
    solver->AddBoundary(0, periodic, periodic);


    writer = std::make_unique<VtkWriter>(config.outputDir);
}

void Simulation::WriteInitialState() const {
    if (writer) writer->Write(*layer, /*step*/0, /*time*/0.0);
}

bool Simulation::ShouldWrite(std::size_t s) const {
    if (config.outputEvery == 0) return false;
    return (s % config.outputEvery) == 0;
}

void Simulation::WriteStepState(std::size_t s, double t) const {
    if (writer && ShouldWrite(s)) {
        writer->Write(*layer, s, t);
    }
}

void Simulation::Run() {
    Initialize();
    WriteInitialState();

    time = 0.0;
    step = 0;

    while (time < config.tEnd) {
        solver->Step(*layer, time, config.tEnd);
        ++step;
        WriteStepState(step, time);
        time += 0.01;
    }
}
