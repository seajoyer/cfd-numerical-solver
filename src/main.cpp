#include "Simulation.hpp"
#include "config/Config.hpp"

int main() {
    Config cfg;
    cfg.n = 200;
    cfg.padding = 2;
    cfg.dim = 1;
    cfg.cfl = 0.5;
    cfg.tEnd = 0.2;
    cfg.outputDir = "data/output";
    cfg.outputEvery = 1;

    Simulation sim(cfg);
    sim.Run();
    return 0;
}