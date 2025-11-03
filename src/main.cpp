#include <stdio.h>
#include <config/ConfigParser.hpp>
#include "Simulation.hpp"
#include "config/Config.hpp"

int main() {
    Config cfg;
    cfg.n = 200;
    cfg.padding = 2;
    cfg.dim = 1;
    cfg.cfl = 0.5;
    cfg.t_end = 0.2;
    cfg.output_dir = "data/output";
    cfg.output_every = 1;

    Simulation sim(cfg);
    sim.Run();
    return 0;

    ConfigParser().Parse("../config.yml", "soda1");
}
