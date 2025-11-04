#include <config/ConfigParser.hpp>
#include "Simulation.hpp"
#include "config/Settings.hpp"

auto main() -> int {
    Settings cfg;
    cfg.N = 200;
    cfg.padding = 2;
    cfg.dim = 1;
    cfg.CFL = 0.5;
    cfg.t_end = 0.2;
    cfg.output_dir = "data/output";
    cfg.output_every_steps = 1;

    Simulation sim(cfg);
    sim.Run();
    return 0;

    ConfigParser().Parse("../config.yml", "soda1");
}
