#include <exception>
#include <iostream>

#include "Simulation.hpp"
#include "config/ConfigParser.hpp"

auto main(int argc, char* argv[]) -> int {
    try {
        std::string config_path = "../config.yaml";
        std::cout << "Loading configuration from: " << config_path << '\n';

        auto parser = ConfigParser();
        parser.Parse(config_path);
        const Settings& settings = parser.GetSettings();

        int test_num = 1;
        InitialConditions initial_conditions = parser.GetSODTest(test_num);
        std::cout << "Using initial conditions: sod" << test_num << '\n';

        Simulation sim(settings, initial_conditions);
        sim.Run();

        return 0;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }
}
