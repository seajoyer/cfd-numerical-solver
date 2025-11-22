#include <exception>
#include <iostream>

#include "Simulation.hpp"
#include "config/ConfigParser.hpp"

auto main(int argc, char* argv[]) -> int {
    try {
        ConfigParser parser;

        auto result = parser.Parse("../config.yaml", argc, argv);

        if (!result.has_value()) {
            return 0;  // --help was requested, exit gracefully
        }
        if (!result.value()) {
            return 1;  // Parsing error occurred 
        }

        const Settings& settings = parser.GetSettings();
        std::cout << "Configuration loaded from: " << parser.GetConfigPath() << "\n\n";

        int test_num = settings.sod_test_num;
        InitialConditions ic = parser.GetSODTest(test_num);
        std::cout << "Using initial conditions: sod" << test_num << "\n";
        std::cout << "Grid cells: " << settings.N << ", CFL: " << settings.cfl << "\n\n";

        Simulation sim(settings, ic);
        sim.Run();

        return 0;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }
}
