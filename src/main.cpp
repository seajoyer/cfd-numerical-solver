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

        // Check if user wants to run all cases
        const std::string& case_name = settings.simulation_case;
        
        if (case_name == "all") {
            // Run all available initial conditions
            auto case_names = parser.GetAllCaseNames();
            
            if (case_names.empty()) {
                std::cerr << "Error: No initial conditions defined in configuration file\n";
                return 1;
            }
            
            std::cout << "Running all " << case_names.size() << " simulation cases:\n";
            for (const auto& name : case_names) {
                std::cout << "  - " << name << "\n";
            }
            std::cout << "\n";
            
            for (const auto& name : case_names) {
                std::cout << "========================================\n";
                std::cout << "Starting simulation case: " << name << "\n";
                std::cout << "========================================\n\n";
                
                // Create modified settings with case-specific output directory
                Settings case_settings = settings;
                case_settings.simulation_case = name;
                case_settings.output_dir = settings.output_dir + "/" + name;
                case_settings.t_end = parser.GetCaseEndTime(name);

                InitialConditions ic = parser.GetInitialCondition(name);
                case_settings.x0 = ic.x0;
                
                Simulation sim(case_settings, ic);
                sim.Run();
                
                std::cout << "\n";
            }
            
            std::cout << "========================================\n";
            std::cout << "All simulations completed successfully!\n";
            std::cout << "========================================\n";
            
        } else {
            // Run single case
            if (!parser.HasInitialCondition(case_name)) {
                std::cerr << "Error: Initial condition '" << case_name << "' not found\n";
                std::cerr << "Available cases:\n";
                for (const auto& name : parser.GetAllCaseNames()) {
                    std::cerr << "  - " << name << "\n";
                }
                return 1;
            }
            
            InitialConditions ic = parser.GetInitialCondition(case_name);
            std::cout << "Using initial condition: " << case_name << "\n";
            std::cout << "Grid cells: " << settings.N << ", CFL: " << settings.cfl << "\n\n";

            Simulation sim(settings, ic);
            sim.Run();
        }

        return 0;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }
}
