#include <vtkAutoInit.h>
VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkRenderingContextOpenGL2);

#include <exception>
#include <iostream>

#include "Simulation.hpp"
#include "config/ConfigParser.hpp"
#include "utils/StringUtils.hpp"

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

        const Settings& global_settings = parser.GetSettings();
        const std::vector<std::string>& run_cases = parser.GetRunCases();

        std::cout << "Configuration loaded from: " << parser.GetConfigPath() << "\n\n";

        // Generate timestamp for this run
        std::string timestamp = utils::GetTimestamp();
        std::string run_dir = global_settings.output_dir + "/run_" + timestamp;

        std::cout << "Run directory: " << run_dir << "\n\n";

        // Determine which cases to run
        std::vector<std::string> cases_to_run;

        bool run_all = false;
        for (const auto& case_name : run_cases) {
            if (case_name == "all") {
                run_all = true;
                break;
            }
        }

        if (run_all) {
            // Run all available initial conditions
            cases_to_run = parser.GetAllCaseNames();

            if (cases_to_run.empty()) {
                std::cerr << "Error: No cases defined in configuration file\n";
                return 1;
            }

            std::cout << "Running all " << cases_to_run.size() << " simulation cases:\n";
            for (const auto& name : cases_to_run) {
                std::cout << "  - " << name << "\n";
            }
            std::cout << "\n";
        } else {
            // Run only specified cases
            cases_to_run = run_cases;

            // Validate that all requested cases exist
            for (const auto& case_name : cases_to_run) {
                if (!parser.HasInitialCondition(case_name)) {
                    std::cerr << "Error: Case '" << case_name << "' not found\n";
                    std::cerr << "Available cases:\n";
                    for (const auto& name : parser.GetAllCaseNames()) {
                        std::cerr << "  - " << name << "\n";
                    }
                    return 1;
                }
            }

            std::cout << "Running " << cases_to_run.size() << " simulation case(s):\n";
            for (const auto& name : cases_to_run) {
                std::cout << "  - " << name << "\n";
            }
            std::cout << "\n";
        }

        // Run each case
        for (const auto& case_name : cases_to_run) {
            std::cout << "========================================\n";
            std::cout << "Starting simulation case: " << case_name << "\n";
            std::cout << "========================================\n\n";

            // Get merged settings (global + case-specific overrides)
            Settings case_settings = parser.GetCaseSettings(case_name);

            // Set the case name for output directory
            case_settings.simulation_case = case_name;
            // Use timestamped run directory structure
            case_settings.output_dir = run_dir + "/" + case_name;

            // Get initial conditions
            InitialConditions ic = parser.GetInitialCondition(case_name);

            Simulation sim(case_settings, ic);
            sim.Run();

            std::cout << "\n";
        }

        if (cases_to_run.size() > 1) {
            std::cout << "========================================\n";
            std::cout << "All simulations completed successfully!\n";
            std::cout << "========================================\n";
        }

        std::cout << "Results saved to: " << run_dir << "\n";

        return 0;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }
}
