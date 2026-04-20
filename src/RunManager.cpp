#include "RunManager.hpp"

#include <iostream>
#include <stdexcept>

#include "utils/StringUtils.hpp"

void RunManager::ValidateRuntimeMode() const {
    // MPI mode is now allowed.
    // Additional runtime checks can be added here later if needed.
}

auto RunManager::Run(int argc, char* argv[]) -> int {
    if (!LoadConfiguration(argc, argv)) {
        return 0;
    }

    if (!ResolveCasesToRun()) {
        return 1;
    }

    ValidateRuntimeMode();

    BuildRunDirectory();

    if (is_root_) {
        std::cout << "Configuration loaded from: " << parser_.GetConfigPath() << "\n\n";
    }

    PrintSelectedCases();

    return RunCases();
}

auto RunManager::LoadConfiguration(int argc, char* argv[]) -> bool {
    mpi_ = std::make_unique<MPIContext>();
    is_root_ = mpi_->IsRoot();

    auto result = parser_.Parse("../config.yaml", argc, argv);

    if (!result.has_value()) {
        return false;
    }

    if (!result.value()) {
        throw std::runtime_error("Failed to parse configuration");
    }

    return true;
}

auto RunManager::ResolveCasesToRun() -> bool {
    const std::vector<std::string>& run_cases = parser_.GetRunCases();

    bool run_all = false;
    for (const auto& case_name : run_cases) {
        if (case_name == "all") {
            run_all = true;
            break;
        }
    }

    if (run_all) {
        cases_to_run_ = parser_.GetAllCaseNames();

        if (cases_to_run_.empty()) {
            if (is_root_) {
                std::cerr << "Error: no cases defined in configuration\n";
            }
            return false;
        }

        return true;
    }

    cases_to_run_ = run_cases;

    for (const auto& case_name : cases_to_run_) {
        if (!parser_.HasInitialCondition(case_name)) {
            if (is_root_) {
                std::cerr << "Error: case '" << case_name << "' not found\n";
                std::cerr << "Available cases:\n";
                for (const auto& available_case : parser_.GetAllCaseNames()) {
                    std::cerr << "  - " << available_case << '\n';
                }
            }
            return false;
        }
    }

    return true;
}

void RunManager::BuildRunDirectory() {
    if (cases_to_run_.empty()) {
        throw std::runtime_error("RunManager::BuildRunDirectory: no cases selected");
    }

    if (is_root_) {
        const Settings first_case_settings = parser_.GetCaseSettings(cases_to_run_.front());
        const std::string timestamp = utils::GetTimestamp();
        run_dir_ = first_case_settings.output_dir + "/run_" + timestamp;
        std::cout << "Run directory: " << run_dir_ << "\n\n";
    }

    if (mpi_ && mpi_->Size() > 1) {
        run_dir_ = mpi_->BroadcastString(run_dir_, 0);
    }
}

void RunManager::PrintSelectedCases() const {
    if (!is_root_) {
        return;
    }

    if (cases_to_run_.size() == 1) {
        std::cout << "Running 1 simulation case:\n";
    }
    else {
        std::cout << "Running " << cases_to_run_.size() << " simulation cases:\n";
    }

    for (const auto& case_name : cases_to_run_) {
        std::cout << "  - " << case_name << '\n';
    }

    std::cout << '\n';
}

auto RunManager::RunCases() -> int {
    for (const auto& case_name : cases_to_run_) {
        if (is_root_) {
            std::cout << "========================================\n";
            std::cout << "Starting simulation case: " << case_name << '\n';
            std::cout << "========================================\n\n";
        }

        Settings case_settings = parser_.GetCaseSettings(case_name);
        case_settings.simulation_case = case_name;
        case_settings.output_dir = run_dir_ + "/" + case_name;

        const InitialConditions& initial_conditions =
            parser_.GetInitialCondition(case_name);

        Simulation simulation(case_settings, initial_conditions);
        simulation.Run();
    }

    if (cases_to_run_.size() > 1 && is_root_) {
        std::cout << "========================================\n";
        std::cout << "All simulations completed successfully\n";
        std::cout << "========================================\n";
    }

    if (is_root_) {
        std::cout << "\nResults saved to: " << run_dir_ << '\n';
    }

    return 0;
}
