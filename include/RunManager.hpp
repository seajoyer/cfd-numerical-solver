#ifndef RUNMANAGER_HPP
#define RUNMANAGER_HPP

#include <string>
#include <vector>

#include "Simulation.hpp"
#include "config/ConfigParser.hpp"
#include "config/InitialConditions.hpp"
#include "config/Settings.hpp"
// #include "parallel/MPIContext.hpp"

/**
 * @file RunManager.hpp
 * @brief Top-level manager for configuration loading and case execution
 */

/**
 * @class RunManager
 * @brief Manages full application run and launches simulation cases
 *
 * Responsibilities:
 * - Load configuration
 * - Resolve selected cases
 * - Create common run directory
 * - Launch simulation cases sequentially
 *
 * RunManager does not perform numerical calculations itself.
 * Each individual case is executed by Simulation.
 */
class RunManager {
public:
    RunManager() = default;

    /**
     * @brief Executes full program workflow
     * @param argc Command-line argument count
     * @param argv Command-line argument values
     * @return Process exit code
     */
    auto Run(int argc, char* argv[]) -> int;

private:
    /**
     * @brief Loads configuration from file and command line
     * @return false if execution should stop
     */
    auto LoadConfiguration(int argc, char* argv[]) -> bool;

    /**
     * @brief Resolves final list of cases to run
     * @return false on validation error
     */
    auto ResolveCasesToRun() -> bool;

    /**
     * @brief Creates timestamped directory for current run
     */
    void BuildRunDirectory();

    /**
     * @brief Runs all resolved simulation cases
     * @return Exit code
     */
    auto RunCases() -> int;

    /**
     * @brief Prints selected cases
     */
    void PrintSelectedCases() const;

    void ValidateRuntimeMode() const;

    bool is_root_ = true;
    ConfigParser parser_;
    std::vector<std::string> cases_to_run_;
    std::string run_dir_;
};

#endif  // RUNMANAGER_HPP
