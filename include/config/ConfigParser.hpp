#ifndef CONFIGPARSER_HPP
#define CONFIGPARSER_HPP

#include <yaml-cpp/yaml.h>

#include <map>
#include <optional>
#include <string>
#include <vector>

#include "config/InitialConditions.hpp"
#include "config/Settings.hpp"

/**
 * @file ConfigParser.hpp
 * @brief YAML configuration file parser and command-line argument handler
 */

/**
 * @class ConfigParser
 * @brief Parses configuration from YAML files and command-line arguments
 * 
 * ConfigParser is responsible for:
 * - Loading simulation settings from YAML configuration files
 * - Parsing command-line arguments to override file settings
 * - Managing multiple initial condition definitions
 * - Validating configuration consistency
 * 
 * The parser supports a two-stage initialization:
 * 1. Parse YAML configuration file
 * 2. Override settings with command-line arguments
 * 
 * Usage pattern:
 * @code
 * ConfigParser parser;
 * auto result = parser.Parse("config.yaml", argc, argv);
 * if (!result.has_value()) {
 *     return 0;  // Help was displayed
 * }
 * if (!result.value()) {
 *     return 1;  // Parsing error
 * }
 * const Settings& settings = parser.GetSettings();
 * @endcode
 * 
 * @note Uses yaml-cpp library for YAML parsing
 * @note Uses cxxopts library for command-line parsing
 * @see Settings, InitialConditions
 */
class ConfigParser {
public:
    /**
     * @brief Default constructor
     * 
     * Initializes an empty parser with default settings.
     */
    ConfigParser();

    /**
     * @brief Parse configuration from a YAML file
     * 
     * Loads all settings and initial conditions from the specified YAML file.
     * The file should contain:
     * - config/settings section with solver parameters
     * - config/initial_conditions section with Riemann problem states
     * 
     * Expected YAML structure:
     * @code{.yaml}
     * config:
     *   settings:
     *     solver: godunov
     *     N: 1000
     *     # ... other settings
     *   initial_conditions:
     *     sod1:
     *       rho_L: 1.0
     *       # ... state definition
     * @endcode
     * 
     * @param filename Path to YAML configuration file
     * @return true if parsing succeeded, false on error
     * @note Prints error messages to std::cerr on failure
     */
    auto ParseFile(const std::string& filename) -> bool;

    /**
     * @brief Parse command-line arguments
     * 
     * Parses command-line arguments and overrides current settings.
     * Supports all settings parameters and provides help functionality.
     * 
     * Special return values:
     * - std::nullopt: User requested --help (graceful exit)
     * - false: Parsing error occurred
     * - true: Parsing succeeded
     * 
     * Command-line options include:
     * - Solver options: --solver, --riemann-solver, --reconstruction
     * - Grid options: --N-cells, --cfl, --padding
     * - Physics options: --gamma, --t-end
     * - Initial conditions: --simulation-case, --x0
     * - Output options: --output-dir, --output-format
     * 
     * @param argc Argument count from main()
     * @param argv Argument vector from main()
     * @return Optional bool indicating success/failure/help
     * @note Use --help to display all available options
     */
    auto ParseCommandLine(int argc, char* argv[]) -> std::optional<bool>;

    /**
     * @brief Parse configuration from file and command-line arguments
     * 
     * Convenience method that combines file and command-line parsing:
     * 1. Check for --config argument to override default file path
     * 2. Parse the configuration file
     * 3. Override file settings with command-line arguments
     * 
     * This is the recommended entry point for configuration parsing.
     * 
     * @param default_config Default configuration file path
     * @param argc Argument count from main()
     * @param argv Argument vector from main()
     * @return Optional bool indicating success/failure/help
     * 
     * @see ParseFile(), ParseCommandLine()
     */
    auto Parse(const std::string& default_config, int argc, char* argv[])
        -> std::optional<bool>;

    /**
     * @brief Get parsed settings structure
     * @return Const reference to Settings object
     * @note Must call Parse() successfully before accessing
     */
    [[nodiscard]] auto GetSettings() const -> const Settings&;
    
    /**
     * @brief Get all initial condition definitions
     * @return Map from case name to InitialConditions
     * @note Keys match the simulation_case names in YAML
     */
    [[nodiscard]] auto GetInitialConditions() const 
        -> const std::map<std::string, InitialConditions>&;
    
    /**
     * @brief Get specific initial condition by name
     * @param case_name Name of the simulation case
     * @return Const reference to InitialConditions
     * @throw std::out_of_range if case_name not found
     */
    [[nodiscard]] auto GetInitialCondition(const std::string& case_name) const 
        -> const InitialConditions&;
    
    /**
     * @brief Get end time for a specific case
     * 
     * Each initial condition can specify its own t_end value in the
     * YAML configuration, which overrides the global settings.t_end
     * for that specific case.
     * 
     * @param case_name Name of the simulation case
     * @return End time for the specified case
     * @throw std::out_of_range if case_name not found
     */
    [[nodiscard]] auto GetCaseEndTime(const std::string& case_name) const 
        -> const double;
    
    /**
     * @brief Get list of all available case names
     * @return Sorted vector of case names
     * @note Useful for "all" mode that runs every case
     */
    [[nodiscard]] auto GetAllCaseNames() const -> std::vector<std::string>;
    
    /**
     * @brief Check if a specific case exists
     * @param case_name Name to check
     * @return true if the case is defined in configuration
     */
    [[nodiscard]] auto HasInitialCondition(const std::string& case_name) const 
        -> bool;

    /**
     * @brief Get the path to the loaded configuration file
     * @return Configuration file path as string
     */
    [[nodiscard]] auto GetConfigPath() const -> const std::string&;

private:
    /** @brief Parsed settings structure */
    Settings settings_;
    
    /** @brief Map of initial condition definitions */
    std::map<std::string, InitialConditions> initial_conditions_;
    
    /** @brief Map of case-specific end times */
    std::map<std::string, double> end_times_;
    
    /** @brief Path to the loaded configuration file */
    std::string config_path_;

    /**
     * @brief Load settings from YAML node
     * @param node YAML node containing settings
     * @param settings Output settings structure
     */
    static void LoadSettings(const YAML::Node& node, Settings& settings);
    
    /**
     * @brief Load initial conditions from YAML node
     * @param node YAML node containing initial_conditions
     * @param initial_conditions Output map of initial conditions
     */
    static void LoadInitialConditions(
        const YAML::Node& node,
        std::map<std::string, InitialConditions>& initial_conditions);
    
    /**
     * @brief Load case-specific end times from YAML node
     * @param node YAML node containing initial_conditions with t_end
     * @param end_times Output map of end times
     */
    static void LoadEndTimes(
        const YAML::Node& node,
        std::map<std::string, double>& end_times);
};

#endif  // CONFIGPARSER_HPP
