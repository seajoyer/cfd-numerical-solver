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
     * @brief Get list of cases to run
     * @return Vector of case names to execute
     * @note Can be specific case names or ["all"] for all cases
     */
    [[nodiscard]] auto GetRunCases() const -> const std::vector<std::string>&;
    
    /**
     * @brief Get all initial condition definitions
     * @return Map from case name to InitialConditions
     * @note Keys match the case names in YAML
     */
    [[nodiscard]] auto GetInitialConditions() const 
        -> const std::map<std::string, InitialConditions>&;
    
    /**
     * @brief Get specific initial condition (case to run) by name
     * @param case_name Name of the simulation case
     * @return Const reference to InitialConditions
     * @throw std::out_of_range if case_name not found
     */
    [[nodiscard]] auto GetInitialCondition(const std::string& case_name) const 
        -> const InitialConditions&;
    
    /**
     * @brief Get merged settings for a specific case
     * 
     * Applies overrides in priority order:
     * 1. Global settings (base)
     * 2. Case-specific overrides from YAML
     * 3. CLI overrides (highest priority)
     * 
     * This ensures CLI arguments always override config file values.
     * 
     * @param case_name Name of the simulation case
     * @return Merged settings (global + case + CLI overrides)
     * @throw std::out_of_range if case_name not found
     */
    [[nodiscard]] auto GetCaseSettings(const std::string& case_name) const 
        -> Settings;
    
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

    /**
     * @brief Print list of supported solvers and their reconstructions
     */
    void PrintSolversList() const;

    /**
     * @brief Print list of supported Riemann solvers
     */
    void PrintRiemannSolversList() const;

    /**
     * @brief Print list of supported boundary conditions
     */
    void PrintBoundaryConditionsList() const;

    /**
     * @brief Print list of all available simulation cases
     */
    void PrintCasesList() const;

    /**
     * @brief Print list of supported output formats
     */
    void PrintOutputFormatsList() const;

private:
    /** @brief Parsed global settings structure */
    Settings settings_;
    
    /** @brief List of cases to run (can include "all") */
    std::vector<std::string> run_cases_;
    
    /** @brief Map of initial condition definitions */
    std::map<std::string, InitialConditions> initial_conditions_;
    
    /** @brief CLI overrides (highest priority) */
    CaseSettings cli_overrides_;
    
    /** @brief Path to the loaded configuration file */
    std::string config_path_;

    /**
     * @brief Parse command line for --config and --help only
     */
    auto ParseCommandLineForConfigAndHelp(int argc, char* argv[]) 
        -> std::optional<bool>;

    /**
     * @brief Load settings from YAML node
     */
    static void LoadSettings(const YAML::Node& node, Settings& settings);
    
    /**
     * @brief Load initial conditions and case settings from YAML node
     */
    static void LoadCases(
        const YAML::Node& node,
        std::map<std::string, InitialConditions>& initial_conditions);
    
    /**
     * @brief Load case-specific setting overrides from YAML node
     */
    static void LoadCaseOverrides(const YAML::Node& node, CaseSettings& overrides);
    
    /**
     * @brief Load run_cases list from YAML node
     */
    static void LoadRunCases(const YAML::Node& node, std::vector<std::string>& run_cases);

    /**
     * @brief Parse output formats from YAML node (can be string or list)
     */
    static auto ParseOutputFormats(const YAML::Node& node) -> std::vector<std::string>;

    /**
     * @brief Parse comma-separated format string from CLI
     */
    static auto ParseOutputFormatsString(const std::string& formats_str) -> std::vector<std::string>;
};

#endif  // CONFIGPARSER_HPP
