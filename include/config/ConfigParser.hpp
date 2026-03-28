#ifndef CONFIGPARSER_HPP
#define CONFIGPARSER_HPP

#include <map>
#include <optional>
#include <string>
#include <vector>

#include "config/CliParser.hpp"
#include "config/InitialConditions.hpp"
#include "config/Settings.hpp"
#include "config/YamlConfigParser.hpp"

class ConfigParser {
public:
    ConfigParser() = default;

    /**
     * @brief Parse YAML config and CLI arguments.
     * @param default_config Default config path used when CLI does not override it.
     * @param argc CLI argument count.
     * @param argv CLI argument vector.
     * @return std::nullopt when help was printed or parsing failed, otherwise true.
     */
    [[nodiscard]] std::optional<bool> Parse(const std::string& default_config,
                                            int argc,
                                            char* argv[]);

    [[nodiscard]] const Settings& GetDefaults() const;
    [[nodiscard]] const std::vector<std::string>& GetRunCases() const;
    [[nodiscard]] const std::map<std::string, InitialConditions>& GetInitialConditions() const;
    [[nodiscard]] const InitialConditions& GetInitialCondition(const std::string& case_name) const;
    [[nodiscard]] Settings GetCaseSettings(const std::string& case_name) const;
    [[nodiscard]] std::vector<std::string> GetAllCaseNames() const;
    [[nodiscard]] bool HasInitialCondition(const std::string& case_name) const;
    [[nodiscard]] const std::string& GetConfigPath() const;

private:
    Settings defaults_;
    std::vector<std::string> run_cases_;
    std::map<std::string, InitialConditions> initial_conditions_;
    std::string config_path_;
    CaseSettings cli_overrides_;
};

#endif  // CONFIGPARSER_HPP
