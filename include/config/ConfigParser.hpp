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

    auto Parse(const std::string& default_config, int argc, char* argv[])
        -> std::optional<bool>;

    [[nodiscard]] auto GetSettings() const -> const Settings&;
    [[nodiscard]] auto GetRunCases() const -> const std::vector<std::string>&;
    [[nodiscard]] auto GetInitialConditions() const
        -> const std::map<std::string, InitialConditions>&;
    [[nodiscard]] auto GetInitialCondition(const std::string& case_name) const
        -> const InitialConditions&;
    [[nodiscard]] auto GetCaseSettings(const std::string& case_name) const -> Settings;
    [[nodiscard]] auto GetAllCaseNames() const -> std::vector<std::string>;
    [[nodiscard]] auto HasInitialCondition(const std::string& case_name) const -> bool;
    [[nodiscard]] auto GetConfigPath() const -> const std::string&;

private:
    Settings settings_;
    std::vector<std::string> run_cases_;
    std::map<std::string, InitialConditions> initial_conditions_;
    std::string config_path_;
};

#endif  // CONFIGPARSER_HPP