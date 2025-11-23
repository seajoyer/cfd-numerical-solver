#ifndef CONFIGPARSER_HPP
#define CONFIGPARSER_HPP

#include <yaml-cpp/yaml.h>

#include <map>
#include <optional>
#include <string>
#include <vector>

#include "config/InitialConditions.hpp"
#include "config/Settings.hpp"

class ConfigParser {
   public:
    ConfigParser();

    auto ParseFile(const std::string& filename) -> bool;

    auto ParseCommandLine(int argc, char* argv[]) -> std::optional<bool>;

    auto Parse(const std::string& default_config, int argc, char* argv[]) -> std::optional<bool>;

    [[nodiscard]] auto GetSettings() const -> const Settings&;
    
    // New interface for general-purpose initial conditions
    [[nodiscard]] auto GetInitialConditions() const -> const std::map<std::string, InitialConditions>&;
    [[nodiscard]] auto GetInitialCondition(const std::string& case_name) const -> const InitialConditions&;
    [[nodiscard]] auto GetAllCaseNames() const -> std::vector<std::string>;
    [[nodiscard]] auto HasInitialCondition(const std::string& case_name) const -> bool;

    [[nodiscard]] auto GetConfigPath() const -> const std::string&;

   private:
    Settings settings_;
    std::map<std::string, InitialConditions> initial_conditions_;
    std::string config_path_;

    static void LoadSettings(const YAML::Node& node, Settings& settings);
    static void LoadInitialConditions(const YAML::Node& node,
                                      std::map<std::string, InitialConditions>& initial_conditions);
};

#endif  // CONFIGPARSER_HPP
