#ifndef CONFIGPARSER_HPP
#define CONFIGPARSER_HPP

#include <yaml-cpp/yaml.h>
#include "config/InitialConditions.hpp"
#include "config/Settings.hpp"

class ConfigParser {
public:
    ConfigParser();
    auto Parse(const std::string& filename, const std::string& initial) -> bool;

    [[nodiscard]] auto GetSettings() const -> const Settings&;
    [[nodiscard]] auto GetInitialConditions() const -> const InitialConditions&;

private:
    Settings settings_;
    InitialConditions initial_conditions_;

    static void LoadSettings(const YAML::Node& node, Settings& settings);
    static void LoadInitialConditions(const YAML::Node& node, InitialConditions& conditions);
};

#endif // CONFIGPARSER_HPP
