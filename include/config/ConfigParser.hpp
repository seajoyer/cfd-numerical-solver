#ifndef CONFIGPARSER_HPP
#define CONFIGPARSER_HPP

#include <yaml-cpp/yaml.h>
#include <array>
#include "config/InitialConditions.hpp"
#include "config/Settings.hpp"

class ConfigParser {
public:
    ConfigParser();
    auto Parse(const std::string& filename) -> bool;

    [[nodiscard]] auto GetSettings() const -> const Settings&;
    [[nodiscard]] auto GetSODTests() const -> const std::array<InitialConditions, 5>&;
    [[nodiscard]] auto GetSODTest(int test_num) const -> const InitialConditions&;

private:
    Settings settings_;
    std::array<InitialConditions, 5> sod_tests_;

    static void LoadSettings(const YAML::Node& node, Settings& settings);
    static void LoadInitialConditions(const YAML::Node& node, std::array<InitialConditions, 5>& sod_tests);
};

#endif // CONFIGPARSER_HPP
