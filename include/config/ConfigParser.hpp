#ifndef CONFIGPARSER_HPP
#define CONFIGPARSER_HPP

#include <yaml-cpp/yaml.h>

#include <array>
#include <optional>
#include <string>

#include "config/InitialConditions.hpp"
#include "config/Settings.hpp"

class ConfigParser {
   public:
    ConfigParser();

    auto ParseFile(const std::string& filename) -> bool;

    auto ParseCommandLine(int argc, char* argv[]) -> std::optional<bool>;

    auto Parse(const std::string& default_config, int argc, char* argv[]) -> std::optional<bool>;

    [[nodiscard]] auto GetSettings() const -> const Settings&;
    [[nodiscard]] auto GetSODTests() const -> const std::array<InitialConditions, 5>&;
    [[nodiscard]] auto GetSODTest(int test_num) const -> const InitialConditions&;

    [[nodiscard]] auto GetConfigPath() const -> const std::string&;

   private:
    Settings settings_;
    std::array<InitialConditions, 5> sod_tests_;
    std::string config_path_;

    static void LoadSettings(const YAML::Node& node, Settings& settings);
    static void LoadInitialConditions(const YAML::Node& node,
                                      std::array<InitialConditions, 5>& sod_tests);
};

#endif  // CONFIGPARSER_HPP
