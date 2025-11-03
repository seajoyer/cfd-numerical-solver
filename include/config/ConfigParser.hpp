#pragma once
#include <yaml-cpp/yaml.h>
#include "data/InitialConditions.hpp"
#include "data/GlobalVariables.hpp"

class ConfigParser {
public:
    ConfigParser();
    auto Parse(const std::string& filename, const std::string& initial) -> bool;
    [[nodiscard]] auto GetInitialConditions() const -> const InitialConditions&;

private:
    InitialConditions initial_conditions_;
    GlobalVariables global_variables_;


    static void LoadInitialConditions(const YAML::Node& node, InitialConditions& conditions);
    static void LoadGlobalVariables(const YAML::Node& node, GlobalVariables& globals);
};
