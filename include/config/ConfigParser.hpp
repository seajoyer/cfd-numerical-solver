#pragma once
#include <yaml-cpp/yaml.h>
#include "data/InitialConditions.hpp"
#include "data/GlobalVariables.hpp"

class ConfigParser {
public:
    ConfigParser();
    bool parse(const std::string& filename, const std::string& initial);
    const InitialConditions& getInitialConditions() const;

private:
    InitialConditions initial_conditions;
    GlobalVariables global_variables;


    static void loadInitialConditions(const YAML::Node& node, InitialConditions& conditions);
    static void loadGlobalVariables(const YAML::Node& node, GlobalVariables& globals);
};