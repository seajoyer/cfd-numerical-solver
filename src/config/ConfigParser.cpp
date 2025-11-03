#include "config/ConfigParser.hpp"
#include <iostream>
#include "data/InitialConditions.hpp"

ConfigParser::ConfigParser() = default;

auto ConfigParser::Parse(const std::string& filename, const std::string& initial) -> bool {
    try {
        YAML::Node config = YAML::LoadFile(filename);
        YAML::Node  params = config["config"]["params"];
        LoadInitialConditions(params["initial"][initial], initial_conditions_);
        // std::cout << params["initial"][initial] << std::endl;
        LoadGlobalVariables(params, global_variables_);
        // std::cout << global_variables.c << std::endl;
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Error parsing YAML: " << e.what() << '\n';
        return false;
    }
}

void ConfigParser::LoadInitialConditions(const YAML::Node& node, InitialConditions& conditions) {
    conditions.rho_L = node["rho_L"].as<double>();
    conditions.u_L = node["u_L"].as<double>();
    conditions.P_L = node["P_L"].as<double>();
    conditions.rho_R = node["rho_R"].as<double>();
    conditions.u_R = node["u_R"].as<double>();
    conditions.P_R = node["P_R"].as<double>();
}

void ConfigParser::LoadGlobalVariables(const YAML::Node& node, GlobalVariables& globals) {
    globals.t_final = node["t_final"].as<double>();
    globals.N = node["N"].as<int>();
    globals.padding = node["padding"].as<int>();
    globals.CFL = node["CFL"].as<double>();
    globals.c = node["c"].as<double>();
    globals.gamma = node["gamma"].as<double>();
    globals.dim = node["dim"].as<int>();
    globals.L_x = node["L_x"].as<double>();
    globals.L_y = node["L_y"].as<double>();
    globals.L_z = node["L_z"].as<double>();
}

auto ConfigParser::GetInitialConditions() const -> const InitialConditions& {
    return initial_conditions_;
}
