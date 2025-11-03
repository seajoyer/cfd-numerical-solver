#include "config/ConfigParser.hpp"
#include "data/InitialConditions.hpp"
#include "data/DataLayer.hpp"
#include <iostream>

ConfigParser::ConfigParser() {}

bool ConfigParser::parse(const std::string& filename, const std::string& initial) {
    try {
        YAML::Node config = YAML::LoadFile(filename);
        YAML::Node  params = config["config"]["params"];
        loadInitialConditions(params["initial"][initial], initial_conditions);
        // std::cout << params["initial"][initial] << std::endl;
        loadGlobalVariables(params, global_variables);
        // std::cout << global_variables.c << std::endl;
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Error parsing YAML: " << e.what() << std::endl;
        return false;
    }
}

void ConfigParser::loadInitialConditions(const YAML::Node& node, InitialConditions& conditions) {
    conditions.rho_L = node["rho_L"].as<double>();
    conditions.u_L = node["u_L"].as<double>();
    conditions.P_L = node["P_L"].as<double>();
    conditions.rho_R = node["rho_R"].as<double>();
    conditions.u_R = node["u_R"].as<double>();
    conditions.P_R = node["P_R"].as<double>();
}

void ConfigParser::loadGlobalVariables(const YAML::Node& node, GlobalVariables& globals) {
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

const InitialConditions& ConfigParser::getInitialConditions() const {
    return initial_conditions;
}
