#include "config/ConfigParser.hpp"
#include <iostream>
#include "config/Settings.hpp"
#include "config/InitialConditions.hpp"

ConfigParser::ConfigParser() = default;

auto ConfigParser::Parse(const std::string& filename) -> bool {
    try {
        YAML::Node config_node = YAML::LoadFile(filename);
        YAML::Node settings_node = config_node["config"]["settings"];
        YAML::Node initial_conditions_node = config_node["config"]["initial_conditions"];

        LoadSettings(settings_node, settings_);
        LoadInitialConditions(initial_conditions_node, initial_conditions_);

        return true;

    } catch (const std::exception& e) {
        std::cerr << "Error parsing YAML: " << e.what() << '\n';
        return false;
    }
}

auto ConfigParser::GetSettings() const -> const Settings& {
    return settings_;
}

auto ConfigParser::GetInitialConditions() const -> const InitialConditions& {
    return initial_conditions_;
}

void ConfigParser::LoadSettings(const YAML::Node& node, Settings& settings) {
    settings.solver = node["solver"].as<std::string>();
    settings.boundary = node["boundary"].as<std::string>();

    settings.N = node["N"].as<int>();
    settings.CFL = node["CFL"].as<double>();
    settings.t_end = node["t_end"].as<double>();
    settings.padding = node["padding"].as<int>();
    settings.c = node["c"].as<double>();
    settings.gamma = node["gamma"].as<double>();
    settings.dim = node["dim"].as<int>();
    settings.L_x = node["L_x"].as<double>();
    settings.L_y = node["L_y"].as<double>();
    settings.L_z = node["L_z"].as<double>();

    settings.output_every_steps = node["output_every_steps"].as<std::size_t>();
    settings.output_types = node["output_types"].as<std::vector<std::string>>();
    settings.output_dir = node["output_dir"].as<std::string>();
}

void ConfigParser::LoadInitialConditions(const YAML::Node& node, InitialConditions& conditions) {
    conditions.rho_L = node["rho_L"].as<double>();
    conditions.u_L = node["u_L"].as<double>();
    conditions.P_L = node["P_L"].as<double>();
    conditions.rho_R = node["rho_R"].as<double>();
    conditions.u_R = node["u_R"].as<double>();
    conditions.P_R = node["P_R"].as<double>();
}
