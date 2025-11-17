#include "config/ConfigParser.hpp"

#include <iostream>

#include "config/InitialConditions.hpp"
#include "config/Settings.hpp"

ConfigParser::ConfigParser() = default;

auto ConfigParser::Parse(const std::string& filename) -> bool {
    try {
        YAML::Node config_node = YAML::LoadFile(filename);
        YAML::Node settings_node = config_node["config"]["settings"];
        YAML::Node initial_conditions_node = config_node["config"]["initial_conditions"];

        LoadSettings(settings_node, settings_);
        LoadInitialConditions(initial_conditions_node, sod_tests_);

        return true;

    } catch (const std::exception& e) {
        std::cerr << "Error parsing YAML: " << e.what() << '\n';
        return false;
    }
}

auto ConfigParser::GetSettings() const -> const Settings& { return settings_; }

auto ConfigParser::GetSODTests() const -> const std::array<InitialConditions, 5>& {
    return sod_tests_;
}

auto ConfigParser::GetSODTest(int test_num) const -> const InitialConditions& {
    if (test_num < 1 || test_num > 5) {
        throw std::out_of_range("SOD test index must be between 1 and 5");
    }
    return sod_tests_[test_num - 1];
}

void ConfigParser::LoadSettings(const YAML::Node& node, Settings& settings) {
    settings.solver = node["solver"].as<std::string>();
    settings.riemann_solver = node["riemann_solver"].as<std::string>();
    settings.reconstruction = node["reconstruction"].as<std::string>();
    settings.right_boundary = node["left_boundary"].as<std::string>();
    settings.left_boundary = node["right_boundary"].as<std::string>();

    settings.N = node["N"].as<int>();
    settings.cfl = node["cfl"].as<double>();
    settings.t_end = node["t_end"].as<double>();
    settings.padding = node["padding"].as<int>();
    settings.gamma = node["gamma"].as<double>();
    settings.dim = node["dim"].as<int>();
    settings.L_x = node["L_x"].as<double>();
    settings.L_y = node["L_y"].as<double>();
    settings.L_z = node["L_z"].as<double>();

    settings.Q_user = node["Q_user"].as<double>();

    settings.sod_test_num = node["sod_test_num"].as<int>();
    settings.x0 = node["x0"].as<double>();
    settings.analytical = node["analytical"].as<std::string>() == "true";

    settings.output_every_steps = node["output_every_steps"].as<std::size_t>();
    settings.output_format = node["output_format"].as<std::string>();
    settings.output_dir = node["output_dir"].as<std::string>();
}

void ConfigParser::LoadInitialConditions(const YAML::Node& node,
                                         std::array<InitialConditions, 5>& sod_tests) {
    for (int i = 1; i <= 5; ++i) {
        std::string key = "sod" + std::to_string(i);
        sod_tests[i - 1].rho_L = node[key]["rho_L"].as<double>();
        sod_tests[i - 1].u_L = node[key]["u_L"].as<double>();
        sod_tests[i - 1].P_L = node[key]["P_L"].as<double>();
        sod_tests[i - 1].rho_R = node[key]["rho_R"].as<double>();
        sod_tests[i - 1].u_R = node[key]["u_R"].as<double>();
        sod_tests[i - 1].P_R = node[key]["P_R"].as<double>();
    }
}
