#include "config/ConfigParser.hpp"

#include <cxxopts.hpp>
#include <iostream>

#include "config/InitialConditions.hpp"
#include "config/Settings.hpp"

ConfigParser::ConfigParser() = default;

auto ConfigParser::ParseFile(const std::string& filename) -> bool {
    try {
        YAML::Node config = YAML::LoadFile(filename);
        YAML::Node settings_node = config["config"]["settings"];
        YAML::Node ic_node = config["config"]["initial_conditions"];

        LoadSettings(settings_node, settings_);
        LoadInitialConditions(ic_node, sod_tests_);
        config_path_ = filename;

        return true;
    } catch (const std::exception& e) {
        std::cerr << "Error parsing YAML: " << e.what() << '\n';
        return false;
    }
}

auto ConfigParser::ParseCommandLine(int argc, char* argv[]) -> std::optional<bool> {
    try {
        cxxopts::Options opts("cfd_numerical_solver", "CFD Numerical Solver - Computational Fluid Dynamics");

        // clang-format off
        opts.add_options()
            ("h,help", "Print help")
            ("c,config", "Config file path", cxxopts::value<std::string>())
        ;

        opts.add_options("Solver")
            ("solver", "Solver type (analytical, godunov, godunov-kolgan, godunov-kolgan-rodionov)", cxxopts::value<std::string>())
            ("riemann-solver", "Riemann solver (exact, hll, hllc, acoustic)", cxxopts::value<std::string>())
            ("reconstruction", "Reconstruction scheme (P0, P1)", cxxopts::value<std::string>())
            ("left-boundary", "Left boundary condition", cxxopts::value<std::string>())
            ("right-boundary", "Right boundary condition", cxxopts::value<std::string>())
        ;

        opts.add_options("Grid")
            ("N,N-cells", "Number of cells", cxxopts::value<int>())
            ("dim", "Dimensions (only 1D is supported for now)", cxxopts::value<int>())
            ("Lx", "Domain length X", cxxopts::value<double>())
            ("Ly", "Domain length Y (not yet supported)", cxxopts::value<double>())
            ("Lz", "Domain length Z (not yet supported)", cxxopts::value<double>())
            ("padding", "Ghost cell padding", cxxopts::value<int>())
        ;

        opts.add_options("Physics")
            ("cfl", "CFL number", cxxopts::value<double>())
            ("gamma", "Heat capacity ratio", cxxopts::value<double>())
            ("t-end", "End time", cxxopts::value<double>())
            ("Q-user", "User-defined Q parameter", cxxopts::value<double>())
        ;

        opts.add_options("Initial Conditions")
            ("sod-test", "SOD test number (1-5)", cxxopts::value<int>())
            ("x0", "Discontinuity position", cxxopts::value<double>())
            ("analytical", "Enable analytical solution", cxxopts::value<bool>())
        ;

        opts.add_options("Output")
            ("output-steps", "Output every N steps", cxxopts::value<std::size_t>())
            ("output-time", "Output every N time units", cxxopts::value<double>())
            ("output-format", "Output format (only VTK is supported for now)", cxxopts::value<std::string>())
            ("output-dir", "Output directory", cxxopts::value<std::string>())
        ;
        // clang-format on

        auto result = opts.parse(argc, argv);

        if (result.count("help")) {
            std::cout << opts.help() << '\n';
            return std::nullopt;  // Signal that help was shown
        }

        // Override settings with CLI arguments
        if (result.count("config")) config_path_ = result["config"].as<std::string>();

        // Solver options
        if (result.count("solver")) settings_.solver = result["solver"].as<std::string>();
        if (result.count("riemann-solver")) {
            settings_.riemann_solver = result["riemann-solver"].as<std::string>();
        }
        if (result.count("reconstruction")) {
            settings_.reconstruction = result["reconstruction"].as<std::string>();
        }
        if (result.count("left-boundary")) {
            settings_.left_boundary = result["left-boundary"].as<std::string>();
        }
        if (result.count("right-boundary")) {
            settings_.right_boundary = result["right-boundary"].as<std::string>();
        }

        // Grid options
        if (result.count("N-cells")) settings_.N = result["N-cells"].as<int>();
        if (result.count("dim")) settings_.dim = result["dim"].as<int>();
        if (result.count("Lx")) settings_.L_x = result["Lx"].as<double>();
        if (result.count("Ly")) settings_.L_y = result["Ly"].as<double>();
        if (result.count("Lz")) settings_.L_z = result["Lz"].as<double>();
        if (result.count("padding")) settings_.padding = result["padding"].as<int>();

        // Physics options
        if (result.count("cfl")) settings_.cfl = result["cfl"].as<double>();
        if (result.count("gamma")) settings_.gamma = result["gamma"].as<double>();
        if (result.count("t-end")) settings_.t_end = result["t-end"].as<double>();
        if (result.count("Q-user")) settings_.Q_user = result["Q-user"].as<double>();

        // Initial conditions
        if (result.count("sod-test")) {
            settings_.sod_test_num = result["sod-test"].as<int>();
        }
        if (result.count("x0")) settings_.x0 = result["x0"].as<double>();
        if (result.count("analytical")) {
            settings_.analytical = result["analytical"].as<bool>();
        }

        // Output options
        if (result.count("output-steps")) {
            settings_.output_every_steps = result["output-steps"].as<std::size_t>();
        }
        if (result.count("output-time")) {
            settings_.output_every_time = result["output-time"].as<double>();
        }
        if (result.count("output-format")) {
            settings_.output_format = result["output-format"].as<std::string>();
        }
        if (result.count("output-dir")) {
            settings_.output_dir = result["output-dir"].as<std::string>();
        }

        return true;

    } catch (const cxxopts::exceptions::exception& e) {
        std::cerr << "CLI parsing error: " << e.what() << '\n';
        return false;
    }
}

auto ConfigParser::Parse(const std::string& default_config, int argc, char* argv[])
    -> std::optional<bool> {
    // First, do a preliminary CLI parse to check for --config override
    config_path_ = default_config;

    auto cli_result = ParseCommandLine(argc, argv);
    if (!cli_result.has_value()) {
        return std::nullopt;  // Help was shown
    }
    if (!cli_result.value()) {
        return false;  // CLI parsing error
    }

    // Now parse the config file (possibly overridden path)
    if (!ParseFile(config_path_)) {
        return false;
    }

    // Re-parse CLI to override file settings
    return ParseCommandLine(argc, argv);
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

auto ConfigParser::GetConfigPath() const -> const std::string& { return config_path_; }

void ConfigParser::LoadSettings(const YAML::Node& node, Settings& settings) {
    settings.solver = node["solver"].as<std::string>();
    settings.riemann_solver = node["riemann_solver"].as<std::string>();
    settings.reconstruction = node["reconstruction"].as<std::string>();
    settings.left_boundary = node["left_boundary"].as<std::string>();
    settings.right_boundary = node["right_boundary"].as<std::string>();

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
    settings.output_every_time = node["output_every_time"].as<double>();
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
