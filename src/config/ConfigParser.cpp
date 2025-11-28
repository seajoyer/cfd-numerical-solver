#include "config/ConfigParser.hpp"

#include <cxxopts.hpp>
#include <iostream>
#include <algorithm>

#include "config/InitialConditions.hpp"
#include "config/Settings.hpp"
#include "utils/StringUtils.hpp"

ConfigParser::ConfigParser() = default;

auto ConfigParser::ParseFile(const std::string& filename) -> bool {
    try {
        YAML::Node config = YAML::LoadFile(filename);
        
        // Load global settings from config/global
        if (config["config"]["global"]) {
            YAML::Node global_node = config["config"]["global"];
            LoadSettings(global_node, settings_);
        } else {
            std::cerr << "Warning: No 'global' section found in config. Using defaults.\n";
        }
        
        // Load cases from config/cases
        if (config["config"]["cases"]) {
            YAML::Node cases_node = config["config"]["cases"];
            LoadCases(cases_node, initial_conditions_);
        } else {
            std::cerr << "Error: No 'cases' section found in config.\n";
            return false;
        }
        
        // Load run_cases from config/run_cases
        if (config["config"]["run_cases"]) {
            YAML::Node run_cases_node = config["config"]["run_cases"];
            LoadRunCases(run_cases_node, run_cases_);
        } else {
            // Default: run all cases
            run_cases_ = {"all"};
        }
        
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
            ("s,solver", "Solver type (analytical, godunov, godunov-kolgan, godunov-kolgan-rodionov)", cxxopts::value<std::string>())
            ("riemann-solver", "Riemann solver (exact, hll, hllc, acoustic)", cxxopts::value<std::string>())
            ("reconstruction", "Reconstruction scheme (P0, P1, ENO3, WENO5, etc.)", cxxopts::value<std::string>())
            ("l,left-boundary", "Left boundary condition", cxxopts::value<std::string>())
            ("r,right-boundary", "Right boundary condition", cxxopts::value<std::string>())
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
            ("x0", "Discontinuity position", cxxopts::value<double>())
            ("a,analytical", "Enable analytical solution", cxxopts::value<bool>())
        ;

        opts.add_options("Output")
            ("output-steps", "Output every N steps", cxxopts::value<std::size_t>())
            ("output-time", "Output every N time units", cxxopts::value<double>())
            ("f,output-format", "Output format (only VTK is supported for now)", cxxopts::value<std::string>())
            ("o,output-dir", "Output directory", cxxopts::value<std::string>())
        ;
        // clang-format on

        auto result = opts.parse(argc, argv);

        if (result.count("help")) {
            std::cout << opts.help() << '\n';
            return std::nullopt;  // Signal that help was shown
        }

        // Override settings with CLI arguments (only affects global settings)
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

auto ConfigParser::GetRunCases() const -> const std::vector<std::string>& {
    return run_cases_;
}

auto ConfigParser::GetInitialConditions() const -> const std::map<std::string, InitialConditions>& {
    return initial_conditions_;
}

auto ConfigParser::GetInitialCondition(const std::string& case_name) const -> const InitialConditions& {
    auto it = initial_conditions_.find(case_name);
    if (it == initial_conditions_.end()) {
        throw std::out_of_range("Initial condition '" + case_name + "' not found in configuration");
    }
    return it->second;
}

auto ConfigParser::GetCaseSettings(const std::string& case_name) const -> Settings {
    const InitialConditions& ic = GetInitialCondition(case_name);
    return MergeSettings(settings_, ic.overrides);
}

auto ConfigParser::GetAllCaseNames() const -> std::vector<std::string> {
    std::vector<std::string> names;
    names.reserve(initial_conditions_.size());
    for (const auto& [name, _] : initial_conditions_) {
        names.push_back(name);
    }
    // Sort for consistent ordering
    std::sort(names.begin(), names.end());
    return names;
}

auto ConfigParser::HasInitialCondition(const std::string& case_name) const -> bool {
    return initial_conditions_.find(case_name) != initial_conditions_.end();
}

auto ConfigParser::GetConfigPath() const -> const std::string& { return config_path_; }

void ConfigParser::LoadSettings(const YAML::Node& node, Settings& settings) {
    if (node["solver"]) settings.solver = node["solver"].as<std::string>();
    if (node["solver"]) settings.solver = utils::ToLower(settings.solver);
    
    if (node["riemann_solver"]) settings.riemann_solver = node["riemann_solver"].as<std::string>();
    if (node["riemann_solver"]) settings.riemann_solver = utils::ToLower(settings.riemann_solver);
    
    if (node["reconstruction"]) settings.reconstruction = node["reconstruction"].as<std::string>();
    if (node["reconstruction"]) settings.reconstruction = utils::ToLower(settings.reconstruction);
    
    if (node["left_boundary"]) settings.left_boundary = node["left_boundary"].as<std::string>();
    if (node["right_boundary"]) settings.right_boundary = node["right_boundary"].as<std::string>();

    if (node["N"]) settings.N = node["N"].as<int>();
    if (node["cfl"]) settings.cfl = node["cfl"].as<double>();
    if (node["t_end"]) settings.t_end = node["t_end"].as<double>();
    if (node["step_end"]) settings.step_end = node["step_end"].as<std::size_t>();
    if (node["padding"]) settings.padding = node["padding"].as<int>();
    if (node["gamma"]) settings.gamma = node["gamma"].as<double>();
    if (node["dim"]) settings.dim = node["dim"].as<int>();
    if (node["L_x"]) settings.L_x = node["L_x"].as<double>();
    if (node["L_y"]) settings.L_y = node["L_y"].as<double>();
    if (node["L_z"]) settings.L_z = node["L_z"].as<double>();

    if (node["Q_user"]) settings.Q_user = node["Q_user"].as<double>();

    if (node["x0"]) settings.x0 = node["x0"].as<double>();
    if (node["analytical"]) {
        auto analytical_str = node["analytical"].as<std::string>();
        settings.analytical = (analytical_str == "true" || analytical_str == "True" || analytical_str == "TRUE");
    }

    if (node["log_every_steps"]) settings.log_every_steps = node["log_every_steps"].as<std::size_t>();
    if (node["log_every_time"]) settings.log_every_time = node["log_every_time"].as<double>();
    
    if (node["output_every_steps"]) settings.output_every_steps = node["output_every_steps"].as<std::size_t>();
    if (node["output_every_time"]) settings.output_every_time = node["output_every_time"].as<double>();
    if (node["output_format"]) settings.output_format = node["output_format"].as<std::string>();
    if (node["output_dir"]) settings.output_dir = node["output_dir"].as<std::string>();
}

void ConfigParser::LoadCases(const YAML::Node& node,
                              std::map<std::string, InitialConditions>& initial_conditions) {
    // Iterate through all entries in the cases section
    for (const auto& entry : node) {
        auto case_name = entry.first.as<std::string>();
        const YAML::Node& case_data = entry.second;
        
        InitialConditions ic;
        ic.rho_L = case_data["rho_L"].as<double>();
        ic.u_L = case_data["u_L"].as<double>();
        ic.P_L = case_data["P_L"].as<double>();
        ic.rho_R = case_data["rho_R"].as<double>();
        ic.u_R = case_data["u_R"].as<double>();
        ic.P_R = case_data["P_R"].as<double>();
        
        // Load case-specific overrides
        LoadCaseOverrides(case_data, ic.overrides);

        initial_conditions[case_name] = ic;
    }
}

void ConfigParser::LoadCaseOverrides(const YAML::Node& node, CaseSettings& overrides) {
    // Solver Configuration
    if (node["solver"]) {
        auto solver = node["solver"].as<std::string>();
        overrides.solver = utils::ToLower(solver);
    }
    if (node["riemann_solver"]) {
        auto riemann = node["riemann_solver"].as<std::string>();
        overrides.riemann_solver = utils::ToLower(riemann);
    }
    if (node["reconstruction"]) {
        auto recon = node["reconstruction"].as<std::string>();
        overrides.reconstruction = utils::ToLower(recon);
    }
    if (node["left_boundary"]) {
        overrides.left_boundary = node["left_boundary"].as<std::string>();
    }
    if (node["right_boundary"]) {
        overrides.right_boundary = node["right_boundary"].as<std::string>();
    }
    
    // Grid Configuration
    if (node["N"]) overrides.N = node["N"].as<int>();
    if (node["cfl"]) overrides.cfl = node["cfl"].as<double>();
    if (node["padding"]) overrides.padding = node["padding"].as<int>();
    if (node["gamma"]) overrides.gamma = node["gamma"].as<double>();
    if (node["dim"]) overrides.dim = node["dim"].as<int>();
    if (node["L_x"]) overrides.L_x = node["L_x"].as<double>();
    if (node["L_y"]) overrides.L_y = node["L_y"].as<double>();
    if (node["L_z"]) overrides.L_z = node["L_z"].as<double>();
    
    // Physical Parameters
    if (node["Q_user"]) overrides.Q_user = node["Q_user"].as<double>();
    
    // Initial Conditions
    if (node["x0"]) overrides.x0 = node["x0"].as<double>();
    if (node["analytical"]) {
        auto analytical_str = node["analytical"].as<std::string>();
        overrides.analytical = (analytical_str == "true" || analytical_str == "True" || analytical_str == "TRUE");
    }
    
    // Time Control
    if (node["t_end"]) overrides.t_end = node["t_end"].as<double>();
    if (node["step_end"]) overrides.step_end = node["step_end"].as<std::size_t>();
    
    // Logging Configuration
    if (node["log_every_steps"]) overrides.log_every_steps = node["log_every_steps"].as<std::size_t>();
    if (node["log_every_time"]) overrides.log_every_time = node["log_every_time"].as<double>();
    
    // Output Configuration
    if (node["output_every_steps"]) overrides.output_every_steps = node["output_every_steps"].as<std::size_t>();
    if (node["output_every_time"]) overrides.output_every_time = node["output_every_time"].as<double>();
    if (node["output_format"]) overrides.output_format = node["output_format"].as<std::string>();
    if (node["output_dir"]) overrides.output_dir = node["output_dir"].as<std::string>();
}

void ConfigParser::LoadRunCases(const YAML::Node& node, std::vector<std::string>& run_cases) {
    run_cases.clear();
    
    if (node.IsSequence()) {
        for (const auto& case_node : node) {
            run_cases.push_back(case_node.as<std::string>());
        }
    } else if (node.IsScalar()) {
        // Single case or "all"
        run_cases.push_back(node.as<std::string>());
    } else {
        std::cerr << "Warning: run_cases should be a list or scalar. Defaulting to 'all'.\n";
        run_cases.emplace_back("all");
    }
}
