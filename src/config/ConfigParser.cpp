#include "config/ConfigParser.hpp"

#include <cxxopts.hpp>
#include <iostream>
#include <algorithm>
#include <sstream>

#include "config/InitialConditions.hpp"
#include "config/Settings.hpp"
#include "utils/StringUtils.hpp"

ConfigParser::ConfigParser() = default;

auto ConfigParser::ParseOutputFormats(const YAML::Node& node) -> std::vector<std::string> {
    std::vector<std::string> formats;
    
    if (node.IsSequence()) {
        // Handle list format: [vtk, png]
        for (const auto& format_node : node) {
            std::string fmt = format_node.as<std::string>();
            // Convert to lowercase
            std::transform(fmt.begin(), fmt.end(), fmt.begin(),
                          [](unsigned char c) { return std::tolower(c); });
            formats.push_back(fmt);
        }
    } else if (node.IsScalar()) {
        // Handle single string format: vtk
        std::string fmt = node.as<std::string>();
        std::transform(fmt.begin(), fmt.end(), fmt.begin(),
                      [](unsigned char c) { return std::tolower(c); });
        formats.push_back(fmt);
    }
    
    return formats;
}

auto ConfigParser::ParseOutputFormatsString(const std::string& formats_str) -> std::vector<std::string> {
    std::vector<std::string> formats;
    std::istringstream ss(formats_str);
    std::string format;
    
    while (std::getline(ss, format, ',')) {
        // Trim whitespace
        format.erase(0, format.find_first_not_of(" \t"));
        format.erase(format.find_last_not_of(" \t") + 1);
        
        if (!format.empty()) {
            // Convert to lowercase
            std::transform(format.begin(), format.end(), format.begin(),
                          [](unsigned char c) { return std::tolower(c); });
            formats.push_back(format);
        }
    }
    
    return formats;
}

auto ConfigParser::ParseCommandLineForConfigAndHelp(int argc, char* argv[]) -> std::optional<bool> {
    try {
        cxxopts::Options opts("cfd_numerical_solver", "CFD Numerical Solver - Computational Fluid Dynamics");

        opts.add_options()
            ("h,help", "Print help")
            ("c,config", "Config file path", cxxopts::value<std::string>())
        ;
        
        opts.allow_unrecognised_options();

        auto result = opts.parse(argc, argv);

        if (result.count("help")) {
            cxxopts::Options full_opts("cfd_numerical_solver", "CFD Numerical Solver - Computational Fluid Dynamics");
            
            full_opts.add_options()
                ("h,help", "Print help")
                ("c,config", "Config file path", cxxopts::value<std::string>())
                ("list-solvers", "List all supported solvers and their reconstructions")
                ("list-riemann", "List all supported Riemann solvers")
                ("list-boundaries", "List all supported boundary conditions")
                ("list-cases", "List all available simulation cases")
                ("list-formats", "List all supported output formats")
            ;

            full_opts.add_options("Solver")
                ("s,solver", "Solver type (analytical, godunov, godunov-kolgan, godunov-kolgan-rodionov)", cxxopts::value<std::string>())
                ("riemann-solver", "Riemann solver (exact, hll, hllc, acoustic)", cxxopts::value<std::string>())
                ("reconstruction", "Reconstruction scheme (P0, P1, ENO3, WENO5, etc.)", cxxopts::value<std::string>())
                ("l,left-boundary", "Left boundary condition", cxxopts::value<std::string>())
                ("r,right-boundary", "Right boundary condition", cxxopts::value<std::string>())
            ;

            full_opts.add_options("Grid")
                ("N,N-cells", "Number of cells", cxxopts::value<int>())
                ("d,dim", "Dimensions (only 1D is supported for now)", cxxopts::value<int>())
                ("x,Lx", "Domain length X", cxxopts::value<double>())
                ("y,Ly", "Domain length Y (not yet supported)", cxxopts::value<double>())
                ("z,Lz", "Domain length Z (not yet supported)", cxxopts::value<double>())
                ("p,padding", "Ghost cell padding", cxxopts::value<int>())
            ;

            full_opts.add_options("Physics")
                ("cfl", "CFL number", cxxopts::value<double>())
                ("gamma", "Heat capacity ratio", cxxopts::value<double>())
                ("t-end", "End time", cxxopts::value<double>())
                ("step-end", "Maximum number of steps", cxxopts::value<std::size_t>())
                ("Q-user", "User-defined Q parameter", cxxopts::value<double>())
            ;

            full_opts.add_options("Execution")
                ("run-cases", "Comma-separated list of cases to run (or 'all')", cxxopts::value<std::string>())
            ;

            full_opts.add_options("Initial Conditions")
                ("x0", "Discontinuity position", cxxopts::value<double>())
                ("a,analytical", "Enable analytical solution", cxxopts::value<bool>())
            ;

            full_opts.add_options("Output")
                ("output-steps", "Output every N steps", cxxopts::value<std::size_t>())
                ("output-time", "Output every N time units", cxxopts::value<double>())
                ("output-formats", "Comma-separated output formats", cxxopts::value<std::string>())
                ("o,output-dir", "Output directory", cxxopts::value<std::string>())
            ;
            
            std::cout << full_opts.help() << '\n';
            return std::nullopt;
        }

        // Extract --config if provided
        if (result.count("config")) {
            config_path_ = result["config"].as<std::string>();
        }

        return true;

    } catch (const cxxopts::exceptions::exception& e) {
        std::cerr << "CLI parsing error: " << e.what() << '\n';
        return false;
    }
}

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
            ("list-solvers", "List all supported solvers and their reconstructions")
            ("list-riemann", "List all supported Riemann solvers")
            ("list-boundaries", "List all supported boundary conditions")
            ("list-cases", "List all available simulation cases")
            ("list-formats", "List all supported output formats")
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
            ("step-end", "Maximum number of steps", cxxopts::value<std::size_t>())
            ("Q-user", "User-defined Q parameter", cxxopts::value<double>())
        ;

        opts.add_options("Execution")
            ("run-cases", "Comma-separated list of cases to run (or 'all')", cxxopts::value<std::string>())
        ;

        opts.add_options("Initial Conditions")
            ("x0", "Discontinuity position", cxxopts::value<double>())
            ("a,analytical", "Enable analytical solution", cxxopts::value<bool>())
        ;

        opts.add_options("Output")
            ("output-steps", "Output every N steps", cxxopts::value<std::size_t>())
            ("output-time", "Output every N time units", cxxopts::value<double>())
            ("output-formats", "Comma-separated output formats", cxxopts::value<std::string>())
            ("o,output-dir", "Output directory", cxxopts::value<std::string>())
        ;
        // clang-format on

        auto result = opts.parse(argc, argv);

        if (result.count("help")) {
            std::cout << opts.help() << '\n';
            return std::nullopt;  // Signal that help was shown
        }

        // Handle listing flags
        if (result.count("list-solvers")) {
            PrintSolversList();
            return std::nullopt;
        }

        if (result.count("list-riemann")) {
            PrintRiemannSolversList();
            return std::nullopt;
        }

        if (result.count("list-boundaries")) {
            PrintBoundaryConditionsList();
            return std::nullopt;
        }

        if (result.count("list-cases")) {
            PrintCasesList();
            return std::nullopt;
        }

        if (result.count("list-formats")) {
            PrintOutputFormatsList();
            return std::nullopt;
        }

        // Override settings with CLI arguments
        if (result.count("config")) config_path_ = result["config"].as<std::string>();

        // Solver options
        if (result.count("solver")) {
            cli_overrides_.solver = result["solver"].as<std::string>();
        }
        if (result.count("riemann-solver")) {
            cli_overrides_.riemann_solver = result["riemann-solver"].as<std::string>();
        }
        if (result.count("reconstruction")) {
            cli_overrides_.reconstruction = result["reconstruction"].as<std::string>();
        }
        if (result.count("left-boundary")) {
            cli_overrides_.left_boundary = result["left-boundary"].as<std::string>();
        }
        if (result.count("right-boundary")) {
            cli_overrides_.right_boundary = result["right-boundary"].as<std::string>();
        }

        // Grid options
        if (result.count("N-cells")) cli_overrides_.N = result["N-cells"].as<int>();
        if (result.count("dim")) cli_overrides_.dim = result["dim"].as<int>();
        if (result.count("Lx")) cli_overrides_.L_x = result["Lx"].as<double>();
        if (result.count("Ly")) cli_overrides_.L_y = result["Ly"].as<double>();
        if (result.count("Lz")) cli_overrides_.L_z = result["Lz"].as<double>();
        if (result.count("padding")) cli_overrides_.padding = result["padding"].as<int>();

        // Physics options
        if (result.count("cfl")) cli_overrides_.cfl = result["cfl"].as<double>();
        if (result.count("gamma")) cli_overrides_.gamma = result["gamma"].as<double>();
        if (result.count("t-end")) cli_overrides_.t_end = result["t-end"].as<double>();
        if (result.count("step-end")) cli_overrides_.step_end = result["step-end"].as<std::size_t>();
        if (result.count("Q-user")) cli_overrides_.Q_user = result["Q-user"].as<double>();

        // Execution options
        if (result.count("run-cases")) {
            std::string cases_str = result["run-cases"].as<std::string>();
            run_cases_.clear();
            
            // Parse comma-separated list
            std::istringstream ss(cases_str);
            std::string case_name;
            while (std::getline(ss, case_name, ',')) {
                // Trim whitespace
                case_name.erase(0, case_name.find_first_not_of(" \t"));
                case_name.erase(case_name.find_last_not_of(" \t") + 1);
                if (!case_name.empty()) {
                    run_cases_.push_back(case_name);
                }
            }
        }

        // Initial conditions
        if (result.count("x0")) cli_overrides_.x0 = result["x0"].as<double>();
        if (result.count("analytical")) {
            cli_overrides_.analytical = result["analytical"].as<bool>();
        }

        // Output options
        if (result.count("output-steps")) {
            cli_overrides_.output_every_steps = result["output-steps"].as<std::size_t>();
        }
        if (result.count("output-time")) {
            cli_overrides_.output_every_time = result["output-time"].as<double>();
        }
        if (result.count("output-formats")) {
            std::string formats_str = result["output-formats"].as<std::string>();
            cli_overrides_.output_formats = ParseOutputFormatsString(formats_str);
        }
        if (result.count("output-dir")) {
            cli_overrides_.output_dir = result["output-dir"].as<std::string>();
        }

        return true;

    } catch (const cxxopts::exceptions::exception& e) {
        std::cerr << "CLI parsing error: " << e.what() << '\n';
        return false;
    }
}

auto ConfigParser::Parse(const std::string& default_config, int argc, char* argv[])
    -> std::optional<bool> {
    // First, do a preliminary CLI parse to check for --config override and --help
    config_path_ = default_config;

    auto cli_result = ParseCommandLineForConfigAndHelp(argc, argv);
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

    // Now parse CLI again for list flags and settings overrides
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
    return MergeSettings(settings_, ic.overrides, cli_overrides_);
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
    
    // Handle output_formats (can be string or list)
    if (node["output_formats"]) {
        settings.output_formats = ParseOutputFormats(node["output_formats"]);
    } else if (node["output_format"]) {
        // Backward compatibility: support old single format option
        std::string fmt = node["output_format"].as<std::string>();
        std::transform(fmt.begin(), fmt.end(), fmt.begin(),
                      [](unsigned char c) { return std::tolower(c); });
        settings.output_formats = {fmt};
    }
    
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
    
    // Handle output_formats (can be string or list)
    if (node["output_formats"]) {
        overrides.output_formats = ParseOutputFormats(node["output_formats"]);
    } else if (node["output_format"]) {
        std::string fmt = node["output_format"].as<std::string>();
        std::transform(fmt.begin(), fmt.end(), fmt.begin(),
                      [](unsigned char c) { return std::tolower(c); });
        overrides.output_formats = std::vector<std::string>{fmt};
    }
    
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

void ConfigParser::PrintSolversList() const {
    std::cout << "Supported Solvers and Compatible Reconstructions:\n";
    std::cout << "=================================================\n\n";
    
    std::cout << "  godunov:\n";
    std::cout << "      Compatible reconstructions: P0\n";
    std::cout << "      Description: First-order Godunov method\n";
    std::cout << "      Reconstruction: P0 (piecewise constant, first-order)\n\n";
    
    std::cout << "  godunov-kolgan:\n";
    std::cout << "      Compatible reconstructions: P0, P1, ENO, WENO\n";
    std::cout << "      Description: Godunov method with higher-order reconstruction\n";
    std::cout << "      Reconstructions:\n";
    std::cout << "        - P0 (piecewise constant, first-order)\n";
    std::cout << "        - P1 (piecewise linear with slope limiting, second-order)\n";
    std::cout << "        - ENO<N> (Essentially Non-Oscillatory, N = order, e.g., ENO3, ENO5)\n";
    std::cout << "        - WENO<N> (Weighted ENO, N = 3, 5, â€¦; e.g., WENO3, WENO5, â€¦)\n\n";
    
    std::cout << "  godunov-kolgan-rodionov:\n";
    std::cout << "      Compatible reconstructions: P1, ENO, WENO\n";
    std::cout << "      Description: Second-order MUSCL-Hancock scheme\n";
    std::cout << "      Reconstructions:\n";
    std::cout << "        - P1 (piecewise linear with slope limiting, MUSCL-Hancock)\n";
    std::cout << "        - ENO<N> (Essentially Non-Oscillatory, N = order, e.g., ENO3, ENO5)\n";
    std::cout << "        - WENO<N> (Weighted ENO, N = 3, 5, â€¦; e.g., WENO3, WENO5, â€¦)\n\n";
    
    std::cout << "  analytical:\n";
    std::cout << "      Compatible reconstructions: N/A\n";
    std::cout << "      Description: Exact Riemann solver for verification\n";
    std::cout << "      Note: No reconstruction needed (exact solution)\n\n";
    
    std::cout << "Note: Solver and reconstruction names are case-insensitive\n";
}

void ConfigParser::PrintRiemannSolversList() const {
    std::cout << "Supported Riemann Solvers:\n";
    std::cout << "==========================\n\n";
    
    std::cout << "  exact:\n";
    std::cout << "      Exact ideal gas Riemann solver (iterative)\n";
    std::cout << "      Most accurate but computationally expensive\n\n";
    
    std::cout << "  hll:\n";
    std::cout << "      HLL (Harten-Lax-van Leer) approximate solver\n";
    std::cout << "      Robust but diffusive, two-wave approximation\n\n";
    
    std::cout << "  hllc:\n";
    std::cout << "      HLLC (HLL-Contact) approximate solver\n";
    std::cout << "      Restores contact discontinuity, less diffusive than HLL\n\n";
    
    std::cout << "  acoustic:\n";
    std::cout << "      Linearized acoustic Riemann solver\n";
    std::cout << "      Fast but only accurate for weak perturbations\n\n";
    
    std::cout << "Note: Riemann solver names are case-insensitive\n";
}

void ConfigParser::PrintBoundaryConditionsList() const {
    std::cout << "Supported Boundary Conditions:\n";
    std::cout << "==============================\n\n";
    
    std::cout << "  free_stream:\n";
    std::cout << "      Far-field boundary with specified external state\n\n";
    
    std::cout << "  inlet:\n";
    std::cout << "      Inflow boundary with prescribed values\n\n";
    
    std::cout << "  outlet:\n";
    std::cout << "      Outflow boundary (zero-gradient extrapolation)\n\n";
    
    std::cout << "  reflective:\n";
    std::cout << "      Reflecting wall (inverts normal velocity)\n\n";
    
    std::cout << "  non_reflective:\n";
    std::cout << "      Non-reflecting characteristic boundary\n\n";
    
    std::cout << "  periodic:\n";
    std::cout << "      Periodic boundary condition\n\n";
    
    std::cout << "  symmetry:\n";
    std::cout << "      Symmetry plane\n\n";
    
    std::cout << "  wall:\n";
    std::cout << "      Solid wall (zero normal velocity)\n\n";
}

void ConfigParser::PrintCasesList() const {
    if (initial_conditions_.empty()) {
        std::cout << "No cases loaded. Please load a configuration file first.\n";
        return;
    }
    
    std::cout << "Available Simulation Cases:\n";
    std::cout << "===========================\n\n";
    
    auto case_names = GetAllCaseNames();
    for (const auto& case_name : case_names) {
        const auto& ic = initial_conditions_.at(case_name);
        
        std::cout << "  " << case_name << ":\n";
        std::cout << "      Left state:   rho_L = " << ic.rho_L 
                  << ",  u_L = " << ic.u_L 
                  << ",  P_L = " << ic.P_L << "\n";
        std::cout << "      Right state:  rho_R = " << ic.rho_R 
                  << ",  u_R = " << ic.u_R 
                  << ",  P_R = " << ic.P_R << "\n";
        
        // Show case-specific overrides if any
        bool has_overrides = false;
        std::ostringstream overrides_str;
        
        if (ic.overrides.x0) {
            overrides_str << "x0=" << *ic.overrides.x0 << ", ";
            has_overrides = true;
        }
        if (ic.overrides.t_end) {
            overrides_str << "t_end=" << *ic.overrides.t_end << ", ";
            has_overrides = true;
        }
        if (ic.overrides.reconstruction) {
            overrides_str << "reconstruction=" << *ic.overrides.reconstruction << ", ";
            has_overrides = true;
        }
        if (ic.overrides.solver) {
            overrides_str << "solver=" << *ic.overrides.solver << ", ";
            has_overrides = true;
        }
        if (ic.overrides.output_formats) {
            overrides_str << "output_formats=[";
            for (size_t i = 0; i < ic.overrides.output_formats->size(); ++i) {
                if (i > 0) overrides_str << ",";
                overrides_str << (*ic.overrides.output_formats)[i];
            }
            overrides_str << "], ";
            has_overrides = true;
        }
        
        if (has_overrides) {
            std::string overrides = overrides_str.str();
            // Remove trailing ", "
            if (overrides.size() >= 2) {
                overrides = overrides.substr(0, overrides.size() - 2);
            }
            std::cout << "      Overrides:    " << overrides << "\n";
        }
        
        std::cout << "\n";
    }
    
    std::cout << "Total: " << case_names.size() << " case(s)\n";
}

void ConfigParser::PrintOutputFormatsList() const {
    std::cout << "Supported Output Formats:\n";
    std::cout << "=========================\n\n";
    
    std::cout << "  vtk:\n";
    std::cout << "      VTK structured grid format\n";
    std::cout << "      Compatible with ParaView, VisIt, and other visualization tools\n";
    std::cout << "      Contains all simulation fields (density, velocity, pressure, etc.)\n\n";
    
    std::cout << "  png or png<width>x<height>:\n";
    std::cout << "      PNG image with 4 subplot panels:\n";
    std::cout << "        - Top-left:     Density vs X\n";
    std::cout << "        - Top-right:    Velocity vs X\n";
    std::cout << "        - Bottom-left:  Pressure vs X\n";
    std::cout << "        - Bottom-right: Specific Internal Energy vs X\n";
    std::cout << "      Numerical solution: red line\n";
    std::cout << "      Analytical solution: black line (if enabled)\n";
    std::cout << "      Default resolution: 1200x900\n";
    std::cout << "      Custom resolution examples: png1920x1080, png800x600, png3840x2160\n";
    std::cout << "      Font sizes and line widths scale automatically with resolution\n\n";
    
    std::cout << "  gif or gif<width>x<height>:\n";
    std::cout << "      Animated GIF with same 4 subplot panels as PNG:\n";
    std::cout << "        - Top-left:     Density vs X\n";
    std::cout << "        - Top-right:    Velocity vs X\n";
    std::cout << "        - Bottom-left:  Pressure vs X\n";
    std::cout << "        - Bottom-right: Specific Internal Energy vs X\n";
    std::cout << "      Numerical solution: red line\n";
    std::cout << "      Analytical solution: black line (if enabled)\n";
    std::cout << "      Default resolution: 1200x900\n";
    std::cout << "      Custom resolution examples: gif1920x1080, gif800x600\n";
    std::cout << "      Font sizes and line widths scale automatically with resolution\n";
    std::cout << "      Frame delay: 10 centiseconds (100ms) per frame\n\n";
    
    std::cout << "      Memory note: GIF stores all frames in memory until finalization.\n";
    std::cout << "      For long simulations, consider reducing output frequency or resolution.\n\n";
    
    std::cout << "Multiple formats can be specified:\n";
    std::cout << "  YAML:  output_formats: [vtk, png1920x1080, gif]\n";
    std::cout << "  CLI:   --output-formats vtk,png1920x1080,gif\n\n";
    
    std::cout << "Note: Format names are case-insensitive\n";
}
