#include "config/YamlConfigParser.hpp"

#include <algorithm>
#include <stdexcept>

#include "utils/StringUtils.hpp"

namespace {
    template <typename T>
    void AssignIfPresent(const YAML::Node& node, const char* key, T& target) {
        if (node[key]) {
            target = node[key].as<T>();
        }
    }

    template <typename T>
    void AssignOptionalIfPresent(const YAML::Node& node, const char* key, std::optional<T>& target) {
        if (node[key]) {
            target = node[key].as<T>();
        }
    }

    auto ReadLowerString(const YAML::Node& node, const char* key) -> std::optional<std::string> {
        if (!node[key]) {
            return std::nullopt;
        }
        return utils::ToLower(node[key].as<std::string>());
    }

    auto ParseBoundaryState(const YAML::Node& node) -> BoundaryStateSettings {
        BoundaryStateSettings state;
        AssignIfPresent(node, "rho", state.rho);
        AssignIfPresent(node, "u", state.u);
        AssignIfPresent(node, "v", state.v);
        AssignIfPresent(node, "w", state.w);
        AssignIfPresent(node, "p", state.p);
        return state;
    }
} // namespace

auto YamlConfigParser::ParseFile(const std::string& filename) -> ParsedYamlConfig {
    const YAML::Node root = YAML::LoadFile(filename);

    ParsedYamlConfig result;

    if (root["run"]) {
        ParseRunCases(root["run"], result.run_cases);
    }
    else {
        result.run_cases = {"all"};
    }

    if (root["defaults"]) {
        ParseDefaults(root["defaults"], result.settings);
    }

    if (!root["cases"]) {
        throw std::runtime_error("Missing 'cases' section in YAML config");
    }

    ParseCases(root["cases"], result.initial_conditions, result.settings);
    return result;
}

void YamlConfigParser::ParseRunCases(const YAML::Node& run_node,
                                     std::vector<std::string>& run_cases) {
    run_cases.clear();

    if (!run_node["cases"]) {
        run_cases.emplace_back("all");
        return;
    }

    const YAML::Node cases_node = run_node["cases"];
    if (!cases_node.IsSequence()) {
        throw std::runtime_error("'run.cases' must be a sequence");
    }

    for (const auto& item : cases_node) {
        run_cases.push_back(item.as<std::string>());
    }
}

void YamlConfigParser::ParseDefaults(const YAML::Node& defaults_node, Settings& settings) {
    if (defaults_node["mesh"]) {
        ParseMesh(defaults_node["mesh"], settings);
    }
    if (defaults_node["physics"]) {
        ParsePhysics(defaults_node["physics"], settings);
    }
    if (defaults_node["chemistry"]) {
        ParseChemistry(defaults_node["chemistry"], settings);
    }
    if (defaults_node["numerics"]) {
        ParseNumerics(defaults_node["numerics"], settings);
    }
    if (defaults_node["boundary_conditions"]) {
        ParseBoundaryConditions(defaults_node["boundary_conditions"], settings);
    }
    if (defaults_node["parallel"]) {
        ParseParallel(defaults_node["parallel"], settings);
    }
    if (defaults_node["stopping"]) {
        ParseStopping(defaults_node["stopping"], settings);
    }
    if (defaults_node["logging"]) {
        ParseLogging(defaults_node["logging"], settings);
    }
    if (defaults_node["output"]) {
        ParseOutput(defaults_node["output"], settings);
    }
    if (defaults_node["immersed_boundaries"]) {
        ParseImmersedBoundaries(defaults_node["immersed_boundaries"], settings);
    }
}

void YamlConfigParser::ParseCases(
    const YAML::Node& cases_node,
    std::map<std::string, InitialConditions>& initial_conditions,
    const Settings& defaults) {
    for (const auto& entry : cases_node) {
        const std::string case_name = entry.first.as<std::string>();
        const YAML::Node case_node = entry.second;

        InitialConditions ic;
        ic.overrides = CaseSettings{};

        ApplyCaseOverrides(case_node, ic);

        const Settings effective_settings = MergeSettings(defaults, ic.overrides, CaseSettings{});

        if (!case_node["initial_condition"]) {
            throw std::runtime_error("Case '" + case_name + "' does not contain 'initial_condition'");
        }

        const YAML::Node ic_node = case_node["initial_condition"];
        if (!ic_node["type"]) {
            throw std::runtime_error("Case '" + case_name + "' initial_condition must contain 'type'");
        }

        const std::string initial_condition_type = utils::ToLower(ic_node["type"].as<std::string>());
        ic.ic_type = initial_condition_type;

        if (initial_condition_type == "structured_regions") {
            ParseStructuredInitialCondition(ic_node, ic, effective_settings);
        }
        else if (initial_condition_type == "taylor_green") {
            // No region arrays are needed here.
            // Taylor-Green will be initialized later directly into
            // PressureVelocityState inside Simulation.
        }
        else if (initial_condition_type == "region_markers") {
            throw std::runtime_error("initial_condition.type=region_markers is not implemented yet");
        }
        else {
            throw std::runtime_error("Unsupported initial_condition.type: " + initial_condition_type);
        }

        initial_conditions[case_name] = ic;
    }
}

void YamlConfigParser::ParseMesh(const YAML::Node& node, Settings& settings) {
    if (node["type"]) {
        const std::string mesh_type = utils::ToLower(node["type"].as<std::string>());
        if (mesh_type != "structured") {
            throw std::runtime_error("Only mesh.type=structured is supported at the moment");
        }
    }

    AssignIfPresent(node, "dim", settings.dim);

    if (node["cells"]) {
        const YAML::Node cells = node["cells"];
        AssignIfPresent(cells, "x", settings.Nx);
        AssignIfPresent(cells, "y", settings.Ny);
        AssignIfPresent(cells, "z", settings.Nz);
    }

    if (node["domain"]) {
        const YAML::Node domain = node["domain"];
        AssignIfPresent(domain, "x", settings.L_x);
        AssignIfPresent(domain, "y", settings.L_y);
        AssignIfPresent(domain, "z", settings.L_z);
    }

    if (node["source"]) {
        throw std::runtime_error("Unstructured mesh source is not implemented yet");
    }
}

void YamlConfigParser::ParsePhysics(const YAML::Node& node, Settings& settings) {
    if (node["eos"]) {
        settings.EOS = utils::ToLower(node["eos"].as<std::string>());
    }

    if (node["parameters"]) {
        const YAML::Node parameters = node["parameters"];
        AssignIfPresent(parameters, "gamma", settings.gamma);
        AssignIfPresent(parameters, "Q_user", settings.Q_user);

        AssignIfPresent(parameters, "hg_rho0", settings.hg_rho0);
        AssignIfPresent(parameters, "hg_C", settings.hg_C);
        AssignIfPresent(parameters, "hg_S", settings.hg_S);
        AssignIfPresent(parameters, "hg_gamma_s", settings.hg_gamma_s);
        AssignIfPresent(parameters, "hg_c_v", settings.hg_c_v);
    }
}

void YamlConfigParser::ParseChemistry(const YAML::Node& node, Settings& settings) {
    AssignIfPresent(node, "enabled", settings.chemistry_enabled);
    AssignIfPresent(node, "z_freq", settings.chemistry_z_freq);
    AssignIfPresent(node, "activation_energy", settings.chemistry_activation_energy);
    AssignIfPresent(node, "gas_constant", settings.chemistry_gas_constant);
    AssignIfPresent(node, "heat_release", settings.chemistry_heat_release);
}

void YamlConfigParser::ParseNumerics(const YAML::Node& node, Settings& settings) {
    if (node["method"]) {
        settings.solver = utils::ToLower(node["method"].as<std::string>());
    }

    if (!node["parameters"]) {
        return;
    }

    const YAML::Node parameters = node["parameters"];

    if (parameters["time_integrator"]) {
        settings.time_integrator = utils::ToLower(parameters["time_integrator"].as<std::string>());
    }
    if (parameters["reconstruction"]) {
        settings.reconstruction = utils::ToLower(parameters["reconstruction"].as<std::string>());
    }
    if (parameters["riemann_solver"]) {
        settings.riemann_solver = utils::ToLower(parameters["riemann_solver"].as<std::string>());
    }
    if (parameters["transport_model"]) {
        settings.mader_transport = utils::ToLower(parameters["transport_model"].as<std::string>());
    }

    AssignIfPresent(parameters, "cfl", settings.cfl);
    AssignIfPresent(parameters, "global_limiter", settings.global_limiter);
    AssignIfPresent(parameters, "vacuum_fix_limiter", settings.vacuum_fix_limiter);
    AssignIfPresent(parameters, "viscosity", settings.viscosity);
    AssignIfPresent(parameters, "diffusion", settings.diffusion);

    // pressure-velocity solvers
    AssignIfPresent(parameters, "n_pressure_correctors", settings.n_pressure_correctors);
    AssignIfPresent(parameters, "n_outer_correctors", settings.n_outer_correctors);
    AssignIfPresent(parameters, "alpha_u", settings.alpha_u);
    AssignIfPresent(parameters, "alpha_p", settings.alpha_p);
    AssignIfPresent(parameters, "kinematic_viscosity", settings.kinematic_viscosity);
}

void YamlConfigParser::ParseBoundaryConditions(const YAML::Node& node, Settings& settings) {
    if (node["default"]) {
        const std::string default_bc = utils::ToLower(node["default"].as<std::string>());
        settings.left_boundary = default_bc;
        settings.right_boundary = default_bc;
        settings.bottom_boundary = default_bc;
        settings.top_boundary = default_bc;
        settings.back_boundary = default_bc;
        settings.front_boundary = default_bc;
    }

    if (node["x_min"]) settings.left_boundary = utils::ToLower(node["x_min"].as<std::string>());
    if (node["x_max"]) settings.right_boundary = utils::ToLower(node["x_max"].as<std::string>());
    if (node["y_min"]) settings.bottom_boundary = utils::ToLower(node["y_min"].as<std::string>());
    if (node["y_max"]) settings.top_boundary = utils::ToLower(node["y_max"].as<std::string>());
    if (node["z_min"]) settings.back_boundary = utils::ToLower(node["z_min"].as<std::string>());
    if (node["z_max"]) settings.front_boundary = utils::ToLower(node["z_max"].as<std::string>());

    if (node["states"]) {
        const YAML::Node states = node["states"];

        if (states["x_min"]) settings.boundary_states.x_min = ParseBoundaryState(states["x_min"]);
        if (states["x_max"]) settings.boundary_states.x_max = ParseBoundaryState(states["x_max"]);
        if (states["y_min"]) settings.boundary_states.y_min = ParseBoundaryState(states["y_min"]);
        if (states["y_max"]) settings.boundary_states.y_max = ParseBoundaryState(states["y_max"]);
        if (states["z_min"]) settings.boundary_states.z_min = ParseBoundaryState(states["z_min"]);
        if (states["z_max"]) settings.boundary_states.z_max = ParseBoundaryState(states["z_max"]);
    }
}

void YamlConfigParser::ParseParallel(const YAML::Node& node, Settings& settings) {
    AssignIfPresent(node, "mpi", settings.mpi_enabled);
}

void YamlConfigParser::ParseStopping(const YAML::Node& node, Settings& settings) {
    AssignIfPresent(node, "t_end", settings.t_end);
    AssignIfPresent(node, "max_steps", settings.step_end);
}

void YamlConfigParser::ParseLogging(const YAML::Node& node, Settings& settings) {
    AssignIfPresent(node, "every_steps", settings.log_every_steps);
    AssignIfPresent(node, "every_time", settings.log_every_time);
}

void YamlConfigParser::ParseOutput(const YAML::Node& node, Settings& settings) {
    if (node["format"]) {
        settings.output_formats = {utils::ToLower(node["format"].as<std::string>())};
    }

    AssignIfPresent(node, "directory", settings.output_dir);
    AssignIfPresent(node, "every_steps", settings.output_every_steps);
    AssignIfPresent(node, "every_time", settings.output_every_time);
}

void YamlConfigParser::ParseImmersedBoundaries(const YAML::Node& node, Settings& settings) {
    AssignIfPresent(node, "enabled", settings.immersed_enabled);

    if (node["objects"]) {
        settings.immersed_objects = ParseImmersedObjects(node["objects"]);
    }
}

void YamlConfigParser::ApplyCaseOverrides(const YAML::Node& case_node, InitialConditions& ic) {
    CaseSettings& overrides = ic.overrides;

    if (case_node["mesh"]) {
        const YAML::Node mesh_node = case_node["mesh"];

        if (mesh_node["type"]) {
            const std::string mesh_type = utils::ToLower(mesh_node["type"].as<std::string>());
            if (mesh_type != "structured") {
                throw std::runtime_error("Only mesh.type=structured is supported at the moment");
            }
        }

        AssignOptionalIfPresent(mesh_node, "dim", overrides.dim);

        if (mesh_node["cells"]) {
            const YAML::Node cells = mesh_node["cells"];
            AssignOptionalIfPresent(cells, "x", overrides.Nx);
            AssignOptionalIfPresent(cells, "y", overrides.Ny);
            AssignOptionalIfPresent(cells, "z", overrides.Nz);
        }

        if (mesh_node["domain"]) {
            const YAML::Node domain = mesh_node["domain"];
            AssignOptionalIfPresent(domain, "x", overrides.L_x);
            AssignOptionalIfPresent(domain, "y", overrides.L_y);
            AssignOptionalIfPresent(domain, "z", overrides.L_z);
        }

        if (mesh_node["source"]) {
            throw std::runtime_error("Unstructured mesh source is not implemented yet");
        }
    }

    if (case_node["physics"]) {
        const YAML::Node physics_node = case_node["physics"];

        if (physics_node["eos"]) {
            overrides.EOS = utils::ToLower(physics_node["eos"].as<std::string>());
        }

        if (physics_node["parameters"]) {
            const YAML::Node parameters = physics_node["parameters"];
            AssignOptionalIfPresent(parameters, "gamma", overrides.gamma);
            AssignOptionalIfPresent(parameters, "Q_user", overrides.Q_user);
            
            AssignOptionalIfPresent(parameters, "hg_rho0", overrides.hg_rho0);
            AssignOptionalIfPresent(parameters, "hg_C", overrides.hg_C);
            AssignOptionalIfPresent(parameters, "hg_S", overrides.hg_S);
            AssignOptionalIfPresent(parameters, "hg_gamma_s", overrides.hg_gamma_s);
            AssignOptionalIfPresent(parameters, "hg_c_v", overrides.hg_c_v);
        }
    }

    if (case_node["chemistry"]) {
        const YAML::Node chemistry_node = case_node["chemistry"];
        AssignOptionalIfPresent(chemistry_node, "enabled", overrides.chemistry_enabled);
        AssignOptionalIfPresent(chemistry_node, "z_freq", overrides.chemistry_z_freq);
        AssignOptionalIfPresent(chemistry_node, "activation_energy", overrides.chemistry_activation_energy);
        AssignOptionalIfPresent(chemistry_node, "gas_constant", overrides.chemistry_gas_constant);
        AssignOptionalIfPresent(chemistry_node, "heat_release", overrides.chemistry_heat_release);
    }

    if (case_node["numerics"]) {
        const YAML::Node numerics_node = case_node["numerics"];

        if (numerics_node["method"]) {
            overrides.solver = utils::ToLower(numerics_node["method"].as<std::string>());
        }

        if (numerics_node["parameters"]) {
            const YAML::Node parameters = numerics_node["parameters"];

            if (parameters["time_integrator"]) {
                overrides.time_integrator = utils::ToLower(parameters["time_integrator"].as<std::string>());
            }
            if (parameters["reconstruction"]) {
                overrides.reconstruction = utils::ToLower(parameters["reconstruction"].as<std::string>());
            }
            if (parameters["riemann_solver"]) {
                overrides.riemann_solver = utils::ToLower(parameters["riemann_solver"].as<std::string>());
            }
            if (parameters["transport_model"]) {
                overrides.mader_transport = utils::ToLower(parameters["transport_model"].as<std::string>());
            }

            AssignOptionalIfPresent(parameters, "cfl", overrides.cfl);
            AssignOptionalIfPresent(parameters, "global_limiter", overrides.global_limiter);
            AssignOptionalIfPresent(parameters, "vacuum_fix_limiter", overrides.vacuum_fix_limiter);
            AssignOptionalIfPresent(parameters, "viscosity", overrides.viscosity);
            AssignOptionalIfPresent(parameters, "diffusion", overrides.diffusion);

            // pressure-velocity solvers
            AssignOptionalIfPresent(parameters, "n_pressure_correctors", overrides.n_pressure_correctors);
            AssignOptionalIfPresent(parameters, "n_outer_correctors", overrides.n_outer_correctors);
            AssignOptionalIfPresent(parameters, "alpha_u", overrides.alpha_u);
            AssignOptionalIfPresent(parameters, "alpha_p", overrides.alpha_p);
            AssignOptionalIfPresent(parameters, "kinematic_viscosity", overrides.kinematic_viscosity);
        }
    }

    if (case_node["boundary_conditions"]) {
        const YAML::Node bc_node = case_node["boundary_conditions"];

        if (bc_node["default"]) {
            const std::string default_bc = utils::ToLower(bc_node["default"].as<std::string>());
            overrides.left_boundary = default_bc;
            overrides.right_boundary = default_bc;
            overrides.bottom_boundary = default_bc;
            overrides.top_boundary = default_bc;
            overrides.back_boundary = default_bc;
            overrides.front_boundary = default_bc;
        }

        if (bc_node["x_min"]) overrides.left_boundary = utils::ToLower(bc_node["x_min"].as<std::string>());
        if (bc_node["x_max"]) overrides.right_boundary = utils::ToLower(bc_node["x_max"].as<std::string>());
        if (bc_node["y_min"]) overrides.bottom_boundary = utils::ToLower(bc_node["y_min"].as<std::string>());
        if (bc_node["y_max"]) overrides.top_boundary = utils::ToLower(bc_node["y_max"].as<std::string>());
        if (bc_node["z_min"]) overrides.back_boundary = utils::ToLower(bc_node["z_min"].as<std::string>());
        if (bc_node["z_max"]) overrides.front_boundary = utils::ToLower(bc_node["z_max"].as<std::string>());

        if (bc_node["states"]) {
            BoundaryStatesSettings states_override;

            const YAML::Node states = bc_node["states"];
            if (states["x_min"]) states_override.x_min = ParseBoundaryState(states["x_min"]);
            if (states["x_max"]) states_override.x_max = ParseBoundaryState(states["x_max"]);
            if (states["y_min"]) states_override.y_min = ParseBoundaryState(states["y_min"]);
            if (states["y_max"]) states_override.y_max = ParseBoundaryState(states["y_max"]);
            if (states["z_min"]) states_override.z_min = ParseBoundaryState(states["z_min"]);
            if (states["z_max"]) states_override.z_max = ParseBoundaryState(states["z_max"]);

            overrides.boundary_states = states_override;
        }
    }

    if (case_node["parallel"]) {
        const YAML::Node parallel_node = case_node["parallel"];
        AssignOptionalIfPresent(parallel_node, "mpi", overrides.mpi_enabled);
    }

    if (case_node["stopping"]) {
        const YAML::Node stopping_node = case_node["stopping"];
        AssignOptionalIfPresent(stopping_node, "t_end", overrides.t_end);
        AssignOptionalIfPresent(stopping_node, "max_steps", overrides.step_end);
    }

    if (case_node["logging"]) {
        const YAML::Node logging_node = case_node["logging"];
        AssignOptionalIfPresent(logging_node, "every_steps", overrides.log_every_steps);
        AssignOptionalIfPresent(logging_node, "every_time", overrides.log_every_time);
    }

    if (case_node["output"]) {
        const YAML::Node output_node = case_node["output"];

        if (output_node["format"]) {
            overrides.output_formats = std::vector<std::string>{
                utils::ToLower(output_node["format"].as<std::string>())
            };
        }

        AssignOptionalIfPresent(output_node, "directory", overrides.output_dir);
        AssignOptionalIfPresent(output_node, "every_steps", overrides.output_every_steps);
        AssignOptionalIfPresent(output_node, "every_time", overrides.output_every_time);
    }

    if (case_node["immersed_boundaries"]) {
        const YAML::Node immersed_node = case_node["immersed_boundaries"];
        AssignOptionalIfPresent(immersed_node, "enabled", overrides.immersed_enabled);

        if (immersed_node["objects"]) {
            overrides.immersed_objects = ParseImmersedObjects(immersed_node["objects"]);
        }
    }
}

void YamlConfigParser::ParseStructuredInitialCondition(
    const YAML::Node& ic_node,
    InitialConditions& ic,
    const Settings& effective_settings) {
    ValidateStructuredShape(ic_node, effective_settings.dim);

    const YAML::Node interfaces = ic_node["interfaces"];
    ic.interfaces_x = interfaces["x"] ? ReadVectorDouble(interfaces["x"]) : std::vector<double>{};
    ic.interfaces_y = interfaces["y"] ? ReadVectorDouble(interfaces["y"]) : std::vector<double>{};
    ic.interfaces_z = interfaces["z"] ? ReadVectorDouble(interfaces["z"]) : std::vector<double>{};

    if (!ic.interfaces_x.empty()) {
        ic.overrides.x0 = ic.interfaces_x.front();
    }
    if (!ic.interfaces_y.empty()) {
        ic.overrides.y0 = ic.interfaces_y.front();
    }
    if (!ic.interfaces_z.empty()) {
        ic.overrides.z0 = ic.interfaces_z.front();
    }

    if (effective_settings.dim == 1) {
        ParseStructured1D(ic_node, ic);
        return;
    }
    if (effective_settings.dim == 2) {
        ParseStructured2D(ic_node, ic);
        return;
    }
    if (effective_settings.dim == 3) {
        ParseStructured3D(ic_node, ic);
        return;
    }

    throw std::runtime_error("Unsupported dimension in structured initial condition");
}

void YamlConfigParser::ParseStructured1D(const YAML::Node& ic_node, InitialConditions& ic) {
    const auto rho_1d = ReadVectorDouble(ic_node["rho"]);
    const auto u_1d = ReadVectorDouble(ic_node["u"]);
    const auto v_1d = ReadVectorDouble(ic_node["v"]);
    const auto w_1d = ReadVectorDouble(ic_node["w"]);
    const auto p_1d = ReadVectorDouble(ic_node["p"]);

    const std::size_t nx = rho_1d.size();

    if (u_1d.size() != nx || v_1d.size() != nx || w_1d.size() != nx || p_1d.size() != nx) {
        throw std::runtime_error("1D initial-condition arrays must have identical size");
    }

    if (nx != ic.RegionCountX() || ic.RegionCountY() != 1 || ic.RegionCountZ() != 1) {
        throw std::runtime_error("1D initial-condition shape does not match interfaces");
    }

    auto lift_1d = [](const std::vector<double>& src) -> Field3DValues {
        Field3DValues dst;
        dst.values.resize(src.size());
        for (std::size_t ix = 0; ix < src.size(); ++ix) {
            dst.values[ix].resize(1);
            dst.values[ix][0].resize(1);
            dst.values[ix][0][0] = src[ix];
        }
        return dst;
    };

    ic.rho = lift_1d(rho_1d);
    ic.u = lift_1d(u_1d);
    ic.v = lift_1d(v_1d);
    ic.w = lift_1d(w_1d);
    ic.p = lift_1d(p_1d);

    if (ic_node["reactant_mass_fraction"]) {
        const auto lambda_1d = ReadVectorDouble(ic_node["reactant_mass_fraction"]);
        if (lambda_1d.size() != nx) {
            throw std::runtime_error("1D reactant_mass_fraction size must match state size");
        }
        ic.reactant_mass_fraction = lift_1d(lambda_1d);
    }
}

void YamlConfigParser::ParseStructured2D(const YAML::Node& ic_node, InitialConditions& ic) {
    const std::size_t nx = ic.RegionCountX();
    const std::size_t ny = ic.RegionCountY();

    if (nx >= 1 && ny == 1) {
        ParseStructured1D(ic_node, ic);
        return;
    }

    if (nx == 1 && ny >= 1) {
        const auto rho_1d = ReadVectorDouble(ic_node["rho"]);
        const auto u_1d = ReadVectorDouble(ic_node["u"]);
        const auto v_1d = ReadVectorDouble(ic_node["v"]);
        const auto w_1d = ReadVectorDouble(ic_node["w"]);
        const auto p_1d = ReadVectorDouble(ic_node["p"]);

        if (rho_1d.size() != ny || u_1d.size() != ny || v_1d.size() != ny ||
            w_1d.size() != ny || p_1d.size() != ny) {
            throw std::runtime_error("2D y-only initial-condition arrays must match y-region count");
        }

        auto lift_y_only = [&](const std::vector<double>& src) -> Field3DValues {
            Field3DValues dst;
            dst.values.resize(1);
            dst.values[0].resize(src.size());
            for (std::size_t iy = 0; iy < src.size(); ++iy) {
                dst.values[0][iy].resize(1);
                dst.values[0][iy][0] = src[iy];
            }
            return dst;
        };

        ic.rho = lift_y_only(rho_1d);
        ic.u = lift_y_only(u_1d);
        ic.v = lift_y_only(v_1d);
        ic.w = lift_y_only(w_1d);
        ic.p = lift_y_only(p_1d);

        if (ic_node["reactant_mass_fraction"]) {
            const auto lambda_1d = ReadVectorDouble(ic_node["reactant_mass_fraction"]);
            if (lambda_1d.size() != ny) {
                throw std::runtime_error("2D y-only reactant_mass_fraction size must match y-region count");
            }
            ic.reactant_mass_fraction = lift_y_only(lambda_1d);
        }
        return;
    }

    const auto rho_2d = ReadMatrixDouble(ic_node["rho"]);
    const auto u_2d = ReadMatrixDouble(ic_node["u"]);
    const auto v_2d = ReadMatrixDouble(ic_node["v"]);
    const auto w_2d = ReadMatrixDouble(ic_node["w"]);
    const auto p_2d = ReadMatrixDouble(ic_node["p"]);

    auto lift_2d = [&](const std::vector<std::vector<double>>& src) -> Field3DValues {
        if (src.size() != nx) {
            throw std::runtime_error("2D initial-condition x-size does not match interfaces");
        }

        Field3DValues dst;
        dst.values.resize(nx);

        for (std::size_t ix = 0; ix < nx; ++ix) {
            if (src[ix].size() != ny) {
                throw std::runtime_error("2D initial-condition y-size does not match interfaces");
            }

            dst.values[ix].resize(ny);
            for (std::size_t iy = 0; iy < ny; ++iy) {
                dst.values[ix][iy].resize(1);
                dst.values[ix][iy][0] = src[ix][iy];
            }
        }

        return dst;
    };

    ic.rho = lift_2d(rho_2d);
    ic.u = lift_2d(u_2d);
    ic.v = lift_2d(v_2d);
    ic.w = lift_2d(w_2d);
    ic.p = lift_2d(p_2d);

    if (ic_node["reactant_mass_fraction"]) {
        const auto lambda_2d = ReadMatrixDouble(ic_node["reactant_mass_fraction"]);
        ic.reactant_mass_fraction = lift_2d(lambda_2d);
    }
}

void YamlConfigParser::ParseStructured3D(const YAML::Node& ic_node, InitialConditions& ic) {
    const std::size_t nx = ic.RegionCountX();
    const std::size_t ny = ic.RegionCountY();
    const std::size_t nz = ic.RegionCountZ();

    const auto rho_3d = ReadTensorDouble(ic_node["rho"]);
    const auto u_3d = ReadTensorDouble(ic_node["u"]);
    const auto v_3d = ReadTensorDouble(ic_node["v"]);
    const auto w_3d = ReadTensorDouble(ic_node["w"]);
    const auto p_3d = ReadTensorDouble(ic_node["p"]);

    auto lift_3d = [&](const std::vector<std::vector<std::vector<double>>>& src) -> Field3DValues {
        if (src.size() != nx) {
            throw std::runtime_error("3D initial-condition x-size does not match interfaces");
        }

        Field3DValues dst;
        dst.values.resize(nx);

        for (std::size_t ix = 0; ix < nx; ++ix) {
            if (src[ix].size() != ny) {
                throw std::runtime_error("3D initial-condition y-size does not match interfaces");
            }

            dst.values[ix].resize(ny);

            for (std::size_t iy = 0; iy < ny; ++iy) {
                if (src[ix][iy].size() != nz) {
                    throw std::runtime_error("3D initial-condition z-size does not match interfaces");
                }

                dst.values[ix][iy] = src[ix][iy];
            }
        }

        return dst;
    };

    ic.rho = lift_3d(rho_3d);
    ic.u = lift_3d(u_3d);
    ic.v = lift_3d(v_3d);
    ic.w = lift_3d(w_3d);
    ic.p = lift_3d(p_3d);

    if (ic_node["reactant_mass_fraction"]) {
        const auto lambda_3d = ReadTensorDouble(ic_node["reactant_mass_fraction"]);
        ic.reactant_mass_fraction = lift_3d(lambda_3d);
    }
}

void YamlConfigParser::ValidateStructuredShape(const YAML::Node& ic_node, int dim) {
    if (!ic_node["interfaces"]) {
        throw std::runtime_error("Missing 'interfaces' in structured initial condition");
    }

    if (!ic_node["rho"] || !ic_node["u"] || !ic_node["v"] || !ic_node["w"] || !ic_node["p"]) {
        throw std::runtime_error("structured initial condition requires rho/u/v/w/p");
    }

    if (dim < 1 || dim > 3) {
        throw std::runtime_error("Only dim=1,2,3 are supported");
    }
}

auto YamlConfigParser::ParseImmersedObjects(const YAML::Node& node)
    -> std::vector<ImmersedObjectSettings> {
    std::vector<ImmersedObjectSettings> objects;

    if (!node.IsSequence()) {
        throw std::runtime_error("'immersed_boundaries.objects' must be a sequence");
    }

    for (const auto& item : node) {
        ImmersedObjectSettings object;
        AssignIfPresent(item, "type", object.type);
        AssignIfPresent(item, "cx", object.cx);
        AssignIfPresent(item, "cy", object.cy);
        AssignIfPresent(item, "cz", object.cz);
        AssignIfPresent(item, "radius", object.radius);
        AssignIfPresent(item, "size_x", object.size_x);
        AssignIfPresent(item, "size_y", object.size_y);
        AssignIfPresent(item, "size_z", object.size_z);
        objects.push_back(object);
    }

    return objects;
}

auto YamlConfigParser::ReadVectorDouble(const YAML::Node& node) -> std::vector<double> {
    if (!node.IsSequence()) {
        throw std::runtime_error("Expected 1D sequence");
    }

    std::vector<double> values;
    values.reserve(node.size());

    for (const auto& item : node) {
        values.push_back(item.as<double>());
    }

    return values;
}

auto YamlConfigParser::ReadMatrixDouble(const YAML::Node& node)
    -> std::vector<std::vector<double>> {
    if (!node.IsSequence()) {
        throw std::runtime_error("Expected 2D sequence");
    }

    std::vector<std::vector<double>> values;
    values.reserve(node.size());

    for (const auto& row : node) {
        values.push_back(ReadVectorDouble(row));
    }

    return values;
}

auto YamlConfigParser::ReadTensorDouble(const YAML::Node& node)
    -> std::vector<std::vector<std::vector<double>>> {
    if (!node.IsSequence()) {
        throw std::runtime_error("Expected 3D sequence");
    }

    std::vector<std::vector<std::vector<double>>> values;
    values.reserve(node.size());

    for (const auto& plane : node) {
        values.push_back(ReadMatrixDouble(plane));
    }

    return values;
}

auto YamlConfigParser::HasKey(const YAML::Node& node, const char* key) -> bool {
    return static_cast<bool>(node[key]);
}
