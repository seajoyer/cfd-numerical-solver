#ifndef SETTINGS_HPP
#define SETTINGS_HPP

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

struct ImmersedObjectSettings {
    std::string type;

    double cx = 0.0;
    double cy = 0.0;
    double cz = 0.0;

    double radius = 0.0;

    double size_x = 0.0;
    double size_y = 0.0;
    double size_z = 0.0;
};



struct BoundaryStateSettings {
    double rho = 0.0;
    double u = 0.0;
    double v = 0.0;
    double w = 0.0;
    double p = 0.0;
};

struct BoundaryStatesSettings {
    std::optional<BoundaryStateSettings> x_min;
    std::optional<BoundaryStateSettings> x_max;
    std::optional<BoundaryStateSettings> y_min;
    std::optional<BoundaryStateSettings> y_max;
    std::optional<BoundaryStateSettings> z_min;
    std::optional<BoundaryStateSettings> z_max;
};

/**
 * @struct CaseSettings
 * @brief Optional per-case overrides applied over global settings
 */
struct CaseSettings {
    // Solver / numerics
    std::optional<std::string> solver;
    std::optional<std::string> time_integrator;
    std::optional<std::string> riemann_solver;
    std::optional<std::string> reconstruction;
    std::optional<std::string> EOS;
    std::optional<std::string> mader_transport;

    // Boundary conditions
    std::optional<std::string> left_boundary;
    std::optional<std::string> right_boundary;
    std::optional<std::string> bottom_boundary;
    std::optional<std::string> top_boundary;
    std::optional<std::string> back_boundary;
    std::optional<std::string> front_boundary;

    // Grid / geometry
    std::optional<int> Nx;
    std::optional<int> Ny;
    std::optional<int> Nz;
    std::optional<int> dim;
    std::optional<double> L_x;
    std::optional<double> L_y;
    std::optional<double> L_z;

    // Physical / model parameters
    std::optional<double> gamma;
    std::optional<double> Q_user;

    std::optional<bool> chemistry_enabled;
    std::optional<double> chemistry_z_freq;
    std::optional<double> chemistry_activation_energy;
    std::optional<double> chemistry_gas_constant;
    std::optional<double> chemistry_heat_release;

    std::optional<double> cfl;
    std::optional<bool> global_limiter;
    std::optional<bool> vacuum_fix_limiter;
    std::optional<bool> viscosity;
    std::optional<bool> diffusion;
    std::optional<bool> analytical;

    // IC interfaces
    std::optional<double> x0;
    std::optional<double> y0;
    std::optional<double> z0;

    // Parallel / immersed
    std::optional<bool> mpi_enabled;
    std::optional<bool> immersed_enabled;
    std::optional<std::vector<ImmersedObjectSettings>> immersed_objects;

    // Run control
    std::optional<double> t_end;
    std::optional<std::size_t> step_end;

    // Logging
    std::optional<std::size_t> log_every_steps;
    std::optional<double> log_every_time;

    // Output
    std::optional<std::size_t> output_every_steps;
    std::optional<double> output_every_time;
    std::optional<std::vector<std::string>> output_formats;
    std::optional<std::string> output_dir;

    std::optional<BoundaryStatesSettings> boundary_states;
};

/**
 * @struct Settings
 * @brief Runtime settings for one simulation case
 */
struct Settings {
    // ==================== Solver / numerics ====================
    std::string solver = "godunov";
    std::string riemann_solver = "exact";
    std::string reconstruction = "p0";
    std::string time_integrator = "euler";

    std::optional<std::string> EOS;
    std::optional<std::string> mader_transport;

    // ==================== Boundary conditions ====================
    std::string left_boundary = "free_stream";
    std::string right_boundary = "free_stream";
    std::string bottom_boundary = "free_stream";
    std::string top_boundary = "free_stream";
    std::string back_boundary = "free_stream";
    std::string front_boundary = "free_stream";

    BoundaryStatesSettings boundary_states;

    // ==================== Grid ====================
    int Nx = 200;
    int Ny = 1;
    int Nz = 1;
    int dim = 1;

    double L_x = 1.0;
    double L_y = 1.0;
    double L_z = 1.0;

    // ==================== Physical / model parameters ====================
    double gamma = 1.4;
    double Q_user = 1.0;

    bool chemistry_enabled = false;
    double chemistry_z_freq = 0.0;
    double chemistry_activation_energy = 0.0;
    double chemistry_gas_constant = 1.0;
    double chemistry_heat_release = 0.0;

    double cfl = 0.5;
    bool global_limiter = false;
    bool vacuum_fix_limiter = false;
    bool viscosity = false;
    bool diffusion = false;

    // ==================== Initial-condition related ====================
    std::string simulation_case = "case";
    double x0 = 0.5;
    double y0 = 0.5;
    double z0 = 0.5;
    bool analytical = false;

    // ==================== Run control ====================
    double t_end = 0.0;
    std::size_t step_end = 0;

    // ==================== Logging ====================
    std::size_t log_every_steps = 1;
    double log_every_time = 0.0;

    // ==================== Output ====================
    std::size_t output_every_steps = 1;
    double output_every_time = 0.0;
    std::vector<std::string> output_formats = {"vtk"};
    std::string output_dir = "data/output";

    // ==================== Parallel / immersed ====================
    bool mpi_enabled = false;
    bool immersed_enabled = false;
    std::vector<ImmersedObjectSettings> immersed_objects;

    [[nodiscard]] auto GetNx() const -> int {
        return Nx;
    }

    [[nodiscard]] auto GetNy() const -> int {
        return Ny;
    }

    [[nodiscard]] auto GetNz() const -> int {
        return Nz;
    }

    [[nodiscard]] auto HasOutputFormat(const std::string& format) const -> bool {
        for (const auto& fmt : output_formats) {
            if (fmt == format) {
                return true;
            }
        }
        return false;
    }
};

/**
 * @brief Merges case-specific overrides into global settings
 */
inline auto MergeSettings(const Settings& global,
                          const CaseSettings& case_overrides,
                          const CaseSettings& cli_overrides) -> Settings {
    Settings merged = global;

#define APPLY_OVERRIDE(field)                  \
    if (case_overrides.field) {                \
        merged.field = *case_overrides.field;  \
    }                                          \
    if (cli_overrides.field) {                 \
        merged.field = *cli_overrides.field;   \
    }

    APPLY_OVERRIDE(solver)
    APPLY_OVERRIDE(time_integrator)
    APPLY_OVERRIDE(riemann_solver)
    APPLY_OVERRIDE(reconstruction)
    APPLY_OVERRIDE(EOS)
    APPLY_OVERRIDE(mader_transport)

    APPLY_OVERRIDE(left_boundary)
    APPLY_OVERRIDE(right_boundary)
    APPLY_OVERRIDE(bottom_boundary)
    APPLY_OVERRIDE(top_boundary)
    APPLY_OVERRIDE(back_boundary)
    APPLY_OVERRIDE(front_boundary)

    APPLY_OVERRIDE(boundary_states)

    APPLY_OVERRIDE(Nx)
    APPLY_OVERRIDE(Ny)
    APPLY_OVERRIDE(Nz)
    APPLY_OVERRIDE(dim)
    APPLY_OVERRIDE(L_x)
    APPLY_OVERRIDE(L_y)
    APPLY_OVERRIDE(L_z)

    APPLY_OVERRIDE(gamma)
    APPLY_OVERRIDE(Q_user)

    APPLY_OVERRIDE(chemistry_enabled)
    APPLY_OVERRIDE(chemistry_z_freq)
    APPLY_OVERRIDE(chemistry_activation_energy)
    APPLY_OVERRIDE(chemistry_gas_constant)
    APPLY_OVERRIDE(chemistry_heat_release)

    APPLY_OVERRIDE(cfl)
    APPLY_OVERRIDE(global_limiter)
    APPLY_OVERRIDE(vacuum_fix_limiter)
    APPLY_OVERRIDE(viscosity)
    APPLY_OVERRIDE(diffusion)

    APPLY_OVERRIDE(x0)
    APPLY_OVERRIDE(y0)
    APPLY_OVERRIDE(z0)
    APPLY_OVERRIDE(analytical)

    APPLY_OVERRIDE(mpi_enabled)
    APPLY_OVERRIDE(immersed_enabled)
    APPLY_OVERRIDE(immersed_objects)

    APPLY_OVERRIDE(t_end)
    APPLY_OVERRIDE(step_end)

    APPLY_OVERRIDE(log_every_steps)
    APPLY_OVERRIDE(log_every_time)

    APPLY_OVERRIDE(output_every_steps)
    APPLY_OVERRIDE(output_every_time)
    APPLY_OVERRIDE(output_formats)
    APPLY_OVERRIDE(output_dir)

#undef APPLY_OVERRIDE

    return merged;
}

#endif  // SETTINGS_HPP
