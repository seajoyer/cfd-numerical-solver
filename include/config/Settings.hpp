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

/**
 * @struct CaseSettings
 * @brief Optional overrides for case-specific settings
 */
struct CaseSettings {
    // Solver Configuration
    std::optional<std::string> solver;
    std::optional<std::string> time_integrator;
    std::optional<std::string> riemann_solver;
    std::optional<std::string> reconstruction;
    std::optional<std::string> left_boundary;
    std::optional<std::string> right_boundary;
    std::optional<std::string> bottom_boundary; // 2D: y-min boundary
    std::optional<std::string> top_boundary; // 2D: y-max boundary
    std::optional<std::string> back_boundary; // 3D: z-min boundary
    std::optional<std::string> front_boundary; // 3D: z-max boundary

    // Grid Configuration
    std::optional<int> N;
    std::optional<int> Nx; // 3D: cells in x
    std::optional<int> Ny; // 3D: cells in y
    std::optional<int> Nz; // 3D: cells in z
    std::optional<double> cfl;
    std::optional<int> padding;
    std::optional<double> gamma;
    std::optional<int> dim;
    std::optional<double> L_x;
    std::optional<double> L_y;
    std::optional<double> L_z;

    // Physical Parameters
    std::optional<double> Q_user;

    // Initial Conditions
    std::optional<double> x0;
    std::optional<double> y0;
    std::optional<double> z0;
    std::optional<bool> analytical;

    std::optional<bool> immersed_enabled;
    std::optional<bool> mpi_enabled;

    std::optional<std::vector<ImmersedObjectSettings>> immersed_objects;

    // Time Control
    std::optional<double> t_end;
    std::optional<std::size_t> step_end;

    // Logging Configuration
    std::optional<std::size_t> log_every_steps;
    std::optional<double> log_every_time;

    // Output Configuration
    std::optional<std::size_t> output_every_steps;
    std::optional<double> output_every_time;
    std::optional<std::vector<std::string>> output_formats;
    std::optional<std::string> output_dir;
};

/**
 * @struct Settings
 * @brief Container for all simulation configuration parameters
 */
struct Settings {
    // ==================== Solver Configuration ====================
    std::string solver = "godunov";
    std::string riemann_solver = "exact";
    std::string reconstruction = "p0";
    std::string time_integrator = "euler";

    std::string left_boundary = "free_stream";
    std::string right_boundary = "free_stream";
    std::string bottom_boundary = "free_stream";
    std::string top_boundary = "free_stream";
    std::string back_boundary = "free_stream";
    std::string front_boundary = "free_stream";

    // ==================== Grid Configuration ====================
    int N = 200; ///< Default cells (used when Nx/Ny/Nz not specified)
    int Nx = 1; ///< Cells in x-direction
    int Ny = 1; ///< Cells in y-direction
    int Nz = 1; ///< Cells in z-direction
    double cfl = 0.5;
    int padding = 2;
    double gamma = 1.4;
    int dim = 1;
    double L_x = 1.0;
    double L_y = 1.0;
    double L_z = 1.0;

    // ==================== Physical Parameters ====================
    double Q_user = 2.0;

    // ========================== Limiter =========================
    bool global_limiter = false;
    bool vacuum_fix_limiter = false;
    bool viscosity = false;
    bool diffusion = false;

    // ==================== Initial Conditions ====================
    std::string simulation_case = "sod1";
    double x0 = 0.5; ///< Discontinuity x-position
    double y0 = 0.5; ///< Discontinuity y-position
    double z0 = 0.5; ///< Discontinuity z-position
    bool analytical = false;

    // ==================== Time Control ====================
    double t_end = 0.0;
    std::size_t step_end = 0;

    // ==================== Logging Configuration ====================
    std::size_t log_every_steps = 1;
    double log_every_time = 0.0;

    // ==================== Output Configuration ====================
    std::size_t output_every_steps = 1;
    double output_every_time = 0.0;
    std::vector<std::string> output_formats = {"vtk"};
    std::string output_dir = "data/output";

    bool immersed_enabled = false;
    bool mpi_enabled = false;

    std::vector<ImmersedObjectSettings> immersed_objects;

    /**
     * @brief Returns effective Nx (uses N if Nx not set)
     */
    [[nodiscard]] auto GetNx() const -> int {
        return Nx > 0 ? Nx : N;
    }

    /**
     * @brief Returns effective Ny (uses N if Ny not set)
     */
    [[nodiscard]] auto GetNy() const -> int {
        return Ny > 0 ? Ny : N;
    }

    /**
     * @brief Returns effective Nz (uses N if Nz not set)
     */
    [[nodiscard]] auto GetNz() const -> int {
        return Nz > 0 ? Nz : N;
    }

    [[nodiscard]] auto HasOutputFormat(const std::string& format) const -> bool {
        for (const auto& fmt : output_formats) {
            if (fmt == format) return true;
        }
        return false;
    }
};

/**
 * @brief Merges case-specific and CLI overrides into global settings
 */
inline auto MergeSettings(const Settings& global,
                          const CaseSettings& case_overrides,
                          const CaseSettings& cli_overrides) -> Settings {
    Settings merged = global;

    // Helper macro for applying overrides in order
#define APPLY_OVERRIDE(field) \
        if (case_overrides.field) merged.field = *case_overrides.field; \
        if (cli_overrides.field) merged.field = *cli_overrides.field;

    APPLY_OVERRIDE(solver)
    APPLY_OVERRIDE(time_integrator)
    APPLY_OVERRIDE(riemann_solver)
    APPLY_OVERRIDE(reconstruction)
    APPLY_OVERRIDE(left_boundary)
    APPLY_OVERRIDE(right_boundary)
    APPLY_OVERRIDE(bottom_boundary)
    APPLY_OVERRIDE(top_boundary)
    APPLY_OVERRIDE(back_boundary)
    APPLY_OVERRIDE(front_boundary)
    APPLY_OVERRIDE(N)
    APPLY_OVERRIDE(Nx)
    APPLY_OVERRIDE(Ny)
    APPLY_OVERRIDE(Nz)
    APPLY_OVERRIDE(cfl)
    APPLY_OVERRIDE(padding)
    APPLY_OVERRIDE(gamma)
    APPLY_OVERRIDE(dim)
    APPLY_OVERRIDE(L_x)
    APPLY_OVERRIDE(L_y)
    APPLY_OVERRIDE(L_z)
    APPLY_OVERRIDE(Q_user)
    APPLY_OVERRIDE(x0)
    APPLY_OVERRIDE(y0)
    APPLY_OVERRIDE(z0)
    APPLY_OVERRIDE(analytical)
    APPLY_OVERRIDE(t_end)
    APPLY_OVERRIDE(step_end)
    APPLY_OVERRIDE(log_every_steps)
    APPLY_OVERRIDE(log_every_time)
    APPLY_OVERRIDE(output_every_steps)
    APPLY_OVERRIDE(output_every_time)
    APPLY_OVERRIDE(output_formats)
    APPLY_OVERRIDE(output_dir)
    APPLY_OVERRIDE(immersed_enabled)
    APPLY_OVERRIDE(mpi_enabled)
    APPLY_OVERRIDE(immersed_objects)

#undef APPLY_OVERRIDE

    return merged;
}

#endif  // SETTINGS_HPP
