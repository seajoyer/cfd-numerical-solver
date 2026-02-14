#ifndef SETTINGS_HPP
#define SETTINGS_HPP

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

/**
 * @struct CaseSettings
 * @brief Optional overrides for case-specific settings
 */
struct CaseSettings {
    // Solver Configuration
    std::optional<std::string> solver;
    std::optional<std::string> riemann_solver;
    std::optional<std::string> reconstruction;
    std::optional<std::string> left_boundary;
    std::optional<std::string> right_boundary;
    std::optional<std::string> bottom_boundary;  // 2D: y-min boundary
    std::optional<std::string> top_boundary;      // 2D: y-max boundary
    
    // Grid Configuration
    std::optional<int> N;
    std::optional<int> Nx;    // 2D: cells in x
    std::optional<int> Ny;    // 2D: cells in y
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
    std::optional<double> y0;   // 2D: discontinuity y-position
    std::optional<bool> analytical;
    
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
 * @brief Container for all simulation configuration parameters (1D and 2D)
 */
struct Settings {
    // ==================== Solver Configuration ====================
    std::string solver = "godunov";
    std::string riemann_solver = "exact";
    std::string reconstruction = "p0";
    std::string time_integrator = "euler";

    std::string left_boundary = "free_stream";
    std::string right_boundary = "free_stream";
    std::string bottom_boundary = "free_stream";  // 2D: y-min
    std::string top_boundary = "free_stream";      // 2D: y-max

    // ==================== Grid Configuration ====================
    int N = 200;     ///< Default cells (used when Nx/Ny not specified)
    int Nx = 0;      ///< Cells in x-direction (0 = use N)
    int Ny = 0;      ///< Cells in y-direction (0 = use N)
    double cfl = 0.5;
    int padding = 2;
    double gamma = 1.4;
    int dim = 1;
    double L_x = 10.0;
    double L_y = 10.0;
    double L_z = 10.0;

    // ==================== Physical Parameters ====================
    double Q_user = 2.0;

    // ========================== Limiter =========================
    bool global_limiter = false;
    bool vacuum_fix_limiter = false;
    bool viscosity = false;
    bool diffusion = false;

    // ==================== Initial Conditions ====================
    std::string simulation_case = "sod1";
    double x0 = 0.5;    ///< Discontinuity x-position
    double y0 = 0.5;    ///< Discontinuity y-position (2D)
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

    /**
     * @brief Returns effective Nx (uses N if Nx not set)
     */
    [[nodiscard]] auto GetNx() const -> int { return (Nx > 0) ? Nx : N; }

    /**
     * @brief Returns effective Ny (uses N if Ny not set)
     */
    [[nodiscard]] auto GetNy() const -> int { return (Ny > 0) ? Ny : N; }

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
    APPLY_OVERRIDE(riemann_solver)
    APPLY_OVERRIDE(reconstruction)
    APPLY_OVERRIDE(left_boundary)
    APPLY_OVERRIDE(right_boundary)
    APPLY_OVERRIDE(bottom_boundary)
    APPLY_OVERRIDE(top_boundary)
    APPLY_OVERRIDE(N)
    APPLY_OVERRIDE(Nx)
    APPLY_OVERRIDE(Ny)
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
    APPLY_OVERRIDE(analytical)
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
