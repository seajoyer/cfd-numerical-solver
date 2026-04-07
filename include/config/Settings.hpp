#ifndef SETTINGS_HPP
#define SETTINGS_HPP

#include <cstddef>
#include <map>
#include <optional>
#include <string>
#include <vector>

/**
 * @brief Runtime mesh source type.
 */
enum class MeshSourceType {
    StructuredCartesian,
    GmshFile,
    GmshGeo,
    DelaunayGeo
};

/**
 * @brief Structured Cartesian mesh settings.
 */
struct StructuredMeshSettings {
    int nx = 200;
    int ny = 1;
    int nz = 1;

    double x_min = 0.0;
    double x_max = 1.0;
    double y_min = 0.0;
    double y_max = 1.0;
    double z_min = 0.0;
    double z_max = 1.0;
};

/**
 * @brief External mesh loaded from Gmsh file.
 */
struct GmshFileMeshSettings {
    std::string file_path;
};

/**
 * @brief Geo file passed for in-house Delaunay mesh generation.
 */
struct DelaunayGeoMeshSettings {
    std::string file_path;
};

/**
 * @brief Geo file passed for mesh generation.
 */
struct GmshGeoMeshSettings {
    std::string file_path;
};

/**
 * @brief Generic mesh settings.
 */
struct MeshSettings {
    int dim = 1;
    MeshSourceType source_type = MeshSourceType::StructuredCartesian;

    std::optional<StructuredMeshSettings> structured;
    std::optional<GmshFileMeshSettings> gmsh_file;
    std::optional<GmshGeoMeshSettings> gmsh_geo;
    std::optional<DelaunayGeoMeshSettings> delaunay_geo;
};


/**
 * @brief One primitive boundary state.
 */
struct BoundaryStateSettings {
    double rho = 0.0;
    double u = 0.0;
    double v = 0.0;
    double w = 0.0;
    double p = 0.0;
};

/**
 * @brief One boundary-condition configuration for one boundary tag.
 */
struct BoundaryConditionSettings {
    std::string type;
    std::optional<BoundaryStateSettings> state;
};

/**
 * @brief Boundary setup indexed by boundary tag.
 */
struct BoundarySettings {
    std::map<int, BoundaryConditionSettings> by_tag;
};

/**
 * @brief Immersed object settings reserved for future use.
 */
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
 * @brief Optional per-case overrides applied over global settings.
 */
struct CaseSettings {
    // Solver / numerics
    std::optional<std::string> solver;
    std::optional<std::string> time_integrator;
    std::optional<std::string> riemann_solver;
    std::optional<std::string> reconstruction;
    std::optional<std::string> eos;
    std::optional<std::string> transport_model;

    // Mesh
    std::optional<MeshSettings> mesh;

    // Boundary conditions
    std::optional<BoundarySettings> boundary;

    // Physical / model parameters
    std::optional<double> gamma;
    std::optional<double> Q_user;

    // Numerical flags
    std::optional<double> cfl;
    std::optional<bool> global_limiter;
    std::optional<bool> vacuum_fix_limiter;
    std::optional<bool> viscosity;
    std::optional<bool> diffusion;

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

    // Misc
    std::optional<bool> analytical;
    std::optional<bool> mpi_enabled;
    std::optional<bool> immersed_enabled;
    std::optional<std::vector<ImmersedObjectSettings>> immersed_objects;
};

/**
 * @brief Runtime settings for one simulation case.
 */
struct Settings {
    // ==================== Solver / numerics ====================
    std::string solver = "godunov";
    std::string riemann_solver = "hllc";
    std::string reconstruction = "p0";
    std::string time_integrator = "euler";

    std::optional<std::string> eos;
    std::optional<std::string> transport_model;

    // ==================== Mesh ====================
    MeshSettings mesh;

    // ==================== Boundary conditions ====================
    BoundarySettings boundary;

    // ==================== Physical / model parameters ====================
    double gamma = 1.4;
    double Q_user = 1.0;

    // ==================== Numerical flags ====================
    double cfl = 0.5;
    bool global_limiter = false;
    bool vacuum_fix_limiter = false;
    bool viscosity = false;
    bool diffusion = false;

    // ==================== Case metadata ====================
    std::string simulation_case = "case";
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

    [[nodiscard]] bool HasOutputFormat(const std::string& format) const {
        for (const std::string& fmt : output_formats) {
            if (fmt == format) {
                return true;
            }
        }
        return false;
    }
};

/**
 * @brief Merge case-specific and CLI-specific overrides into one runtime settings object.
 */
inline Settings MergeSettings(const Settings& global,
                              const CaseSettings& case_overrides,
                              const CaseSettings& cli_overrides) {
    Settings merged = global;

#define APPLY_OVERRIDE(field)                 \
    if (case_overrides.field) {               \
        merged.field = *case_overrides.field; \
    }                                         \
    if (cli_overrides.field) {                \
        merged.field = *cli_overrides.field;  \
    }

    APPLY_OVERRIDE(solver)
    APPLY_OVERRIDE(time_integrator)
    APPLY_OVERRIDE(riemann_solver)
    APPLY_OVERRIDE(reconstruction)
    APPLY_OVERRIDE(eos)
    APPLY_OVERRIDE(transport_model)

    APPLY_OVERRIDE(mesh)
    APPLY_OVERRIDE(boundary)

    APPLY_OVERRIDE(gamma)
    APPLY_OVERRIDE(Q_user)

    APPLY_OVERRIDE(cfl)
    APPLY_OVERRIDE(global_limiter)
    APPLY_OVERRIDE(vacuum_fix_limiter)
    APPLY_OVERRIDE(viscosity)
    APPLY_OVERRIDE(diffusion)

    APPLY_OVERRIDE(t_end)
    APPLY_OVERRIDE(step_end)

    APPLY_OVERRIDE(log_every_steps)
    APPLY_OVERRIDE(log_every_time)

    APPLY_OVERRIDE(output_every_steps)
    APPLY_OVERRIDE(output_every_time)
    APPLY_OVERRIDE(output_formats)
    APPLY_OVERRIDE(output_dir)

    APPLY_OVERRIDE(analytical)
    APPLY_OVERRIDE(mpi_enabled)
    APPLY_OVERRIDE(immersed_enabled)
    APPLY_OVERRIDE(immersed_objects)

#undef APPLY_OVERRIDE

    return merged;
}

#endif  // SETTINGS_HPP
