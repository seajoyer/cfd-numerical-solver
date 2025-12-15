#ifndef SETTINGS_HPP
#define SETTINGS_HPP

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

/**
 * @file Settings.hpp
 * @brief Global configuration parameters for CFD simulation
 */

/**
 * @struct CaseSettings
 * @brief Optional overrides for case-specific settings
 * 
 * Each field is optional. When set, it overrides the corresponding
 * global setting for that specific simulation case.
 */
struct CaseSettings {
    // Solver Configuration
    std::optional<std::string> solver;
    std::optional<std::string> riemann_solver;
    std::optional<std::string> reconstruction;
    std::optional<std::string> left_boundary;
    std::optional<std::string> right_boundary;
    
    // Grid Configuration
    std::optional<int> N;
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
 * @brief Container for all simulation configuration parameters
 * 
 * This structure holds all parameters that control the simulation behavior,
 * including solver selection, grid configuration, physical parameters,
 * output settings, and initial conditions.
 * 
 * Parameters are typically loaded from a YAML configuration file via
 * ConfigParser, but can also be set programmatically or overridden via
 * command-line arguments.
 * 
 * @see ConfigParser
 */
struct Settings {
    // ==================== Solver Configuration ====================
    
    /** @brief Solver type identifier
     * 
     * Supported values:
     * - "godunov": First-order Godunov method with P0 reconstruction
     * - "godunov-kolgan": Godunov method with P1 reconstruction
     * - "godunov-kolgan-rodionov": Second-order MUSCL-Hancock scheme
     * - "analytical": Exact Riemann solver (for verification only)
     * 
     * @note Case-insensitive
     */
    std::string solver = "godunov";
    
    /** @brief Riemann solver type
     * 
     * Supported values:
     * - "exact": Exact ideal gas Riemann solver (iterative)
     * - "hll": HLL approximate Riemann solver
     * - "hllc": HLLC approximate Riemann solver (with contact wave)
     * - "acoustic": Linearized acoustic Riemann solver
     * 
     * @note Case-insensitive
     */
    std::string riemann_solver = "exact";
    
    /** @brief Spatial reconstruction scheme
     * 
     * Supported values:
     * - "p0": Piecewise constant (first-order)
     * - "p1": Piecewise linear with slope limiting (second-order)
     * - "eno3", "eno5", etc.: ENO reconstruction of specified order
     * - "weno3", "weno5": WENO reconstruction of specified order
     * 
     * @note Case-insensitive
     */
    std::string reconstruction = "p0";

    /** @brief Time integration scheme identifier
     *
     * Supported values:
     * - "euler"      : First-order Forward Euler scheme
     * - "ssprk2"     : Second-order SSP Runge–Kutta (TVD RK2)
     * - "ssprk3"     : Third-order SSP Runge–Kutta (SSPRK(3,3))
     * - "maccormack" : Two-step MacCormack predictor–corrector scheme
     *
     * @note Case-insensitive
     * @note Used by solver factory to construct the appropriate TimeIntegrator
     */
    std::string time_integrator = "euler";

    /** @brief Left boundary condition type
     * 
     * Supported values:
     * - "free_stream": Far-field boundary with specified external state
     * - "inlet": Inflow boundary with prescribed values
     * - "outlet": Outflow boundary (zero-gradient extrapolation)
     * - "reflective": Reflecting wall (inverts normal velocity)
     * - "non_reflective": Non-reflecting characteristic boundary
     * - "periodic": Periodic boundary
     * - "symmetry": Symmetry plane
     * - "wall": Solid wall (zero normal velocity)
     * 
     * @see BoundaryCondition, BoundaryFactory
     */
    std::string left_boundary = "free_stream";
    
    /** @brief Right boundary condition type
     * @see left_boundary for supported values
     */
    std::string right_boundary = "free_stream";

    // ==================== Grid Configuration ====================
    
    /** @brief Number of physical grid cells (excluding ghost cells)
     * @note Must be positive
     */
    int N = 200;
    
    /** @brief CFL number for time step calculation
     * @note Should be in range (0, 1] for stability
     * @note Typical values: 0.4-0.9 for explicit schemes
     */
    double cfl = 0.5;
    
    /** @brief Number of ghost cells on each boundary
     * @note Must be non-negative
     * @note Should be at least equal to reconstruction stencil size
     */
    int padding = 2;
    
    /** @brief Ratio of specific heat
     * @note For ideal gas EOS: γ = cp/cv
     * @note Common values: 1.4 (air), 1.67 (monatomic), 1.33 (polyatomic)
     */
    double gamma = 1.4;
    
    /** @brief Spatial dimension of the problem
     * @note Currently only 1D is fully supported
     * @note Valid values: 1, 2, 3
     */
    int dim = 1;
    
    /** @brief Domain length in x-direction */
    double L_x = 10.0;
    
    /** @brief Domain length in y-direction (for future 2D/3D support) */
    double L_y = 10.0;
    
    /** @brief Domain length in z-direction (for future 3D support) */
    double L_z = 10.0;

    // ==================== Physical Parameters ====================
    
    /** @brief User-defined parameter Q for exact Riemann solver
     * @note Controls initial pressure ratio threshold in PVRS guess
     * @note Typical value: 2.0
     */
    double Q_user = 2.0;

    // ========================== Limiter =========================

    /** @brief Enable global post-step limiter (simple minmod-like energy smoothing). */
    bool global_limiter = false;

    /** @brief Enable vacuum-fix limiter that caps specific internal energy in near-vacuum regions. */
    bool vacuum_fix_limiter = false;

    /** @brief Add Von Neumann Richtmyer artificial viscosity. */
    bool viscosity = false;

    /** @brief Add Georgievskiy-Levin-Sutyrin diffusion / anti-diffusion. */
    bool diffusion = false;

    // ==================== Initial Conditions ====================
    
    /** @brief Name of the simulation case to run
     * 
     * Should match a key in the initial_conditions section of config.yaml,
     * or use "all" to run all available cases sequentially.
     * 
     * Examples: "sod1", "sod2", "blast_wave", "all"
     * 
     * @see InitialConditions
     */
    std::string simulation_case = "sod1";
    
    /** @brief Position of initial discontinuity (x-coordinate)
     * @note For Riemann problems, separates left and right states
     */
    double x0 = 0.5;
    
    /** @brief Enable analytical solution computation
     * 
     * When true, runs exact Riemann solver in parallel with numerical
     * solution for comparison and verification purposes.
     * 
     * @note Only applicable for Riemann problem initial conditions
     */
    bool analytical = false;

    // ==================== Time Control ====================
    
    /** @brief Final simulation time
     * @note Set to 0 to disable time-based stopping criterion
     */
    double t_end = 0.0;
    
    /** @brief Maximum number of time steps
     * @note Set to 0 to disable step-based stopping criterion
     */
    std::size_t step_end = 0;

    // ==================== Logging Configuration ====================
    
    /** @brief Print progress log every N steps
     * @note Set to 0 to disable step-based logging
     */
    std::size_t log_every_steps = 1;
    
    /** @brief Print progress log every T time units
     * @note Set to 0.0 to disable time-based logging
     */
    double log_every_time = 0.0;

    // ==================== Output Configuration ====================
    
    /** @brief Write output file every N steps
     * @note Set to 0 to disable step-based output
     */
    std::size_t output_every_steps = 1;
    
    /** @brief Write output file every T time units
     * @note Set to 0.0 to disable time-based output
     */
    double output_every_time = 0.0;
    
    /** @brief Output file formats
     *
     * Supported values:
     * - "vtk": VTK structured grid format
     * - "png": PNG image with 4 subplot panels
     *
     * Multiple formats can be specified simultaneously.
     *
     * @see VTKWriter, PNGWriter
     */
    std::vector<std::string> output_formats = {"vtk"};
    
    /** @brief Base directory for output files
     * 
     * Output structure:
     * ```
     * output_dir/run_<timestamp>/
     *   └── <case_name>/
     *       ├── png/
     *       │   └── step_NNNN.png
     *       └── vtk/
     *           ├── analytical/
     *           │   └── step_NNNN.vtk
     *           └── solver__R_reconstruction__N_size__CFL_value/
     *               └── solver__...__step_NNNN.vtk
     * ```
     */
    std::string output_dir = "data/output";

    /**
     * @brief Check if a specific output format is enabled
     * @param format Format to check (e.g., "vtk", "png")
     * @return true if the format is in output_formats
     */
    [[nodiscard]] auto HasOutputFormat(const std::string& format) const -> bool {
        for (const auto& fmt : output_formats) {
            if (fmt == format) return true;
        }
        return false;
    }
};

/**
 * @brief Merges case-specific and CLI overrides into global settings
 * 
 * Creates a new Settings object by applying overrides in priority order:
 * 1. Start with global settings
 * 2. Apply case-specific overrides (from YAML)
 * 3. Apply CLI overrides (highest priority - always win)
 * 
 * This ensures CLI arguments always take precedence over config file values.
 * 
 * @param global Global settings
 * @param case_overrides Case-specific overrides from YAML
 * @param cli_overrides CLI overrides (highest priority)
 * @return Merged settings with proper precedence
 */
inline auto MergeSettings(const Settings& global, 
                          const CaseSettings& case_overrides,
                          const CaseSettings& cli_overrides) -> Settings {
    Settings merged = global;
    
    // Apply case-specific overrides first (lower priority)
    // Solver Configuration
    if (case_overrides.solver) merged.solver = *case_overrides.solver;
    if (case_overrides.riemann_solver) merged.riemann_solver = *case_overrides.riemann_solver;
    if (case_overrides.reconstruction) merged.reconstruction = *case_overrides.reconstruction;
    if (case_overrides.left_boundary) merged.left_boundary = *case_overrides.left_boundary;
    if (case_overrides.right_boundary) merged.right_boundary = *case_overrides.right_boundary;
    
    // Grid Configuration
    if (case_overrides.N) merged.N = *case_overrides.N;
    if (case_overrides.cfl) merged.cfl = *case_overrides.cfl;
    if (case_overrides.padding) merged.padding = *case_overrides.padding;
    if (case_overrides.gamma) merged.gamma = *case_overrides.gamma;
    if (case_overrides.dim) merged.dim = *case_overrides.dim;
    if (case_overrides.L_x) merged.L_x = *case_overrides.L_x;
    if (case_overrides.L_y) merged.L_y = *case_overrides.L_y;
    if (case_overrides.L_z) merged.L_z = *case_overrides.L_z;
    
    // Physical Parameters
    if (case_overrides.Q_user) merged.Q_user = *case_overrides.Q_user;
    
    // Initial Conditions
    if (case_overrides.x0) merged.x0 = *case_overrides.x0;
    if (case_overrides.analytical) merged.analytical = *case_overrides.analytical;
    
    // Time Control
    if (case_overrides.t_end) merged.t_end = *case_overrides.t_end;
    if (case_overrides.step_end) merged.step_end = *case_overrides.step_end;
    
    // Logging Configuration
    if (case_overrides.log_every_steps) merged.log_every_steps = *case_overrides.log_every_steps;
    if (case_overrides.log_every_time) merged.log_every_time = *case_overrides.log_every_time;
    
    // Output Configuration
    if (case_overrides.output_every_steps) merged.output_every_steps = *case_overrides.output_every_steps;
    if (case_overrides.output_every_time) merged.output_every_time = *case_overrides.output_every_time;
    if (case_overrides.output_formats) merged.output_formats = *case_overrides.output_formats;
    if (case_overrides.output_dir) merged.output_dir = *case_overrides.output_dir;
    
    // Apply CLI overrides last (highest priority - always win)
    // Solver Configuration
    if (cli_overrides.solver) merged.solver = *cli_overrides.solver;
    if (cli_overrides.riemann_solver) merged.riemann_solver = *cli_overrides.riemann_solver;
    if (cli_overrides.reconstruction) merged.reconstruction = *cli_overrides.reconstruction;
    if (cli_overrides.left_boundary) merged.left_boundary = *cli_overrides.left_boundary;
    if (cli_overrides.right_boundary) merged.right_boundary = *cli_overrides.right_boundary;
    
    // Grid Configuration
    if (cli_overrides.N) merged.N = *cli_overrides.N;
    if (cli_overrides.cfl) merged.cfl = *cli_overrides.cfl;
    if (cli_overrides.padding) merged.padding = *cli_overrides.padding;
    if (cli_overrides.gamma) merged.gamma = *cli_overrides.gamma;
    if (cli_overrides.dim) merged.dim = *cli_overrides.dim;
    if (cli_overrides.L_x) merged.L_x = *cli_overrides.L_x;
    if (cli_overrides.L_y) merged.L_y = *cli_overrides.L_y;
    if (cli_overrides.L_z) merged.L_z = *cli_overrides.L_z;
    
    // Physical Parameters
    if (cli_overrides.Q_user) merged.Q_user = *cli_overrides.Q_user;
    
    // Initial Conditions
    if (cli_overrides.x0) merged.x0 = *cli_overrides.x0;
    if (cli_overrides.analytical) merged.analytical = *cli_overrides.analytical;
    
    // Time Control
    if (cli_overrides.t_end) merged.t_end = *cli_overrides.t_end;
    if (cli_overrides.step_end) merged.step_end = *cli_overrides.step_end;
    
    // Logging Configuration
    if (cli_overrides.log_every_steps) merged.log_every_steps = *cli_overrides.log_every_steps;
    if (cli_overrides.log_every_time) merged.log_every_time = *cli_overrides.log_every_time;
    
    // Output Configuration
    if (cli_overrides.output_every_steps) merged.output_every_steps = *cli_overrides.output_every_steps;
    if (cli_overrides.output_every_time) merged.output_every_time = *cli_overrides.output_every_time;
    if (cli_overrides.output_formats) merged.output_formats = *cli_overrides.output_formats;
    if (cli_overrides.output_dir) merged.output_dir = *cli_overrides.output_dir;
    
    return merged;
}

#endif  // SETTINGS_HPP
