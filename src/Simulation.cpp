#include "Simulation.hpp"

#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <chrono>

#include "bc/BoundaryFactory.hpp"
#include "utils/StringUtils.hpp"
#include "output/WriterFactory.hpp"
#include "solver/AnalyticalSolver.hpp"
#include "solver/SolverFactory.hpp"

Simulation::Simulation(Settings settings, const InitialConditions& initial_conditions)
    : settings_(std::move(settings)), initial_conditions_(initial_conditions) {}

auto Simulation::CreateSolver() -> std::unique_ptr<Solver> {
    return SolverFactory::Create(settings_);
}

auto Simulation::CreateBoundaryCondition(const std::string& boundary_type, double rho_inf,
                                         double u_inf, double p_inf)
    -> std::shared_ptr<BoundaryCondition> {
    return BoundaryFactory::Create(boundary_type, rho_inf, u_inf, p_inf);
}

void Simulation::CreateWriters() {
    // Base output directory for this case is already set in settings_.output_dir
    case_output_dir_ = settings_.output_dir;
    
    // Create writers for each enabled format
    for (const auto& format : settings_.output_formats) {
        // Convert to lowercase for comparison
        std::string format_lower = format;
        std::transform(format_lower.begin(), format_lower.end(), format_lower.begin(),
                       [](unsigned char c) { return std::tolower(c); });

        if (format_lower == "vtk") {
            // VTK format: create subdirectory structure
            // vtk/solver__R_recon__N_size__CFL_value/
            std::ostringstream subdir_oss;
            subdir_oss << case_output_dir_ << "/vtk/"
                       << settings_.solver << "__R_" << settings_.reconstruction 
                       << "__N_" << settings_.N 
                       << "__CFL_" << utils::DoubleWithoutDot(settings_.cfl);
            std::string vtk_numerical_dir = subdir_oss.str();
            
            vtk_writer_ = WriterFactory::Create("vtk", vtk_numerical_dir, false);
            
            // Create analytical VTK writer if analytical is enabled
            if (settings_.analytical) {
                std::string vtk_analytical_dir = case_output_dir_ + "/vtk/analytical";
                vtk_analytical_writer_ = WriterFactory::Create("vtk", vtk_analytical_dir, true);
            }
        } else if (format_lower.substr(0, 3) == "png") {
            // PNG format: create single directory
            std::string png_dir = case_output_dir_ + "/png";
            png_writer_ = WriterFactory::Create(format, png_dir, false);
        } else if (format_lower.substr(0, 3) == "gif") {
            // GIF format: write directly to case directory
            // GIF file will be named: solver__R_recon__N_size__CFL_value.gif
            gif_writer_ = WriterFactory::Create(format, case_output_dir_, false);
        }
    }
}

void Simulation::ApplyInitialConditions(DataLayer& layer) {
    if (!layer_) {
        throw std::runtime_error(
            "DataLayer must be initialized before applying initial conditions");
    }
    const int n_ghost = layer.GetNGhostCells();
    const int n = layer.GetN();
    const int total = layer.GetTotalSize();

    const double dx = settings_.L_x / static_cast<double>(n);
    const double x_mid = settings_.x0;
    const double gamma = settings_.gamma;

    // Coordinates
    for (int i = 0; i < total; ++i) {
        const int cell_index = i - n_ghost;
        layer.xb(i) = cell_index * dx;
        layer.xc(i) = (cell_index + 0.5) * dx;
    }

    const int core_start = layer.GetCoreStart(0);
    const int core_end = layer.GetCoreEndExclusive(0);

    for (int i = core_start; i < core_end; ++i) {
        const double x = layer.xc(i);

        double rho, u, P;
        if (x < x_mid) {
            rho = initial_conditions_.rho_L;
            u = initial_conditions_.u_L;
            P = initial_conditions_.P_L;
        } else {
            rho = initial_conditions_.rho_R;
            u = initial_conditions_.u_R;
            P = initial_conditions_.P_R;
        }

        layer.rho(i) = rho;
        layer.u(i) = u;
        layer.P(i) = P;

        layer.m(i) = rho * dx;
        layer.p(i) = rho * u;
        layer.V(i) = 1.0 / rho;

        const double kinetic = 0.5 * rho * u * u;
        const double Eint = P / (gamma - 1.0);
        const double Etot = Eint + kinetic;
        const double eint = Eint / rho;

        layer.U(i) = eint;  // specific internal
        layer.e(i) = Etot;  // total density
    }
}

void Simulation::Initialize() {
    std::cout << "Initializing simulation..." << '\n';

    // Initialize data layer
    layer_ = std::make_unique<DataLayer>(settings_.N, settings_.padding, settings_.dim);

    // Apply initial conditions
    ApplyInitialConditions(*layer_);

    // Initialize solver based on settings
    solver_ = CreateSolver();
    solver_->SetCfl(settings_.cfl);

    // Set up boundary conditions based on settings
    auto left_bc =
        CreateBoundaryCondition(settings_.left_boundary, initial_conditions_.rho_L,
                                initial_conditions_.u_L, initial_conditions_.P_L);
    auto right_bc =
        CreateBoundaryCondition(settings_.right_boundary, initial_conditions_.rho_R,
                                initial_conditions_.u_R, initial_conditions_.P_R);
    solver_->AddBoundary(0, left_bc, right_bc);

    // Initialize analytical solver if enabled
    if (settings_.analytical) {
        analytical_settings_ = settings_;
        analytical_settings_.solver = "analytical";

        analytical_layer_ = std::make_unique<DataLayer>(
            analytical_settings_.N, analytical_settings_.padding,
            analytical_settings_.dim);

        ApplyInitialConditions(*analytical_layer_);

        analytical_solver_ = std::make_unique<AnalyticalSolver>(analytical_settings_);
    }

    // Create output writers
    CreateWriters();

    // Print configuration summary
    std::cout << '\n';
    std::cout << "Simulation initialized:" << '\n';
    std::cout << ">>> Case:                " << settings_.simulation_case << '\n';
    std::cout << ">>> Solver:              " << settings_.solver << '\n';
    std::cout << ">>> Riemann Solver:      " << settings_.riemann_solver << '\n';
    std::cout << ">>> Reconstruction:      " << settings_.reconstruction << '\n';
    std::cout << ">>> Boundary left:       " << settings_.left_boundary << '\n';
    std::cout << ">>> Boundary right:      " << settings_.right_boundary << '\n';
    std::cout << ">>> Analytical:          " << (settings_.analytical ? "enabled" : "disabled") << "\n";
    std::cout << ">>> Domain length (L_x): " << settings_.L_x << '\n';
    std::cout << ">>> Grid size (N):       " << settings_.N << '\n';
    std::cout << ">>> Ghost cells:         " << settings_.padding << '\n';
    std::cout << ">>> Dimension:           " << settings_.dim << '\n';
    std::cout << ">>> x0:                  " << settings_.x0 << "\n";
    std::cout << ">>> CFL:                 " << settings_.cfl << '\n';
    std::cout << ">>> Q_user:              " << settings_.Q_user << '\n';
    std::cout << ">>> Gamma:               " << settings_.gamma << '\n';
    std::cout << ">>> Max steps:           " << settings_.step_end << '\n';
    std::cout << ">>> Max time:            " << settings_.t_end << '\n';
    std::cout << ">>> Output every steps:  " << settings_.output_every_steps << '\n';
    std::cout << ">>> Output every time:   " << settings_.output_every_time << '\n';
    
    // Print output formats
    std::cout << ">>> Output formats:      [";
    for (size_t i = 0; i < settings_.output_formats.size(); ++i) {
        if (i > 0) std::cout << ", ";
        std::cout << settings_.output_formats[i];
    }
    std::cout << "]\n";
    
    std::cout << ">>> Output directory:    " << settings_.output_dir << '\n';
    std::cout << ">>> Initial conditions:" << '\n';
    std::cout << std::fixed << std::setprecision(4);
    std::cout << ">>>     Left:  rho = " << std::setw(7) << initial_conditions_.rho_L
              << ",    u = " << std::setw(7) << initial_conditions_.u_L
              << ",    P = " << std::setw(7) << initial_conditions_.P_L << '\n';
    std::cout << ">>>     Right: rho = " << std::setw(7) << initial_conditions_.rho_R
              << ",    u = " << std::setw(7) << initial_conditions_.u_R
              << ",    P = " << std::setw(7) << initial_conditions_.P_R << "\n\n";
}

auto Simulation::GetDataLayer() -> DataLayer& {
    if (!layer_) {
        throw std::runtime_error("DataLayer is not initialized.");
    }
    return *layer_;
}

auto Simulation::GetCurrentStep() const -> std::size_t { return step_cur_; }

auto Simulation::GetCurrentTime() const -> double { return t_cur_; }

auto Simulation::ShouldWrite() const -> bool {
    if (settings_.output_every_time == 0.0 && settings_.output_every_steps == 0) {
        return false;
    }

    const bool time_ok = t_cur_ >= settings_.t_end || settings_.log_every_time == 0.0 ||
                         std::floor((t_cur_ - dt_) / settings_.log_every_time) <
                         std::floor(t_cur_ / settings_.log_every_time);

    const bool step_ok = settings_.output_every_steps == 0 ||
                         step_cur_ % settings_.output_every_steps == 0;

    return (time_ok && step_ok) || t_cur_ >= settings_.t_end;
}

auto Simulation::ShouldLog() const -> bool {
    if (settings_.log_every_time == 0.0 && settings_.log_every_steps == 0) {
        return false;
    }

    const bool time_ok = settings_.log_every_time == 0.0 ||
                         std::floor((t_cur_ - dt_) / settings_.log_every_time) <
                         std::floor(t_cur_ / settings_.log_every_time);

    const bool step_ok =
        settings_.log_every_steps == 0 || step_cur_ % settings_.log_every_steps == 0;

    return (time_ok && step_ok) || t_cur_ >= settings_.t_end;
}

auto Simulation::ShouldRun() const -> bool {
    // If both stopping criteria are disabled (both set to 0), don't run
    if (settings_.t_end == 0.0 && settings_.step_end == 0) {
        return false;
    }

    // Check time-based stopping criterion (if enabled)
    const bool time_not_exceeded = 
        (settings_.t_end == 0.0) || (t_cur_ < settings_.t_end);
    
    // Check step-based stopping criterion (if enabled)
    const bool steps_not_exceeded = 
        (settings_.step_end == 0) || (step_cur_ < settings_.step_end);

    // Continue running if BOTH enabled criteria are satisfied
    return time_not_exceeded && steps_not_exceeded;
}

void Simulation::WriteInitialState() const {
    std::cout << "Writing the initial state..." << '\n';

    // Write VTK files
    if (vtk_writer_) {
        vtk_writer_->Write(*layer_, settings_, 0, 0.0);
    }
    if (vtk_analytical_writer_ && analytical_layer_) {
        vtk_analytical_writer_->Write(*analytical_layer_, settings_, 0, 0.0);
    }
    
    // Write PNG file (with analytical comparison if available)
    if (png_writer_) {
        const DataLayer* analytical_ptr = analytical_layer_ ? analytical_layer_.get() : nullptr;
        png_writer_->Write(*layer_, analytical_ptr, settings_, 0, 0.0);
    }

    // Write GIF frame (with analytical comparison if available)
    if (gif_writer_) {
        const DataLayer* analytical_ptr = analytical_layer_ ? analytical_layer_.get() : nullptr;
        gif_writer_->Write(*layer_, analytical_ptr, settings_, 0, 0.0);
    }
}

void Simulation::WriteStepState(double t_cur, std::size_t step_cur) const {
    if (!ShouldWrite()) {
        return;
    }
    
    // Write VTK files
    if (vtk_writer_) {
        vtk_writer_->Write(*layer_, settings_, step_cur, t_cur);
    }
    if (vtk_analytical_writer_ && analytical_layer_) {
        vtk_analytical_writer_->Write(*analytical_layer_, settings_, step_cur, t_cur);
    }
    
    // Write PNG file (with analytical comparison if available)
    if (png_writer_) {
        const DataLayer* analytical_ptr = analytical_layer_ ? analytical_layer_.get() : nullptr;
        png_writer_->Write(*layer_, analytical_ptr, settings_, step_cur, t_cur);
    }

    // Write GIF frame (with analytical comparison if available)
    if (gif_writer_) {
        const DataLayer* analytical_ptr = analytical_layer_ ? analytical_layer_.get() : nullptr;
        gif_writer_->Write(*layer_, analytical_ptr, settings_, step_cur, t_cur);
    }
}

void Simulation::PrintLog() const {
    if (ShouldLog()) {
        double progress = t_cur_ / settings_.t_end * 100.0;
        int percent = static_cast<int>(progress);
        std::cout << '\r';
        std::cout << ">>> [PROGRESS]: Step " << step_cur_ << ", " << percent
            << "% processed, time: " << t_cur_ << " of " << settings_.t_end;
        std::cout.flush();
    }
}

void Simulation::FinalizeWriters() {
    // Finalize any writers that require it
    if (gif_writer_ && gif_writer_->RequiresFinalization()) {
        std::string gif_path = gif_writer_->Finalize(settings_);
    }
    
    // VTK and PNG writers don't require finalization
    // but we check anyway for future extensibility
    if (vtk_writer_ && vtk_writer_->RequiresFinalization()) {
        vtk_writer_->Finalize(settings_);
    }
    if (vtk_analytical_writer_ && vtk_analytical_writer_->RequiresFinalization()) {
        vtk_analytical_writer_->Finalize(settings_);
    }
    if (png_writer_ && png_writer_->RequiresFinalization()) {
        png_writer_->Finalize(settings_);
    }
}

void Simulation::Run() {
    Initialize();
    WriteInitialState();

    t_cur_ = 0.0;
    step_cur_ = 0;

    std::cout << "\nStarting simulation..." << '\n';

    std::chrono::duration<double> runtime{0};
    std::chrono::duration<double> wall_time{0};

    auto start_wall = std::chrono::high_resolution_clock::now();
    while (ShouldRun()) {
        // Advance one time step and get the actual dt used
        auto start = std::chrono::high_resolution_clock::now();
        dt_ = solver_->Step(*layer_, t_cur_);
        auto end = std::chrono::high_resolution_clock::now();
        runtime += end - start;

        if (settings_.solver == "analytical") {
            t_cur_ += dt_;
        }
        ++step_cur_;

        if (settings_.analytical && analytical_solver_ && analytical_layer_) {
            // Force this step size
            analytical_solver_->SetDt(dt_);
            double analytical_t = t_cur_;
            analytical_solver_->Step(*analytical_layer_, analytical_t);
        }

        WriteStepState(t_cur_, step_cur_);

        PrintLog();
    }
    auto end_wall = std::chrono::high_resolution_clock::now();
    wall_time = end_wall - start_wall;

    std::cout << '\n';
    std::cout << "\nSimulation completed!" << '\n';
    std::cout << ">>> Final time:  " << t_cur_ << '\n';
    std::cout << ">>> Total steps: " << step_cur_ << '\n';
    std::cout << ">>> Wall time: " << wall_time.count() << "s\n";
    std::cout << ">>> Computation time: " << runtime.count() << "s\n";

    // Finalize writers (e.g., encode GIF from accumulated frames)
    FinalizeWriters();
}
