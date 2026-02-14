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

auto Simulation::CreateBoundaryCondition2D(const std::string& boundary_type,
                                           double rho_inf, double u_inf,
                                           double v_inf, double p_inf)
    -> std::shared_ptr<BoundaryCondition> {
    return BoundaryFactory::Create2D(boundary_type, rho_inf, u_inf, v_inf, p_inf);
}

void Simulation::CreateWriters() {
    case_output_dir_ = settings_.output_dir;
    
    for (const auto& format : settings_.output_formats) {
        std::string format_lower = format;
        std::transform(format_lower.begin(), format_lower.end(), format_lower.begin(),
                       [](unsigned char c) { return std::tolower(c); });

        if (format_lower == "vtk") {
            std::ostringstream subdir_oss;
            subdir_oss << case_output_dir_ << "/vtk/"
                       << settings_.solver << "__R_" << settings_.reconstruction 
                       << "__N_" << settings_.GetNx();
            if (settings_.dim >= 2) {
                subdir_oss << "x" << settings_.GetNy();
            }
            subdir_oss << "__CFL_" << utils::DoubleWithoutDot(settings_.cfl);
            std::string vtk_numerical_dir = subdir_oss.str();
            
            vtk_writer_ = WriterFactory::Create("vtk", vtk_numerical_dir, false);
            
            if (settings_.analytical && settings_.dim == 1) {
                std::string vtk_analytical_dir = case_output_dir_ + "/vtk/analytical";
                vtk_analytical_writer_ = WriterFactory::Create("vtk", vtk_analytical_dir, true);
            }
        } else if (format_lower.substr(0, 3) == "png") {
            std::string png_dir = case_output_dir_ + "/png";
            png_writer_ = WriterFactory::Create(format, png_dir, false);
        } else if (format_lower.substr(0, 3) == "gif") {
            gif_writer_ = WriterFactory::Create(format, case_output_dir_, false);
        }
    }
}

// ============================================================================
// Initial Conditions (1D and 2D)
// ============================================================================

void Simulation::ApplyInitialConditions(DataLayer& layer) {
    if (!layer_) {
        throw std::runtime_error("DataLayer must be initialized before applying initial conditions");
    }

    if (settings_.dim >= 2) {
        ApplyInitialConditions2D(layer);
        return;
    }

    // --- Original 1D logic ---
    const int n_ghost = layer.GetNGhostCells();
    const int n = layer.GetN();
    const int total = layer.GetTotalSize();

    const double dx = settings_.L_x / static_cast<double>(n);
    const double x_mid = settings_.x0;
    const double gamma = settings_.gamma;

    for (int i = 0; i < total; ++i) {
        const int cell_index = i - n_ghost;
        layer.xb(i) = cell_index * dx;
        layer.xc(i) = (cell_index + 0.5) * dx;
    }

    const int core_start = layer.GetCoreStart(0);
    const int core_end = layer.GetCoreEndExclusive(0);

    for (int i = core_start; i < core_end; ++i) {
        const double x = layer.xc(i);
        double rho, u_val, P;
        if (x < x_mid) {
            rho = initial_conditions_.rho_L;
            u_val = initial_conditions_.u_L;
            P = initial_conditions_.P_L;
        } else {
            rho = initial_conditions_.rho_R;
            u_val = initial_conditions_.u_R;
            P = initial_conditions_.P_R;
        }

        layer.rho(i) = rho;
        layer.u(i) = u_val;
        layer.P(i) = P;
        layer.m(i) = rho * dx;
        layer.p(i) = rho * u_val;
        layer.V(i) = 1.0 / rho;

        const double kinetic = 0.5 * rho * u_val * u_val;
        const double Eint = P / (gamma - 1.0);
        const double Etot = Eint + kinetic;
        const double eint = Eint / rho;
        layer.U(i) = eint;
        layer.e(i) = Etot;
    }
}

void Simulation::ApplyInitialConditions2D(DataLayer& layer) {
    const int pad = layer.GetPadding();
    const int nx = layer.GetNx();
    const int ny = layer.GetNy();
    const int tx = layer.GetTotalSize(0);
    const int ty = layer.GetTotalSize(1);
    const double dx = settings_.L_x / static_cast<double>(nx);
    const double dy = settings_.L_y / static_cast<double>(ny);
    const double gamma = settings_.gamma;
    const double x_mid = settings_.x0;
    const double y_mid = settings_.y0;

    // Set up coordinate arrays
    for (int i = 0; i < tx; ++i) {
        layer.xb(i) = (i - pad) * dx;
        layer.xc(i) = (i - pad + 0.5) * dx;
    }
    for (int j = 0; j < ty; ++j) {
        layer.yb(j) = (j - pad) * dy;
        layer.yc(j) = (j - pad + 0.5) * dy;
    }

    const int cs_x = layer.GetCoreStart(0);
    const int ce_x = layer.GetCoreEndExclusive(0);
    const int cs_y = layer.GetCoreStart(1);
    const int ce_y = layer.GetCoreEndExclusive(1);

    for (int i = cs_x; i < ce_x; ++i) {
        const double x = layer.xc(i);
        for (int j = cs_y; j < ce_y; ++j) {
            const double y = layer.yc(j);

            double rho, u_val, v_val, P;

            if (initial_conditions_.IsQuadrant()) {
                // Four-quadrant IC
                if (x >= x_mid && y >= y_mid) {
                    // Q1: top-right
                    rho = initial_conditions_.rho_Q1;
                    u_val = initial_conditions_.u_Q1;
                    v_val = initial_conditions_.v_Q1;
                    P = initial_conditions_.P_Q1;
                } else if (x < x_mid && y >= y_mid) {
                    // Q2: top-left
                    rho = initial_conditions_.rho_Q2;
                    u_val = initial_conditions_.u_Q2;
                    v_val = initial_conditions_.v_Q2;
                    P = initial_conditions_.P_Q2;
                } else if (x < x_mid && y < y_mid) {
                    // Q3: bottom-left
                    rho = initial_conditions_.rho_Q3;
                    u_val = initial_conditions_.u_Q3;
                    v_val = initial_conditions_.v_Q3;
                    P = initial_conditions_.P_Q3;
                } else {
                    // Q4: bottom-right
                    rho = initial_conditions_.rho_Q4;
                    u_val = initial_conditions_.u_Q4;
                    v_val = initial_conditions_.v_Q4;
                    P = initial_conditions_.P_Q4;
                }
            } else if (initial_conditions_.ic_type == "y_riemann") {
                // Discontinuity along y
                if (y < y_mid) {
                    rho = initial_conditions_.rho_L;
                    u_val = initial_conditions_.u_L;
                    v_val = initial_conditions_.v_L;
                    P = initial_conditions_.P_L;
                } else {
                    rho = initial_conditions_.rho_R;
                    u_val = initial_conditions_.u_R;
                    v_val = initial_conditions_.v_R;
                    P = initial_conditions_.P_R;
                }
            } else {
                // Default: x_riemann (discontinuity along x)
                if (x < x_mid) {
                    rho = initial_conditions_.rho_L;
                    u_val = initial_conditions_.u_L;
                    v_val = initial_conditions_.v_L;
                    P = initial_conditions_.P_L;
                } else {
                    rho = initial_conditions_.rho_R;
                    u_val = initial_conditions_.u_R;
                    v_val = initial_conditions_.v_R;
                    P = initial_conditions_.P_R;
                }
            }

            layer.rho(i, j) = rho;
            layer.u(i, j) = u_val;
            layer.v(i, j) = v_val;
            layer.P(i, j) = P;
            layer.p(i, j) = rho * u_val;
            layer.q(i, j) = rho * v_val;
            layer.V(i, j) = 1.0 / rho;
            layer.m(i, j) = rho * dx * dy;

            const double kinetic = 0.5 * rho * (u_val * u_val + v_val * v_val);
            const double Eint = P / (gamma - 1.0);
            const double Etot = Eint + kinetic;
            const double eint = Eint / rho;
            layer.U(i, j) = eint;
            layer.e(i, j) = Etot;
        }
    }
}

// ============================================================================
// Initialize
// ============================================================================

void Simulation::Initialize() {
    std::cout << "Initializing simulation..." << '\n';

    // Initialize data layer
    if (settings_.dim >= 2) {
        layer_ = std::make_unique<DataLayer>(settings_.GetNx(), settings_.GetNy(),
                                             settings_.padding, settings_.dim);
    } else {
        layer_ = std::make_unique<DataLayer>(settings_.N, settings_.padding, settings_.dim);
    }

    ApplyInitialConditions(*layer_);

    solver_ = CreateSolver();
    solver_->SetCfl(settings_.cfl);

    // Set up boundary conditions
    if (settings_.dim >= 2) {
        // X-axis boundaries
        auto left_bc = CreateBoundaryCondition2D(settings_.left_boundary,
            initial_conditions_.rho_L, initial_conditions_.u_L,
            initial_conditions_.v_L, initial_conditions_.P_L);
        auto right_bc = CreateBoundaryCondition2D(settings_.right_boundary,
            initial_conditions_.rho_R, initial_conditions_.u_R,
            initial_conditions_.v_R, initial_conditions_.P_R);
        solver_->AddBoundary(0, left_bc, right_bc);

        // Y-axis boundaries
        auto bottom_bc = CreateBoundaryCondition2D(settings_.bottom_boundary,
            initial_conditions_.rho_L, initial_conditions_.u_L,
            initial_conditions_.v_L, initial_conditions_.P_L);
        auto top_bc = CreateBoundaryCondition2D(settings_.top_boundary,
            initial_conditions_.rho_R, initial_conditions_.u_R,
            initial_conditions_.v_R, initial_conditions_.P_R);
        solver_->AddBoundary(1, bottom_bc, top_bc);
    } else {
        auto left_bc = CreateBoundaryCondition(settings_.left_boundary,
            initial_conditions_.rho_L, initial_conditions_.u_L, initial_conditions_.P_L);
        auto right_bc = CreateBoundaryCondition(settings_.right_boundary,
            initial_conditions_.rho_R, initial_conditions_.u_R, initial_conditions_.P_R);
        solver_->AddBoundary(0, left_bc, right_bc);
    }

    // Analytical solver (1D only)
    if (settings_.analytical && settings_.dim == 1) {
        analytical_settings_ = settings_;
        analytical_settings_.solver = "analytical";
        analytical_layer_ = std::make_unique<DataLayer>(
            analytical_settings_.N, analytical_settings_.padding, analytical_settings_.dim);
        ApplyInitialConditions(*analytical_layer_);
        analytical_solver_ = std::make_unique<AnalyticalSolver>(analytical_settings_);
    }

    CreateWriters();

    // Print configuration
    std::cout << '\n';
    std::cout << "Simulation initialized:" << '\n';
    std::cout << ">>> Case:                " << settings_.simulation_case << '\n';
    std::cout << ">>> Solver:              " << settings_.solver << '\n';
    std::cout << ">>> Riemann Solver:      " << settings_.riemann_solver << '\n';
    std::cout << ">>> Reconstruction:      " << settings_.reconstruction << '\n';
    std::cout << ">>> Dimension:           " << settings_.dim << '\n';
    if (settings_.dim >= 2) {
        std::cout << ">>> Grid size (Nx x Ny): " << settings_.GetNx() << " x " << settings_.GetNy() << '\n';
        std::cout << ">>> Domain (Lx x Ly):    " << settings_.L_x << " x " << settings_.L_y << '\n';
        std::cout << ">>> Boundary left:       " << settings_.left_boundary << '\n';
        std::cout << ">>> Boundary right:      " << settings_.right_boundary << '\n';
        std::cout << ">>> Boundary bottom:     " << settings_.bottom_boundary << '\n';
        std::cout << ">>> Boundary top:        " << settings_.top_boundary << '\n';
        std::cout << ">>> x0, y0:              " << settings_.x0 << ", " << settings_.y0 << '\n';
    } else {
        std::cout << ">>> Grid size (N):       " << settings_.N << '\n';
        std::cout << ">>> Domain length (L_x): " << settings_.L_x << '\n';
        std::cout << ">>> Boundary left:       " << settings_.left_boundary << '\n';
        std::cout << ">>> Boundary right:      " << settings_.right_boundary << '\n';
        std::cout << ">>> x0:                  " << settings_.x0 << '\n';
    }
    std::cout << ">>> Ghost cells:         " << settings_.padding << '\n';
    std::cout << ">>> CFL:                 " << settings_.cfl << '\n';
    std::cout << ">>> Gamma:               " << settings_.gamma << '\n';
    std::cout << ">>> Max time:            " << settings_.t_end << '\n';
    std::cout << ">>> Output formats:      [";
    for (size_t i = 0; i < settings_.output_formats.size(); ++i) {
        if (i > 0) std::cout << ", ";
        std::cout << settings_.output_formats[i];
    }
    std::cout << "]\n";
    std::cout << ">>> Output directory:    " << settings_.output_dir << "\n\n";
}

// ============================================================================
// Rest of Simulation (unchanged except Run loop is identical)
// ============================================================================

auto Simulation::GetDataLayer() -> DataLayer& {
    if (!layer_) throw std::runtime_error("DataLayer is not initialized.");
    return *layer_;
}

auto Simulation::GetCurrentStep() const -> std::size_t { return step_cur_; }
auto Simulation::GetCurrentTime() const -> double { return t_cur_; }

auto Simulation::ShouldWrite() const -> bool {
    if (settings_.output_every_time == 0.0 && settings_.output_every_steps == 0) return false;
    const bool time_ok = t_cur_ >= settings_.t_end || settings_.log_every_time == 0.0 ||
                         std::floor((t_cur_ - dt_) / settings_.log_every_time) <
                         std::floor(t_cur_ / settings_.log_every_time);
    const bool step_ok = settings_.output_every_steps == 0 ||
                         step_cur_ % settings_.output_every_steps == 0;
    return (time_ok && step_ok) || t_cur_ >= settings_.t_end;
}

auto Simulation::ShouldLog() const -> bool {
    if (settings_.log_every_time == 0.0 && settings_.log_every_steps == 0) return false;
    const bool time_ok = settings_.log_every_time == 0.0 ||
                         std::floor((t_cur_ - dt_) / settings_.log_every_time) <
                         std::floor(t_cur_ / settings_.log_every_time);
    const bool step_ok = settings_.log_every_steps == 0 ||
                         step_cur_ % settings_.log_every_steps == 0;
    return (time_ok && step_ok) || t_cur_ >= settings_.t_end;
}

auto Simulation::ShouldRun() const -> bool {
    if (settings_.t_end == 0.0 && settings_.step_end == 0) return false;
    const bool time_not_exceeded = (settings_.t_end == 0.0) || (t_cur_ < settings_.t_end);
    const bool steps_not_exceeded = (settings_.step_end == 0) || (step_cur_ < settings_.step_end);
    return time_not_exceeded && steps_not_exceeded;
}

void Simulation::WriteInitialState() const {
    std::cout << "Writing the initial state..." << '\n';
    if (vtk_writer_) vtk_writer_->Write(*layer_, settings_, 0, 0.0);
    if (vtk_analytical_writer_ && analytical_layer_)
        vtk_analytical_writer_->Write(*analytical_layer_, settings_, 0, 0.0);
    if (png_writer_) {
        const DataLayer* ap = analytical_layer_ ? analytical_layer_.get() : nullptr;
        png_writer_->Write(*layer_, ap, settings_, 0, 0.0);
    }
    if (gif_writer_) {
        const DataLayer* ap = analytical_layer_ ? analytical_layer_.get() : nullptr;
        gif_writer_->Write(*layer_, ap, settings_, 0, 0.0);
    }
}

void Simulation::WriteStepState(double t_cur, std::size_t step_cur) const {
    if (!ShouldWrite()) return;
    if (vtk_writer_) vtk_writer_->Write(*layer_, settings_, step_cur, t_cur);
    if (vtk_analytical_writer_ && analytical_layer_)
        vtk_analytical_writer_->Write(*analytical_layer_, settings_, step_cur, t_cur);
    if (png_writer_) {
        const DataLayer* ap = analytical_layer_ ? analytical_layer_.get() : nullptr;
        png_writer_->Write(*layer_, ap, settings_, step_cur, t_cur);
    }
    if (gif_writer_) {
        const DataLayer* ap = analytical_layer_ ? analytical_layer_.get() : nullptr;
        gif_writer_->Write(*layer_, ap, settings_, step_cur, t_cur);
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
    if (gif_writer_ && gif_writer_->RequiresFinalization())
        gif_writer_->Finalize(settings_);
    if (vtk_writer_ && vtk_writer_->RequiresFinalization())
        vtk_writer_->Finalize(settings_);
    if (vtk_analytical_writer_ && vtk_analytical_writer_->RequiresFinalization())
        vtk_analytical_writer_->Finalize(settings_);
    if (png_writer_ && png_writer_->RequiresFinalization())
        png_writer_->Finalize(settings_);
}

void Simulation::Run() {
    Initialize();
    WriteInitialState();

    t_cur_ = 0.0;
    step_cur_ = 0;

    std::cout << "\nStarting simulation..." << '\n';

    std::chrono::duration<double> runtime{0};
    auto start_wall = std::chrono::high_resolution_clock::now();

    while (ShouldRun()) {
        auto start = std::chrono::high_resolution_clock::now();
        dt_ = solver_->Step(*layer_, t_cur_);
        auto end = std::chrono::high_resolution_clock::now();
        runtime += end - start;

        if (settings_.solver == "analytical") t_cur_ += dt_;
        ++step_cur_;

        if (settings_.analytical && analytical_solver_ && analytical_layer_) {
            analytical_solver_->SetDt(dt_);
            double analytical_t = t_cur_;
            analytical_solver_->Step(*analytical_layer_, analytical_t);
        }

        WriteStepState(t_cur_, step_cur_);
        PrintLog();
    }

    auto end_wall = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> wall_time = end_wall - start_wall;

    std::cout << '\n';
    std::cout << "\nSimulation completed!" << '\n';
    std::cout << ">>> Final time:  " << t_cur_ << '\n';
    std::cout << ">>> Total steps: " << step_cur_ << '\n';
    std::cout << ">>> Wall time: " << wall_time.count() << "s\n";
    std::cout << ">>> Computation time: " << runtime.count() << "s\n";

    FinalizeWriters();
}
