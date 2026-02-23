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
#include "solver/SolverFactory.hpp"

Simulation::Simulation(Settings settings, const InitialConditions& initial_conditions)
    : settings_(std::move(settings)), initial_conditions_(initial_conditions),
      boundary_manager_(std::make_shared<BoundaryManager>()) {}

auto Simulation::CreateSolver() -> std::unique_ptr<Solver> {
    return SolverFactory::Create(settings_, boundary_manager_);
}

auto Simulation::CreateBoundaryCondition(const std::string& boundary_type,
                                         const FarfieldConservative& far_field_U)
    -> std::shared_ptr<BoundaryCondition> {
    return BoundaryFactory::Create(boundary_type, far_field_U);
}

void Simulation::CreateWriters() {
    case_output_dir_ = settings_.output_dir;

    const int nx = settings_.GetNx();
    const int ny = settings_.GetNy();
    const int nz = settings_.GetNz();

    for (const auto& format : settings_.output_formats) {
        std::string f = format;
        std::transform(f.begin(), f.end(), f.begin(),
                       [](unsigned char c) {
                           return static_cast<char>(std::tolower(c));
                       });

        if (f == "vtk") {
            std::ostringstream subdir;
            subdir << case_output_dir_ << "/vtk/"
                << "dim_" << settings_.dim
                << "__" << settings_.solver
                << "__R_" << settings_.reconstruction
                << "__grid_" << nx << "x" << ny << "x" << nz
                << "__CFL_" << utils::DoubleWithoutDot(settings_.cfl);

            vtk_writer_ = WriterFactory::Create("vtk", subdir.str(), false);

            if (settings_.analytical) {
                std::string vtk_analytical_dir = case_output_dir_ + "/vtk/analytical";
                vtk_analytical_writer_ = WriterFactory::Create("vtk", vtk_analytical_dir, true);
            }
        }
    }
}

void Simulation::ApplyInitialConditions(DataLayer& layer) {
    const int dim = layer.GetDim();
    const int pad = layer.GetPadding();

    const int sx = layer.GetSx();
    const int sy = layer.GetSy();
    const int sz = layer.GetSz();

    const int nx = layer.GetNx(); // core
    const int ny = layer.GetNy(); // уже 1 если dim<2
    const int nz = layer.GetNz(); // уже 1 если dim<3

    const double gamma = settings_.gamma;

    // --- Uniform grid spacing for now ---
    const double dx = settings_.L_x / static_cast<double>(nx);
    const double dy = (dim >= 2) ? (settings_.L_y / static_cast<double>(ny)) : 1.0;
    const double dz = (dim >= 3) ? (settings_.L_z / static_cast<double>(nz)) : 1.0;

    // --- Build boundary coordinates (size s? + 1) and centers (size s?) ---
    // X
    {
        auto& xb = layer.Xb();
        auto& xc = layer.Xc();
        for (int i = 0; i <= sx; ++i) {
            const int ii = i - pad;
            xb(static_cast<std::size_t>(i)) = static_cast<double>(ii) * dx;
        }
        for (int i = 0; i < sx; ++i) {
            const int ii = i - pad;
            xc(static_cast<std::size_t>(i)) = (static_cast<double>(ii) + 0.5) * dx;
        }
    }

    // Y
    {
        auto& yb = layer.Yb();
        auto& yc = layer.Yc();
        for (int j = 0; j <= sy; ++j) {
            const int jj = j - pad;
            yb(static_cast<std::size_t>(j)) = static_cast<double>(jj) * dy;
        }
        for (int j = 0; j < sy; ++j) {
            const int jj = j - pad;
            yc(static_cast<std::size_t>(j)) = (static_cast<double>(jj) + 0.5) * dy;
        }
    }

    // Z
    {
        auto& zb = layer.Zb();
        auto& zc = layer.Zc();
        for (int k = 0; k <= sz; ++k) {
            const int kk = k - pad;
            zb(static_cast<std::size_t>(k)) = static_cast<double>(kk) * dz;
        }
        for (int k = 0; k < sz; ++k) {
            const int kk = k - pad;
            zc(static_cast<std::size_t>(k)) = (static_cast<double>(kk) + 0.5) * dz;
        }
    }

    // Fill dx/dy/dz arrays + inverses from boundary coords (no allocations)
    layer.UpdateMetricsFromCoordinates();

    // --- Core region indices ---
    const int i0 = layer.GetCoreStartX();
    const int i1 = layer.GetCoreEndExclusiveX();
    const int j0 = layer.GetCoreStartY();
    const int j1 = layer.GetCoreEndExclusiveY();
    const int k0 = layer.GetCoreStartZ();
    const int k1 = layer.GetCoreEndExclusiveZ();

    const std::string& ic_type = initial_conditions_.ic_type;

    auto& U = layer.U();

    auto write_conservative = [&](int i, int j, int k,
                                  double rho, double u, double v, double w, double P) {
        // conservative
        const double rhoU = rho * u;
        const double rhoV = rho * v;
        const double rhoW = rho * w;
        const double kinetic = 0.5 * rho * (u * u + v * v + w * w);
        const double E = P / (gamma - 1.0) + kinetic;

        U(DataLayer::k_rho, i, j, k) = rho;
        U(DataLayer::k_rhoU, i, j, k) = rhoU;
        U(DataLayer::k_rhoV, i, j, k) = rhoV;
        U(DataLayer::k_rhoW, i, j, k) = rhoW;
        U(DataLayer::k_E, i, j, k) = E;
    };

    // --- Apply IC on core cells ---
    for (int k = k0; k < k1; ++k) {
        const double z = layer.Zc()(static_cast<std::size_t>(k));
        for (int j = j0; j < j1; ++j) {
            const double y = layer.Yc()(static_cast<std::size_t>(j));
            for (int i = i0; i < i1; ++i) {
                const double x = layer.Xc()(static_cast<std::size_t>(i));

                double rho = 0.0, u = 0.0, v = 0.0, w = 0.0, P = 0.0;

                if (ic_type == "x_riemann") {
                    if (x < settings_.x0) {
                        rho = initial_conditions_.rho_L;
                        u = initial_conditions_.u_L;
                        v = initial_conditions_.v_L;
                        w = 0.0;
                        P = initial_conditions_.P_L;
                    }
                    else {
                        rho = initial_conditions_.rho_R;
                        u = initial_conditions_.u_R;
                        v = initial_conditions_.v_R;
                        w = 0.0;
                        P = initial_conditions_.P_R;
                    }
                }
                else if (ic_type == "y_riemann") {
                    if (dim < 2) throw std::runtime_error("ic_type=y_riemann requires dim>=2");
                    if (y < settings_.y0) {
                        rho = initial_conditions_.rho_L;
                        u = initial_conditions_.u_L;
                        v = initial_conditions_.v_L;
                        w = 0.0;
                        P = initial_conditions_.P_L;
                    }
                    else {
                        rho = initial_conditions_.rho_R;
                        u = initial_conditions_.u_R;
                        v = initial_conditions_.v_R;
                        w = 0.0;
                        P = initial_conditions_.P_R;
                    }
                }
                else if (ic_type == "quadrant") {
                    if (dim != 2) throw std::runtime_error("ic_type=quadrant currently supported only for dim==2");
                    const bool right = x >= settings_.x0;
                    const bool top = y >= settings_.y0;

                    if (right && top) {
                        rho = initial_conditions_.rho_Q1;
                        u = initial_conditions_.u_Q1;
                        v = initial_conditions_.v_Q1;
                        w = 0.0;
                        P = initial_conditions_.P_Q1;
                    }
                    else if (!right && top) {
                        rho = initial_conditions_.rho_Q2;
                        u = initial_conditions_.u_Q2;
                        v = initial_conditions_.v_Q2;
                        w = 0.0;
                        P = initial_conditions_.P_Q2;
                    }
                    else if (!right && !top) {
                        rho = initial_conditions_.rho_Q3;
                        u = initial_conditions_.u_Q3;
                        v = initial_conditions_.v_Q3;
                        w = 0.0;
                        P = initial_conditions_.P_Q3;
                    }
                    else {
                        rho = initial_conditions_.rho_Q4;
                        u = initial_conditions_.u_Q4;
                        v = initial_conditions_.v_Q4;
                        w = 0.0;
                        P = initial_conditions_.P_Q4;
                    }
                }
                else if (ic_type == "octant" || ic_type == "octants") {
                    if (dim != 3) throw std::runtime_error("ic_type=octant requires dim==3");
                    const bool xp = x >= settings_.x0;
                    const bool yp = y >= settings_.y0;
                    const bool zp = z >= settings_.z0;

                    if (!zp) {
                        if (xp && yp) {
                            rho = initial_conditions_.rho_Q1;
                            u = initial_conditions_.u_Q1;
                            v = initial_conditions_.v_Q1;
                            w = initial_conditions_.w_Q1;
                            P = initial_conditions_.P_Q1;
                        }
                        else if (!xp && yp) {
                            rho = initial_conditions_.rho_Q2;
                            u = initial_conditions_.u_Q2;
                            v = initial_conditions_.v_Q2;
                            w = initial_conditions_.w_Q2;
                            P = initial_conditions_.P_Q2;
                        }
                        else if (!xp && !yp) {
                            rho = initial_conditions_.rho_Q3;
                            u = initial_conditions_.u_Q3;
                            v = initial_conditions_.v_Q3;
                            w = initial_conditions_.w_Q3;
                            P = initial_conditions_.P_Q3;
                        }
                        else {
                            rho = initial_conditions_.rho_Q4;
                            u = initial_conditions_.u_Q4;
                            v = initial_conditions_.v_Q4;
                            w = initial_conditions_.w_Q4;
                            P = initial_conditions_.P_Q4;
                        }
                    }
                    else {
                        if (xp && yp) {
                            rho = initial_conditions_.rho_Q5;
                            u = initial_conditions_.u_Q5;
                            v = initial_conditions_.v_Q5;
                            w = initial_conditions_.w_Q5;
                            P = initial_conditions_.P_Q5;
                        }
                        else if (!xp && yp) {
                            rho = initial_conditions_.rho_Q6;
                            u = initial_conditions_.u_Q6;
                            v = initial_conditions_.v_Q6;
                            w = initial_conditions_.w_Q6;
                            P = initial_conditions_.P_Q6;
                        }
                        else if (!xp && !yp) {
                            rho = initial_conditions_.rho_Q7;
                            u = initial_conditions_.u_Q7;
                            v = initial_conditions_.v_Q7;
                            w = initial_conditions_.w_Q7;
                            P = initial_conditions_.P_Q7;
                        }
                        else {
                            rho = initial_conditions_.rho_Q8;
                            u = initial_conditions_.u_Q8;
                            v = initial_conditions_.v_Q8;
                            w = initial_conditions_.w_Q8;
                            P = initial_conditions_.P_Q8;
                        }
                    }
                }
                else {
                    throw std::runtime_error("Unknown ic_type: " + ic_type);
                }

                // For 1D/2D safety: enforce unused velocity components = 0
                if (dim < 2) v = 0.0;
                if (dim < 3) w = 0.0;

                write_conservative(i, j, k, rho, u, v, w, P);
            }
        }
    }
}

void Simulation::Initialize() {
    std::cout << "Initializing simulation..." << '\n';

    layer_ = std::make_unique<DataLayer>(settings_.Nx, settings_.Ny, settings_.Nz,
                                         settings_.padding, settings_.dim);

    // single-domain: all 6 faces are global boundaries
    layer_->SetAllGlobalBoundaries(true);

    ApplyInitialConditions(*layer_);

    FarfieldConservative far_field_U{};

    // X boundaries always exist
    auto left_bc = CreateBoundaryCondition(settings_.left_boundary, far_field_U);
    auto right_bc = CreateBoundaryCondition(settings_.right_boundary, far_field_U);
    boundary_manager_->Set(Axis::X, left_bc, right_bc);

    // Y boundaries only if dim >= 2
    if (settings_.dim >= 2) {
        auto bottom_bc = CreateBoundaryCondition(settings_.bottom_boundary, far_field_U);
        auto top_bc = CreateBoundaryCondition(settings_.top_boundary, far_field_U);
        boundary_manager_->Set(Axis::Y, bottom_bc, top_bc);
    }

    // Z boundaries only if dim >= 3
    if (settings_.dim >= 3) {
        auto back_bc = CreateBoundaryCondition(settings_.back_boundary, far_field_U);
        auto front_bc = CreateBoundaryCondition(settings_.front_boundary, far_field_U);
        boundary_manager_->Set(Axis::Z, back_bc, front_bc);
    }

    solver_ = CreateSolver();
    solver_->SetCfl(settings_.cfl);

    if (settings_.analytical) {
        analytical_settings_ = settings_;
        analytical_settings_.solver = "analytical";

        analytical_layer_ = std::make_unique<DataLayer>(analytical_settings_.Nx, analytical_settings_.Ny,
                                                        analytical_settings_.Nz, analytical_settings_.padding,
                                                        analytical_settings_.dim);

        analytical_layer_->SetAllGlobalBoundaries(true);

        ApplyInitialConditions(*analytical_layer_);
    }

    CreateWriters();
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

        if (settings_.solver == "analytical") {
            t_cur_ += dt_;
        }
        ++step_cur_;

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

auto Simulation::GetDataLayer() -> DataLayer& {
    if (!layer_) {
        throw std::runtime_error("DataLayer is not initialized.");
    }
    return *layer_;
}

auto Simulation::GetCurrentStep() const -> std::size_t {
    return step_cur_;
}

auto Simulation::GetCurrentTime() const -> double {
    return t_cur_;
}

auto Simulation::ShouldWrite() const -> bool {
    if (settings_.output_every_time == 0.0 && settings_.output_every_steps == 0) {
        return false;
    }

    const bool time_ok = t_cur_ >= settings_.t_end || settings_.output_every_time == 0.0 ||
        std::floor((t_cur_ - dt_) / settings_.output_every_time) <
        std::floor(t_cur_ / settings_.output_every_time);

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
    if (settings_.t_end == 0.0 && settings_.step_end == 0) {
        return false;
    }

    const bool time_not_exceeded =
        settings_.t_end == 0.0 || t_cur_ < settings_.t_end;

    const bool steps_not_exceeded =
        settings_.step_end == 0 || step_cur_ < settings_.step_end;

    return time_not_exceeded && steps_not_exceeded;
}

void Simulation::WriteInitialState() const {
    std::cout << "Writing the initial state..." << '\n';
    if (vtk_writer_) {
        vtk_writer_->Write(*layer_, settings_, 0, 0.0);
    }
    if (vtk_analytical_writer_ && analytical_layer_) {
        vtk_analytical_writer_->Write(*analytical_layer_, settings_, 0, 0.0);
    }
}

void Simulation::WriteStepState(double t_cur, std::size_t step_cur) const {
    if (!ShouldWrite()) {
        return;
    }
    if (vtk_writer_) {
        vtk_writer_->Write(*layer_, settings_, step_cur, t_cur);
    }
    if (vtk_analytical_writer_ && analytical_layer_) {
        vtk_analytical_writer_->Write(*analytical_layer_, settings_, step_cur, t_cur);
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
    if (gif_writer_ && gif_writer_->RequiresFinalization()) {
        gif_writer_->Finalize(settings_);
    }
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
