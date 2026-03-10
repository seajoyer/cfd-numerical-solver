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
#include "data/geometry/GeometryFactory.hpp"

Simulation::Simulation(Settings settings, const InitialConditions& initial_conditions)
    : settings_(std::move(settings)), initial_conditions_(initial_conditions),
      boundary_manager_(std::make_shared<BoundaryManager>(nullptr)) {}

auto Simulation::CreateSolver() -> std::unique_ptr<Solver> {
    return SolverFactory::Create(settings_, *mesh_, boundary_manager_, mpi_context_.get());
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

    const int rank = mpi_context_ ? mpi_context_->Rank() : 0;
    const int size = mpi_context_ ? mpi_context_->Size() : 1;

    for (const auto& format : settings_.output_formats) {
        std::string f = format;
        std::transform(f.begin(), f.end(), f.begin(),
                       [](unsigned char c) {
                           return static_cast<char>(std::tolower(c));
                       });

        if (f == "vtk") {
            std::ostringstream subdir;
            subdir << case_output_dir_ << "/vtk/"
                << settings_.solver
                << "__R_" << settings_.reconstruction
                << "__N_" << nx << "x" << ny << "x" << nz
                << "__CFL_" << utils::DoubleWithoutDot(settings_.cfl);

            vtk_writer_ = WriterFactory::Create("vtk", subdir.str(), false, rank, size);

            if (settings_.analytical) {
                std::string vtk_analytical_dir = case_output_dir_ + "/vtk/analytical";
                vtk_analytical_writer_ = WriterFactory::Create("vtk", vtk_analytical_dir, true, rank, size);
            }
        }
    }
}

void Simulation::ApplyInitialConditions(DataLayer& layer, Mesh& mesh) {
    const int dim = mesh.GetDim();
    const int pad = mesh.GetPadding();

    const int sx = mesh.GetSx();
    const int sy = mesh.GetSy();
    const int sz = mesh.GetSz();

    const double gamma = settings_.gamma;

    const int global_nx = mesh.GetGlobalNx();
    const int global_ny = mesh.GetGlobalNy();
    const int global_nz = mesh.GetGlobalNz();

    const double dx = settings_.L_x / static_cast<double>(global_nx);
    const double dy = dim >= 2 ? settings_.L_y / static_cast<double>(global_ny) : 1.0;
    const double dz = dim >= 3 ? settings_.L_z / static_cast<double>(global_nz) : 1.0;

    {
        auto& xb = mesh.Xb();
        auto& xc = mesh.Xc();

        const int offset_x = mesh.GetOffsetX();

        for (int i = 0; i <= sx; ++i) {
            const int ig = offset_x + (i - pad);
            xb(static_cast<std::size_t>(i)) = static_cast<double>(ig) * dx;
        }
        for (int i = 0; i < sx; ++i) {
            const int ig = offset_x + (i - pad);
            xc(static_cast<std::size_t>(i)) = (static_cast<double>(ig) + 0.5) * dx;
        }
    }

    {
        auto& yb = mesh.Yb();
        auto& yc = mesh.Yc();

        const int offset_y = mesh.GetOffsetY();

        for (int j = 0; j <= sy; ++j) {
            const int jg = offset_y + (j - pad);
            yb(static_cast<std::size_t>(j)) = static_cast<double>(jg) * dy;
        }
        for (int j = 0; j < sy; ++j) {
            const int jg = offset_y + (j - pad);
            yc(static_cast<std::size_t>(j)) = (static_cast<double>(jg) + 0.5) * dy;
        }
    }

    {
        auto& zb = mesh.Zb();
        auto& zc = mesh.Zc();

        const int offset_z = mesh.GetOffsetZ();

        for (int k = 0; k <= sz; ++k) {
            const int kg = offset_z + (k - pad);
            zb(static_cast<std::size_t>(k)) = static_cast<double>(kg) * dz;
        }
        for (int k = 0; k < sz; ++k) {
            const int kg = offset_z + (k - pad);
            zc(static_cast<std::size_t>(k)) = (static_cast<double>(kg) + 0.5) * dz;
        }
    }

    mesh.UpdateMetricsFromCoordinates();
    mesh.SetAllCellsFluid();

    const int i0 = mesh.GetCoreStartX();
    const int i1 = mesh.GetCoreEndExclusiveX();
    const int j0 = mesh.GetCoreStartY();
    const int j1 = mesh.GetCoreEndExclusiveY();
    const int k0 = mesh.GetCoreStartZ();
    const int k1 = mesh.GetCoreEndExclusiveZ();

    const std::string& ic_type = initial_conditions_.ic_type;

    auto& U = layer.U();

    auto write_conservative = [&](int i, int j, int k,
                                  double rho, double u, double v, double w, double P) {
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

    for (int k = k0; k < k1; ++k) {
        const double z = mesh.Zc()(static_cast<std::size_t>(k));
        for (int j = j0; j < j1; ++j) {
            const double y = mesh.Yc()(static_cast<std::size_t>(j));
            for (int i = i0; i < i1; ++i) {
                const double x = mesh.Xc()(static_cast<std::size_t>(i));

                double rho = 0.0;
                double u = 0.0;
                double v = 0.0;
                double w = 0.0;
                double P = 0.0;

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
                    if (dim < 2) {
                        throw std::runtime_error("ic_type=y_riemann requires dim>=2");
                    }
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
                    if (dim != 2) {
                        throw std::runtime_error("ic_type=quadrant currently supported only for dim==2");
                    }

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
                    if (dim != 3) {
                        throw std::runtime_error("ic_type=octant requires dim==3");
                    }

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

                if (dim < 2) {
                    v = 0.0;
                }
                if (dim < 3) {
                    w = 0.0;
                }

                write_conservative(i, j, k, rho, u, v, w, P);
            }
        }
    }
}

void Simulation::Initialize() {
    std::cout << "Initializing simulation..." << '\n';

    if (settings_.mpi_enabled) {
        mpi_context_ = std::make_unique<MPIContext>(MPI_COMM_WORLD, false);
        decomposition_ = std::make_unique<DomainDecomposition>(settings_, *mpi_context_);

        auto halo_exchange = std::make_shared<HaloExchange>(*mpi_context_);
        boundary_manager_ = std::make_shared<BoundaryManager>(halo_exchange);

        mesh_ = std::make_unique<Mesh>(
            decomposition_->LocalNx(),
            decomposition_->LocalNy(),
            decomposition_->LocalNz(),
            settings_.padding,
            settings_.dim
        );

        decomposition_->ApplyToMesh(*mesh_);
    }
    else {
        mesh_ = std::make_unique<Mesh>(
            settings_.GetNx(),
            settings_.GetNy(),
            settings_.GetNz(),
            settings_.padding,
            settings_.dim
        );

        mesh_->SetAllGlobalBoundaries(true);
        mesh_->SetGlobalDecomposition(
            mesh_->GetNx(),
            mesh_->GetNy(),
            mesh_->GetNz(),
            0, 0, 0
        );
    }

    layer_ = std::make_unique<DataLayer>(
        mesh_->GetSx(),
        mesh_->GetSy(),
        mesh_->GetSz()
    );

    ApplyInitialConditions(*layer_, *mesh_);

    if (settings_.immersed_enabled) {
        for (const auto& object : settings_.immersed_objects) {
            mesh_->AddPrimitive(CreateGeometryPrimitive(object));
        }

        mesh_->BuildCellTypesFromPrimitives();
        mesh_->BuildImmersedFaces();
    }

    FarfieldConservative far_field_U{};

    auto left_bc = CreateBoundaryCondition(settings_.left_boundary, far_field_U);
    auto right_bc = CreateBoundaryCondition(settings_.right_boundary, far_field_U);
    boundary_manager_->Set(Axis::X, left_bc, right_bc);

    if (settings_.dim >= 2) {
        auto bottom_bc = CreateBoundaryCondition(settings_.bottom_boundary, far_field_U);
        auto top_bc = CreateBoundaryCondition(settings_.top_boundary, far_field_U);
        boundary_manager_->Set(Axis::Y, bottom_bc, top_bc);
    }

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

        if (settings_.mpi_enabled) {
            analytical_mesh_ = std::make_unique<Mesh>(
                decomposition_->LocalNx(),
                decomposition_->LocalNy(),
                decomposition_->LocalNz(),
                analytical_settings_.padding,
                analytical_settings_.dim
            );

            decomposition_->ApplyToMesh(*analytical_mesh_);
        }
        else {
            analytical_mesh_ = std::make_unique<Mesh>(
                analytical_settings_.GetNx(),
                analytical_settings_.GetNy(),
                analytical_settings_.GetNz(),
                analytical_settings_.padding,
                analytical_settings_.dim
            );

            analytical_mesh_->SetAllGlobalBoundaries(true);
            analytical_mesh_->SetGlobalDecomposition(
                analytical_mesh_->GetNx(),
                analytical_mesh_->GetNy(),
                analytical_mesh_->GetNz(),
                0, 0, 0
            );
        }

        analytical_layer_ = std::make_unique<DataLayer>(
            analytical_mesh_->GetSx(),
            analytical_mesh_->GetSy(),
            analytical_mesh_->GetSz()
        );

        ApplyInitialConditions(*analytical_layer_, *analytical_mesh_);

        if (settings_.immersed_enabled) {
            for (const auto& object : settings_.immersed_objects) {
                analytical_mesh_->AddPrimitive(CreateGeometryPrimitive(object));
            }

            analytical_mesh_->BuildCellTypesFromPrimitives();
            analytical_mesh_->BuildImmersedFaces();
        }
    }

    CreateWriters();
}


void Simulation::Run() {
    Initialize();
    WriteInitialState();

    t_cur_ = 0.0;
    step_cur_ = 0;

    const bool is_root = !mpi_context_ || mpi_context_->IsRoot();

    if (is_root) {
        std::cout << "\nStarting simulation..." << '\n';
    }

    std::chrono::duration<double> runtime{0};
    const auto start_wall = std::chrono::high_resolution_clock::now();

    while (ShouldRun()) {
        const auto start = std::chrono::high_resolution_clock::now();
        dt_ = solver_->Step(*layer_, t_cur_);
        const auto end = std::chrono::high_resolution_clock::now();
        runtime += end - start;

        if (settings_.solver == "analytical") {
            t_cur_ += dt_;
        }

        ++step_cur_;

        WriteStepState(t_cur_, step_cur_);

        if (is_root) {
            PrintLog();
        }
    }

    if (mpi_context_) {
        mpi_context_->Barrier();
    }

    const auto end_wall = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> wall_time_local = end_wall - start_wall;
    const double runtime_local = runtime.count();

    double wall_time = wall_time_local.count();
    double computation_time = runtime_local;

    if (mpi_context_) {
        wall_time = mpi_context_->GlobalMax(wall_time);
        computation_time = mpi_context_->GlobalMax(computation_time);
    }

    if (is_root) {
        std::cout << '\n';
        std::cout << "\nSimulation completed!" << '\n';
        std::cout << ">>> Final time:  " << t_cur_ << '\n';
        std::cout << ">>> Total steps: " << step_cur_ << '\n';
        std::cout << ">>> Wall time: " << wall_time << "s\n";
        std::cout << ">>> Computation time: " << computation_time << "s\n";
    }

    FinalizeWriters();

    if (mpi_context_) {
        mpi_context_->Barrier();
    }
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
        vtk_writer_->Write(*layer_, *mesh_, settings_, 0, 0.0);
    }
    if (vtk_analytical_writer_ && analytical_layer_) {
        vtk_analytical_writer_->Write(*analytical_layer_, *analytical_mesh_, settings_, 0, 0.0);
    }
}

void Simulation::WriteStepState(double t_cur, std::size_t step_cur) const {
    if (!ShouldWrite()) {
        return;
    }
    if (vtk_writer_) {
        vtk_writer_->Write(*layer_, *mesh_, settings_, step_cur, t_cur);
    }
    if (vtk_analytical_writer_ && analytical_layer_) {
        vtk_analytical_writer_->Write(*analytical_layer_, *analytical_mesh_, settings_, step_cur, t_cur);
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
