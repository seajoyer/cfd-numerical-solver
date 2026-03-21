#include "Simulation.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <utility>

#include "bc/BoundaryFactory.hpp"
#include "data/geometry/GeometryFactory.hpp"
#include "output/WriterFactory.hpp"
#include "solver/SolverFactory.hpp"
#include "utils/StringUtils.hpp"

Simulation::Simulation(Settings settings, const InitialConditions& initial_conditions)
    : settings_(std::move(settings)),
      initial_conditions_(initial_conditions),
      boundary_manager_(std::make_shared<BoundaryManager>(nullptr)) {}

auto Simulation::CreateSolver() -> std::unique_ptr<Solver> {
    return SolverFactory::Create(settings_, *mesh_, boundary_manager_, mpi_context_.get());
}

void Simulation::ApplyInitialConditions(DataLayer& layer, Mesh& mesh) {
    const int dim = mesh.GetDim();
    const double gamma = settings_.gamma;

    const int i0 = mesh.GetCoreStartX();
    const int i1 = mesh.GetCoreEndExclusiveX();
    const int j0 = mesh.GetCoreStartY();
    const int j1 = mesh.GetCoreEndExclusiveY();
    const int k0 = mesh.GetCoreStartZ();
    const int k1 = mesh.GetCoreEndExclusiveZ();

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

    auto region_index = [](double coord, const std::vector<double>& interfaces) -> std::size_t {
        return static_cast<std::size_t>(
            std::upper_bound(interfaces.begin(), interfaces.end(), coord) - interfaces.begin()
        );
    };

    const std::size_t expected_nx = initial_conditions_.RegionCountX();
    const std::size_t expected_ny = initial_conditions_.RegionCountY();
    const std::size_t expected_nz = initial_conditions_.RegionCountZ();

    if (initial_conditions_.rho.Nx() != expected_nx ||
        initial_conditions_.u.Nx() != expected_nx ||
        initial_conditions_.v.Nx() != expected_nx ||
        initial_conditions_.w.Nx() != expected_nx ||
        initial_conditions_.p.Nx() != expected_nx) {
        throw std::runtime_error("Initial condition x-shape does not match interface count");
    }

    if (initial_conditions_.rho.Ny() != expected_ny ||
        initial_conditions_.u.Ny() != expected_ny ||
        initial_conditions_.v.Ny() != expected_ny ||
        initial_conditions_.w.Ny() != expected_ny ||
        initial_conditions_.p.Ny() != expected_ny) {
        throw std::runtime_error("Initial condition y-shape does not match interface count");
    }

    if (initial_conditions_.rho.Nz() != expected_nz ||
        initial_conditions_.u.Nz() != expected_nz ||
        initial_conditions_.v.Nz() != expected_nz ||
        initial_conditions_.w.Nz() != expected_nz ||
        initial_conditions_.p.Nz() != expected_nz) {
        throw std::runtime_error("Initial condition z-shape does not match interface count");
    }

    for (int k = k0; k < k1; ++k) {
        const double z = mesh.Zc()(static_cast<std::size_t>(k));
        const std::size_t iz = dim >= 3
                                   ? region_index(z, initial_conditions_.interfaces_z)
                                   : 0;

        for (int j = j0; j < j1; ++j) {
            const double y = mesh.Yc()(static_cast<std::size_t>(j));
            const std::size_t iy = dim >= 2
                                       ? region_index(y, initial_conditions_.interfaces_y)
                                       : 0;

            for (int i = i0; i < i1; ++i) {
                const double x = mesh.Xc()(static_cast<std::size_t>(i));
                const std::size_t ix = region_index(x, initial_conditions_.interfaces_x);

                double rho = initial_conditions_.rho.At(ix, iy, iz);
                double u = initial_conditions_.u.At(ix, iy, iz);
                double v = initial_conditions_.v.At(ix, iy, iz);
                double w = initial_conditions_.w.At(ix, iy, iz);
                double P = initial_conditions_.p.At(ix, iy, iz);

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

void Simulation::InitializeParallel() {
    if (!settings_.mpi_enabled) {
        return;
    }

    mpi_context_ = std::make_unique<MPIContext>(MPI_COMM_WORLD, false);
    decomposition_ = std::make_unique<DomainDecomposition>(settings_, *mpi_context_);

    auto halo_exchange = std::make_shared<HaloExchange>(*mpi_context_);
    boundary_manager_ = std::make_shared<BoundaryManager>(halo_exchange);
}

auto Simulation::DeterminePadding() const -> int {
    int padding = 1;

    const std::string reconstruction = utils::ToLower(settings_.reconstruction);

    if (reconstruction == "p0") {
        padding = 1;
    }
    else if (reconstruction == "p1") {
        padding = 2;
    }
    else if (reconstruction == "eno3") {
        padding = 2;
    }
    else if (reconstruction == "weno5") {
        padding = 3;
    }
    else {
        throw std::runtime_error("Unknown reconstruction for automatic padding: " + reconstruction);
    }

    if (settings_.viscosity) {
        padding = std::max(padding, 2);
    }

    return padding;
}

void Simulation::ValidateConfiguration() const {
    if (settings_.dim < 1 || settings_.dim > 3) {
        throw std::runtime_error("dim must be 1, 2, or 3");
    }

    if (settings_.GetNx() <= 0) {
        throw std::runtime_error("Nx must be positive");
    }
    if (settings_.dim >= 2 && settings_.GetNy() <= 0) {
        throw std::runtime_error("Ny must be positive for dim >= 2");
    }
    if (settings_.dim >= 3 && settings_.GetNz() <= 0) {
        throw std::runtime_error("Nz must be positive for dim >= 3");
    }

    if (settings_.L_x <= 0.0) {
        throw std::runtime_error("L_x must be positive");
    }
    if (settings_.dim >= 2 && settings_.L_y <= 0.0) {
        throw std::runtime_error("L_y must be positive for dim >= 2");
    }
    if (settings_.dim >= 3 && settings_.L_z <= 0.0) {
        throw std::runtime_error("L_z must be positive for dim >= 3");
    }

    if (settings_.gamma <= 1.0) {
        throw std::runtime_error("gamma must be greater than 1");
    }

    if (settings_.cfl <= 0.0) {
        throw std::runtime_error("cfl must be positive");
    }

    if (settings_.t_end == 0.0 && settings_.step_end == 0) {
        throw std::runtime_error("Both t_end and step_end are zero; simulation would never run");
    }

    if (!IsKnownSolver(settings_.solver)) {
        throw std::runtime_error("Unknown solver type: " + settings_.solver);
    }

    if (!IsKnownTimeIntegrator(settings_.time_integrator)) {
        throw std::runtime_error("Unknown time integrator type: " + settings_.time_integrator);
    }

    if (!IsKnownReconstruction(settings_.reconstruction)) {
        throw std::runtime_error("Unknown reconstruction type: " + settings_.reconstruction);
    }

    if (!IsKnownRiemannSolver(settings_.riemann_solver)) {
        throw std::runtime_error("Unknown Riemann solver type: " + settings_.riemann_solver);
    }

    if (!IsKnownBoundaryCondition(settings_.left_boundary) ||
        !IsKnownBoundaryCondition(settings_.right_boundary)) {
        throw std::runtime_error("Unknown x-boundary condition");
    }

    if (settings_.dim >= 2) {
        if (!IsKnownBoundaryCondition(settings_.bottom_boundary) ||
            !IsKnownBoundaryCondition(settings_.top_boundary)) {
            throw std::runtime_error("Unknown y-boundary condition");
        }
    }

    if (settings_.dim >= 3) {
        if (!IsKnownBoundaryCondition(settings_.back_boundary) ||
            !IsKnownBoundaryCondition(settings_.front_boundary)) {
            throw std::runtime_error("Unknown z-boundary condition");
        }
    }

    for (const auto& format : settings_.output_formats) {
        if (!IsKnownOutputFormat(format)) {
            throw std::runtime_error("Unsupported output format: " + format);
        }
    }

    const std::string solver = utils::ToLower(settings_.solver);
    const std::string reconstruction = utils::ToLower(settings_.reconstruction);
    const std::string time_integrator = utils::ToLower(settings_.time_integrator);

    if (solver == "godunov") {
        if (reconstruction != "p0") {
            throw std::runtime_error("godunov requires reconstruction p0");
        }
    }

    if (solver == "godunov-kolgan-rodionov") {
        if (reconstruction == "p0") {
            throw std::runtime_error("godunov-kolgan-rodionov requires at least p1");
        }
    }

    if (solver == "mader") {
        if (time_integrator != "mader") {
            throw std::runtime_error("mader solver requires mader time_integrator");
        }
    }

    if (settings_.immersed_enabled) {
        for (const auto& object : settings_.immersed_objects) {
            const std::string type = utils::ToLower(object.type);

            if (type != "circle" && type != "rectangle") {
                throw std::runtime_error("Unknown immersed object type: " + object.type);
            }

            if (settings_.dim == 1) {
                throw std::runtime_error("Immersed objects are not supported for dim = 1");
            }

            if (type == "circle") {
                if (object.radius <= 0.0) {
                    throw std::runtime_error("Circle immersed object must have positive radius");
                }
            }

            if (type == "rectangle") {
                if (object.size_x <= 0.0 || object.size_y <= 0.0) {
                    throw std::runtime_error("Rectangle immersed object must have positive size_x and size_y");
                }
            }
        }
    }
}

void Simulation::InitializeMesh() {
    const int padding = DeterminePadding();

    if (settings_.mpi_enabled) {
        mesh_ = std::make_unique<Mesh>(
            decomposition_->LocalNx(),
            decomposition_->LocalNy(),
            decomposition_->LocalNz(),
            padding,
            settings_.dim
        );

        decomposition_->ApplyToMesh(*mesh_);
        return;
    }

    mesh_ = std::make_unique<Mesh>(
        settings_.GetNx(),
        settings_.GetNy(),
        settings_.GetNz(),
        padding,
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

void Simulation::InitializeDataLayer() {
    layer_ = std::make_unique<DataLayer>(
        mesh_->GetSx(),
        mesh_->GetSy(),
        mesh_->GetSz()
    );

    ApplyInitialConditions(*layer_, *mesh_);
}

void Simulation::InitializeGeometry() {
    if (!settings_.immersed_enabled) {
        return;
    }

    for (const auto& object : settings_.immersed_objects) {
        mesh_->AddPrimitive(CreateGeometryPrimitive(object));
    }

    mesh_->BuildCellTypesFromPrimitives();
    mesh_->BuildImmersedFaces();
}

void Simulation::InitializeBoundaryConditions() {
    FarfieldConservative far_field_U{};

    auto left_bc = BoundaryFactory::Create(settings_.left_boundary, far_field_U);
    auto right_bc = BoundaryFactory::Create(settings_.right_boundary, far_field_U);
    boundary_manager_->Set(Axis::X, left_bc, right_bc);

    if (settings_.dim >= 2) {
        auto bottom_bc = BoundaryFactory::Create(settings_.bottom_boundary, far_field_U);
        auto top_bc = BoundaryFactory::Create(settings_.top_boundary, far_field_U);
        boundary_manager_->Set(Axis::Y, bottom_bc, top_bc);
    }

    if (settings_.dim >= 3) {
        auto back_bc = BoundaryFactory::Create(settings_.back_boundary, far_field_U);
        auto front_bc = BoundaryFactory::Create(settings_.front_boundary, far_field_U);
        boundary_manager_->Set(Axis::Z, back_bc, front_bc);
    }
}

void Simulation::InitializeCoordinates() {
    const int dim = mesh_->GetDim();
    const int pad = mesh_->GetPadding();

    const int sx = mesh_->GetSx();
    const int sy = mesh_->GetSy();
    const int sz = mesh_->GetSz();

    const int global_nx = mesh_->GetGlobalNx();
    const int global_ny = mesh_->GetGlobalNy();
    const int global_nz = mesh_->GetGlobalNz();

    const double dx = settings_.L_x / static_cast<double>(global_nx);
    const double dy = dim >= 2 ? settings_.L_y / static_cast<double>(global_ny) : 1.0;
    const double dz = dim >= 3 ? settings_.L_z / static_cast<double>(global_nz) : 1.0;

    {
        auto& xb = mesh_->Xb();
        auto& xc = mesh_->Xc();
        const int offset_x = mesh_->GetOffsetX();

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
        auto& yb = mesh_->Yb();
        auto& yc = mesh_->Yc();
        const int offset_y = mesh_->GetOffsetY();

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
        auto& zb = mesh_->Zb();
        auto& zc = mesh_->Zc();
        const int offset_z = mesh_->GetOffsetZ();

        for (int k = 0; k <= sz; ++k) {
            const int kg = offset_z + (k - pad);
            zb(static_cast<std::size_t>(k)) = static_cast<double>(kg) * dz;
        }
        for (int k = 0; k < sz; ++k) {
            const int kg = offset_z + (k - pad);
            zc(static_cast<std::size_t>(k)) = (static_cast<double>(kg) + 0.5) * dz;
        }
    }

    mesh_->UpdateMetricsFromCoordinates();
    mesh_->SetAllCellsFluid();
}

void Simulation::InitializeSolver() {
    solver_ = CreateSolver();
    solver_->SetCfl(settings_.cfl);
}

void Simulation::InitializeWriter() {
    case_output_dir_ = settings_.output_dir;

    const int nx = settings_.GetNx();
    const int ny = settings_.GetNy();
    const int nz = settings_.GetNz();

    const int rank = mpi_context_ ? mpi_context_->Rank() : 0;
    const int size = mpi_context_ ? mpi_context_->Size() : 1;

    if (!settings_.HasOutputFormat("vtk")) {
        return;
    }

    std::ostringstream subdir;
    subdir << case_output_dir_ << "/vtk/"
        << settings_.solver
        << "__R_" << settings_.reconstruction
        << "__N_" << nx << "x" << ny << "x" << nz
        << "__CFL_" << utils::DoubleWithoutDot(settings_.cfl);

    vtk_writer_ = WriterFactory::Create("vtk", subdir.str(), false, rank, size);
}

auto Simulation::IsKnownSolver(const std::string& solver) const -> bool {
    const std::string s = utils::ToLower(solver);
    return s == "godunov" ||
        s == "godunov-kolgan" ||
        s == "godunov-kolgan-rodionov" ||
        s == "flic" ||
        s == "mader";
}

auto Simulation::IsKnownTimeIntegrator(const std::string& time_integrator) const -> bool {
    const std::string t = utils::ToLower(time_integrator);
    return t == "euler" ||
        t == "ssprk2" ||
        t == "ssprk3" ||
        t == "maccormack" ||
        t == "mader";
}

auto Simulation::IsKnownReconstruction(const std::string& reconstruction) const -> bool {
    const std::string r = utils::ToLower(reconstruction);
    return r == "p0" ||
        r == "p1" ||
        r == "eno3" ||
        r == "weno5";
}

auto Simulation::IsKnownRiemannSolver(const std::string& riemann_solver) const -> bool {
    const std::string r = utils::ToLower(riemann_solver);
    return r == "exact" ||
        r == "hll" ||
        r == "hllc" ||
        r == "acoustic" ||
        r == "roe" ||
        r == "rusanov" ||
        r == "osher";
}

auto Simulation::IsKnownBoundaryCondition(const std::string& bc) const -> bool {
    const std::string b = utils::ToLower(bc);
    return b == "outlet" ||
        b == "reflective" ||
        b == "free_stream" ||
        b == "inlet" ||
        b == "periodic" ||
        b == "symmetry" ||
        b == "wall" ||
        b == "non_reflective";
}

auto Simulation::IsKnownOutputFormat(const std::string& format) const -> bool {
    return utils::ToLower(format) == "vtk";
}

void Simulation::Initialize() {
    std::cout << "Initializing simulation...\n";

    ValidateConfiguration();

    InitializeParallel();
    InitializeMesh();
    InitializeCoordinates();
    InitializeDataLayer();
    InitializeGeometry();
    InitializeBoundaryConditions();
    InitializeSolver();
    InitializeWriter();
}

void Simulation::Run() {
    Initialize();
    WriteInitialState();

    t_cur_ = 0.0;
    step_cur_ = 0;

    const bool is_root = !mpi_context_ || mpi_context_->IsRoot();

    if (is_root) {
        std::cout << "\nStarting simulation...\n";
    }
    if (mpi_context_) {
        mpi_context_->Barrier();
    }

    std::chrono::duration<double> runtime{0};
    const auto start_wall = std::chrono::high_resolution_clock::now();

    while (ShouldRun()) {
        const auto start = std::chrono::high_resolution_clock::now();
        dt_ = solver_->Step(*layer_, t_cur_);
        const auto end = std::chrono::high_resolution_clock::now();
        runtime += end - start;

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
        std::cout << "\nSimulation completed!\n";
        std::cout << ">>> Final time:  " << t_cur_ << '\n';
        std::cout << ">>> Total steps: " << step_cur_ << '\n';
        std::cout << ">>> Wall time: " << wall_time << "s\n";
        std::cout << ">>> Computation time: " << computation_time << "s\n";
    }

    FinalizeWriter();

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

    const bool time_ok = t_cur_ >= settings_.t_end ||
        settings_.output_every_time == 0.0 ||
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

    const bool step_ok = settings_.log_every_steps == 0 ||
        step_cur_ % settings_.log_every_steps == 0;

    return (time_ok && step_ok) || t_cur_ >= settings_.t_end;
}

auto Simulation::ShouldRun() const -> bool {
    if (settings_.t_end == 0.0 && settings_.step_end == 0) {
        return false;
    }

    const bool time_not_exceeded = settings_.t_end == 0.0 || t_cur_ < settings_.t_end;
    const bool steps_not_exceeded = settings_.step_end == 0 || step_cur_ < settings_.step_end;

    return time_not_exceeded && steps_not_exceeded;
}

void Simulation::WriteInitialState() const {
    if (!mpi_context_ || mpi_context_->IsRoot()) {
        std::cout << "Writing the initial state...\n";
    }

    if (vtk_writer_) {
        vtk_writer_->Write(*layer_, *mesh_, settings_, 0, 0.0);
    }
}

void Simulation::WriteStepState(double t_cur, std::size_t step_cur) const {
    if (!ShouldWrite()) {
        return;
    }

    if (vtk_writer_) {
        vtk_writer_->Write(*layer_, *mesh_, settings_, step_cur, t_cur);
    }
}

void Simulation::PrintLog() const {
    if (!ShouldLog()) {
        return;
    }

    double progress = 0.0;
    if (settings_.t_end > 0.0) {
        progress = t_cur_ / settings_.t_end * 100.0;
    }

    const int percent = static_cast<int>(progress);

    std::cout << "\r \r";
    std::cout << ">>> [PROGRESS]: Step " << step_cur_
        << ", " << percent
        << "% processed, time: " << t_cur_
        << " of " << settings_.t_end;
    std::cout.flush();
}

void Simulation::FinalizeWriter() {
    if (vtk_writer_ && vtk_writer_->RequiresFinalization()) {
        vtk_writer_->Finalize(settings_);
    }
}
