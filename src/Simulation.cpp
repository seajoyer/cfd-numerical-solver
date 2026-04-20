#include "Simulation.hpp"

#include <chrono>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <utility>

#include "bc/BoundaryFactory.hpp"
#include "bc/BoundaryManager.hpp"
#include "config/InitialConditionInitializer.hpp"
#include "data/DataLayer.hpp"
#include "data/Workspace.hpp"
#include "geometry/DelaunayMeshBuilder.hpp"
#include "geometry/Face.hpp"
#include "geometry/GmshMeshBuilder.hpp"
#include "geometry/Mesh.hpp"
#include "geometry/StructuredMeshBuilder.hpp"
#include "output/StepWriter.hpp"
#include "output/WriterFactory.hpp"
#include "parallel/DomainDecomposition.hpp"
#include "parallel/HaloExchange.hpp"
#include "parallel/MPIContext.hpp"
#include "solver/Solver.hpp"
#include "solver/SolverFactory.hpp"
#include "utils/StringUtils.hpp"
#include "output/VTKRecomposer.hpp"
#include "parallel/MeshDistribution.hpp"

Simulation::Simulation(Settings settings, InitialConditions initial_conditions)
    : settings_(std::move(settings)),
      initial_conditions_(std::move(initial_conditions)) {}

Simulation::~Simulation() = default;

void Simulation::Initialize() {
    ValidateConfiguration();
    BuildMesh();
    ValidateBoundaryCoverage();
    AllocateState();
    InitializeFields();
    InitializeBoundaryConditions();
    InitializeSolver();
    InitializeWriter();
}

void Simulation::Run() {
    Initialize();

    t_cur_ = 0.0;
    step_cur_ = 0;
    dt_ = 0.0;

    WriteInitialState();

    if (IsRootRank()) {
        std::cout << "\nStarting simulation...\n";
    }

    std::chrono::duration<double> runtime{0.0};
    const auto start_wall = std::chrono::high_resolution_clock::now();

    while (ShouldRun()) {
        const auto start = std::chrono::high_resolution_clock::now();

        dt_ = solver_->Step(*layer_, t_cur_);

        const auto end = std::chrono::high_resolution_clock::now();

        runtime += end - start;
        ++step_cur_;

        WriteStepState();
        PrintLog();
    }

    const auto end_wall = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> wall_time = end_wall - start_wall;

    if (IsRootRank()) {
        std::cout << "\n\nSimulation completed!\n";
        std::cout << ">>> Final time:  " << t_cur_ << '\n';
        std::cout << ">>> Total steps: " << step_cur_ << '\n';
        std::cout << ">>> Wall time: " << wall_time.count() << "s\n";
        std::cout << ">>> Computation time: " << runtime.count() << "s\n";
    }

    FinalizeWriter();
}

DataLayer& Simulation::GetDataLayer() {
    if (!layer_) {
        throw std::runtime_error("Simulation: DataLayer is not initialized");
    }
    return *layer_;
}

std::size_t Simulation::GetCurrentStep() const {
    return step_cur_;
}

double Simulation::GetCurrentTime() const {
    return t_cur_;
}

bool Simulation::IsParallelRun() const {
    return settings_.mpi_enabled && mpi_context_ && mpi_context_->Size() > 1;
}

bool Simulation::IsRootRank() const {
    if (!mpi_context_) {
        return true;
    }
    return mpi_context_->IsRoot();
}

void Simulation::ValidateConfiguration() const {
    if (settings_.mesh.dim < 1 || settings_.mesh.dim > 3) {
        throw std::runtime_error("Simulation: mesh.dim must be 1, 2, or 3");
    }

    if (settings_.gamma <= 1.0) {
        throw std::runtime_error("Simulation: gamma must be greater than 1");
    }

    if (settings_.cfl <= 0.0) {
        throw std::runtime_error("Simulation: cfl must be positive");
    }

    if (settings_.t_end == 0.0 && settings_.step_end == 0) {
        throw std::runtime_error(
            "Simulation: both t_end and step_end are zero; simulation would never run"
        );
    }

    if (!IsKnownSolver(settings_.solver)) {
        throw std::runtime_error("Simulation: unknown solver type: " + settings_.solver);
    }

    if (!IsKnownTimeIntegrator(settings_.time_integrator)) {
        throw std::runtime_error(
            "Simulation: unknown time integrator type: " + settings_.time_integrator
        );
    }

    if (!IsKnownReconstruction(settings_.reconstruction)) {
        throw std::runtime_error(
            "Simulation: unknown reconstruction type: " + settings_.reconstruction
        );
    }

    if (!IsKnownRiemannSolver(settings_.riemann_solver)) {
        throw std::runtime_error(
            "Simulation: unknown Riemann solver type: " + settings_.riemann_solver
        );
    }

    for (const std::string& format : settings_.output_formats) {
        if (!IsKnownOutputFormat(format)) {
            throw std::runtime_error("Simulation: unsupported output format: " + format);
        }
    }

    const std::string solver = utils::ToLower(settings_.solver);
    const std::string reconstruction = utils::ToLower(settings_.reconstruction);
    const std::string time_integrator = utils::ToLower(settings_.time_integrator);

    if (solver == "godunov" && reconstruction != "p0") {
        throw std::runtime_error("Simulation: godunov requires reconstruction p0");
    }

    if (solver == "godunov-kolgan-rodionov" && reconstruction == "p0") {
        throw std::runtime_error(
            "Simulation: godunov-kolgan-rodionov requires at least p1"
        );
    }

    if (solver == "mader" && time_integrator != "mader") {
        throw std::runtime_error(
            "Simulation: mader solver requires mader time_integrator"
        );
    }

    if (solver == "mader") {
        const bool has_lambda_structured =
            initial_conditions_.structured_regions.has_value() &&
            initial_conditions_.structured_regions->reactant_mass_fraction.has_value();

        const bool has_lambda_constant =
            initial_conditions_.constant.has_value() &&
            initial_conditions_.constant->reactant_mass_fraction.has_value();

        if (!has_lambda_structured && !has_lambda_constant) {
            throw std::runtime_error(
                "Simulation: mader solver requires reactant_mass_fraction in initial condition"
            );
        }
    }
}

void Simulation::ValidateBoundaryCoverage() const {
    if (!mesh_) {
        throw std::runtime_error("Simulation::ValidateBoundaryCoverage: mesh is not initialized");
    }

    for (const Face& face : mesh_->Faces()) {
        if (!face.IsPhysicalBoundary()) {
            continue;
        }

        if (face.boundary_tag < 0) {
            throw std::runtime_error(
                "Simulation::ValidateBoundaryCoverage: boundary face has invalid boundary_tag"
            );
        }

        if (settings_.boundary.by_tag.find(face.boundary_tag) == settings_.boundary.by_tag.end()) {
            throw std::runtime_error(
                "Simulation::ValidateBoundaryCoverage: no boundary condition configured for boundary tag " +
                std::to_string(face.boundary_tag)
            );
        }
    }
}

void Simulation::BuildMesh() {
    mpi_context_ = std::make_unique<MPIContext>();

    if (!settings_.mpi_enabled || mpi_context_->Size() == 1) {
        mesh_ = CreateMesh();
        mesh_->Validate();
        return;
    }

    MeshDistribution::LocalPartition local_partition;

    if (mpi_context_->IsRoot()) {
        std::shared_ptr<Mesh> global_mesh = CreateMesh();
        local_partition =
            MeshDistribution::DistributeFromRoot(global_mesh.get(), *mpi_context_);
    }
    else {
        local_partition =
            MeshDistribution::DistributeFromRoot(nullptr, *mpi_context_);
    }

    mesh_ = std::make_shared<Mesh>(std::move(local_partition.mesh));
    halo_exchange_ = std::make_unique<HaloExchange>(*mpi_context_, std::move(local_partition.halos));

    mesh_->Validate();
}

void Simulation::AllocateState() {
    layer_ = std::make_unique<DataLayer>();
    layer_->ResizeFrom(*mesh_);

    workspace_ = std::make_unique<Workspace>();
    workspace_->ResizeFrom(*mesh_);
}

void Simulation::InitializeFields() {
    InitialConditionInitializer initializer(settings_, initial_conditions_);
    initializer.Apply(*mesh_, *layer_);

    if (halo_exchange_) {
        halo_exchange_->Synchronize(*layer_);
    }
}

void Simulation::InitializeBoundaryConditions() {
    boundary_manager_ = std::make_shared<BoundaryManager>();

    for (const auto& [tag, boundary_settings] : settings_.boundary.by_tag) {
        boundary_manager_->Register(tag, BoundaryFactory::Create(boundary_settings, settings_));
    }
}

void Simulation::InitializeSolver() {
    solver_ = CreateSolver();
    solver_->SetCfl(settings_.cfl);
}

void Simulation::InitializeWriter() {
    if (!settings_.HasOutputFormat("vtk")) {
        return;
    }

    vtk_writer_ = WriterFactory::Create("vtk", settings_.output_dir);
}

std::shared_ptr<Mesh> Simulation::CreateMesh() const {
    if (settings_.mesh.source_type == MeshSourceType::StructuredCartesian) {
        if (!settings_.mesh.structured.has_value()) {
            throw std::runtime_error("Simulation: structured mesh settings are missing");
        }

        const StructuredMeshSettings& s = *settings_.mesh.structured;

        if (settings_.mesh.dim == 2) {
            return std::make_shared<Mesh>(
                StructuredMeshBuilder::BuildUniformCartesian2D(
                    s.nx, s.ny,
                    s.x_min, s.x_max,
                    s.y_min, s.y_max
                )
            );
        }

        if (settings_.mesh.dim == 3) {
            return std::make_shared<Mesh>(
                StructuredMeshBuilder::BuildUniformCartesian3D(
                    s.nx, s.ny, s.nz,
                    s.x_min, s.x_max,
                    s.y_min, s.y_max,
                    s.z_min, s.z_max
                )
            );
        }

        throw std::runtime_error(
            "Simulation: only dim=2 and dim=3 are currently supported by StructuredMeshBuilder"
        );
    }

    if (settings_.mesh.source_type == MeshSourceType::GmshFile) {
        if (!settings_.mesh.gmsh_file.has_value()) {
            throw std::runtime_error("Simulation: gmsh file settings are missing");
        }

        return std::make_shared<Mesh>(
            GmshMeshBuilder::BuildFromFile(
                settings_.mesh.gmsh_file->file_path,
                settings_.mesh.dim
            )
        );
    }

    if (settings_.mesh.source_type == MeshSourceType::GmshGeo) {
        if (!settings_.mesh.gmsh_geo.has_value()) {
            throw std::runtime_error("Simulation: gmsh geo settings are missing");
        }

        return std::make_shared<Mesh>(
            GmshMeshBuilder::BuildFromGeoFile(
                settings_.mesh.gmsh_geo->file_path,
                settings_.mesh.dim
            )
        );
    }

    if (settings_.mesh.source_type == MeshSourceType::DelaunayGeo) {
        if (!settings_.mesh.delaunay_geo.has_value()) {
            throw std::runtime_error("Simulation: delaunay geo settings are missing");
        }

        return std::make_shared<Mesh>(
            DelaunayMeshBuilder::BuildFromGeoFile(
                settings_.mesh.delaunay_geo->file_path,
                settings_.mesh.dim
            )
        );
    }

    throw std::runtime_error("Simulation: unsupported mesh source type");
}

std::unique_ptr<Solver> Simulation::CreateSolver() {
    return SolverFactory::Create(settings_, mesh_, boundary_manager_, mpi_context_.get(), halo_exchange_.get());
}

bool Simulation::IsKnownSolver(const std::string& solver) const {
    const std::string s = utils::ToLower(solver);
    return s == "godunov" ||
        s == "godunov-kolgan" ||
        s == "godunov-kolgan-rodionov" ||
        s == "flic" ||
        s == "mader";
}

bool Simulation::IsKnownTimeIntegrator(const std::string& time_integrator) const {
    const std::string t = utils::ToLower(time_integrator);
    return t == "euler" ||
        t == "ssprk2" ||
        t == "ssprk3" ||
        t == "maccormack" ||
        t == "mader";
}

bool Simulation::IsKnownReconstruction(const std::string& reconstruction) const {
    const std::string r = utils::ToLower(reconstruction);
    return r == "p0" ||
        r == "p1" ||
        r == "eno3" ||
        r == "weno5";
}

bool Simulation::IsKnownRiemannSolver(const std::string& riemann_solver) const {
    const std::string r = utils::ToLower(riemann_solver);
    return r == "exact" ||
        r == "hll" ||
        r == "hllc" ||
        r == "acoustic" ||
        r == "roe" ||
        r == "rusanov" ||
        r == "osher";
}

bool Simulation::IsKnownOutputFormat(const std::string& format) const {
    return utils::ToLower(format) == "vtk";
}

bool Simulation::ShouldWrite() const {
    if (settings_.output_every_time == 0.0 && settings_.output_every_steps == 0) {
        return false;
    }

    const bool time_ok =
        t_cur_ >= settings_.t_end ||
        settings_.output_every_time == 0.0 ||
        std::floor((t_cur_ - dt_) / settings_.output_every_time) <
        std::floor(t_cur_ / settings_.output_every_time);

    const bool step_ok =
        settings_.output_every_steps == 0 ||
        step_cur_ % settings_.output_every_steps == 0;

    return (time_ok && step_ok) || t_cur_ >= settings_.t_end;
}

bool Simulation::ShouldLog() const {
    if (settings_.log_every_time == 0.0 && settings_.log_every_steps == 0) {
        return false;
    }

    const bool time_ok =
        settings_.log_every_time == 0.0 ||
        std::floor((t_cur_ - dt_) / settings_.log_every_time) <
        std::floor(t_cur_ / settings_.log_every_time);

    const bool step_ok =
        settings_.log_every_steps == 0 ||
        step_cur_ % settings_.log_every_steps == 0;

    return (time_ok && step_ok) || t_cur_ >= settings_.t_end;
}

bool Simulation::ShouldRun() const {
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
    if (vtk_writer_) {
        vtk_writer_->Write(*layer_, *mesh_, settings_, 0, 0.0);
    }
}

void Simulation::WriteStepState() const {
    if (!ShouldWrite()) {
        return;
    }

    if (vtk_writer_) {
        vtk_writer_->Write(*layer_, *mesh_, settings_, step_cur_, t_cur_);
    }
}

void Simulation::PrintLog() const {
    if (!ShouldLog()) {
        return;
    }

    if (!IsRootRank()) {
        return;
    }

    double progress = 0.0;
    if (settings_.t_end > 0.0) {
        progress = t_cur_ / settings_.t_end * 100.0;
    }

    const int percent = static_cast<int>(progress);

    std::cout << "\r";
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
