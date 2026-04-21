#include "Simulation.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <filesystem>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <utility>

#include "bc/BoundaryFactory.hpp"
#include "data/geometry/GeometryFactory.hpp"
#include "output/WriterFactory.hpp"
#include "solver/SolverFactory.hpp"
#include "utils/StringUtils.hpp"
#include "solver/PressureVelocitySolver.hpp"

Simulation::Simulation(Settings settings, const InitialConditions& initial_conditions)
    : settings_(std::move(settings)),
      initial_conditions_(initial_conditions),
      boundary_manager_(std::make_shared<BoundaryManager>(nullptr)) {}

auto Simulation::CreateSolver() -> std::unique_ptr<Solver> {
    return SolverFactory::Create(settings_, *mesh_, boundary_manager_, mpi_context_.get(), eos_);
}

auto Simulation::PrimitiveToFarfieldConservative(const BoundaryStateSettings& s) -> FarfieldConservative {
    FarfieldConservative out;
    out.rho = s.rho;
    out.rhoU = s.rho * s.u;
    out.rhoV = s.rho * s.v;
    out.rhoW = s.rho * s.w;

    const double kinetic = 0.5 * (s.u * s.u + s.v * s.v + s.w * s.w);

    const double I = eos_->ComputeInternalEnergy(s.rho, s.p, 1.0);

    out.E = s.rho * (I + kinetic);
    return out;
}

auto Simulation::IsPressureVelocitySolver() const -> bool {
    const std::string solver = utils::ToLower(settings_.solver);
    return solver == "simple" || solver == "piso" || solver == "pimple";
}

void Simulation::ApplyInitialConditions(DataLayer& layer, Mesh& mesh) {
    const int dim = mesh.GetDim();

    const int i0 = mesh.GetCoreStartX();
    const int i1 = mesh.GetCoreEndExclusiveX();
    const int j0 = mesh.GetCoreStartY();
    const int j1 = mesh.GetCoreEndExclusiveY();
    const int k0 = mesh.GetCoreStartZ();
    const int k1 = mesh.GetCoreEndExclusiveZ();

    auto& U = layer.U();
    auto* lambda_ptr = initial_conditions_.reactant_mass_fraction.has_value()
                           ? &layer.ReactantMassFraction()
                           : nullptr;

    auto write_conservative = [&](int i, int j, int k,
                                  double rho, double u, double v, double w, double P, double lambda_val) {
        const double rhoU = rho * u;
        const double rhoV = rho * v;
        const double rhoW = rho * w;

        const double kinetic_specific = 0.5 * (u * u + v * v + w * w);

        const double I = eos_->ComputeInternalEnergy(rho, P, lambda_val);
        const double E = rho * (I + kinetic_specific);

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

    if (initial_conditions_.reactant_mass_fraction.has_value()) {
        const auto& lambda_ic = *initial_conditions_.reactant_mass_fraction;

        if (lambda_ic.Nx() != expected_nx) {
            throw std::runtime_error("Reactant mass fraction x-shape does not match interface count");
        }
        if (lambda_ic.Ny() != expected_ny) {
            throw std::runtime_error("Reactant mass fraction y-shape does not match interface count");
        }
        if (lambda_ic.Nz() != expected_nz) {
            throw std::runtime_error("Reactant mass fraction z-shape does not match interface count");
        }
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

                double lambda_val = 1.0;
                if (lambda_ptr) {
                    lambda_val = initial_conditions_.reactant_mass_fraction->At(ix, iy, iz);
                    (*lambda_ptr)(i, j, k) = lambda_val;
                }

                write_conservative(i, j, k, rho, u, v, w, P, lambda_val);
            }
        }
    }
}

void Simulation::InitializeParallel() {
    if (!settings_.mpi_enabled) {
        is_root_ = true;
        return;
    }

    mpi_context_ = std::make_unique<MPIContext>(MPI_COMM_WORLD, false);
    decomposition_ = std::make_unique<DomainDecomposition>(settings_, *mpi_context_);

    is_root_ = !mpi_context_ || mpi_context_->IsRoot();

    auto halo_exchange = std::make_shared<HaloExchange>(decomposition_->CartComm(), mpi_context_->Size(),
                                                        utils::ToLower(settings_.solver) == "mader");
    boundary_manager_ = std::make_shared<BoundaryManager>(halo_exchange);
}

auto Simulation::DeterminePadding() const -> int {
    const std::string solver = utils::ToLower(settings_.solver);

    if (solver == "simple" || solver == "piso" || solver == "pimple") {
        return 1;
    }

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

    for (const auto& format : settings_.output_formats) {
        if (!IsKnownOutputFormat(format)) {
            throw std::runtime_error("Unsupported output format: " + format);
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

    const std::string solver = utils::ToLower(settings_.solver);
    const bool is_pressure_velocity =
        solver == "simple" || solver == "piso" || solver == "pimple";

    auto validate_periodic_pair = [](const std::string& left_name,
                                     const std::string& right_name,
                                     const std::string& axis_name) -> void {
        const bool left_periodic = utils::ToLower(left_name) == "periodic";
        const bool right_periodic = utils::ToLower(right_name) == "periodic";

        if (left_periodic != right_periodic) {
            throw std::runtime_error(
                "Periodic boundary on axis " + axis_name +
                " must be specified on both sides"
            );
        }
    };

    auto require_boundary_state = [&](const std::string& bc_name,
                                      const std::optional<BoundaryStateSettings>& state,
                                      const std::string& side_name) {
        const std::string bc = utils::ToLower(bc_name);
        if ((bc == "inlet" || bc == "free_stream") && !state.has_value()) {
            throw std::runtime_error(
                "Boundary '" + side_name + "' uses " + bc +
                " but no boundary state is provided in boundary_conditions.states"
            );
        }
    };

    if (is_pressure_velocity) {
        if (settings_.dim != 2) {
            throw std::runtime_error("simple/piso/pimple currently support only dim = 2");
        }

        if (settings_.density <= 0.0) {
            throw std::runtime_error("density must be positive for pressure-velocity solvers");
        }

        if (settings_.kinematic_viscosity <= 0.0) {
            throw std::runtime_error("kinematic_viscosity must be positive for pressure-velocity solvers");
        }

        if (settings_.n_pressure_correctors <= 0) {
            throw std::runtime_error("n_pressure_correctors must be > 0");
        }

        if (settings_.pressure_max_iterations <= 0) {
            throw std::runtime_error("pressure_max_iterations must be > 0");
        }

        if (settings_.momentum_max_iterations <= 0) {
            throw std::runtime_error("momentum_max_iterations must be > 0");
        }

        if (settings_.pressure_tolerance <= 0.0) {
            throw std::runtime_error("pressure_tolerance must be > 0");
        }

        if (settings_.steady_tolerance <= 0.0) {
            throw std::runtime_error("steady_tolerance must be > 0");
        }

        if (solver == "simple") {
            if (settings_.velocity_relaxation <= 0.0 || settings_.velocity_relaxation > 1.0) {
                throw std::runtime_error("velocity_relaxation for SIMPLE must be in (0, 1]");
            }
            if (settings_.pressure_relaxation <= 0.0 || settings_.pressure_relaxation > 1.0) {
                throw std::runtime_error("pressure_relaxation for SIMPLE must be in (0, 1]");
            }
        }

        if (solver == "pimple") {
            if (settings_.n_outer_correctors <= 0) {
                throw std::runtime_error("n_outer_correctors must be > 0 for PIMPLE");
            }
            if (settings_.velocity_relaxation <= 0.0 || settings_.velocity_relaxation > 1.0) {
                throw std::runtime_error("velocity_relaxation for PIMPLE must be in (0, 1]");
            }
            if (settings_.pressure_relaxation <= 0.0 || settings_.pressure_relaxation > 1.0) {
                throw std::runtime_error("pressure_relaxation for PIMPLE must be in (0, 1]");
            }
        }

        if (!IsKnownBoundaryCondition(settings_.left_boundary) ||
            !IsKnownBoundaryCondition(settings_.right_boundary)) {
            throw std::runtime_error("Unknown x-boundary condition");
        }

        if (!IsKnownBoundaryCondition(settings_.bottom_boundary) ||
            !IsKnownBoundaryCondition(settings_.top_boundary)) {
            throw std::runtime_error("Unknown y-boundary condition");
        }

        validate_periodic_pair(settings_.left_boundary, settings_.right_boundary, "X");
        validate_periodic_pair(settings_.bottom_boundary, settings_.top_boundary, "Y");

        require_boundary_state(settings_.left_boundary, settings_.boundary_states.x_min, "x_min");
        require_boundary_state(settings_.right_boundary, settings_.boundary_states.x_max, "x_max");
        require_boundary_state(settings_.bottom_boundary, settings_.boundary_states.y_min, "y_min");
        require_boundary_state(settings_.top_boundary, settings_.boundary_states.y_max, "y_max");

        return;
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

    if (utils::ToLower(settings_.solver) == "mader" && !initial_conditions_.reactant_mass_fraction.has_value()) {
        throw std::runtime_error("Mader solver requires initial_condition.reactant_mass_fraction");
    }

    require_boundary_state(settings_.left_boundary, settings_.boundary_states.x_min, "x_min");
    require_boundary_state(settings_.right_boundary, settings_.boundary_states.x_max, "x_max");

    if (settings_.dim >= 2) {
        require_boundary_state(settings_.bottom_boundary, settings_.boundary_states.y_min, "y_min");
        require_boundary_state(settings_.top_boundary, settings_.boundary_states.y_max, "y_max");
    }

    if (settings_.dim >= 3) {
        require_boundary_state(settings_.back_boundary, settings_.boundary_states.z_min, "z_min");
        require_boundary_state(settings_.front_boundary, settings_.boundary_states.z_max, "z_max");
    }

    validate_periodic_pair(settings_.left_boundary, settings_.right_boundary, "X");

    if (settings_.dim >= 2) {
        validate_periodic_pair(settings_.bottom_boundary, settings_.top_boundary, "Y");
    }

    if (settings_.dim >= 3) {
        validate_periodic_pair(settings_.back_boundary, settings_.front_boundary, "Z");
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

    if (!IsPressureVelocitySolver()) {
        ApplyInitialConditions(*layer_, *mesh_);
    }
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
    FarfieldConservative x_min_state{};
    FarfieldConservative x_max_state{};
    FarfieldConservative y_min_state{};
    FarfieldConservative y_max_state{};
    FarfieldConservative z_min_state{};
    FarfieldConservative z_max_state{};

    if (settings_.boundary_states.x_min) {
        x_min_state = PrimitiveToFarfieldConservative(*settings_.boundary_states.x_min);
    }
    if (settings_.boundary_states.x_max) {
        x_max_state = PrimitiveToFarfieldConservative(*settings_.boundary_states.x_max);
    }
    if (settings_.boundary_states.y_min) {
        y_min_state = PrimitiveToFarfieldConservative(*settings_.boundary_states.y_min);
    }
    if (settings_.boundary_states.y_max) {
        y_max_state = PrimitiveToFarfieldConservative(*settings_.boundary_states.y_max);
    }
    if (settings_.boundary_states.z_min) {
        z_min_state = PrimitiveToFarfieldConservative(*settings_.boundary_states.z_min);
    }
    if (settings_.boundary_states.z_max) {
        z_max_state = PrimitiveToFarfieldConservative(*settings_.boundary_states.z_max);
    }

    int mpi_size = 1;
    if (mpi_context_) {
        mpi_size = mpi_context_->Size();
    }

    const std::string solver = utils::ToLower(settings_.solver);
    const bool is_pressure_velocity =
        solver == "simple" || solver == "piso" || solver == "pimple";

    auto create_bc = [&](const std::string& bc_name,
                         const std::optional<BoundaryStateSettings>& primitive_state,
                         const FarfieldConservative& farfield_state) -> std::shared_ptr<BoundaryCondition> {
        const std::string type = utils::ToLower(bc_name);

        if (is_pressure_velocity) {
            if (type == "inlet") {
                if (!primitive_state.has_value()) {
                    throw std::runtime_error(
                        "Inlet boundary requires primitive boundary state for pressure-velocity solvers");
                }
                return BoundaryFactory::Create(type, *primitive_state, settings_, mpi_size);
            }

            return BoundaryFactory::Create(type);
        }

        return BoundaryFactory::Create(type, farfield_state, settings_, mpi_size);
    };

    auto left_bc = create_bc(settings_.left_boundary,
                             settings_.boundary_states.x_min,
                             x_min_state);
    auto right_bc = create_bc(settings_.right_boundary,
                              settings_.boundary_states.x_max,
                              x_max_state);
    boundary_manager_->Set(Axis::X, left_bc, right_bc);

    if (settings_.dim >= 2) {
        auto bottom_bc = create_bc(settings_.bottom_boundary,
                                   settings_.boundary_states.y_min,
                                   y_min_state);
        auto top_bc = create_bc(settings_.top_boundary,
                                settings_.boundary_states.y_max,
                                y_max_state);
        boundary_manager_->Set(Axis::Y, bottom_bc, top_bc);
    }

    if (settings_.dim >= 3) {
        auto back_bc = create_bc(settings_.back_boundary,
                                 settings_.boundary_states.z_min,
                                 z_min_state);
        auto front_bc = create_bc(settings_.front_boundary,
                                  settings_.boundary_states.z_max,
                                  z_max_state);
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

void Simulation::InitializePressureVelocityInitialConditions() {
    auto* pv_solver = dynamic_cast<PressureVelocitySolver*>(solver_.get());
    if (!pv_solver) {
        throw std::runtime_error(
            "InitializePressureVelocityInitialConditions: solver is not PressureVelocitySolver");
    }

    auto& state = pv_solver->GetState();
    state.ResizeFrom(*mesh_);
    state.ZeroAll();

    if (mesh_->GetDim() != 2) {
        throw std::runtime_error(
            "Pressure-velocity initial conditions are currently supported only for dim = 2");
    }
    // CHECK: TG_INIT
    if (initial_conditions_.ic_type != "taylor_green") {
        // For now we only support Taylor-Green initialization
        // in the pressure-velocity branch.
        state.CopyCurrentToOld();
        boundary_manager_->ApplyPhysicalBc(state, *mesh_);
        return;
    }

    const int i0 = mesh_->GetCoreStartX();
    const int i1 = mesh_->GetCoreEndExclusiveX();
    const int j0 = mesh_->GetCoreStartY();
    const int j1 = mesh_->GetCoreEndExclusiveY();

    const auto& xc = mesh_->Xc();
    const auto& yc = mesh_->Yc();
    const auto& xb = mesh_->Xb();
    const auto& yb = mesh_->Yb();

    auto& p = state.Pressure();
    auto& ux = state.Ux();
    auto& vy = state.Vy();
    auto& wz = state.Wz();

    for (int j = j0; j < j1; ++j) {
        for (int i = i0; i < i1; ++i) {
            const double x = xc(static_cast<std::size_t>(i));
            const double y = yc(static_cast<std::size_t>(j));

            p(i, j, 0) =
                -0.25 * (std::cos(2.0 * M_PI * x) + std::cos(2.0 * M_PI * y));
        }
    }

    for (int j = j0; j < j1; ++j) {
        const double y = yc(static_cast<std::size_t>(j));
        for (int i = i0; i <= i1; ++i) {
            const double x = xb(static_cast<std::size_t>(i));
            ux(i, j, 0) = -std::cos(M_PI * x) * std::sin(M_PI * y);
        }
    }

    for (int j = j0; j <= j1; ++j) {
        const double y = yb(static_cast<std::size_t>(j));
        for (int i = i0; i < i1; ++i) {
            const double x = xc(static_cast<std::size_t>(i));
            vy(i, j, 0) = std::sin(M_PI * x) * std::cos(M_PI * y);
        }
    }

    wz.fill(0.0);

    state.CopyCurrentToOld();
    boundary_manager_->ApplyPhysicalBc(state, *mesh_);
    pv_solver->ApplyImmersedVelocityConstraints();
}

void Simulation::InitializeSolver() {
    solver_ = CreateSolver();
    solver_->SetCfl(settings_.cfl);

    if (IsPressureVelocitySolver()) {
        InitializePressureVelocityInitialConditions();
    }
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

    const std::string vtk_dir = subdir.str();

    if (size > 1) {
        const std::string rank_dir = vtk_dir + "/" + std::format("rank_{:04d}", rank);
        std::filesystem::create_directories(rank_dir);
    }
    else {
        std::filesystem::create_directories(vtk_dir);
    }

    if (mpi_context_) {
        mpi_context_->Barrier();
    }

    vtk_writer_ = WriterFactory::Create("vtk", vtk_dir, eos_, false, rank, size);
}


void Simulation::InitializeEos() {
    eos_ = std::make_shared<EOS>();

    const std::string eos_name = settings_.EOS.has_value() ? utils::ToLower(*settings_.EOS) : "ideal_gas";

    if (eos_name == "hugoniot_gruneisen") {
        eos_->SetType(EosType::HugoniotGruneisen);
        HugoniotGruneisenEosParameters hg_params;
        if (settings_.hg_rho0) hg_params.rho0 = settings_.hg_rho0;
        if (settings_.hg_C) hg_params.C = settings_.hg_C;
        if (settings_.hg_S) hg_params.S = settings_.hg_S;
        if (settings_.hg_gamma_s) hg_params.gamma_s = settings_.hg_gamma_s;
        if (settings_.hg_c_v) hg_params.c_v = settings_.hg_c_v;
        eos_->SetHugoniotGruneisenParameters(hg_params);
    }
    else {
        eos_->SetType(EosType::IdealGas);
        IdealGasEosParameters ig_params;
        ig_params.gamma = settings_.gamma;
        eos_->SetIdealGasParameters(ig_params);
    }
}

auto Simulation::IsKnownSolver(const std::string& solver) const -> bool {
    const std::string s = utils::ToLower(solver);
    return s == "godunov" ||
        s == "godunov-kolgan" ||
        s == "godunov-kolgan-rodionov" ||
        s == "flic" ||
        s == "mader" ||
        s == "simple" ||
        s == "piso" ||
        s == "pimple";
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
    ValidateConfiguration();

    InitializeParallel();
    if (is_root_) {
        std::cout << "Initializing simulation...\n";
    }
    InitializeEos();
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

    t_cur_ = 0.0;
    step_cur_ = 0;

    WriteInitialState();

    if (is_root_) {
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

        if (is_root_) {
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

    if (is_root_) {
        std::cout << "\n\nSimulation completed!\n";
        std::cout << ">>> Final time:  " << t_cur_ << '\n';
        std::cout << ">>> Total steps: " << step_cur_ << '\n';
        std::cout << ">>> Wall time: " << wall_time << "s\n";
        std::cout << ">>> Computation time: " << computation_time << "s\n";
    }

    if (mpi_context_) {
        mpi_context_->Barrier();
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
    if (is_root_) {
        std::cout << "Writing the initial state...\n";
    }

    if (!vtk_writer_) {
        return;
    }

    if (const auto* pv_solver = dynamic_cast<const PressureVelocitySolver*>(solver_.get())) {
        vtk_writer_->Write(pv_solver->GetState(), *mesh_, settings_, 0, 0.0);
        return;
    }

    vtk_writer_->Write(*layer_, *mesh_, settings_, 0, 0.0);
}

void Simulation::WriteStepState(const double t_cur, const std::size_t step_cur) const {
    if (!ShouldWrite()) {
        return;
    }

    if (!vtk_writer_) {
        return;
    }

    if (const auto* pv_solver = dynamic_cast<const PressureVelocitySolver*>(solver_.get())) {
        vtk_writer_->Write(pv_solver->GetState(), *mesh_, settings_, step_cur, t_cur);
        return;
    }

    vtk_writer_->Write(*layer_, *mesh_, settings_, step_cur, t_cur);
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

    std::cout << "\r" << ' ' * 200 << "\r";
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
