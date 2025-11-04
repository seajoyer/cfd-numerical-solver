#include "Simulation.hpp"

#include <iostream>
#include <stdexcept>
#include <utility>

#include "bc/BoundaryFactory.hpp"
#include "output/VTKWriter.hpp"
#include "solver/SolverFactory.hpp"

Simulation::Simulation(Settings settings, InitialConditions initial_conditions)
    : settings_(std::move(settings)), initial_conditions_(initial_conditions) {}

auto Simulation::CreateSolver() -> std::unique_ptr<Solver> {
    return SolverFactory::Create(settings_.solver, settings_.dim);
}

auto Simulation::CreateBoundaryCondition(const std::string& boundary_type, double rho_inf,
                                         double u_inf, double p_inf)
    -> std::shared_ptr<BoundaryCondition> {
    return BoundaryFactory::Create(boundary_type, rho_inf, u_inf, p_inf);
}

void Simulation::ApplyInitialConditions() {
    if (!layer_) {
        throw std::runtime_error(
            "DataLayer must be initialized before applying initial conditions");
    }

    // Get grid parameters
    const int n_ghost = layer_->GetNGhostCells();
    const int n = layer_->GetN();
    const int total_size = layer_->GetTotalSize();

    // Calculate cell positions
    const double dx = settings_.L_x / static_cast<double>(n);
    const double x_mid = settings_.L_x / 2.0;

    // Initialize grid coordinates
    for (int i = 0; i < total_size; ++i) {
        const int cell_index = i - n_ghost;
        layer_->xb(i) = cell_index * dx;
        layer_->xc(i) = (cell_index + 0.5) * dx;
    }

    // Apply left and right states (Riemann problem setup)
    const int core_start = layer_->GetCoreStart(0);
    const int core_end = layer_->GetCoreEndExclusive(0);

    for (int i = core_start; i < core_end; ++i) {
        const double x = layer_->xc(i);

        if (x < x_mid) {
            // Left state
            layer_->rho(i) = initial_conditions_.rho_L;
            layer_->u(i) = initial_conditions_.u_L;
            layer_->P(i) = initial_conditions_.P_L;
        } else {
            // Right state
            layer_->rho(i) = initial_conditions_.rho_R;
            layer_->u(i) = initial_conditions_.u_R;
            layer_->P(i) = initial_conditions_.P_R;
        }

        // Calculate conservative variables and other derived quantities
        const double rho = layer_->rho(i);
        const double u = layer_->u(i);
        const double P = layer_->P(i);
        const double gamma = settings_.gamma;

        layer_->m(i) = rho * u;                                // momentum
        layer_->e(i) = P / (gamma - 1.0) + 0.5 * rho * u * u;  // total energy
        layer_->p(i) = rho * u;  // momentum (alternative representation)
    }
}

void Simulation::Initialize() {
    std::cout << "Initializing simulation..." << '\n';

    // Initialize data layer
    layer_ = std::make_unique<DataLayer>(settings_.N, settings_.padding, settings_.dim);

    // Apply initial conditions
    ApplyInitialConditions();

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

    // Initialize output writer
    writer_ = std::make_unique<VTKWriter>(settings_.output_dir);

    std::cout << '\n';
    std::cout << "Simulation initialized:" << '\n';
    std::cout << "  Solver: " << settings_.solver << '\n';
    std::cout << "  Boundary left:  " << settings_.left_boundary << '\n';
    std::cout << "  Boundary right: " << settings_.right_boundary << '\n';
    std::cout << "  Grid size (N): " << settings_.N << '\n';
    std::cout << "  Dimension: " << settings_.dim << '\n';
    std::cout << "  Domain length (L_x): " << settings_.L_x << '\n';
    std::cout << "  CFL: " << settings_.cfl << '\n';
    std::cout << "  Gamma: " << settings_.gamma << '\n';
    std::cout << "  End time: " << settings_.t_end << '\n';
    std::cout << "  Output directory: " << settings_.output_dir << '\n';
    std::cout << "  Initial conditions:" << '\n';
    std::cout << "    Left:  rho=" << initial_conditions_.rho_L
              << ", u=" << initial_conditions_.u_L << ", P=" << initial_conditions_.P_L
              << '\n';
    std::cout << "    Right: rho=" << initial_conditions_.rho_R
              << ", u=" << initial_conditions_.u_R << ", P=" << initial_conditions_.P_R
              << "\n\n";
}

auto Simulation::GetDataLayer() -> DataLayer& {
    if (!layer_) {
        throw std::runtime_error("DataLayer is not initialized.");
    }
    return *layer_;
}

auto Simulation::GetCurrentStep() const -> std::size_t { return step_; }

auto Simulation::GetCurrentTime() const -> double { return t_cur_; }

auto Simulation::ShouldWrite(std::size_t step) const -> bool {
    if (settings_.output_every_steps == 0) return false;
    return (step % settings_.output_every_steps) == 0;
}

void Simulation::WriteInitialState() const {
    if (writer_) {
        writer_->Write(*layer_, 0, 0.0);
    }
}

void Simulation::WriteStepState(std::size_t step, double t_cur) const {
    if (writer_ && ShouldWrite(step)) {
        writer_->Write(*layer_, step, t_cur);
    }
}

void Simulation::Run() {
    Initialize();
    WriteInitialState();

    t_cur_ = 0.0;
    step_ = 0;

    std::cout << "\nStarting simulation..." << '\n';

    while (t_cur_ < settings_.t_end) {
        // Advance one time step and get the actual dt used
        double dt = solver_->Step(*layer_, t_cur_);

        t_cur_ += dt;
        ++step_;

        // Write output if needed
        WriteStepState(step_, t_cur_);

        // Progress output
        if (0 == step_ % settings_.output_every_steps) {
            std::cout << "Step " << step_ << ", time = " << t_cur_ << '\n';
        }
    }

    std::cout << "Simulation completed!" << '\n';
    std::cout << "  Final time: " << t_cur_ << '\n';
    std::cout << "  Total steps: " << step_ << '\n';
}
