#include "Simulation.hpp"

#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <utility>

#include "bc/BoundaryFactory.hpp"
#include "output/WriterFactory.hpp"
#include "solver/SolverFactory.hpp"

Simulation::Simulation(Settings settings, const InitialConditions& initial_conditions,
                       const bool log_progress) :
    settings_(std::move(settings)),
    initial_conditions_(initial_conditions),
    log_progress_(log_progress) {
}

auto Simulation::CreateSolver() -> std::unique_ptr<Solver> {
    return SolverFactory::Create(settings_);
}

auto Simulation::CreateBoundaryCondition(const std::string& boundary_type, double rho_inf,
                                         double u_inf, double p_inf)
    -> std::shared_ptr<BoundaryCondition> {
    return BoundaryFactory::Create(boundary_type, rho_inf, u_inf, p_inf);
}

auto Simulation::CreateWriter(const std::string& output_format,
                              const std::string& output_dir)
    -> std::unique_ptr<StepWriter> {
    return WriterFactory::Create(output_format, output_dir);
}

void Simulation::ApplyInitialConditions(DataLayer& layer) {
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

        layer.U(i) = eint; // specific internal
        layer.e(i) = Etot; // total density
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

    // Initialize output_rodionov writer
    writer_ = CreateWriter(settings_.output_format, settings_.output_dir);

    if (settings_.analytical) {
        Settings analytical_settings = settings_;
        analytical_settings.solver = "analytical";
        analytical_settings.output_dir += "_analytical";

        analytical_layer_ = std::make_unique<DataLayer>(
            analytical_settings.N,
            analytical_settings.padding,
            analytical_settings.dim);

        ApplyInitialConditions(*analytical_layer_);

        analytical_solver_ = SolverFactory::Create(analytical_settings);

        analytical_writer_ = CreateWriter(analytical_settings.output_format,
                                          analytical_settings.output_dir);
    }

    std::cout << '\n';
    std::cout << "Simulation initialized:" << '\n';
    std::cout << ">>> Solver:           " << settings_.solver << '\n';
    std::cout << ">>> Riemann Solver:   " << settings_.riemann_solver << '\n';
    std::cout << ">>> Reconstruction:   " << settings_.reconstruction << '\n';
    std::cout << ">>> Boundary left:    " << settings_.left_boundary << '\n';
    std::cout << ">>> Boundary right:   " << settings_.right_boundary << '\n';
    std::cout << ">>> Grid size (N):    " << settings_.N << '\n';
    std::cout << ">>> Dimension:        " << settings_.dim << '\n';
    std::cout << ">>> Domain length:    " << settings_.L_x << '\n';
    std::cout << ">>> CFL:              " << settings_.cfl << '\n';
    std::cout << ">>> Gamma:            " << settings_.gamma << '\n';
    std::cout << ">>> End time:         " << settings_.t_end << '\n';
    std::cout << ">>> Output directory: " << settings_.output_dir << '\n';
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

auto Simulation::GetCurrentStep() const -> std::size_t { return step_; }

auto Simulation::GetCurrentTime() const -> double { return t_cur_; }

auto Simulation::ShouldWrite(std::size_t step) const -> bool {
    if (settings_.output_every_steps == 0)
        return false;
    return (step % settings_.output_every_steps == 0 or t_cur_ >= settings_.t_end);
}

void Simulation::WriteInitialState() const {
    std::cout << "\nWriting the initial state..." << '\n';;
    if (writer_) {
        writer_->Write(*layer_, 0, 0.0);
    }
}

void Simulation::WriteStepState(std::size_t step, double t_cur) const {
    if (writer_ && ShouldWrite(step)) {
        writer_->Write(*layer_, step, t_cur);
    }
}

void Simulation::WriteAnalyticalStepState(std::size_t step, double t_cur) const {
    if (analytical_writer_ && analytical_layer_ && ShouldWrite(step)) {
        analytical_writer_->Write(*analytical_layer_, step, t_cur);
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

        if (settings_.solver == "analytical") {
            t_cur_ += dt;
        }

        ++step_;

        // Write output_rodionov if needed
        WriteStepState(step_, t_cur_);

        if (settings_.analytical) {
            // Force this step size
            if (auto* as = dynamic_cast<AnalyticalSolver*>(analytical_solver_.get())) {
                as->SetDt(dt);
            }
            const double dt_a = analytical_solver_->Step(*analytical_layer_, t_cur_);
            (void)dt_a;
            WriteAnalyticalStepState(step_, t_cur_);
        }

        // Progress output_rodionov
        if (log_progress_ && 0 == step_ % settings_.output_every_steps or settings_.t_end
            <= t_cur_) {
            double progress = (t_cur_ / settings_.t_end) * 100.0;
            int percent = static_cast<int>(progress);
            std::cout << '\r' << ">>> [PROGRESS]: Step " << step_ << ", " << percent
                << "% processed, time: " << t_cur_ << " of " << settings_.t_end <<
                std::flush;
            // if (t_cur_ > 0.0) {
            //     double est_total_steps = step_ * (settings_.t_end / t_cur_);
            //     std::cout << ", total steps: " << est_total_steps;
            // } else {
            //     std::cout << ", total steps: unknown";
            // }
        }
    }

    std::cout << '\n';
    std::cout << "Simulation completed!" << '\n';
    std::cout << ">>> Final time:  " << t_cur_ << '\n';
    std::cout << ">>> Total steps: " << step_ << '\n';
}