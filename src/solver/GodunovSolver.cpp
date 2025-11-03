#include "solver/GodunovSolver.hpp"
#include "data/DataLayer.hpp"

GodunovSolver::GodunovSolver()
    : bc_manager_(1) {
}

void GodunovSolver::SetCfl(double value) {
    cfl_ = value;
}

void GodunovSolver::AddBoundary(int axis,
                                std::shared_ptr<BoundaryCondition> min,
                                std::shared_ptr<BoundaryCondition> max) {
    bc_manager_.Set(axis, std::move(min), std::move(max));
}

void GodunovSolver::Step(DataLayer &layer, double &time, double t_end) {
    bc_manager_.ApplyAll(layer);


    double dt = 0.0; // TODO: рассчитать по CFL

    // Вычислить потоки и обновить layer (TODO)


    (void) t_end;
    time += dt;
}
