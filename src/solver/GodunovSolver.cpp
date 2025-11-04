#include "solver/GodunovSolver.hpp"
#include <cstddef>
#include "data/DataLayer.hpp"

GodunovSolver::GodunovSolver() : bc_manager_(1) {}

GodunovSolver::GodunovSolver(int dim) : bc_manager_(dim) {}

void GodunovSolver::SetCfl(double cfl) { cfl_ = cfl; }

void GodunovSolver::AddBoundary(int axis,
                                std::shared_ptr<BoundaryCondition> left,
                                std::shared_ptr<BoundaryCondition> right) {
    bc_manager_.Set(axis, std::move(left), std::move(right));
}

// TODO: implement support of multiple dimensions

auto GodunovSolver::Step(DataLayer &layer, double &t_cur) -> double {
    bc_manager_.ApplyAll(layer);

    double dt = 0.0; // TODO(dmitry): рассчитать по CFL, на каждом шаге или 1 раз при инициализации солвера, хз

    // Вычислить потоки и обновить layer (TODO)

    return dt;
}
