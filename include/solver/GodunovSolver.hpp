#ifndef GODUNOVSOLVER_HPP
#define GODUNOVSOLVER_HPP

#include <cstddef>
#include <memory>
#include "solver/Solver.hpp"
#include "bc/BoundaryManager.hpp"


/**
 * @class GodunovSolver
 * @brief First-order Godunov finite-volume solver for hyperbolic conservation laws.
 *
 * Implements the classical Godunov method for solving 1D Euler equations
 * using piecewise constant reconstruction and an approximate Riemann solver.
 *
 * The class inherits from Solver and defines specific details of the
 * numerical scheme while preserving the interface contract.
 *
 * @note This solver is currently implemented for 1D problems,
 *       but the design allows extension to 2D/3D schemes by adding
 *       loops over spatial dimensions and generalizing flux computation.
 */
class GodunovSolver : public Solver {
public:
    GodunovSolver();

    explicit GodunovSolver(int dim);

    auto Step(DataLayer &layer, double &t_cur) -> double override;

    void SetCfl(double cfl) override;

    void AddBoundary(int axis,
                     std::shared_ptr<BoundaryCondition> left,
                     std::shared_ptr<BoundaryCondition> right) override;

private:
    double cfl_ = 0.5;
    BoundaryManager bc_manager_;
};

#endif // GODUNOVSOLVER_HPP
