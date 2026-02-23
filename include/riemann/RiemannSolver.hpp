#ifndef RIEMANNSOLVER_HPP
#define RIEMANNSOLVER_HPP

#include "data/Variables.hpp"
#include "riemann/RiemannHelpers.hpp"

/**
 * @class RiemannSolver
 * @brief Abstract base class for axis-aligned Riemann solvers (Euler, ideal gas).
 *
 * Computes numerical flux at one interface given left/right primitive states
 * and the interface normal axis.
 */
class RiemannSolver {
public:
    virtual ~RiemannSolver() = default;

    /**
     * @brief Compute numerical flux at one interface.
     *
     * @param left  Left primitive state at interface.
     * @param right Right primitive state at interface.
     * @param gamma Ratio of specific heats.
     * @param axis  Interface normal axis.
     * @return FluxCell in conservative ordering (mass, mom_x, mom_y, mom_z, energy).
     */
    [[nodiscard]] virtual auto ComputeFlux(const PrimitiveCell& left,
                                           const PrimitiveCell& right,
                                           double gamma,
                                           Axis axis) const -> FluxCell = 0;
};

#endif  // RIEMANNSOLVER_HPP