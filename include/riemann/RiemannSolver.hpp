#ifndef RIEMANNSOLVER_HPP
#define RIEMANNSOLVER_HPP

#include "../data/Variables.hpp"

/**
 * @class RiemannSolver
 * @brief Abstract base class for 1D Riemann solvers for the Euler equations.
 *
 * Provides a common interface for computing numerical fluxes at cell interfaces
 * given left and right primitive states. Concrete implementations correspond
 * to different approximate or exact Riemann solvers.
 *
 * Typical usage:
 * - Reconstruction provides left/right Primitive states at each interface.
 * - Solver calls ComputeFlux() for each interface.
 * - Returned Flux is used in the finite-volume update of conservative variables.
 */
class RiemannSolver {
   public:
    /**
     * @brief Virtual destructor for safe polymorphic deletion.
     */
    virtual ~RiemannSolver() = default;

    /**
     * @brief Computes the numerical flux at an interface.
     *
     * @param left Left primitive state at the interface.
     * @param right Right primitive state at the interface.
     * @param gamma Ratio of specific heats.
     * @return Numerical flux for the Euler equations.
     */
    [[nodiscard]] virtual auto ComputeFlux(const Primitive& left, const Primitive& right,
                                           double gamma) const -> Flux = 0;
};

#endif  // RIEMANNSOLVER_HPP
