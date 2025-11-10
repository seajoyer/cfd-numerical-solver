#ifndef HLLCRIEMANNSOLVER_HPP
#define HLLCRIEMANNSOLVER_HPP

#include "RiemannSolver.hpp"
#include "solver/EOS.hpp"

/**
 * @class HLLCRiemannSolver
 * @brief HLLC (HLL-Contact) approximate Riemann solver.
 *
 * Extends HLL by restoring the contact wave, providing better resolution
 * of contact discontinuities and shear layers while remaining robust.
 */
class HLLCRiemannSolver : public RiemannSolver {
public:
    /**
     * @brief Default constructor.
     */
    HLLCRiemannSolver() {}

    /**
     * @brief Computes the HLLC numerical flux.
     *
     * @param left Left primitive state.
     * @param right Right primitive state.
     * @param gamma Ratio of specific heats.
     * @return HLLC flux at the interface.
     */
    Flux ComputeFlux(const Primitive &left,
                     const Primitive &right,
                     double gamma) const override;
};

#endif  // HLLCRIEMANNSOLVER_HPP
