#ifndef P0RECONSTRUCTION_HPP
#define P0RECONSTRUCTION_HPP

#include <vector>
#include "solver/Variables.hpp"
#include "Reconstruction.hpp"


/**
 * @class P0Reconstruction
 * @brief Piecewise-constant (first-order) reconstruction scheme.
 *
 * Implements a zero-order reconstruction where the interface states
 * are taken directly from the adjacent cell-centered values without
 * any slope or limiter. This corresponds to the classical first-order
 * Godunov method.
 */
class P0Reconstruction : public Reconstruction {
public:
    /**
     * @brief Default constructor.
     */
    P0Reconstruction() = default;

    /**
     * @brief Performs piecewise-constant reconstruction.
     *
     * For each interface i, assigns:
     * - leftStates[i]  from the cell on the left side,
     * - rightStates[i] from the cell on the right side.
     *
     * @param cellStates Cell-centered primitive variables.
     * @param leftStates Output array of left states at interfaces.
     * @param rightStates Output array of right states at interfaces.
     */
    void Reconstruct(const std::vector<Primitive> &cell_states,
                     std::vector<Primitive> &left_states,
                     std::vector<Primitive> &right_states) const override;
};

#endif  // P0RECONSTRUCTION_HPP
