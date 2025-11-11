#ifndef RECONSTRUCTION_HPP
#define RECONSTRUCTION_HPP

#include <vector>
#include "../solver/Variables.hpp"

/**
 * @class Reconstruction
 * @brief Abstract base class for spatial reconstruction schemes.
 *
 * Provides an interface for reconstructing left and right primitive
 * states at cell interfaces from cell-centered primitive values.
 *
 * The reconstruction is performed in 1D along a given indexing, and
 * is typically used before evaluating Riemann fluxes.
 */
class Reconstruction {
public:
    /**
     * @brief Virtual destructor for safe polymorphic deletion.
     */
    virtual ~Reconstruction() = default;

    /**
     * @brief Reconstructs left/right states at cell interfaces.
     *
     * Given an array of cell-centered primitive states, this method
     * computes the corresponding left and right states at each
     * interface between cells.
     *
     * The expected sizes:
     * - cellStates.size() = number of cells (including ghosts if used).
     * - leftStates.size() and rightStates.size() must be preallocated
     *   by the caller to the number of interfaces to be filled.
     *
     * @param cellStates Cell-centered primitive variables.
     * @param leftStates Output array of left states at interfaces.
     * @param rightStates Output array of right states at interfaces.
     */
    virtual void Reconstruct(const std::vector<Primitive> &cell_states,
                             std::vector<Primitive> &left_states,
                             std::vector<Primitive> &right_states) const = 0;
};

#endif  // RECONSTRUCTION_HPP
