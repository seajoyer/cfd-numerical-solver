#ifndef RIEMANNSOLVER_HPP
#define RIEMANNSOLVER_HPP

#include "data/Variables.hpp"

/**
 * @class RiemannSolver
 * @brief Abstract numerical flux solver for one face of a generic mesh.
 *
 * Computes conservative numerical flux through one face from left/right
 * primitive states and the face unit normal.
 *
 * State orientation:
 * - left  : owner-side state
 * - right : neighbor-side or exterior boundary state
 *
 * Returned flux orientation:
 * - Flux is directed along the given face normal.
 * - For internal faces, mesh face normal is oriented from owner to neighbor.
 * - For boundary faces, mesh face normal is outward from owner.
 */
class RiemannSolver {
public:
    virtual ~RiemannSolver() = default;

    /**
     * @brief Compute numerical flux through one face.
     *
     * @param left Owner-side primitive state.
     * @param right Neighbor-side or exterior primitive state.
     * @param gamma Ratio of specific heats.
     * @param normal Face unit normal.
     * @return Conservative numerical flux through the face.
     */
    [[nodiscard]] virtual ConservativeCell ComputeFlux(const PrimitiveCell& left,
                                                       const PrimitiveCell& right,
                                                       double gamma,
                                                       const FaceNormal& normal) const = 0;
};

#endif  // RIEMANNSOLVER_HPP
