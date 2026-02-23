#ifndef RECONSTRUCTION_HPP
#define RECONSTRUCTION_HPP

#include "data/Variables.hpp"  // Axis, PrimitiveCell

/**
 * @class Reconstruction
 * @brief Interface for axis-aligned reconstruction schemes on cell-centered primitives W.
 *
 * Input:
 *  - W(var,i,j,k) in Workspace, var=(rho,u,v,w,P), shape (5,sx,sy,sz).
 *
 * Reconstruction is performed per face (interface) and returns:
 *  - WL: state on the left side of the face
 *  - WR: state on the right side of the face
 *
 * Face indexing convention:
 *  - Provide the index of the LEFT cell adjacent to the face.
 *  - The RIGHT cell is at +1 along the chosen axis.
 *
 * Contract:
 *  - No allocations in hot path.
 *  - Must not modify W.
 */
class Reconstruction {
public:
    virtual ~Reconstruction() = default;

    virtual void ReconstructFace(const xt::xtensor<double, 4>& W,
                                 Axis axis,
                                 int i, int j, int k,
                                 PrimitiveCell& WL,
                                 PrimitiveCell& WR) const = 0;
};

#endif  // RECONSTRUCTION_HPP