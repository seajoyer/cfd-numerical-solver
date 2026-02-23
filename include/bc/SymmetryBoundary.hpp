#ifndef SYMMETRYBOUNDARY_HPP
#define SYMMETRYBOUNDARY_HPP

#include "bc/BoundaryCondition.hpp"
#include "data/Variables.hpp"

/**
 * @class SymmetryBoundary
 * @brief Symmetry (mirror) boundary for conservative Euler state U.
 *
 * Ghost cells are filled by mirroring interior cells across the boundary.
 * The normal momentum component is inverted; tangential components and scalars are copied.
 *
 * For axis:
 *  - X: rhoU changes sign
 *  - Y: rhoV changes sign
 *  - Z: rhoW changes sign
 */
class SymmetryBoundary final : public BoundaryCondition {
public:
    void Apply(DataLayer& layer, Axis axis, Side side) const override;
};

#endif  // SYMMETRYBOUNDARY_HPP