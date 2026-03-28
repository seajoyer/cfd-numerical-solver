#ifndef SYMMETRYBOUNDARY_HPP
#define SYMMETRYBOUNDARY_HPP

#include "bc/BoundaryCondition.hpp"

/**
 * @class SymmetryBoundary
 * @brief Symmetry boundary.
 *
 * Equivalent to reflective slip boundary for Euler equations.
 */
class SymmetryBoundary final : public BoundaryCondition {
public:
    [[nodiscard]] PrimitiveCell BuildExteriorState(const DataLayer& layer,
                                                   const Mesh& mesh,
                                                   const Face& face,
                                                   const PrimitiveCell& interior_state) const override;
};

#endif  // SYMMETRYBOUNDARY_HPP
