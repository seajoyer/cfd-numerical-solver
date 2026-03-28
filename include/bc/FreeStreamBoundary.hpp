#ifndef FREESTREAMBOUNDARY_HPP
#define FREESTREAMBOUNDARY_HPP

#include "bc/BoundaryCondition.hpp"

/**
 * @class FreeStreamBoundary
 * @brief Prescribed far-field boundary state.
 */
class FreeStreamBoundary final : public BoundaryCondition {
public:
    explicit FreeStreamBoundary(PrimitiveCell farfield_state);

    [[nodiscard]] PrimitiveCell BuildExteriorState(const DataLayer& layer,
                                                   const Mesh& mesh,
                                                   const Face& face,
                                                   const PrimitiveCell& interior_state) const override;

private:
    PrimitiveCell farfield_state_;
};

#endif  // FREESTREAMBOUNDARY_HPP
