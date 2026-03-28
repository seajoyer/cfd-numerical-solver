#ifndef NONREFLECTIVEBOUNDARY_HPP
#define NONREFLECTIVEBOUNDARY_HPP

#include "bc/BoundaryCondition.hpp"

/**
 * @class NonReflectiveBoundary
 * @brief Simple far-field non-reflective placeholder.
 *
 * Current version returns prescribed far-field exterior state.
 * More advanced characteristic treatment can be added later.
 */
class NonReflectiveBoundary final : public BoundaryCondition {
public:
    NonReflectiveBoundary(PrimitiveCell farfield_state, double gamma);

    [[nodiscard]] PrimitiveCell BuildExteriorState(const DataLayer& layer,
                                                   const Mesh& mesh,
                                                   const Face& face,
                                                   const PrimitiveCell& interior_state) const override;

private:
    PrimitiveCell farfield_state_;
    double gamma_ = 1.4;
};

#endif  // NONREFLECTIVEBOUNDARY_HPP
