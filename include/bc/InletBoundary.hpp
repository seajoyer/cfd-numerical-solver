#ifndef INLETBOUNDARY_HPP
#define INLETBOUNDARY_HPP

#include "bc/BoundaryCondition.hpp"

/**
 * @class InletBoundary
 * @brief Prescribed primitive inflow state.
 *
 * Current implementation always returns prescribed exterior state.
 * More advanced inflow/outflow switching can be added later if needed.
 */
class InletBoundary final : public BoundaryCondition {
public:
    explicit InletBoundary(PrimitiveCell inflow_state);

    [[nodiscard]] PrimitiveCell BuildExteriorState(const DataLayer& layer,
                                                   const Mesh& mesh,
                                                   const Face& face,
                                                   const PrimitiveCell& interior_state) const override;

private:
    PrimitiveCell inflow_state_;
};

#endif  // INLETBOUNDARY_HPP
