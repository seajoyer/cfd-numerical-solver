#ifndef FREESTREAMBOUNDARY_HPP
#define FREESTREAMBOUNDARY_HPP

#include "bc/BoundaryCondition.hpp"
#include "bc/BoundaryFactory.hpp"
#include "data/Variables.hpp"

/**
 * @class FreeStreamBoundary
 * @brief Conditional freestream boundary: inflow uses prescribed state, outflow uses zero-gradient.
 *
 * The inflow/outflow decision is made from the local normal velocity component
 * at the nearest interior core cell (computed from conservative state U).
 *
 * Works only with conservative state U(var,i,j,k).
 */
class FreeStreamBoundary final : public BoundaryCondition {
public:
    /**
     * @brief Constructs boundary condition with a prescribed freestream conservative state.
     * @param freestream_U Conservative state imposed during inflow.
     */
    explicit FreeStreamBoundary(const FarfieldConservative& freestream_U);

    void Apply(DataLayer& layer, const Mesh& mesh, Axis axis, Side side) const override;

private:
    FarfieldConservative freestream_U_;
};

#endif  // FREESTREAMBOUNDARY_HPP