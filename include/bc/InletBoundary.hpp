#ifndef INLETBOUNDARY_HPP
#define INLETBOUNDARY_HPP

#include "bc/BoundaryCondition.hpp"
#include "bc/BoundaryFactory.hpp"
#include "data/Variables.hpp"

/**
 * @class InletBoundary
 * @brief Conditional inlet boundary condition with fixed external conservative state.
 *
 * If the flow at the boundary is directed into the domain (based on the normal
 * velocity component at the nearest interior core cell), ghost cells are set
 * to the prescribed inflow conservative state (rho, rhoU, rhoV, rhoW, E).
 * Otherwise, behaves as an outlet (zero-gradient): copies nearest interior layer.
 *
 * Works only with conservative state U(var,i,j,k).
 */
class InletBoundary final : public BoundaryCondition {
public:
    /**
     * @brief Constructs inlet boundary with prescribed inflow conservative state.
     * @param inflow_U Conservative state to impose during inflow.
     */
    explicit InletBoundary(const FarfieldConservative& inflow_U);

    void Apply(DataLayer& layer, Axis axis, Side side) const override;

private:
    FarfieldConservative inflow_U_;
};

#endif  // INLETBOUNDARY_HPP