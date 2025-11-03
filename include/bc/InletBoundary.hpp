#ifndef INLETBOUNDARY_HPP
#define INLETBOUNDARY_HPP

#include "bc/BoundaryCondition.hpp"

/**
 * @class InletBoundary
 * @brief Implements inlet (inflow) boundary conditions.
 *
 * This boundary condition specifies fixed external (freestream) values
 * for density, velocity, and pressure on the side where the flow enters
 * the computational domain.
 *
 * The direction of flow is determined automatically:
 *  - On the left boundary (Side::Min), inflow occurs when u > 0.
 *  - On the right boundary (Side::Max), inflow occurs when u < 0.
 *
 * If the flow is directed *into* the domain, the specified inlet
 * state (ρ, u, P) is applied in the ghost cells.
 * Otherwise, the boundary behaves as an outlet (zero-gradient condition).
 *
 *
 * @see OutletBoundary
 * @see FreeStreamBoundary
 */
class InletBoundary : public BoundaryCondition {
public:
    /**
     * @brief Constructs an inlet boundary with prescribed external conditions.
     *
     * Initializes the boundary with constant inflow values for density, velocity,
     * and pressure that will be imposed on ghost cells whenever the flow direction
     * indicates inflow into the domain.
     *
     * @param rhoIn  External density (ρ) at the inlet.
     * @param uIn    External velocity (u) at the inlet.
     * @param pIn    External pressure (P) at the inlet.
     */
    InletBoundary(double rho_in, double u_in, double p_in)
        : rho_in_(rho_in), u_in_(u_in), p_in_(p_in) {
    }

    /**
     * @brief Applies inlet boundary condition along a specific axis.
     *
     * Updates ghost cells with user-defined external values when the flow
     * direction indicates inflow; otherwise applies zero-gradient copying.
     *
     * @param layer Reference to the DataLayer being modified.
     * @param axis Axis index (0 for 1D).
     * @param side Which side to apply: Side::Min (left) or Side::Max (right).
     *
     * @note For 1D problems, only axis = 0 is used. In higher dimensions,
     *       the same logic would be applied per direction.
     */
    void Apply(DataLayer &layer, int axis, Side side) const override;

private:
    double rho_in_;
    double u_in_;
    double p_in_;
};

#endif // INLETBOUNDARY_HPP
