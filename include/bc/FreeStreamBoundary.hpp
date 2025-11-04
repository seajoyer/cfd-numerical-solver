#ifndef FREESTREAMBOUNDARY_HPP
#define FREESTREAMBOUNDARY_HPP

#include "bc/BoundaryCondition.hpp"

/**
 * @class FreeStreamBoundary
 * @brief Implements free-stream (open flow) boundary conditions.
 *
 * This boundary condition represents the interface between the computational
 * domain and an external undisturbed flow region. It automatically distinguishes
 * between inflow and outflow, applying either external (freestream) values
 * or a zero-gradient condition accordingly.
 *
 * Behavior in 1D:
 *  - If the flow enters the domain (inflow): ghost cells are filled with
 *    prescribed freestream values (ρ∞, u∞, P∞).
 *  - If the flow exits the domain (outflow): ghost cells copy values from
 *    adjacent interior cells (zero-gradient).
 *
 * Physically, this condition minimizes artificial reflections at boundaries
 * and is suitable for external flow simulations (e.g., airfoil, jet outflow).
 *
 *
 * @see InletBoundary
 * @see OutletBoundary
 * @see NonReflectingBoundary
 */
class FreeStreamBoundary : public BoundaryCondition {
   public:
    /**
     * @brief Constructs a free-stream boundary with prescribed external flow parameters.
     *
     * Initializes the boundary with freestream values that are applied in the ghost
     * cells whenever the flow direction indicates inflow. When the flow leaves the
     * domain, these values are not enforced, maintaining natural outflow.
     *
     * @param rhoInf External density (ρ∞) in the freestream.
     * @param uInf   External velocity (u∞) in the freestream.
     * @param pInf   External pressure (P∞) in the freestream.
     */
    FreeStreamBoundary(double rho_inf, double u_inf, double p_inf);

    /**
     * @brief Applies free-stream boundary condition along the specified axis.
     *
     * Detects local flow direction and applies either:
     *  - Inflow: sets ghost cells to external (freestream) values.
     *  - Outflow: applies zero-gradient copying from inner cells.
     *
     * @param layer Reference to the DataLayer being modified.
     * @param axis Axis index (0 for 1D problems).
     * @param side Which boundary to process (Side::Min or Side::Max).
     *
     * @note For 1D implementation, only axis = 0 is used. The same logic
     *       generalizes easily to 2D and 3D.
     */
    void Apply(DataLayer& layer, int axis, Side side) const override;

   private:
    double rho_inf_;
    double u_inf_;
    double p_inf_;
};

#endif  // FREESTREAMBOUNDARY_HPP
