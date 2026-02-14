#ifndef INLETBOUNDARY_HPP
#define INLETBOUNDARY_HPP

#include "bc/BoundaryCondition.hpp"

/**
 * @class InletBoundary
 * @brief Implements inlet (inflow) boundary conditions for 1D and 2D.
 *
 * This boundary condition specifies fixed external (freestream) values
 * for density, velocity, and pressure on the side where the flow enters
 * the computational domain.
 *
 * The direction of flow is determined automatically:
 *  - On the left boundary (Side::kLeft), inflow occurs when normal velocity > 0.
 *  - On the right boundary (Side::kRight), inflow occurs when normal velocity < 0.
 *
 * If the flow is directed *into* the domain, the specified inlet
 * state is applied in the ghost cells.
 * Otherwise, the boundary behaves as an outlet (zero-gradient condition).
 *
 * In 2D:
 *  - For x-axis boundaries: normal velocity is u, tangential is v.
 *  - For y-axis boundaries: normal velocity is v, tangential is u.
 *
 * @see OutletBoundary
 * @see FreeStreamBoundary
 */
class InletBoundary : public BoundaryCondition {
   public:
    /**
     * @brief Constructs a 1D inlet boundary with prescribed external conditions.
     *
     * @param rho_in  External density at the inlet.
     * @param u_in    External velocity (u) at the inlet.
     * @param p_in    External pressure (P) at the inlet.
     */
    InletBoundary(double rho_in, double u_in, double p_in);

    /**
     * @brief Constructs a 2D inlet boundary with prescribed external conditions.
     *
     * @param rho_in  External density at the inlet.
     * @param u_in    External x-velocity at the inlet.
     * @param v_in    External y-velocity at the inlet.
     * @param p_in    External pressure at the inlet.
     */
    InletBoundary(double rho_in, double u_in, double v_in, double p_in);

    /**
     * @brief Applies inlet boundary condition along a specific axis.
     *
     * For 2D, checks the normal velocity component for the given axis
     * to determine inflow/outflow. On inflow: applies prescribed values.
     * On outflow: zero-gradient (copies from nearest core cell).
     *
     * @param layer Reference to the DataLayer being modified.
     * @param axis Axis index (0=x, 1=y).
     * @param side Which side to apply: Side::kLeft or Side::kRight.
     */
    void Apply(DataLayer& layer, int axis, Side side) const override;

   private:
    double rho_in_;
    double u_in_;
    double v_in_;
    double p_in_;
};

#endif  // INLETBOUNDARY_HPP
