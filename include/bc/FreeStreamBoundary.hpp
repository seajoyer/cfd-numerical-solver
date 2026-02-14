#ifndef FREESTREAMBOUNDARY_HPP
#define FREESTREAMBOUNDARY_HPP

#include "bc/BoundaryCondition.hpp"

/**
 * @class FreeStreamBoundary
 * @brief Far-field boundary condition for 1D and 2D.
 *
 * Supports both 1D and 2D data layouts.
 * For 2D, fills ghost columns (axis=0) or ghost rows (axis=1).
 */
class FreeStreamBoundary : public BoundaryCondition {
public:
    FreeStreamBoundary(double rho_inf, double u_inf, double p_inf);
    FreeStreamBoundary(double rho_inf, double u_inf, double v_inf, double p_inf);

    void Apply(DataLayer& layer, int axis, Side side) const override;

private:
    double rho_inf_;
    double u_inf_;
    double p_inf_;
    double v_inf_;

    void Apply2D(DataLayer& layer, int axis, Side side) const;
};

#endif  // FREESTREAMBOUNDARY_HPP
