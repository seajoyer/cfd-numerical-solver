#ifndef INITIALCONDITIONS_HPP
#define INITIALCONDITIONS_HPP

#include "config/Settings.hpp"

#include <string>

/**
 * @file InitialConditions.hpp
 * @brief Data structure for Riemann problem initial conditions (1D and 2D)
 *
 * For 1D: standard left/right Riemann problem along x.
 * For 2D: supports several IC types:
 *   - "x_riemann":  discontinuity along x (vertical interface)
 *   - "y_riemann":  discontinuity along y (horizontal interface)
 *   - "quadrant":   four-quadrant 2D Riemann problem (Lax-Liu type)
 *
 * Quadrant layout (when ic_type == "quadrant"):
 *   Q2 (x<x0, y>=y0) | Q1 (x>=x0, y>=y0)
 *   -----------------+------------------
 *   Q3 (x<x0, y<y0)  | Q4 (x>=x0, y<y0)
 *
 * For non-quadrant types, only left/right states are used.
 */
struct InitialConditions {
    // Left state (x < x0 for x_riemann; y < y0 for y_riemann)
    double rho_L = 1.0;
    double u_L = 0.0;
    double v_L = 0.0;   ///< y-velocity for left state (2D)
    double P_L = 1.0;
    
    // Right state
    double rho_R = 0.125;
    double u_R = 0.0;
    double v_R = 0.0;   ///< y-velocity for right state (2D)
    double P_R = 0.1;

    // ==================== 2D Quadrant IC ====================
    // Quadrant 1: x >= x0, y >= y0 (top-right)
    double rho_Q1 = 0.0, u_Q1 = 0.0, v_Q1 = 0.0, P_Q1 = 0.0;
    // Quadrant 2: x < x0, y >= y0 (top-left)
    double rho_Q2 = 0.0, u_Q2 = 0.0, v_Q2 = 0.0, P_Q2 = 0.0;
    // Quadrant 3: x < x0, y < y0 (bottom-left)
    double rho_Q3 = 0.0, u_Q3 = 0.0, v_Q3 = 0.0, P_Q3 = 0.0;
    // Quadrant 4: x >= x0, y < y0 (bottom-right)
    double rho_Q4 = 0.0, u_Q4 = 0.0, v_Q4 = 0.0, P_Q4 = 0.0;

    /**
     * @brief Initial condition type for 2D problems.
     *
     * Values:
     *   "x_riemann" - discontinuity along x (default, also used for 1D)
     *   "y_riemann" - discontinuity along y
     *   "quadrant"  - four-quadrant 2D Riemann problem
     */
    std::string ic_type = "x_riemann";

    /** @brief Case-specific setting overrides */
    CaseSettings overrides;

    /**
     * @brief Returns true if this is a quadrant-type IC
     */
    [[nodiscard]] auto IsQuadrant() const -> bool {
        return ic_type == "quadrant";
    }
};

#endif  // INITIALCONDITIONS_HPP
