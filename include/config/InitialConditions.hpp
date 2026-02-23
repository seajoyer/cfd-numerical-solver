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
    double u_L = 0.0;   ///< x-velocity for left state
    double v_L = 0.0;   ///< y-velocity for left state
    double w_L = 0.0;   ///< z-velocity for left state
    double P_L = 1.0;
    
    // Right state
    double rho_R = 0.125;
    double u_R = 0.0;   ///< z-velocity for right state
    double v_R = 0.0;   ///< y-velocity for right state
    double w_R = 0.0;   ///< z-velocity for right state
    double P_R = 0.1;

    // ==================== 2D Quadrant IC ====================
    // Quadrant 1: x >= x0, y >= y0, z < z0
    double rho_Q1 = 0.0, u_Q1 = 0.0, v_Q1 = 0.0, w_Q1 = 0.0, P_Q1 = 0.0;
    // Quadrant 2: x < x0, y >= y0, z < z0
    double rho_Q2 = 0.0, u_Q2 = 0.0, v_Q2 = 0.0, w_Q2 = 0.0, P_Q2 = 0.0;
    // Quadrant 3: x < x0, y < y0, z < z0
    double rho_Q3 = 0.0, u_Q3 = 0.0, v_Q3 = 0.0, w_Q3 = 0.0, P_Q3 = 0.0;
    // Quadrant 4: x >= x0, y < y0, z < z0
    double rho_Q4 = 0.0, u_Q4 = 0.0, v_Q4 = 0.0, w_Q4 = 0.0, P_Q4 = 0.0;
    // Quadrant 5: x >= x0, y >= y0, z >= z0
    double rho_Q5 = 0.0, u_Q5 = 0.0, v_Q5 = 0.0, w_Q5 = 0.0, P_Q5 = 0.0;
    // Quadrant 6: x < x0, y >= y0, z >= z0
    double rho_Q6 = 0.0, u_Q6 = 0.0, v_Q6 = 0.0, w_Q6 = 0.0, P_Q6 = 0.0;
    // Quadrant 7: x < x0, y < y0, z >= z0
    double rho_Q7 = 0.0, u_Q7 = 0.0, v_Q7 = 0.0, w_Q7 = 0.0, P_Q7 = 0.0;
    // Quadrant 8: x >= x0, y < y0, z >= z0
    double rho_Q8 = 0.0, u_Q8 = 0.0, v_Q8 = 0.0, w_Q8 = 0.0, P_Q8 = 0.0;

    /**
     * @brief Initial condition type for 3D problems.
     *
     * Values:
     *   "x_riemann" - discontinuity along x (default)
     *   "y_riemann" - discontinuity along y
     *   "z_riemann" - discontinuity along z
     *   "quadrant"  - four-quadrant 2D Riemann problem
     *   "octant"    - eight-octant 3D Riemann problem
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

    /**
     * @brief Returns true if this is a Octant-type IC
     */
    [[nodiscard]] auto IsOctant() const -> bool {
        return ic_type == "octant";
    }
};

#endif  // INITIALCONDITIONS_HPP
