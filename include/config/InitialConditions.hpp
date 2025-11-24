#ifndef INITIALCONDITIONS_HPP
#define INITIALCONDITIONS_HPP

/**
 * @file InitialConditions.hpp
 * @brief Data structure for Riemann problem initial conditions
 */

/**
 * @struct InitialConditions
 * @brief Defines left and right states for a 1D Riemann problem
 * 
 * This structure contains the primitive variables (density, velocity, pressure)
 * for both the left and right states of a Riemann problem, along with the
 * position of the discontinuity.
 * 
 * The Riemann problem consists of two constant states separated by a
 * discontinuity at position x0:
 * 
 * - Left state (x < x0): (ρ_L, u_L, P_L)
 * - Right state (x ≥ x0): (ρ_R, u_R, P_R)
 * 
 * Initial conditions are typically loaded from the YAML configuration file
 * by ConfigParser and can represent various flow scenarios:
 * - Shock tubes (Sod problems)
 * - Blast waves
 * - Shock-shock interactions
 * - Rarefaction waves
 * 
 * @note All quantities are in SI units or dimensionless
 * @see Settings, ConfigParser
 * 
 * Example YAML configuration:
 * @code{.yaml}
 * initial_conditions:
 *   sod1:
 *     rho_L: 1.0
 *     u_L: 0.0
 *     P_L: 1.0
 *     rho_R: 0.125
 *     u_R: 0.0
 *     P_R: 0.1
 *     x0: 0.5
 *     t_end: 0.25
 * @endcode
 */
struct InitialConditions {
    /** @brief Left state density (ρ_L)
     * @note Must be positive for physical validity
     */
    double rho_L;
    
    /** @brief Left state velocity (u_L) */
    double u_L;
    
    /** @brief Left state pressure (P_L)
     * @note Must be positive for physical validity
     */
    double P_L;
    
    /** @brief Right state density (ρ_R)
     * @note Must be positive for physical validity
     */
    double rho_R;
    
    /** @brief Right state velocity (u_R) */
    double u_R;
    
    /** @brief Right state pressure (P_R)
     * @note Must be positive for physical validity
     */
    double P_R;
    
    /** @brief Position of the discontinuity
     * @note Typically in range [0, L_x]
     * @note Default value: 0.5 for centered discontinuity in unit domain
     */
    double x0;
};

#endif  // INITIALCONDITIONS_HPP
