#ifndef EOS_HPP
#define EOS_HPP

#include "data/Variables.hpp"

/**
 * @class EOS
 * @brief Ideal gas equation of state utilities.
 *
 * Provides static helper functions to convert between primitive and
 * conservative variables and to compute thermodynamic quantities for
 * a calorically perfect ideal gas with constant heat capacity ratio gamma.
 *
 * All functions are stateless: the value of gamma is passed explicitly.
 */
class EOS {
public:
    /**
     * @brief Converts primitive to conservative variables.
     *
     * @param w Primitive variables (rho, u, P).
     * @param gamma Ratio of specific heats.
     * @return Conservative variables (rho, rhoU, E).
     */
    static Conservative PrimToCons(const Primitive &w, double gamma);

    /**
     * @brief Converts conservative to primitive variables.
     *
     * @param u Conservative variables (rho, rhoU, E).
     * @param gamma Ratio of specific heats.
     * @return Primitive variables (rho, u, P).
     */
    static Primitive ConsToPrim(const Conservative &u, double gamma);

    /**
     * @brief Computes thermodynamic pressure from conservative variables.
     *
     * Uses the ideal gas relation:
     * \f$ P = (\gamma - 1)\left(E - \tfrac{1}{2}\rho u^2\right) \f$.
     *
     * @param u Conservative variables (rho, rhoU, E).
     * @param gamma Ratio of specific heats.
     * @return Pressure P.
     */
    static double Pressure(const Conservative &u, double gamma);

    /**
     * @brief Computes the speed of sound from primitive variables.
     *
     * Uses the ideal gas relation:
     * \f$ a = \sqrt{\gamma P / \rho} \f$.
     *
     * @param w Primitive variables (rho, u, P).
     * @param gamma Ratio of specific heats.
     * @return Speed of sound.
     */
    static double SoundSpeed(const Primitive &w, double gamma);
};

#endif  // EOS_HPP
