#ifndef EOS_HPP
#define EOS_HPP

#include "data/Variables.hpp"

/**
 * @class EOS
 * @brief Ideal gas equation of state utilities (1D/2D).
 *
 * Updated for 2D: kinetic energy includes both u and v components.
 */
class EOS {
   public:
    static auto PrimToCons(const Primitive& w, double gamma) -> Conservative;
    static auto ConsToPrim(const Conservative& u, double gamma) -> Primitive;
    static auto Pressure(const Conservative& u, double gamma) -> double;
    static auto SoundSpeed(const Primitive& w, double gamma) -> double;
};

#endif  // EOS_HPP
