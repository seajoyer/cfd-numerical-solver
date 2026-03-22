#ifndef NONREFLECTIVEBOUNDARY_HPP
#define NONREFLECTIVEBOUNDARY_HPP

#include "bc/BoundaryCondition.hpp"
#include "bc/BoundaryFactory.hpp"

/**
 * @class NonReflectiveBoundary
 * @brief Approximate characteristic open boundary for Euler equations.
 *
 * Uses a simple 1D characteristic treatment along the boundary normal:
 * - supersonic outflow: extrapolate interior state
 * - supersonic inflow: impose far-field state
 * - subsonic: combine outgoing characteristic from interior with incoming
 *   characteristic from far-field
 *
 * Tangential velocities and entropy are taken:
 * - from interior for outflow
 * - from far-field for inflow
 */
class NonReflectiveBoundary final : public BoundaryCondition {
public:
    NonReflectiveBoundary(const FarfieldConservative& farfield_U, double gamma);

    void Apply(DataLayer& layer, const Mesh& mesh, Axis axis, Side side) const override;

private:
    FarfieldConservative farfield_U_;
    double gamma_;
};

#endif  // NONREFLECTIVEBOUNDARY_HPP