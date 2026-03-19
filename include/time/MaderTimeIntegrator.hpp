#ifndef MADERTIMEINTEGRATOR_HPP
#define MADERTIMEINTEGRATOR_HPP

#include "time/TimeIntegrator.hpp"

class DataLayer;
class Mesh;
class Workspace;
class SpatialOperator;

/**
 * @class MaderTimeIntegrator
 * @brief Time integrator for phase-split Mader 2DE scheme.
 *
 * One full time step consists of 5 consecutive phases executed in fixed order.
 */
class MaderTimeIntegrator final : public TimeIntegrator {
public:
    MaderTimeIntegrator() = default;

    void Advance(DataLayer& layer,
                 const Mesh& mesh,
                 Workspace& workspace,
                 double dt,
                 double gamma,
                 const SpatialOperator& op) const override;
};

#endif  // MADERTIMEINTEGRATOR_HPP
