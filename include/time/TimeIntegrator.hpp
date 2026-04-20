#ifndef TIMEINTEGRATOR_HPP
#define TIMEINTEGRATOR_HPP

class DataLayer;
class Mesh;
class Workspace;
class SpatialOperator;
class StateSynchronizer;

/**
 * @class TimeIntegrator
 * @brief Abstract explicit time integration scheme for semi-discrete FV systems.
 *
 * Integrators advance conservative cell-centered state U stored in DataLayer.
 * Boundary conditions and face flux construction are handled inside SpatialOperator.
 */
class TimeIntegrator {
public:
    virtual ~TimeIntegrator() = default;

    /**
     * @brief Advance conservative solution by one time step.
     *
     * @param layer Conservative state storage updated in-place.
     * @param mesh Mesh with cells, faces, geometry, and connectivity.
     * @param workspace Reusable scratch buffers.
     * @param dt Time step size.
     * @param gamma Ratio of specific heats.
     * @param op Spatial operator providing RHS evaluations.
     * @param halo_exchange Optional halo exchange used between RK stages.
     */
    virtual void Advance(DataLayer& layer,
                         const Mesh& mesh,
                         Workspace& workspace,
                         double dt,
                         double gamma,
                         const SpatialOperator& op,
                         const StateSynchronizer* halo_exchange) const = 0;

    /**
     * @brief Set positivity floors used after conservative update.
     *
     * @param rho_min Density floor.
     * @param p_min Pressure floor.
     */
    virtual void SetPositivityThresholds(double rho_min, double p_min) {
        rho_min_ = rho_min;
        p_min_ = p_min;
    }

protected:
    double rho_min_ = 1e-10;
    double p_min_ = 1e-10;
};

#endif  // TIMEINTEGRATOR_HPP
