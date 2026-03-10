#ifndef ARTIFICIALVISCOSITY_HPP
#define ARTIFICIALVISCOSITY_HPP

#include "data/Variables.hpp"

/**
 * @class ArtificialViscosity
 * @brief Interface for artificial viscosity models (2D/3D-ready).
 *
 * Artificial viscosity adds a dissipative contribution to the conservative RHS:
 *   rhs += L_visc(U)
 *
 * SpatialOperator is responsible for:
 *  - halo + physical BC on conservative U
 *  - converting U -> W into Workspace.W()
 *
 * Implementations should:
 *  - read primitive field W(var,i,j,k) = (rho,u,v,w,P)
 *  - use Mesh metrics (inv_dx, inv_dy, inv_dz)
 *  - avoid allocations in hot path (cache buffers internally if needed)
 */
class DataLayer;
class Mesh;

class ArtificialViscosity {
public:
    virtual ~ArtificialViscosity() = default;

    /**
     * @brief Add viscosity contribution to RHS (in-place accumulation).
     *
     * @param layer Conservative state.
     * @param mesh  Structured mesh with geometry, metrics, and cell types.
     * @param W     Primitive field W(var,i,j,k) = (rho,u,v,w,P), including ghosts.
     * @param gamma Ratio of specific heats.
     * @param dt    Local time step (optional for some models).
     * @param rhs   RHS buffer in conservative ordering (accumulated in-place).
     */
    virtual void AddToRhs(const DataLayer& layer,
                          const Mesh& mesh,
                          const xt::xtensor<double, 4>& W,
                          double gamma,
                          double dt,
                          xt::xtensor<double, 4>& rhs) const = 0;

    /** @brief Minimal required padding (ghost cells) for this viscosity stencil. */
    [[nodiscard]] virtual int GetRequiredPadding() const { return 1; }
};

#endif  // ARTIFICIALVISCOSITY_HPP