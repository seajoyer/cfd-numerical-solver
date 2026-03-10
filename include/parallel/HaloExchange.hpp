#ifndef HALOEXCHANGE_HPP
#define HALOEXCHANGE_HPP

#include <vector>

#include "data/DataLayer.hpp"
#include "data/Variables.hpp"
#include "bc/BoundaryCondition.hpp"
#include "data/Mesh.hpp"
#include "parallel/MPIContext.hpp"

/**
 * @class HaloExchange
 * @brief MPI halo exchange for conservative state U on structured Cartesian subdomains.
 *
 * Exchange is performed for ghost layers of width ng along each active axis.
 * Only conservative state U(var,i,j,k) is exchanged.
 */
class HaloExchange final {
public:
    explicit HaloExchange(const MPIContext& mpi);

    /**
     * @brief Exchange halo layers for conservative state U.
     * @param layer Local conservative state storage.
     * @param mesh Local mesh with neighbor metadata.
     */
    void Exchange(DataLayer& layer, const Mesh& mesh) const;

private:
    const MPIContext& mpi_;

    void ExchangeX(DataLayer& layer, const Mesh& mesh) const;
    void ExchangeY(DataLayer& layer, const Mesh& mesh) const;
    void ExchangeZ(DataLayer& layer, const Mesh& mesh) const;

    [[nodiscard]] std::vector<double> PackX(const DataLayer& layer, const Mesh& mesh, int i_begin) const;
    [[nodiscard]] std::vector<double> PackY(const DataLayer& layer, const Mesh& mesh, int j_begin) const;
    [[nodiscard]] std::vector<double> PackZ(const DataLayer& layer, const Mesh& mesh, int k_begin) const;

    void UnpackX(DataLayer& layer, const Mesh& mesh, int i_begin, const std::vector<double>& buffer) const;
    void UnpackY(DataLayer& layer, const Mesh& mesh, int j_begin, const std::vector<double>& buffer) const;
    void UnpackZ(DataLayer& layer, const Mesh& mesh, int k_begin, const std::vector<double>& buffer) const;
};

#endif  // HALOEXCHANGE_HPP