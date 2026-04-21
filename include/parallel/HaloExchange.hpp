#ifndef HALOEXCHANGE_HPP
#define HALOEXCHANGE_HPP

#include <cstddef>
#include <vector>

#include <mpi.h>

#include "parallel/DomainDecomposition.hpp"
#include "parallel/StateSynchronizer.hpp"

class DataLayer;
class MPIContext;
// CHECK: HALO_EXCHANGE
/**
 * @class HaloExchange
 * @brief Exchanges ghost-cell conservative state between neighboring MPI ranks.
 *
 * Data packet sent per cell:
 * - U[5] conservative variables
 * - lambda reactant mass fraction
 *
 * Communication pattern:
 * - post all Irecv
 * - pack send buffers
 * - post all Isend
 * - wait all
 * - unpack received ghost states
 */
class HaloExchange final : public StateSynchronizer {
public:
    struct CellStatePacket final {
        double U[5];
        double lambda = 0.0;
    };

    HaloExchange(const MPIContext& mpi,
                 std::vector<DomainDecomposition::NeighborHalo> halos);

    ~HaloExchange();

    HaloExchange(const HaloExchange&) = delete;
    HaloExchange& operator=(const HaloExchange&) = delete;

    HaloExchange(HaloExchange&&) = delete;
    HaloExchange& operator=(HaloExchange&&) = delete;

    /**
     * @brief Exchange halo values and write them into local ghost cells.
     */
    void Synchronize(DataLayer& layer) const override;

private:
    const MPIContext* mpi_ = nullptr;
    std::vector<DomainDecomposition::NeighborHalo> halos_;

    MPI_Datatype packet_type_ = MPI_DATATYPE_NULL;

    void CreatePacketType();
    void DestroyPacketType();

    [[nodiscard]] std::vector<CellStatePacket> PackSendBuffer(
        const DataLayer& layer,
        const std::vector<std::size_t>& send_local_ids
    ) const;

    void UnpackRecvBuffer(
        DataLayer& layer,
        const std::vector<std::size_t>& recv_local_ids,
        const std::vector<CellStatePacket>& recv_buffer
    ) const;
};

#endif  // HALOEXCHANGE_HPP
