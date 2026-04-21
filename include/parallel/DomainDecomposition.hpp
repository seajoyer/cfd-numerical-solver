#ifndef DOMAINDECOMPOSITION_HPP
#define DOMAINDECOMPOSITION_HPP

#include <cstddef>
#include <unordered_map>
#include <vector>

#include "geometry/Mesh.hpp"

class MPIContext;
// CHECK: DECOMPOSITION
/**
 * @class DomainDecomposition
 * @brief Builds MPI-local mesh partition and halo metadata from a global mesh.
 *
 * Current implementation:
 * - recursive coordinate bisection (RCB)
 * - owned cells first, ghost cells second
 * - MPI faces are materialized explicitly as FaceKind::MPIBoundary
 * - face owner always points to owned cell, neighbor points to local ghost cell
 */
class DomainDecomposition final {
public:
    struct NeighborHalo final {
        int remote_rank = -1;

        /**
         * @brief Local owned cell ids that must be sent to remote rank.
         */
        std::vector<std::size_t> send_local_ids;

        /**
         * @brief Local ghost cell ids that must be received from remote rank.
         */
        std::vector<std::size_t> recv_local_ids;
    };

    struct Result final {
        Mesh local_mesh;

        std::size_t n_owned_cells = 0;
        std::size_t n_ghost_cells = 0;

        /**
         * @brief part[global_cell_id] = rank
         */
        std::vector<int> global_part;

        /**
         * @brief local_cell_id -> global/original cell id
         */
        std::vector<std::size_t> local_to_global_cell;

        /**
         * @brief global/original cell id -> local_cell_id
         */
        std::unordered_map<std::size_t, std::size_t> global_to_local_cell;

        /**
         * @brief Local owned global/original ids in local ordering.
         */
        std::vector<std::size_t> owned_global_ids;

        /**
         * @brief Local ghost global/original ids in local ordering.
         */
        std::vector<std::size_t> ghost_global_ids;

        /**
         * @brief Neighbor halo metadata for future halo exchange.
         */
        std::vector<NeighborHalo> halos;
    };

    /**
     * @brief Build RCB partition and extract MPI-local mesh for current rank.
     *
     * Assumptions:
     * - every rank has access to the same global mesh
     * - mesh cell id is stable and unique in the global mesh
     */
    [[nodiscard]] static Result BuildRCB(const Mesh& global_mesh,
                                         const MPIContext& mpi);

    [[nodiscard]] static std::vector<int> BuildRCBPartition(const Mesh& global_mesh,
                                                            int nproc);

    [[nodiscard]] static Result BuildLocalResultForRank(const Mesh& global_mesh,
                                                        const std::vector<int>& global_part,
                                                        int rank);
};

#endif  // DOMAINDECOMPOSITION_HPP
