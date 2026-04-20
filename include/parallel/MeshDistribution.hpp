#ifndef MESHDISTRIBUTION_HPP
#define MESHDISTRIBUTION_HPP

#include <memory>
#include <vector>

#include "geometry/Mesh.hpp"
#include "parallel/DomainDecomposition.hpp"

class MPIContext;

class MeshDistribution final {
public:
    struct LocalPartition final {
        Mesh mesh;
        std::vector<DomainDecomposition::NeighborHalo> halos;
    };

    /**
     * @brief Root owns global mesh, decomposes it, and distributes local meshes to all ranks.
     *
     * On root:
     * - global_mesh must be non-null
     * On non-root:
     * - global_mesh must be null
     */
    [[nodiscard]] static LocalPartition DistributeFromRoot(const Mesh* global_mesh,
                                                           const MPIContext& mpi);
};

#endif  // MESHDISTRIBUTION_HPP
