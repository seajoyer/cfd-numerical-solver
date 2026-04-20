#include "parallel/MeshDistribution.hpp"

#include <stdexcept>
#include <vector>

#include <mpi.h>

#include "parallel/DomainDecomposition.hpp"
#include "parallel/HaloSerialization.hpp"
#include "parallel/MPIContext.hpp"
#include "parallel/MeshSerialization.hpp"

namespace {
    void SendBuffer(const std::vector<char>& buffer,
                    int dst_rank,
                    int tag_size,
                    int tag_payload,
                    MPI_Comm comm) {
        const int size = static_cast<int>(buffer.size());
        MPI_Send(&size, 1, MPI_INT, dst_rank, tag_size, comm);

        if (size > 0) {
            MPI_Send(buffer.data(), size, MPI_CHAR, dst_rank, tag_payload, comm);
        }
    }

    std::vector<char> RecvBuffer(int src_rank,
                                 int tag_size,
                                 int tag_payload,
                                 MPI_Comm comm) {
        int size = 0;
        MPI_Recv(&size, 1, MPI_INT, src_rank, tag_size, comm, MPI_STATUS_IGNORE);

        if (size < 0) {
            throw std::runtime_error("MeshDistribution: received negative buffer size");
        }

        std::vector<char> buffer(static_cast<std::size_t>(size));
        if (size > 0) {
            MPI_Recv(buffer.data(), size, MPI_CHAR, src_rank, tag_payload, comm, MPI_STATUS_IGNORE);
        }

        return buffer;
    }
}

MeshDistribution::LocalPartition MeshDistribution::DistributeFromRoot(const Mesh* global_mesh,
                                                                      const MPIContext& mpi) {
    constexpr int k_mesh_size_tag = 4001;
    constexpr int k_mesh_payload_tag = 4002;
    constexpr int k_halo_size_tag = 4003;
    constexpr int k_halo_payload_tag = 4004;

    LocalPartition result;

    if (mpi.IsRoot()) {
        if (!global_mesh) {
            throw std::runtime_error("MeshDistribution: root rank requires global mesh");
        }

        const std::vector<int> part =
            DomainDecomposition::BuildRCBPartition(*global_mesh, mpi.Size());

        for (int rank = 0; rank < mpi.Size(); ++rank) {
            DomainDecomposition::Result local_result =
                DomainDecomposition::BuildLocalResultForRank(*global_mesh, part, rank);

            std::vector<char> mesh_buffer =
                MeshSerialization::Pack(local_result.local_mesh);

            std::vector<char> halo_buffer =
                HaloSerialization::Pack(local_result.halos);

            if (rank == mpi.Rank()) {
                result.mesh = std::move(local_result.local_mesh);
                result.halos = std::move(local_result.halos);
            }
            else {
                SendBuffer(mesh_buffer, rank, k_mesh_size_tag, k_mesh_payload_tag, mpi.Comm());
                SendBuffer(halo_buffer, rank, k_halo_size_tag, k_halo_payload_tag, mpi.Comm());
            }
        }
    }
    else {
        if (global_mesh != nullptr) {
            throw std::runtime_error("MeshDistribution: non-root rank must not pass global mesh");
        }

        const std::vector<char> mesh_buffer =
            RecvBuffer(0, k_mesh_size_tag, k_mesh_payload_tag, mpi.Comm());

        const std::vector<char> halo_buffer =
            RecvBuffer(0, k_halo_size_tag, k_halo_payload_tag, mpi.Comm());

        result.mesh = MeshSerialization::Unpack(mesh_buffer);
        result.halos = HaloSerialization::Unpack(halo_buffer);
    }

    result.mesh.Validate();
    return result;
}
