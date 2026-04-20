#include "parallel/HaloExchange.hpp"

#include <cstddef>
#include <stdexcept>
#include <vector>

#include <mpi.h>

#include "data/DataLayer.hpp"
#include "parallel/MPIContext.hpp"

HaloExchange::HaloExchange(const MPIContext& mpi,
                           std::vector<DomainDecomposition::NeighborHalo> halos)
    : mpi_(&mpi),
      halos_(std::move(halos)) {
    if (!mpi_) {
        throw std::runtime_error("HaloExchange: mpi context is null");
    }

    CreatePacketType();
}

HaloExchange::~HaloExchange() {
    DestroyPacketType();
}

void HaloExchange::CreatePacketType() {
    if (packet_type_ != MPI_DATATYPE_NULL) {
        return;
    }

    CellStatePacket sample{};

    int block_lengths[2] = {5, 1};
    MPI_Aint displacements[2] = {};
    MPI_Datatype types[2] = {MPI_DOUBLE, MPI_DOUBLE};

    MPI_Aint base_address = 0;
    MPI_Aint u_address = 0;
    MPI_Aint lambda_address = 0;

    MPI_Get_address(&sample, &base_address);
    MPI_Get_address(&sample.U[0], &u_address);
    MPI_Get_address(&sample.lambda, &lambda_address);

    displacements[0] = u_address - base_address;
    displacements[1] = lambda_address - base_address;

    MPI_Type_create_struct(2, block_lengths, displacements, types, &packet_type_);
    MPI_Type_commit(&packet_type_);
}

void HaloExchange::DestroyPacketType() {
    if (packet_type_ != MPI_DATATYPE_NULL) {
        MPI_Type_free(&packet_type_);
        packet_type_ = MPI_DATATYPE_NULL;
    }
}

std::vector<HaloExchange::CellStatePacket> HaloExchange::PackSendBuffer(
    const DataLayer& layer,
    const std::vector<std::size_t>& send_local_ids
) const {
    const auto& U = layer.U();
    const auto& lambda = layer.ReactantMassFraction();

    std::vector<CellStatePacket> buffer(send_local_ids.size());

    for (std::size_t i = 0; i < send_local_ids.size(); ++i) {
        const std::size_t cell_id = send_local_ids[i];

        buffer[i].U[0] = U(cell_id, DataLayer::k_rho);
        buffer[i].U[1] = U(cell_id, DataLayer::k_rhoU);
        buffer[i].U[2] = U(cell_id, DataLayer::k_rhoV);
        buffer[i].U[3] = U(cell_id, DataLayer::k_rhoW);
        buffer[i].U[4] = U(cell_id, DataLayer::k_E);
        buffer[i].lambda = lambda(cell_id);
    }

    return buffer;
}

void HaloExchange::UnpackRecvBuffer(
    DataLayer& layer,
    const std::vector<std::size_t>& recv_local_ids,
    const std::vector<CellStatePacket>& recv_buffer
) const {
    if (recv_local_ids.size() != recv_buffer.size()) {
        throw std::runtime_error("HaloExchange::UnpackRecvBuffer: buffer size mismatch");
    }

    auto& U = layer.U();
    auto& lambda = layer.ReactantMassFraction();

    for (std::size_t i = 0; i < recv_local_ids.size(); ++i) {
        const std::size_t cell_id = recv_local_ids[i];

        U(cell_id, DataLayer::k_rho) = recv_buffer[i].U[0];
        U(cell_id, DataLayer::k_rhoU) = recv_buffer[i].U[1];
        U(cell_id, DataLayer::k_rhoV) = recv_buffer[i].U[2];
        U(cell_id, DataLayer::k_rhoW) = recv_buffer[i].U[3];
        U(cell_id, DataLayer::k_E) = recv_buffer[i].U[4];
        lambda(cell_id) = recv_buffer[i].lambda;
    }
}

void HaloExchange::Synchronize(DataLayer& layer) const {
    if (!mpi_) {
        throw std::runtime_error("HaloExchange::Exchange: mpi context is null");
    }

    if (packet_type_ == MPI_DATATYPE_NULL) {
        throw std::runtime_error("HaloExchange::Exchange: packet MPI datatype is not initialized");
    }

    if (halos_.empty()) {
        return;
    }

    std::vector<std::vector<CellStatePacket>> send_buffers;
    std::vector<std::vector<CellStatePacket>> recv_buffers;
    std::vector<MPI_Request> requests;

    send_buffers.reserve(halos_.size());
    recv_buffers.reserve(halos_.size());
    requests.reserve(2 * halos_.size());

    constexpr int k_halo_tag = 1001;

    for (const DomainDecomposition::NeighborHalo& halo : halos_) {
        recv_buffers.emplace_back(halo.recv_local_ids.size());

        MPI_Request recv_request = MPI_REQUEST_NULL;
        MPI_Irecv(
            recv_buffers.back().data(),
            static_cast<int>(recv_buffers.back().size()),
            packet_type_,
            halo.remote_rank,
            k_halo_tag,
            mpi_->Comm(),
            &recv_request
        );
        requests.push_back(recv_request);
    }

    for (const DomainDecomposition::NeighborHalo& halo : halos_) {
        send_buffers.push_back(PackSendBuffer(layer, halo.send_local_ids));

        MPI_Request send_request = MPI_REQUEST_NULL;
        MPI_Isend(
            send_buffers.back().data(),
            static_cast<int>(send_buffers.back().size()),
            packet_type_,
            halo.remote_rank,
            k_halo_tag,
            mpi_->Comm(),
            &send_request
        );
        requests.push_back(send_request);
    }

    MPI_Waitall(
        static_cast<int>(requests.size()),
        requests.data(),
        MPI_STATUSES_IGNORE
    );

    for (std::size_t i = 0; i < halos_.size(); ++i) {
        UnpackRecvBuffer(layer, halos_[i].recv_local_ids, recv_buffers[i]);
    }
}
