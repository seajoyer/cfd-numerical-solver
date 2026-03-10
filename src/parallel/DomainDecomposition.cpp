#include "parallel/DomainDecomposition.hpp"

#include <stdexcept>

DomainDecomposition::DomainDecomposition(const Settings& settings, const MPIContext& mpi)
    : dim_(settings.dim),
      global_nx_(settings.GetNx()),
      global_ny_(settings.dim >= 2 ? settings.GetNy() : 1),
      global_nz_(settings.dim >= 3 ? settings.GetNz() : 1) {
    if (dim_ < 1 || dim_ > 3) {
        throw std::invalid_argument("DomainDecomposition: dim must be 1..3");
    }

    BuildCartesianTopology(mpi);
    ComputeLocalSizesAndOffsets();
    ComputeNeighbors();
    ComputeGlobalBoundaryFlags();
}

DomainDecomposition::~DomainDecomposition() {
    if (cart_comm_ != MPI_COMM_NULL) {
        MPI_Comm_free(&cart_comm_);
        cart_comm_ = MPI_COMM_NULL;
    }
}

int DomainDecomposition::Dim() const {
    return dim_;
}

int DomainDecomposition::GlobalNx() const {
    return global_nx_;
}

int DomainDecomposition::GlobalNy() const {
    return global_ny_;
}

int DomainDecomposition::GlobalNz() const {
    return global_nz_;
}

int DomainDecomposition::LocalNx() const {
    return local_nx_;
}

int DomainDecomposition::LocalNy() const {
    return local_ny_;
}

int DomainDecomposition::LocalNz() const {
    return local_nz_;
}

int DomainDecomposition::OffsetX() const {
    return offset_x_;
}

int DomainDecomposition::OffsetY() const {
    return offset_y_;
}

int DomainDecomposition::OffsetZ() const {
    return offset_z_;
}

int DomainDecomposition::ProcCountX() const {
    return proc_dims_[0];
}

int DomainDecomposition::ProcCountY() const {
    return proc_dims_[1];
}

int DomainDecomposition::ProcCountZ() const {
    return proc_dims_[2];
}

int DomainDecomposition::CoordX() const {
    return coords_[0];
}

int DomainDecomposition::CoordY() const {
    return coords_[1];
}

int DomainDecomposition::CoordZ() const {
    return coords_[2];
}

int DomainDecomposition::NeighborRank(const Axis axis, const Side side) const {
    const auto a = static_cast<int>(static_cast<std::uint8_t>(axis));
    const auto s = static_cast<int>(static_cast<std::uint8_t>(side));
    return neighbor_rank_[a][s];
}

bool DomainDecomposition::IsGlobalBoundary(const Axis axis, const Side side) const {
    const auto a = static_cast<int>(static_cast<std::uint8_t>(axis));
    const auto s = static_cast<int>(static_cast<std::uint8_t>(side));
    return is_global_boundary_[a][s];
}

MPI_Comm DomainDecomposition::CartComm() const {
    return cart_comm_;
}

int DomainDecomposition::CartRank() const {
    return cart_rank_;
}

void DomainDecomposition::ApplyToMesh(Mesh& mesh) const {
    mesh.SetGlobalDecomposition(global_nx_, global_ny_, global_nz_,
                                offset_x_, offset_y_, offset_z_);

    mesh.SetNeighborRank(Axis::X, Side::Left, neighbor_rank_[0][0]);
    mesh.SetNeighborRank(Axis::X, Side::Right, neighbor_rank_[0][1]);
    mesh.SetGlobalBoundary(Axis::X, Side::Left, is_global_boundary_[0][0]);
    mesh.SetGlobalBoundary(Axis::X, Side::Right, is_global_boundary_[0][1]);

    if (dim_ >= 2) {
        mesh.SetNeighborRank(Axis::Y, Side::Left, neighbor_rank_[1][0]);
        mesh.SetNeighborRank(Axis::Y, Side::Right, neighbor_rank_[1][1]);
        mesh.SetGlobalBoundary(Axis::Y, Side::Left, is_global_boundary_[1][0]);
        mesh.SetGlobalBoundary(Axis::Y, Side::Right, is_global_boundary_[1][1]);
    }

    if (dim_ >= 3) {
        mesh.SetNeighborRank(Axis::Z, Side::Left, neighbor_rank_[2][0]);
        mesh.SetNeighborRank(Axis::Z, Side::Right, neighbor_rank_[2][1]);
        mesh.SetGlobalBoundary(Axis::Z, Side::Left, is_global_boundary_[2][0]);
        mesh.SetGlobalBoundary(Axis::Z, Side::Right, is_global_boundary_[2][1]);
    }
}

void DomainDecomposition::BuildCartesianTopology(const MPIContext& mpi) {
    int dims[3] = {0, 0, 0};
    int periods[3] = {0, 0, 0};

    dims[0] = 0;
    dims[1] = dim_ >= 2 ? 0 : 1;
    dims[2] = dim_ >= 3 ? 0 : 1;

    const int nnodes = mpi.Size();
    MPI_Dims_create(nnodes, 3, dims);

    proc_dims_[0] = dims[0];
    proc_dims_[1] = dim_ >= 2 ? dims[1] : 1;
    proc_dims_[2] = dim_ >= 3 ? dims[2] : 1;

    MPI_Cart_create(mpi.Comm(), 3, proc_dims_, periods, 0, &cart_comm_);
    if (cart_comm_ == MPI_COMM_NULL) {
        throw std::runtime_error("DomainDecomposition: MPI_Cart_create failed");
    }

    MPI_Comm_rank(cart_comm_, &cart_rank_);
    MPI_Cart_coords(cart_comm_, cart_rank_, 3, coords_);
}

void DomainDecomposition::ComputeLocalSizesAndOffsets() {
    ComputeBalancedPartition(global_nx_, proc_dims_[0], coords_[0], local_nx_, offset_x_);

    if (dim_ >= 2) {
        ComputeBalancedPartition(global_ny_, proc_dims_[1], coords_[1], local_ny_, offset_y_);
    }
    else {
        local_ny_ = 1;
        offset_y_ = 0;
    }

    if (dim_ >= 3) {
        ComputeBalancedPartition(global_nz_, proc_dims_[2], coords_[2], local_nz_, offset_z_);
    }
    else {
        local_nz_ = 1;
        offset_z_ = 0;
    }
}

void DomainDecomposition::ComputeNeighbors() {
    int src = -1;
    int dst = -1;

    MPI_Cart_shift(cart_comm_, 0, 1, &src, &dst);
    neighbor_rank_[0][0] = src;
    neighbor_rank_[0][1] = dst;

    if (dim_ >= 2) {
        MPI_Cart_shift(cart_comm_, 1, 1, &src, &dst);
        neighbor_rank_[1][0] = src;
        neighbor_rank_[1][1] = dst;
    }
    else {
        neighbor_rank_[1][0] = -1;
        neighbor_rank_[1][1] = -1;
    }

    if (dim_ >= 3) {
        MPI_Cart_shift(cart_comm_, 2, 1, &src, &dst);
        neighbor_rank_[2][0] = src;
        neighbor_rank_[2][1] = dst;
    }
    else {
        neighbor_rank_[2][0] = -1;
        neighbor_rank_[2][1] = -1;
    }
}

void DomainDecomposition::ComputeGlobalBoundaryFlags() {
    is_global_boundary_[0][0] = coords_[0] == 0;
    is_global_boundary_[0][1] = coords_[0] == proc_dims_[0] - 1;

    if (dim_ >= 2) {
        is_global_boundary_[1][0] = coords_[1] == 0;
        is_global_boundary_[1][1] = coords_[1] == proc_dims_[1] - 1;
    }
    else {
        is_global_boundary_[1][0] = false;
        is_global_boundary_[1][1] = false;
    }

    if (dim_ >= 3) {
        is_global_boundary_[2][0] = coords_[2] == 0;
        is_global_boundary_[2][1] = coords_[2] == proc_dims_[2] - 1;
    }
    else {
        is_global_boundary_[2][0] = false;
        is_global_boundary_[2][1] = false;
    }
}

void DomainDecomposition::ComputeBalancedPartition(const int global_n,
                                                   const int proc_n,
                                                   const int coord,
                                                   int& local_n,
                                                   int& offset) {
    if (global_n <= 0) {
        throw std::invalid_argument("DomainDecomposition: global_n must be > 0");
    }
    if (proc_n <= 0) {
        throw std::invalid_argument("DomainDecomposition: proc_n must be > 0");
    }
    if (coord < 0 || coord >= proc_n) {
        throw std::invalid_argument("DomainDecomposition: coord is out of range");
    }

    const int base = global_n / proc_n;
    const int rem = global_n % proc_n;

    local_n = base + (coord < rem ? 1 : 0);
    offset = coord * base + (coord < rem ? coord : rem);
}
