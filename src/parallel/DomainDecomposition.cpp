#include "parallel/DomainDecomposition.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include "utils/StringUtils.hpp"

namespace {
    [[nodiscard]] inline double Square(const double x) {
        return x * x;
    }
} // namespace

DomainDecomposition::DomainDecomposition(const Settings& settings, const MPIContext& mpi)
    : dim_(settings.dim),
      global_nx_(settings.GetNx()),
      global_ny_(settings.dim >= 2 ? settings.GetNy() : 1),
      global_nz_(settings.dim >= 3 ? settings.GetNz() : 1) {
    if (dim_ < 1 || dim_ > 3) {
        throw std::invalid_argument("DomainDecomposition: dim must be 1..3");
    }

    if (global_nx_ <= 0) {
        throw std::invalid_argument("DomainDecomposition: global_nx must be > 0");
    }
    if (dim_ >= 2 && global_ny_ <= 0) {
        throw std::invalid_argument("DomainDecomposition: global_ny must be > 0 for dim >= 2");
    }
    if (dim_ >= 3 && global_nz_ <= 0) {
        throw std::invalid_argument("DomainDecomposition: global_nz must be > 0 for dim >= 3");
    }

    periodic_[0] =
        utils::ToLower(settings.left_boundary) == "periodic" &&
        utils::ToLower(settings.right_boundary) == "periodic";

    periodic_[1] =
        dim_ >= 2 &&
        utils::ToLower(settings.bottom_boundary) == "periodic" &&
        utils::ToLower(settings.top_boundary) == "periodic";

    periodic_[2] =
        dim_ >= 3 &&
        utils::ToLower(settings.back_boundary) == "periodic" &&
        utils::ToLower(settings.front_boundary) == "periodic";

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
    ChooseProcessGridFromCells(mpi.Size());

    int periods[3] = {
        periodic_[0] ? 1 : 0,
        periodic_[1] ? 1 : 0,
        periodic_[2] ? 1 : 0
    };

    MPI_Cart_create(mpi.Comm(), 3, proc_dims_, periods, 0, &cart_comm_);
    if (cart_comm_ == MPI_COMM_NULL) {
        throw std::runtime_error("DomainDecomposition: MPI_Cart_create failed");
    }

    MPI_Comm_rank(cart_comm_, &cart_rank_);
    MPI_Cart_coords(cart_comm_, cart_rank_, 3, coords_);
}

void DomainDecomposition::ChooseProcessGridFromCells(const int world_size) {
    if (world_size <= 0) {
        throw std::invalid_argument("DomainDecomposition: world_size must be > 0");
    }

    proc_dims_[0] = 1;
    proc_dims_[1] = 1;
    proc_dims_[2] = 1;

    // 1D
    if (dim_ == 1) {
        if (world_size > global_nx_) {
            throw std::runtime_error(
                "DomainDecomposition: number of MPI ranks exceeds global_nx in 1D");
        }

        proc_dims_[0] = world_size;
        proc_dims_[1] = 1;
        proc_dims_[2] = 1;
        return;
    }

    double best_score = std::numeric_limits<double>::infinity();
    bool found = false;

    // 2D
    if (dim_ == 2) {
        for (int px = 1; px <= world_size; ++px) {
            if (world_size % px != 0) {
                continue;
            }

            const int py = world_size / px;

            if (px > global_nx_ || py > global_ny_) {
                continue;
            }

            const double score = ComputeProcessGridScore(px, py, 1);

            if (!found || score < best_score) {
                best_score = score;
                proc_dims_[0] = px;
                proc_dims_[1] = py;
                proc_dims_[2] = 1;
                found = true;
            }
        }

        if (!found) {
            throw std::runtime_error(
                "DomainDecomposition: cannot build valid 2D process grid "
                "(too many ranks for available cells)");
        }

        return;
    }

    // 3D
    for (int px = 1; px <= world_size; ++px) {
        if (world_size % px != 0) {
            continue;
        }

        const int rem_xy = world_size / px;

        for (int py = 1; py <= rem_xy; ++py) {
            if (rem_xy % py != 0) {
                continue;
            }

            const int pz = rem_xy / py;

            if (px > global_nx_ || py > global_ny_ || pz > global_nz_) {
                continue;
            }

            const double score = ComputeProcessGridScore(px, py, pz);

            if (!found || score < best_score) {
                best_score = score;
                proc_dims_[0] = px;
                proc_dims_[1] = py;
                proc_dims_[2] = pz;
                found = true;
            }
        }
    }

    if (!found) {
        throw std::runtime_error(
            "DomainDecomposition: cannot build valid 3D process grid "
            "(too many ranks for available cells)");
    }
}

double DomainDecomposition::ComputeProcessGridScore(const int px,
                                                    const int py,
                                                    const int pz) const {
    if (px <= 0 || py <= 0 || pz <= 0) {
        return std::numeric_limits<double>::infinity();
    }

    const double lx = static_cast<double>(global_nx_) / static_cast<double>(px);
    const double ly = static_cast<double>(global_ny_) / static_cast<double>(py);
    const double lz = static_cast<double>(global_nz_) / static_cast<double>(pz);

    // Оцениваем, насколько "похожи" локальные размеры по осям.
    // Используем логарифмы отношений, чтобы оценка была симметричной:
    // lx=100, ly=10 и lx=10, ly=100 дают одинаковый штраф.
    if (dim_ == 1) {
        return 0.0;
    }

    if (dim_ == 2) {
        const double rxy = std::log(lx / ly);
        return Square(rxy);
    }

    const double rxy = std::log(lx / ly);
    const double rxz = std::log(lx / lz);
    const double ryz = std::log(ly / lz);

    return Square(rxy) + Square(rxz) + Square(ryz);
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
    if (periodic_[0]) {
        is_global_boundary_[0][0] = false;
        is_global_boundary_[0][1] = false;
    }
    else {
        is_global_boundary_[0][0] = coords_[0] == 0;
        is_global_boundary_[0][1] = coords_[0] == proc_dims_[0] - 1;
    }

    if (dim_ >= 2) {
        if (periodic_[1]) {
            is_global_boundary_[1][0] = false;
            is_global_boundary_[1][1] = false;
        }
        else {
            is_global_boundary_[1][0] = coords_[1] == 0;
            is_global_boundary_[1][1] = coords_[1] == proc_dims_[1] - 1;
        }
    }
    else {
        is_global_boundary_[1][0] = false;
        is_global_boundary_[1][1] = false;
    }

    if (dim_ >= 3) {
        if (periodic_[2]) {
            is_global_boundary_[2][0] = false;
            is_global_boundary_[2][1] = false;
        }
        else {
            is_global_boundary_[2][0] = coords_[2] == 0;
            is_global_boundary_[2][1] = coords_[2] == proc_dims_[2] - 1;
        }
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
    if (proc_n > global_n) {
        throw std::invalid_argument(
            "DomainDecomposition: proc_n must not exceed global_n "
            "(empty subdomains are not allowed)");
    }

    const int base = global_n / proc_n;
    const int rem = global_n % proc_n;

    local_n = base + (coord < rem ? 1 : 0);
    offset = coord * base + (coord < rem ? coord : rem);
}
