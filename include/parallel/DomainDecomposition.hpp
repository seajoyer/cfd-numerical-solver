#ifndef DOMAINDECOMPOSITION_HPP
#define DOMAINDECOMPOSITION_HPP

#include <mpi.h>

#include "bc/BoundaryCondition.hpp"
#include "config/Settings.hpp"
#include "data/Mesh.hpp"
#include "data/Variables.hpp"
#include "parallel/MPIContext.hpp"

/**
 * @class DomainDecomposition
 * @brief Balanced Cartesian MPI domain decomposition for structured grids.
 *
 * The decomposition:
 *  - uses MPI_Dims_create for active dimensions
 *  - builds Cartesian communicator
 *  - computes local physical sizes and global offsets
 *  - computes neighbor ranks for each active axis
 *  - determines which local sides are global physical boundaries
 *
 * Local cell counts are split as evenly as possible:
 *  - base = global_n / proc_n
 *  - first (global_n % proc_n) subdomains get one extra cell
 */
class DomainDecomposition final {
public:
    DomainDecomposition(const Settings& settings, const MPIContext& mpi);
    ~DomainDecomposition();

    DomainDecomposition(const DomainDecomposition&) = delete;
    DomainDecomposition& operator=(const DomainDecomposition&) = delete;

    DomainDecomposition(DomainDecomposition&&) = delete;
    DomainDecomposition& operator=(DomainDecomposition&&) = delete;

    [[nodiscard]] int Dim() const;

    [[nodiscard]] int GlobalNx() const;
    [[nodiscard]] int GlobalNy() const;
    [[nodiscard]] int GlobalNz() const;

    [[nodiscard]] int LocalNx() const;
    [[nodiscard]] int LocalNy() const;
    [[nodiscard]] int LocalNz() const;

    [[nodiscard]] int OffsetX() const;
    [[nodiscard]] int OffsetY() const;
    [[nodiscard]] int OffsetZ() const;

    [[nodiscard]] int ProcCountX() const;
    [[nodiscard]] int ProcCountY() const;
    [[nodiscard]] int ProcCountZ() const;

    [[nodiscard]] int CoordX() const;
    [[nodiscard]] int CoordY() const;
    [[nodiscard]] int CoordZ() const;

    [[nodiscard]] int NeighborRank(Axis axis, Side side) const;
    [[nodiscard]] bool IsGlobalBoundary(Axis axis, Side side) const;

    [[nodiscard]] MPI_Comm CartComm() const;
    [[nodiscard]] int CartRank() const;

    /**
     * @brief Apply decomposition metadata to mesh.
     * @details Sets global sizes, offsets, neighbor ranks, and global boundary flags.
     */
    void ApplyToMesh(Mesh& mesh) const;

private:
    int dim_ = 1;

    int global_nx_ = 1;
    int global_ny_ = 1;
    int global_nz_ = 1;

    int local_nx_ = 1;
    int local_ny_ = 1;
    int local_nz_ = 1;

    int offset_x_ = 0;
    int offset_y_ = 0;
    int offset_z_ = 0;

    int proc_dims_[3] = {1, 1, 1};
    int coords_[3] = {0, 0, 0};
    int neighbor_rank_[3][2] = {
        {-1, -1},
        {-1, -1},
        {-1, -1}
    };
    bool is_global_boundary_[3][2] = {
        {true, true},
        {true, true},
        {true, true}
    };

    MPI_Comm cart_comm_ = MPI_COMM_NULL;
    int cart_rank_ = 0;

    void BuildCartesianTopology(const MPIContext& mpi);
    void ComputeLocalSizesAndOffsets();
    void ComputeNeighbors();
    void ComputeGlobalBoundaryFlags();

    static void ComputeBalancedPartition(int global_n, int proc_n, int coord,
                                         int& local_n, int& offset);
};

#endif  // DOMAINDECOMPOSITION_HPP