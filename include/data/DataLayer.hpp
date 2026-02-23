#ifndef DATALAYER_HPP
#define DATALAYER_HPP

#include <cstddef>
#include <stdexcept>
#include <xtensor.hpp>


enum class Axis : std::uint8_t;
enum class Side : std::uint8_t;

/**
 * @brief Storage owner for conservative state and grid metadata.
 *
 * DataLayer owns:
 *  - Conservative state U(var,i,j,k) with var = (rho, rhoU, rhoV, rhoW, E)
 *  - Coordinate arrays (cell boundaries and centers) for each axis
 *  - Cell-size metrics (dx,dy,dz) and inverses (inv_dx, inv_dy, inv_dz)
 *
 * Storage layout:
 *   U(var, i, j, k), var in {rho, rhoU, rhoV, rhoW, E}.
 *
 * Dimension handling (always stored as 3D):
 *  - dim = 1: ny=1, nz=1, sy=1, sz=1
 *  - dim = 2: nz=1, sz=1
 *  - dim = 3: full
 *
 * Ghost cells always exist: sx = nx + 2*ng, etc.
 *
 * Notes on coordinates:
 *  - boundary arrays xb/yb/zb have size (s? + 1)
 *  - center arrays   xc/yc/zc have size (s?)
 *  - cell sizes dx/dy/dz have size (s?)
 */
class DataLayer final {
public:
    static constexpr std::size_t k_nvar = 5;
    static constexpr std::size_t k_rho = 0;
    static constexpr std::size_t k_rhoU = 1;
    static constexpr std::size_t k_rhoV = 2;
    static constexpr std::size_t k_rhoW = 3;
    static constexpr std::size_t k_E = 4;

    DataLayer() = default;

    /**
     * @brief Construct storage for a given grid size, padding, and dimension.
     * @param nx Physical cells in x.
     * @param ny Physical cells in y (ignored for dim=1).
     * @param nz Physical cells in z (ignored for dim<=2).
     * @param padding Number of ghost cells on each side.
     * @param dim Spatial dimension (1..3).
     */
    DataLayer(int nx, int ny, int nz, int padding, int dim);

    // -------------------- sizes / meta --------------------
    [[nodiscard]] int GetDim() const;
    [[nodiscard]] int GetPadding() const;

    [[nodiscard]] int GetNx() const; // physical (core) cells
    [[nodiscard]] int GetNy() const; // physical (core) cells (returns 1 if dim<2)
    [[nodiscard]] int GetNz() const; // physical (core) cells (returns 1 if dim<3)

    [[nodiscard]] int GetSx() const; // total including ghosts
    [[nodiscard]] int GetSy() const;
    [[nodiscard]] int GetSz() const;

    [[nodiscard]] int GetCoreStartX() const;
    [[nodiscard]] int GetCoreStartY() const;
    [[nodiscard]] int GetCoreStartZ() const;

    [[nodiscard]] int GetCoreEndExclusiveX() const;
    [[nodiscard]] int GetCoreEndExclusiveY() const;
    [[nodiscard]] int GetCoreEndExclusiveZ() const;

    [[nodiscard]] int GetCoreNx() const;
    [[nodiscard]] int GetCoreNy() const;
    [[nodiscard]] int GetCoreNz() const;

    // -------------------- state access --------------------
    /** @brief Conservative state array U(var,i,j,k). Shape (5,sx,sy,sz). */
    [[nodiscard]] xt::xtensor<double, 4>& U();
    [[nodiscard]] const xt::xtensor<double, 4>& U() const;

    // -------------------- coordinates --------------------
    /** @brief Boundary coordinates (size sx+1). */
    [[nodiscard]] xt::xtensor<double, 1>& Xb();
    [[nodiscard]] const xt::xtensor<double, 1>& Xb() const;

    /** @brief Center coordinates (size sx). */
    [[nodiscard]] xt::xtensor<double, 1>& Xc();
    [[nodiscard]] const xt::xtensor<double, 1>& Xc() const;

    [[nodiscard]] xt::xtensor<double, 1>& Yb();
    [[nodiscard]] const xt::xtensor<double, 1>& Yb() const;

    [[nodiscard]] xt::xtensor<double, 1>& Yc();
    [[nodiscard]] const xt::xtensor<double, 1>& Yc() const;

    [[nodiscard]] xt::xtensor<double, 1>& Zb();
    [[nodiscard]] const xt::xtensor<double, 1>& Zb() const;

    [[nodiscard]] xt::xtensor<double, 1>& Zc();
    [[nodiscard]] const xt::xtensor<double, 1>& Zc() const;

    // -------------------- metrics --------------------
    /** @brief Cell sizes dx (size sx). */
    [[nodiscard]] const xt::xtensor<double, 1>& Dx() const;
    [[nodiscard]] const xt::xtensor<double, 1>& Dy() const;
    [[nodiscard]] const xt::xtensor<double, 1>& Dz() const;

    /** @brief Inverse cell sizes 1/dx (size sx). */
    [[nodiscard]] const xt::xtensor<double, 1>& InvDx() const;
    [[nodiscard]] const xt::xtensor<double, 1>& InvDy() const;
    [[nodiscard]] const xt::xtensor<double, 1>& InvDz() const;

    /**
     * @brief Recompute centers and cell-size metrics from boundary coordinates.
     * @details Fills xc/yc/zc, dx/dy/dz, inv_dx/inv_dy/inv_dz. No allocations.
     */
    void UpdateMetricsFromCoordinates();

    /**
     * @brief Check if a given (axis, side) is a global external boundary.
     * @details Physical BC should be applied only when this returns true.
     */
    [[nodiscard]] bool IsGlobalBoundary(Axis axis, Side side) const;

    /**
     * @brief Set whether a given (axis, side) is a global external boundary.
     * @details In single-domain: all true. In MPI: internal interfaces are false.
     */
    void SetGlobalBoundary(Axis axis, Side side, bool is_global);

    /** @brief Convenience: set all 6 sides at once. */
    void SetAllGlobalBoundaries(bool is_global);

private:
    int nx_ = 0;
    int ny_ = 1;
    int nz_ = 1;
    int ng_ = 0;
    int dim_ = 1;

    int sx_ = 1;
    int sy_ = 1;
    int sz_ = 1;

    // [axis][side], axis: X=0,Y=1,Z=2 ; side: Left=0,Right=1
    bool is_global_boundary_[3][2] = {
        {true, true},
        {true, true},
        {true, true}
    };

    xt::xtensor<double, 4> U_;

    // boundary coords: size (s? + 1)
    xt::xtensor<double, 1> xb_;
    xt::xtensor<double, 1> yb_;
    xt::xtensor<double, 1> zb_;

    // centers: size (s?)
    xt::xtensor<double, 1> xc_;
    xt::xtensor<double, 1> yc_;
    xt::xtensor<double, 1> zc_;

    // metrics: size (s?)
    xt::xtensor<double, 1> dx_;
    xt::xtensor<double, 1> dy_;
    xt::xtensor<double, 1> dz_;
    xt::xtensor<double, 1> inv_dx_;
    xt::xtensor<double, 1> inv_dy_;
    xt::xtensor<double, 1> inv_dz_;

    void Validate() const;
    void RecomputeSizes();
    void Allocate();
};

#endif  // DATALAYER_HPP
