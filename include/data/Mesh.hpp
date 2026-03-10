#ifndef MESH_HPP
#define MESH_HPP

#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <vector>

#include <xtensor.hpp>

enum class Axis : std::uint8_t;
enum class Side : std::uint8_t;

class GeometryPrimitive;

/**
 * @enum CellType
 * @brief Cell classification for structured Cartesian mesh.
 */
enum class CellType : std::uint8_t {
    Fluid = 0,
    Solid = 1
};

/**
 * @struct ImmersedFaceInfo
 * @brief Geometry data for one immersed face adjacent to a fluid cell.
 */
struct ImmersedFaceInfo final {
    bool is_active = false;
    double normal_x = 0.0;
    double normal_y = 0.0;
    double normal_z = 0.0;
    double distance = 0.0;
};

/**
 * @brief Structured Cartesian mesh owner for geometry, metrics, decomposition, and masks.
 *
 * Mesh owns:
 *  - Local logical grid sizes and padded sizes
 *  - Optional global grid sizes and local offsets for domain decomposition
 *  - Neighbor rank ids for each axis/side
 *  - Coordinate arrays (cell boundaries and centers) for each axis
 *  - Cell-size metrics (dx,dy,dz) and inverses (inv_dx,inv_dy,inv_dz)
 *  - Global external boundary flags for domain sides
 *  - Cell classification mask (fluid/solid)
 *  - Embedded geometry primitives
 *  - Immersed-face metadata for local cells
 */
class Mesh final {
public:
    Mesh() = default;

    /**
     * @brief Construct mesh for a given local grid size, padding, and dimension.
     * @param nx Local physical cells in x.
     * @param ny Local physical cells in y (ignored for dim=1).
     * @param nz Local physical cells in z (ignored for dim<=2).
     * @param padding Number of ghost cells on each side.
     * @param dim Spatial dimension (1..3).
     */
    Mesh(int nx, int ny, int nz, int padding, int dim);

    // -------------------- local sizes / meta --------------------
    [[nodiscard]] int GetDim() const;
    [[nodiscard]] int GetPadding() const;

    [[nodiscard]] int GetNx() const;
    [[nodiscard]] int GetNy() const;
    [[nodiscard]] int GetNz() const;

    [[nodiscard]] int GetSx() const;
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

    // -------------------- global decomposition metadata --------------------
    /**
     * @brief Set global grid sizes and local physical-cell offsets.
     * @details For single-domain runs, global sizes normally equal local sizes and offsets are zero.
     */
    void SetGlobalDecomposition(int global_nx, int global_ny, int global_nz,
                                int offset_x, int offset_y, int offset_z);

    [[nodiscard]] int GetGlobalNx() const;
    [[nodiscard]] int GetGlobalNy() const;
    [[nodiscard]] int GetGlobalNz() const;

    [[nodiscard]] int GetOffsetX() const;
    [[nodiscard]] int GetOffsetY() const;
    [[nodiscard]] int GetOffsetZ() const;

    /**
     * @brief Set neighbor rank for a local domain side.
     * @param axis Spatial axis.
     * @param side Side along the axis.
     * @param rank Neighbor rank id, or -1 if absent.
     */
    void SetNeighborRank(Axis axis, Side side, int rank);

    /**
     * @brief Get neighbor rank for a local domain side.
     * @return Neighbor rank id, or -1 if absent.
     */
    [[nodiscard]] int GetNeighborRank(Axis axis, Side side) const;

    /** @brief Check whether the local side has an MPI neighbor. */
    [[nodiscard]] bool HasNeighborRank(Axis axis, Side side) const;

    // -------------------- coordinates --------------------
    [[nodiscard]] xt::xtensor<double, 1>& Xb();
    [[nodiscard]] const xt::xtensor<double, 1>& Xb() const;

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
    [[nodiscard]] const xt::xtensor<double, 1>& Dx() const;
    [[nodiscard]] const xt::xtensor<double, 1>& Dy() const;
    [[nodiscard]] const xt::xtensor<double, 1>& Dz() const;

    [[nodiscard]] const xt::xtensor<double, 1>& InvDx() const;
    [[nodiscard]] const xt::xtensor<double, 1>& InvDy() const;
    [[nodiscard]] const xt::xtensor<double, 1>& InvDz() const;

    /**
     * @brief Recompute centers and cell-size metrics from boundary coordinates.
     * @details Fills xc/yc/zc, dx/dy/dz, inv_dx/inv_dy/inv_dz. No allocations.
     */
    void UpdateMetricsFromCoordinates();

    // -------------------- global boundaries --------------------
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

    // -------------------- cell classification --------------------
    /** @brief Cell classification array cell_type(i,j,k). Shape (sx,sy,sz). */
    [[nodiscard]] xt::xtensor<std::uint8_t, 3>& CellTypes();
    [[nodiscard]] const xt::xtensor<std::uint8_t, 3>& CellTypes() const;

    /** @brief Set all cells to Fluid. */
    void SetAllCellsFluid();

    /** @brief Set all cells to Solid. */
    void SetAllCellsSolid();

    /** @brief Set one cell classification. */
    void SetCellType(int i, int j, int k, CellType type);

    /** @brief Get one cell classification. */
    [[nodiscard]] CellType GetCellType(int i, int j, int k) const;

    /** @brief Check whether a cell is fluid. */
    [[nodiscard]] bool IsFluidCell(int i, int j, int k) const;

    /** @brief Check whether a cell is solid. */
    [[nodiscard]] bool IsSolidCell(int i, int j, int k) const;

    // -------------------- geometry primitives --------------------
    /**
     * @brief Add one embedded geometry primitive.
     * @param primitive Shared pointer to geometry primitive.
     */
    void AddPrimitive(std::shared_ptr<GeometryPrimitive> primitive);

    /** @brief Remove all embedded geometry primitives. */
    void ClearPrimitives();

    /** @brief Access embedded geometry primitives. */
    [[nodiscard]] const std::vector<std::shared_ptr<GeometryPrimitive>>& Primitives() const;

    /**
     * @brief Rebuild fluid/solid mask from current primitive list.
     * @details A cell is marked Solid if its center lies inside at least one primitive.
     */
    void BuildCellTypesFromPrimitives();

    /**
     * @brief Rebuild immersed-face metadata from current fluid/solid mask and primitive list.
     * @details An immersed face is created when a fluid cell has a solid neighbor along an active axis.
     */
    void BuildImmersedFaces();

    // -------------------- immersed-face metadata --------------------
    /** @brief Clear all immersed-face metadata on the local mesh. */
    void ClearAllImmersedFaces();

    /** @brief Clear one immersed-face entry. */
    void ClearImmersedFace(Axis axis, Side side, int i, int j, int k);

    /**
     * @brief Set one immersed-face entry.
     * @param axis Spatial axis of the face slot.
     * @param side Side of the face relative to cell (i,j,k).
     * @param i,j,k Local cell indices.
     * @param normal_x Boundary normal x-component.
     * @param normal_y Boundary normal y-component.
     * @param normal_z Boundary normal z-component.
     * @param distance Distance from cell center to embedded boundary.
     */
    void SetImmersedFace(Axis axis, Side side, int i, int j, int k,
                         double normal_x, double normal_y, double normal_z,
                         double distance);

    /** @brief Check whether a cell has immersed face on given axis/side. */
    [[nodiscard]] bool HasImmersedFace(Axis axis, Side side, int i, int j, int k) const;

    /** @brief Get immersed-face data for one axis/side/cell slot. */
    [[nodiscard]] ImmersedFaceInfo GetImmersedFaceInfo(Axis axis, Side side, int i, int j, int k) const;

    /** @brief Check whether a cell has at least one immersed face. */
    [[nodiscard]] bool HasAnyImmersedFace(int i, int j, int k) const;

private:
    int nx_ = 0;
    int ny_ = 1;
    int nz_ = 1;
    int ng_ = 0;
    int dim_ = 1;

    int sx_ = 1;
    int sy_ = 1;
    int sz_ = 1;

    int global_nx_ = 0;
    int global_ny_ = 1;
    int global_nz_ = 1;

    int offset_x_ = 0;
    int offset_y_ = 0;
    int offset_z_ = 0;

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

    xt::xtensor<double, 1> xb_;
    xt::xtensor<double, 1> yb_;
    xt::xtensor<double, 1> zb_;

    xt::xtensor<double, 1> xc_;
    xt::xtensor<double, 1> yc_;
    xt::xtensor<double, 1> zc_;

    xt::xtensor<double, 1> dx_;
    xt::xtensor<double, 1> dy_;
    xt::xtensor<double, 1> dz_;
    xt::xtensor<double, 1> inv_dx_;
    xt::xtensor<double, 1> inv_dy_;
    xt::xtensor<double, 1> inv_dz_;

    xt::xtensor<std::uint8_t, 3> cell_type_;

    xt::xtensor<std::uint8_t, 5> immersed_face_active_;
    xt::xtensor<double, 6> immersed_face_normal_;
    xt::xtensor<double, 5> immersed_face_distance_;

    std::vector<std::shared_ptr<GeometryPrimitive>> primitives_;

    void Validate() const;
    void RecomputeSizes();
    void Allocate();

    void ValidateCellIndex(int i, int j, int k) const;
    void ValidateAxisActive(Axis axis) const;

    [[nodiscard]] int AxisIndex(Axis axis) const;
    [[nodiscard]] int SideIndex(Side side) const;

    [[nodiscard]] const GeometryPrimitive* FindClosestPrimitive(double x, double y, double z) const;
};

#endif  // MESH_HPP