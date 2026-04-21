#ifndef MESH_HPP
#define MESH_HPP

#include <cstddef>
#include <stdexcept>
#include <vector>

#include "geometry/Cell.hpp"
#include "geometry/Face.hpp"
#include "geometry/Node.hpp"
// CHECK: MESH_STRUCTS
/**
 * @brief Single-process or MPI-local unstructured mesh owner for cell-centered finite volume methods.
 *
 * Mesh owns:
 *  - nodes
 *  - faces
 *  - cells
 *
 * Geometry and connectivity must already be constructed before use in solver code.
 *
 * Face normal convention:
 *  - interior face: normal points from owner cell to neighbor cell
 *  - physical boundary face: normal points outward from owner cell
 *  - MPI boundary face: normal points from owned cell to ghost neighbor
 */
class Mesh final {
public:
    Mesh() = default;
    explicit Mesh(int dim);

    /** @brief Spatial dimension. */
    [[nodiscard]] int GetDim() const;

    /** @brief Number of nodes. */
    [[nodiscard]] std::size_t GetNodeCount() const;

    /** @brief Number of faces. */
    [[nodiscard]] std::size_t GetFaceCount() const;

    /** @brief Number of cells. */
    [[nodiscard]] std::size_t GetCellCount() const;

    /** @brief Number of owned cells in MPI-local mesh. */
    [[nodiscard]] std::size_t GetOwnedCellCount() const;

    /** @brief Number of ghost cells in MPI-local mesh. */
    [[nodiscard]] std::size_t GetGhostCellCount() const;

    void SetOwnedCellCount(std::size_t count);
    void SetGhostCellCount(std::size_t count);

    /** @brief Mutable node storage. */
    [[nodiscard]] std::vector<Node>& Nodes();

    /** @brief Read-only node storage. */
    [[nodiscard]] const std::vector<Node>& Nodes() const;

    /** @brief Mutable face storage. */
    [[nodiscard]] std::vector<Face>& Faces();

    /** @brief Read-only face storage. */
    [[nodiscard]] const std::vector<Face>& Faces() const;

    /** @brief Mutable cell storage. */
    [[nodiscard]] std::vector<Cell>& Cells();

    /** @brief Read-only cell storage. */
    [[nodiscard]] const std::vector<Cell>& Cells() const;

    /** @brief Access one node by id. */
    [[nodiscard]] Node& GetNode(std::size_t node_id);

    /** @brief Access one node by id. */
    [[nodiscard]] const Node& GetNode(std::size_t node_id) const;

    /** @brief Access one face by id. */
    [[nodiscard]] Face& GetFace(std::size_t face_id);

    /** @brief Access one face by id. */
    [[nodiscard]] const Face& GetFace(std::size_t face_id) const;

    /** @brief Access one cell by local id. */
    [[nodiscard]] Cell& GetCell(std::size_t cell_id);

    /** @brief Access one cell by local id. */
    [[nodiscard]] const Cell& GetCell(std::size_t cell_id) const;

    /** @brief Check whether a face is physical boundary face. */
    [[nodiscard]] bool IsPhysicalBoundaryFace(std::size_t face_id) const;

    /** @brief Check whether a face is interior face. */
    [[nodiscard]] bool IsInternalFace(std::size_t face_id) const;

    /** @brief Check whether a face is MPI boundary face. */
    [[nodiscard]] bool IsMPIBoundaryFace(std::size_t face_id) const;

    /**
     * @brief Remove all mesh entities.
     * @details Dimension is preserved.
     */
    void Clear();

    /**
     * @brief Validate mesh dimension, ids, connectivity, geometry, and boundary tags.
     * @throws std::runtime_error if mesh data is inconsistent.
     */
    void Validate() const;

private:
    int dim_ = 0;

    std::vector<Node> nodes_;
    std::vector<Face> faces_;
    std::vector<Cell> cells_;

    std::size_t owned_cell_count_ = 0;
    std::size_t ghost_cell_count_ = 0;

    void ValidateDimension() const;
    void ValidateNodeIds() const;
    void ValidateFaceIds() const;
    void ValidateCellLocalIds() const;
    void ValidateFaceConnectivity() const;
    void ValidateCellConnectivity() const;
    void ValidateGeometry() const;
};

#endif  // MESH_HPP
