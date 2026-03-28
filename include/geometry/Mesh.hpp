#ifndef MESH_HPP
#define MESH_HPP

#include <cstddef>
#include <stdexcept>
#include <vector>

#include "geometry/Cell.hpp"
#include "geometry/Face.hpp"
#include "geometry/Node.hpp"

/**
 * @brief Single-process unstructured mesh owner for cell-centered finite volume methods.
 *
 * Mesh owns:
 *  - nodes
 *  - faces
 *  - cells
 *
 * Geometry and connectivity must already be constructed before use in solver code.
 * Face normal convention:
 *  - internal face: normal points from owner cell to neighbor cell
 *  - boundary face: normal points outward from owner cell
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

    /** @brief Access one cell by id. */
    [[nodiscard]] Cell& GetCell(std::size_t cell_id);

    /** @brief Access one cell by id. */
    [[nodiscard]] const Cell& GetCell(std::size_t cell_id) const;

    /** @brief Check whether a face is boundary face. */
    [[nodiscard]] bool IsBoundaryFace(std::size_t face_id) const;

    /** @brief Check whether a face is internal face. */
    [[nodiscard]] bool IsInternalFace(std::size_t face_id) const;

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

    void ValidateDimension() const;
    void ValidateNodeIds() const;
    void ValidateFaceIds() const;
    void ValidateCellIds() const;
    void ValidateFaceConnectivity() const;
    void ValidateCellConnectivity() const;
    void ValidateGeometry() const;
};

#endif  // MESH_HPP
