#ifndef FACE_HPP
#define FACE_HPP

#include <cstddef>
#include <vector>

enum class FaceKind {
    Interior,
    PhysicalBoundary,
    MPIBoundary
};

/**
 * @brief One mesh face separating two cells or a cell and a boundary/MPI interface.
 *
 * Normal is oriented from owner cell toward neighbor cell.
 * For physical boundary faces, normal is outward from owner cell.
 * For MPI boundary faces, neighbor_cell_id refers to a local ghost cell.
 */
struct Face final {
    static constexpr std::size_t k_invalid_cell_id = static_cast<std::size_t>(-1);

    std::size_t id = 0;

    /** @brief Indices of nodes forming this face. */
    std::vector<std::size_t> node_ids;

    /** @brief Owner cell id (local index in current mesh container). */
    std::size_t owner_cell_id = k_invalid_cell_id;

    /**
     * @brief Neighbor cell id for interior or MPI face.
     * @details For physical boundary face equals k_invalid_cell_id.
     */
    std::size_t neighbor_cell_id = k_invalid_cell_id;

    /** @brief Face classification. */
    FaceKind kind = FaceKind::Interior;

    /**
     * @brief Remote rank for MPI boundary faces.
     * @details For non-MPI faces equals -1.
     */
    int remote_rank = -1;

    /**
     * @brief Original/global id of remote real cell behind MPI boundary.
     * @details Optional helper field used during decomposition / halo setup.
     */
    std::size_t remote_cell_id = k_invalid_cell_id;

    /** @brief Face center coordinates. */
    double center_x = 0.0;
    double center_y = 0.0;
    double center_z = 0.0;

    /**
     * @brief Face measure.
     * @details In 2D this is edge length, in 3D this is face area.
     */
    double measure = 0.0;

    /**
     * @brief Unit normal directed from owner to neighbor.
     * @details For physical boundary faces this is outward from owner.
     */
    double normal_x = 0.0;
    double normal_y = 0.0;
    double normal_z = 0.0;

    /**
     * @brief Boundary tag for physical boundary faces.
     * @details For interior and MPI faces usually equals -1.
     */
    int boundary_tag = -1;

    /** @brief Check whether this face is a physical boundary face. */
    [[nodiscard]] bool IsPhysicalBoundary() const;

    /** @brief Check whether this face is an interior face. */
    [[nodiscard]] bool IsInternal() const;

    /** @brief Check whether this face is an MPI boundary face. */
    [[nodiscard]] bool IsMPIBoundary() const;
};

#endif  // FACE_HPP
