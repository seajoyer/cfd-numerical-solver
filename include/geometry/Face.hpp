#ifndef FACE_HPP
#define FACE_HPP

#include <cstddef>
#include <vector>

/**
 * @brief One mesh face separating two cells or a cell and a physical boundary.
 *
 * Normal is oriented from owner cell toward neighbor cell.
 * For boundary faces, normal is outward from owner cell.
 */
struct Face final {
    static constexpr std::size_t k_invalid_cell_id = static_cast<std::size_t>(-1);

    std::size_t id = 0;

    /** @brief Indices of nodes forming this face. */
    std::vector<std::size_t> node_ids;

    /** @brief Owner cell id. */
    std::size_t owner_cell_id = k_invalid_cell_id;

    /**
     * @brief Neighbor cell id for internal face.
     * @details For boundary face equals k_invalid_cell_id.
     */
    std::size_t neighbor_cell_id = k_invalid_cell_id;

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
     * @details For boundary faces this is outward from owner.
     */
    double normal_x = 0.0;
    double normal_y = 0.0;
    double normal_z = 0.0;

    /**
     * @brief Boundary tag for boundary faces.
     * @details For internal faces usually equals -1.
     */
    int boundary_tag = -1;

    /** @brief Check whether this face is a boundary face. */
    [[nodiscard]] bool IsBoundary() const;

    /** @brief Check whether this face is an internal face. */
    [[nodiscard]] bool IsInternal() const;
};

#endif  // FACE_HPP
