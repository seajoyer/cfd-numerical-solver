#ifndef CELL_HPP
#define CELL_HPP

#include <cstddef>
#include <vector>

/**
 * @brief One control volume of an unstructured mesh.
 */
struct Cell final {
    std::size_t id = 0;

    /** @brief Indices of nodes belonging to the cell. */
    std::vector<std::size_t> node_ids;

    /** @brief Indices of faces bounding the cell. */
    std::vector<std::size_t> face_ids;

    /** @brief Cell center coordinates. */
    double center_x = 0.0;
    double center_y = 0.0;
    double center_z = 0.0;

    /**
     * @brief Cell measure.
     * @details In 2D this is area, in 3D this is volume.
     */
    double volume = 0.0;
};

#endif  // CELL_HPP
