#ifndef NODE_HPP
#define NODE_HPP

#include <cstddef>

/**
 * @brief Geometric node of an unstructured mesh.
 */
struct Node final {
    std::size_t id = 0;

    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

#endif  // NODE_HPP
