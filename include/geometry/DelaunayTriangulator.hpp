#ifndef DELAUNAYTRIANGULATOR_HPP
#define DELAUNAYTRIANGULATOR_HPP

#include <array>
#include <cstddef>
#include <vector>

/**
 * @brief Pure 2D Bowyer-Watson Delaunay triangulator.
 *
 * Input:
 *  - arbitrary 2D point set
 *
 * Output:
 *  - triangle connectivity
 *
 * This class does not know anything about:
 *  - boundaries
 *  - loops
 *  - physical tags
 *  - generic Mesh conversion
 */
class DelaunayTriangulator final {
public:
    /**
     * @brief One 2D point.
     */
    struct Point2D final {
        std::size_t id = 0;
        double x = 0.0;
        double y = 0.0;
    };

    /**
     * @brief One triangle.
     */
    struct Triangle final {
        std::array<std::size_t, 3> node_ids{};
    };

    /**
     * @brief Build Delaunay triangulation of the given points.
     */
    [[nodiscard]] static std::vector<Triangle> Triangulate(
        const std::vector<Point2D>& points
    );

private:
    /**
     * @brief Undirected edge key.
     */
    struct EdgeKey final {
        std::size_t a = 0;
        std::size_t b = 0;

        EdgeKey() = default;
        EdgeKey(std::size_t node_a, std::size_t node_b);

        [[nodiscard]] bool operator<(const EdgeKey& other) const;
    };

    [[nodiscard]] static double Orient2D(
        const Point2D& a,
        const Point2D& b,
        const Point2D& c
    );

    [[nodiscard]] static bool IsPointInCircumcircle(
        const Point2D& a,
        const Point2D& b,
        const Point2D& c,
        const Point2D& p
    );

    [[nodiscard]] static Triangle MakeCcwTriangle(
        std::size_t a,
        std::size_t b,
        std::size_t c,
        const std::vector<Point2D>& points
    );

    [[nodiscard]] static bool TriangleContainsSuperNode(
        const Triangle& triangle,
        std::size_t super0,
        std::size_t super1,
        std::size_t super2
    );
};

#endif  // DELAUNAYTRIANGULATOR_HPP