#include "geometry/DelaunayTriangulator.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {
    struct EdgeUse final {
        std::size_t count = 0;
        std::size_t first_a = 0;
        std::size_t first_b = 0;
    };
}

DelaunayTriangulator::EdgeKey::EdgeKey(const std::size_t node_a, const std::size_t node_b) {
    if (node_a < node_b) {
        a = node_a;
        b = node_b;
    }
    else {
        a = node_b;
        b = node_a;
    }
}

bool DelaunayTriangulator::EdgeKey::operator<(const EdgeKey& other) const {
    if (a != other.a) {
        return a < other.a;
    }
    return b < other.b;
}

double DelaunayTriangulator::Orient2D(
    const Point2D& a,
    const Point2D& b,
    const Point2D& c
) {
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

bool DelaunayTriangulator::IsPointInCircumcircle(
    const Point2D& a,
    const Point2D& b,
    const Point2D& c,
    const Point2D& p
) {
    const double ax = a.x - p.x;
    const double ay = a.y - p.y;
    const double bx = b.x - p.x;
    const double by = b.y - p.y;
    const double cx = c.x - p.x;
    const double cy = c.y - p.y;

    const double det =
        (ax * ax + ay * ay) * (bx * cy - by * cx) -
        (bx * bx + by * by) * (ax * cy - ay * cx) +
        (cx * cx + cy * cy) * (ax * by - ay * bx);

    const double orient = Orient2D(a, b, c);
    const double eps = 1e-14;

    if (orient > 0.0) {
        return det > eps;
    }

    return det < -eps;
}

DelaunayTriangulator::Triangle DelaunayTriangulator::MakeCcwTriangle(
    const std::size_t a,
    const std::size_t b,
    const std::size_t c,
    const std::vector<Point2D>& points
) {
    Triangle triangle{{a, b, c}};

    if (Orient2D(points[a], points[b], points[c]) < 0.0) {
        std::swap(triangle.node_ids[1], triangle.node_ids[2]);
    }

    return triangle;
}

bool DelaunayTriangulator::TriangleContainsSuperNode(
    const Triangle& triangle,
    const std::size_t super0,
    const std::size_t super1,
    const std::size_t super2
) {
    return triangle.node_ids[0] == super0 ||
        triangle.node_ids[1] == super0 ||
        triangle.node_ids[2] == super0 ||
        triangle.node_ids[0] == super1 ||
        triangle.node_ids[1] == super1 ||
        triangle.node_ids[2] == super1 ||
        triangle.node_ids[0] == super2 ||
        triangle.node_ids[1] == super2 ||
        triangle.node_ids[2] == super2;
}

std::vector<DelaunayTriangulator::Triangle> DelaunayTriangulator::Triangulate(
    const std::vector<Point2D>& input_points
) {
    if (input_points.size() < 3) {
        throw std::runtime_error("DelaunayTriangulator: at least 3 points are required");
    }

    std::vector<Point2D> points = input_points;

    double min_x = points[0].x;
    double max_x = points[0].x;
    double min_y = points[0].y;
    double max_y = points[0].y;

    for (const Point2D& p : points) {
        min_x = std::min(min_x, p.x);
        max_x = std::max(max_x, p.x);
        min_y = std::min(min_y, p.y);
        max_y = std::max(max_y, p.y);
    }

    const double dx = max_x - min_x;
    const double dy = max_y - min_y;
    const double d = std::max(dx, dy);
    const double cx = 0.5 * (min_x + max_x);
    const double cy = 0.5 * (min_y + max_y);

    const std::size_t super0 = points.size();
    const std::size_t super1 = points.size() + 1;
    const std::size_t super2 = points.size() + 2;

    points.push_back(Point2D{super0, cx - 2.0 * d, cy - d});
    points.push_back(Point2D{super1, cx + 2.0 * d, cy - d});
    points.push_back(Point2D{super2, cx, cy + 2.0 * d});

    std::vector<Triangle> triangles;
    triangles.push_back(MakeCcwTriangle(super0, super1, super2, points));

    for (std::size_t point_id = 0; point_id < input_points.size(); ++point_id) {
        std::vector<std::size_t> bad_triangles;
        bad_triangles.reserve(triangles.size());

        for (std::size_t tri_id = 0; tri_id < triangles.size(); ++tri_id) {
            const Triangle& triangle = triangles[tri_id];

            const Point2D& a = points[triangle.node_ids[0]];
            const Point2D& b = points[triangle.node_ids[1]];
            const Point2D& c = points[triangle.node_ids[2]];
            const Point2D& p = points[point_id];

            if (IsPointInCircumcircle(a, b, c, p)) {
                bad_triangles.push_back(tri_id);
            }
        }

        std::map<EdgeKey, EdgeUse> polygon;
        std::vector<bool> is_bad(triangles.size(), false);

        for (const std::size_t tri_id : bad_triangles) {
            is_bad[tri_id] = true;

            const Triangle& triangle = triangles[tri_id];
            const std::array<std::pair<std::size_t, std::size_t>, 3> edges = {
                {
                    {triangle.node_ids[0], triangle.node_ids[1]},
                    {triangle.node_ids[1], triangle.node_ids[2]},
                    {triangle.node_ids[2], triangle.node_ids[0]}
                }
            };

            for (const auto& [ea, eb] : edges) {
                const EdgeKey key(ea, eb);
                auto& use = polygon[key];
                use.count += 1;
                use.first_a = ea;
                use.first_b = eb;
            }
        }

        std::vector<Triangle> survivors;
        survivors.reserve(triangles.size());

        for (std::size_t tri_id = 0; tri_id < triangles.size(); ++tri_id) {
            if (!is_bad[tri_id]) {
                survivors.push_back(triangles[tri_id]);
            }
        }

        triangles = std::move(survivors);

        for (const auto& [edge, use] : polygon) {
            if (use.count == 1) {
                triangles.push_back(
                    MakeCcwTriangle(use.first_a, use.first_b, point_id, points)
                );
            }
        }
    }

    std::vector<Triangle> result;
    result.reserve(triangles.size());

    for (const Triangle& triangle : triangles) {
        if (!TriangleContainsSuperNode(triangle, super0, super1, super2)) {
            result.push_back(triangle);
        }
    }

    return result;
}
