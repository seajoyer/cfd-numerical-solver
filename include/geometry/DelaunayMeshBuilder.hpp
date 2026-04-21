#ifndef DELAUNAYMESHBUILDER_HPP
#define DELAUNAYMESHBUILDER_HPP

#include <array>
#include <cstddef>
#include <map>
#include <string>
#include <vector>

#include "geometry/PointCloudBuilder.hpp"
#include "geometry/Mesh.hpp"


struct Point2D;
struct BoundaryLoop;
// CHECK: BUILD_FACES
/**
 * @brief Builder that constructs a 2D mesh from .geo geometry using a point set
 *        and an in-house Delaunay triangulator.
 *
 * Responsibilities:
 *  - read physical curves from Gmsh geometry
 *  - sample boundary points
 *  - reconstruct boundary loops
 *  - prepare point cloud for triangulation
 *  - triangulate using DelaunayTriangulator
 *  - remove triangles outside domain
 *  - convert result to generic Mesh
 *
 * This class intentionally keeps point-cloud preparation separate from the
 * triangulation kernel so future PointCloudBuilder strategies can be integrated.
 */
class DelaunayMeshBuilder final {
public:
    /**
     * @brief Build 2D mesh from .geo file.
     * @param file_path Path to .geo file.
     * @param dim Spatial dimension. Must be 2.
     */
    [[nodiscard]] static Mesh BuildFromGeoFile(const std::string& file_path, int dim);

private:

    /**
     * @brief One sampled boundary segment with physical tag.
     */
    struct BoundarySegment final {
        std::size_t a = 0;
        std::size_t b = 0;
        int boundary_tag = -1;
    };

    /**
     * @brief One physical curve read from Gmsh model.
     */
    struct GeometricCurve final {
        int curve_tag = -1;
        int physical_tag = -1;
    };

    /**
     * @brief Internal triangle type used after conversion from triangulator.
     */
    struct Triangle final {
        std::array<std::size_t, 3> node_ids{};
    };

    /**
     * @brief Undirected edge key.
     */
    struct EdgeKey final {
        std::size_t a = 0;
        std::size_t b = 0;

        EdgeKey() = default;
        EdgeKey(std::size_t node_a, std::size_t node_b);

        [[nodiscard]] bool operator<(const EdgeKey& other) const;
        [[nodiscard]] bool operator==(const EdgeKey& other) const;
    };

    /**
     * @brief RAII helper for local Gmsh session ownership.
     */
    class GmshSessionGuard final {
    public:
        explicit GmshSessionGuard(bool finalize_on_destroy);
        ~GmshSessionGuard();

        GmshSessionGuard(const GmshSessionGuard&) = delete;
        auto operator=(const GmshSessionGuard&) -> GmshSessionGuard& = delete;

    private:
        bool finalize_on_destroy_ = false;
    };

    static void ValidateInputDimension(int dim);

    static void ReadGeometryFromGeoFile(
        const std::string& file_path,
        std::vector<Point2D>& points,
        std::vector<BoundarySegment>& boundary_segments
    );

    static void ReadPhysicalCurves(std::vector<GeometricCurve>& curves);

    static void SamplePhysicalCurve(
        int curve_tag,
        int physical_tag,
        std::vector<Point2D>& points,
        std::vector<BoundarySegment>& boundary_segments,
        std::map<EdgeKey, int>& existing_boundary_edges
    );

    [[nodiscard]] static std::vector<BoundaryLoop> BuildBoundaryLoops(
        const std::vector<BoundarySegment>& boundary_segments
    );

    [[nodiscard]] static std::vector<Triangle> FilterTrianglesInsideDomain(
        const std::vector<Point2D>& points,
        const std::vector<Triangle>& triangles,
        const std::vector<BoundaryLoop>& loops
    );

    [[nodiscard]] static Mesh BuildMeshFromTriangles(
        const std::vector<Point2D>& points,
        const std::vector<Triangle>& triangles,
        const std::vector<BoundarySegment>& boundary_segments
    );

    [[nodiscard]] static double ComputeLoopSignedArea(
        const BoundaryLoop& loop,
        const std::vector<Point2D>& points
    );

    [[nodiscard]] static bool IsPointOnSegment(
        const Point2D& a,
        const Point2D& b,
        const Point2D& p,
        double eps = 1e-12
    );

    [[nodiscard]] static bool IsPointInsidePolygon(
        const Point2D& p,
        const BoundaryLoop& loop,
        const std::vector<Point2D>& points
    );

    [[nodiscard]] static bool IsPointInsideDomain(
        const Point2D& p,
        const std::vector<BoundaryLoop>& loops,
        const std::vector<Point2D>& points
    );

    [[nodiscard]] static double EdgeLength(
        const Point2D& a,
        const Point2D& b
    );

    [[nodiscard]] static double Orient2D(
        const Point2D& a,
        const Point2D& b,
        const Point2D& c
    );

    [[nodiscard]] static double DistancePointToSegment(
        const Point2D& p,
        const Point2D& a,
        const Point2D& b
    );

    [[nodiscard]] static int FindBoundaryTagForFace(
        const Point2D& a,
        const Point2D& b,
        const std::vector<Point2D>& points,
        const std::vector<BoundarySegment>& boundary_segments
    );

    static void FinalizeFaceNormals(Mesh& mesh);
};

#endif  // DELAUNAYMESHBUILDER_HPP
