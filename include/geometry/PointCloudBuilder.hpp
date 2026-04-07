#ifndef POINTCLOUDBUILDER_HPP
#define POINTCLOUDBUILDER_HPP

#include <cstddef>
#include <vector>


/**
 * @brief One 2D point of a point cloud.
 */
struct Point2D final {
    std::size_t id = 0;
    double x = 0.0;
    double y = 0.0;
    bool is_boundary = false;
};

/**
 * @brief One closed polygonal loop.
 */
struct BoundaryLoop final {
    std::vector<std::size_t> node_ids;
};

/**
 * @brief Utility builder for preparing 2D point clouds before triangulation.
 *
 * Strategy:
 *  - boundary points are provided externally
 *  - optional boundary-layer-like point rows are added near inner boundaries
 *  - remaining interior points are added using a smooth size law:
 *
 *    h(r) = h_min + (h_max - h_min) * (r / R)^alpha
 *
 * where:
 *  - r is distance to inner object
 *  - R is transition radius
 */
class PointCloudBuilder final {
public:
    /**
     * @brief Parameters of near-wall point layers.
     *
     * If n_layers == 0, no near-wall layers are generated.
     *
     * first_layer_height:
     *  - if <= 0, computed automatically from boundary spacing
     *
     * growth:
     *  - geometric growth factor between layers
     *
     * n_layers:
     *  - number of near-wall layers
     */
    struct BoundaryLayerSettings final {
        double first_layer_height = 0.0;
        double growth = 1.12;
        int n_layers = 7;
    };

    /**
     * @brief Parameters of smooth size law in the remaining domain.
     *
     * h_min:
     *  - if <= 0, computed automatically
     *
     * h_max:
     *  - if <= 0, computed automatically
     *
     * transition_radius:
     *  - if <= 0, computed automatically
     *
     * alpha:
     *  - exponent in h(r)
     */
    struct SizeFunctionSettings final {
        double h_min = 0.0;
        double h_max = 0.0;
        double transition_radius = 0.0;
        double alpha = 1.8;
    };

    /**
     * @brief Append interior points using:
     *  - optional near-wall layers
     *  - smooth size-law-based filling in the remaining domain
     *
     * @details The largest loop is treated as outer boundary.
     */
    static void BuildPointCloud(
        std::vector<Point2D>& points,
        const std::vector<BoundaryLoop>& loops,
        const BoundaryLayerSettings& boundary_layer_settings = BoundaryLayerSettings{
            .first_layer_height = 0.0, .growth = 1.12, .n_layers = 7
        },
        const SizeFunctionSettings& size_function_settings = SizeFunctionSettings{
            .h_min = 0.0, .h_max = 0.0, .transition_radius = 0.0, .alpha = 1.8
        }
    );

private:
    struct GridKey final {
        int ix = 0;
        int iy = 0;

        [[nodiscard]] bool operator<(const GridKey& other) const {
            if (ix != other.ix) {
                return ix < other.ix;
            }
            return iy < other.iy;
        }
    };

    [[nodiscard]] static double EdgeLength(
        const Point2D& a,
        const Point2D& b
    );

    [[nodiscard]] static double Orient2D(
        const Point2D& a,
        const Point2D& b,
        const Point2D& c
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

    [[nodiscard]] static std::size_t FindOuterLoopId(
        const std::vector<BoundaryLoop>& loops,
        const std::vector<Point2D>& points
    );

    [[nodiscard]] static double DistancePointToSegment(
        const Point2D& p,
        const Point2D& a,
        const Point2D& b
    );

    [[nodiscard]] static double DistanceToInnerBoundaries(
        const Point2D& p,
        const std::vector<Point2D>& points,
        const std::vector<BoundaryLoop>& loops,
        std::size_t outer_loop_id
    );

    [[nodiscard]] static double ComputeMedianInnerBoundarySpacing(
        const std::vector<Point2D>& points,
        const std::vector<BoundaryLoop>& loops,
        std::size_t outer_loop_id
    );

    [[nodiscard]] static double ComputeLocalLoopSpacing(
        const BoundaryLoop& loop,
        std::size_t local_id,
        const std::vector<Point2D>& points
    );

    [[nodiscard]] static double Clamp(
        double x,
        double lo,
        double hi
    );

    [[nodiscard]] static double ComputeSizeByDistanceLaw(
        double distance_to_body,
        double h_min,
        double h_max,
        double transition_radius,
        double alpha
    );

    static void RenumberPoints(std::vector<Point2D>& points);

    [[nodiscard]] static double AppendInnerBoundaryLayerPoints(
        std::vector<Point2D>& points,
        const std::vector<BoundaryLoop>& loops,
        std::size_t outer_loop_id,
        const BoundaryLayerSettings& settings,
        double reference_spacing
    );

    static void AppendPointsBySizeFunction(
        std::vector<Point2D>& points,
        const std::vector<BoundaryLoop>& loops,
        std::size_t outer_loop_id,
        double exclusion_distance,
        const SizeFunctionSettings& settings,
        double reference_spacing
    );
};

#endif  // POINTCLOUDBUILDER_HPP
