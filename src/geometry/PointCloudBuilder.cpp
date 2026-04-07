#include "geometry/PointCloudBuilder.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <vector>

double PointCloudBuilder::EdgeLength(const Point2D& a, const Point2D& b) {
    const double dx = b.x - a.x;
    const double dy = b.y - a.y;
    return std::sqrt(dx * dx + dy * dy);
}

double PointCloudBuilder::Orient2D(
    const Point2D& a,
    const Point2D& b,
    const Point2D& c
) {
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

double PointCloudBuilder::ComputeLoopSignedArea(
    const BoundaryLoop& loop,
    const std::vector<Point2D>& points
) {
    double twice_area = 0.0;

    for (std::size_t i = 0; i < loop.node_ids.size(); ++i) {
        const Point2D& a = points[loop.node_ids[i]];
        const Point2D& b = points[loop.node_ids[(i + 1) % loop.node_ids.size()]];
        twice_area += a.x * b.y - b.x * a.y;
    }

    return 0.5 * twice_area;
}

bool PointCloudBuilder::IsPointOnSegment(
    const Point2D& a,
    const Point2D& b,
    const Point2D& p,
    const double eps
) {
    if (std::abs(Orient2D(a, b, p)) > eps) {
        return false;
    }

    const double min_x = std::min(a.x, b.x) - eps;
    const double max_x = std::max(a.x, b.x) + eps;
    const double min_y = std::min(a.y, b.y) - eps;
    const double max_y = std::max(a.y, b.y) + eps;

    return p.x >= min_x && p.x <= max_x && p.y >= min_y && p.y <= max_y;
}

bool PointCloudBuilder::IsPointInsidePolygon(
    const Point2D& p,
    const BoundaryLoop& loop,
    const std::vector<Point2D>& points
) {
    bool inside = false;

    for (std::size_t i = 0, j = loop.node_ids.size() - 1; i < loop.node_ids.size(); j = i++) {
        const Point2D& a = points[loop.node_ids[i]];
        const Point2D& b = points[loop.node_ids[j]];

        if (IsPointOnSegment(a, b, p)) {
            return true;
        }

        const bool intersect =
            ((a.y > p.y) != (b.y > p.y)) &&
            (p.x < (b.x - a.x) * (p.y - a.y) / (b.y - a.y) + a.x);

        if (intersect) {
            inside = !inside;
        }
    }

    return inside;
}

std::size_t PointCloudBuilder::FindOuterLoopId(
    const std::vector<BoundaryLoop>& loops,
    const std::vector<Point2D>& points
) {
    if (loops.empty()) {
        throw std::runtime_error("PointCloudBuilder: loops are empty");
    }

    std::size_t outer_id = 0;
    double max_abs_area = std::abs(ComputeLoopSignedArea(loops[0], points));

    for (std::size_t i = 1; i < loops.size(); ++i) {
        const double abs_area = std::abs(ComputeLoopSignedArea(loops[i], points));
        if (abs_area > max_abs_area) {
            max_abs_area = abs_area;
            outer_id = i;
        }
    }

    return outer_id;
}

bool PointCloudBuilder::IsPointInsideDomain(
    const Point2D& p,
    const std::vector<BoundaryLoop>& loops,
    const std::vector<Point2D>& points
) {
    if (loops.empty()) {
        return false;
    }

    const std::size_t outer_id = FindOuterLoopId(loops, points);

    if (!IsPointInsidePolygon(p, loops[outer_id], points)) {
        return false;
    }

    for (std::size_t i = 0; i < loops.size(); ++i) {
        if (i == outer_id) {
            continue;
        }

        if (IsPointInsidePolygon(p, loops[i], points)) {
            return false;
        }
    }

    return true;
}

double PointCloudBuilder::DistancePointToSegment(
    const Point2D& p,
    const Point2D& a,
    const Point2D& b
) {
    const double abx = b.x - a.x;
    const double aby = b.y - a.y;
    const double apx = p.x - a.x;
    const double apy = p.y - a.y;

    const double ab2 = abx * abx + aby * aby;
    if (!(ab2 > 0.0)) {
        return EdgeLength(p, a);
    }

    const double t = Clamp((apx * abx + apy * aby) / ab2, 0.0, 1.0);

    const double qx = a.x + t * abx;
    const double qy = a.y + t * aby;

    const double dx = p.x - qx;
    const double dy = p.y - qy;
    return std::sqrt(dx * dx + dy * dy);
}

double PointCloudBuilder::DistanceToInnerBoundaries(
    const Point2D& p,
    const std::vector<Point2D>& points,
    const std::vector<BoundaryLoop>& loops,
    const std::size_t outer_loop_id
) {
    double best = std::numeric_limits<double>::max();

    for (std::size_t loop_id = 0; loop_id < loops.size(); ++loop_id) {
        if (loop_id == outer_loop_id) {
            continue;
        }

        const BoundaryLoop& loop = loops[loop_id];
        for (std::size_t i = 0; i < loop.node_ids.size(); ++i) {
            const Point2D& a = points[loop.node_ids[i]];
            const Point2D& b = points[loop.node_ids[(i + 1) % loop.node_ids.size()]];
            best = std::min(best, DistancePointToSegment(p, a, b));
        }
    }

    return best;
}

double PointCloudBuilder::ComputeMedianInnerBoundarySpacing(
    const std::vector<Point2D>& points,
    const std::vector<BoundaryLoop>& loops,
    const std::size_t outer_loop_id
) {
    std::vector<double> lengths;

    for (std::size_t loop_id = 0; loop_id < loops.size(); ++loop_id) {
        if (loop_id == outer_loop_id) {
            continue;
        }

        const BoundaryLoop& loop = loops[loop_id];
        for (std::size_t i = 0; i < loop.node_ids.size(); ++i) {
            const Point2D& a = points[loop.node_ids[i]];
            const Point2D& b = points[loop.node_ids[(i + 1) % loop.node_ids.size()]];
            const double len = EdgeLength(a, b);
            if (len > 0.0) {
                lengths.push_back(len);
            }
        }
    }

    if (lengths.empty()) {
        throw std::runtime_error("PointCloudBuilder: no inner-boundary segment lengths found");
    }

    std::sort(lengths.begin(), lengths.end());
    return lengths[lengths.size() / 2];
}

double PointCloudBuilder::ComputeLocalLoopSpacing(
    const BoundaryLoop& loop,
    const std::size_t local_id,
    const std::vector<Point2D>& points
) {
    const std::size_t n = loop.node_ids.size();

    const Point2D& p_prev = points[loop.node_ids[(local_id + n - 1) % n]];
    const Point2D& p_cur = points[loop.node_ids[local_id]];
    const Point2D& p_next = points[loop.node_ids[(local_id + 1) % n]];

    const double l_prev = EdgeLength(p_prev, p_cur);
    const double l_next = EdgeLength(p_cur, p_next);

    return 0.5 * (l_prev + l_next);
}

double PointCloudBuilder::Clamp(
    const double x,
    const double lo,
    const double hi
) {
    return std::max(lo, std::min(hi, x));
}

double PointCloudBuilder::ComputeSizeByDistanceLaw(
    const double distance_to_body,
    const double h_min,
    const double h_max,
    const double transition_radius,
    const double alpha
) {
    if (distance_to_body <= 0.0) {
        return h_min;
    }

    if (!(transition_radius > 0.0)) {
        return h_max;
    }

    const double q = Clamp(distance_to_body / transition_radius, 0.0, 1.0);
    return h_min + (h_max - h_min) * std::pow(q, alpha);
}

void PointCloudBuilder::RenumberPoints(std::vector<Point2D>& points) {
    for (std::size_t i = 0; i < points.size(); ++i) {
        points[i].id = i;
    }
}

double PointCloudBuilder::AppendInnerBoundaryLayerPoints(
    std::vector<Point2D>& points,
    const std::vector<BoundaryLoop>& loops,
    const std::size_t outer_loop_id,
    const BoundaryLayerSettings& settings,
    const double reference_spacing
) {
    if (settings.n_layers <= 0) {
        return 0.0;
    }

    const int n_layers = settings.n_layers;
    const double growth = (settings.growth > 1.0) ? settings.growth : 1.12;
    const double first_layer_height =
        (settings.first_layer_height > 0.0) ? settings.first_layer_height : 0.85 * reference_spacing;

    double max_total_thickness = 0.0;

    for (std::size_t loop_id = 0; loop_id < loops.size(); ++loop_id) {
        if (loop_id == outer_loop_id) {
            continue;
        }

        const BoundaryLoop& loop = loops[loop_id];
        const double signed_area = ComputeLoopSignedArea(loop, points);
        const bool is_ccw = signed_area > 0.0;

        double total_thickness = 0.0;
        for (int layer = 0; layer < n_layers; ++layer) {
            total_thickness += first_layer_height * std::pow(growth, static_cast<double>(layer));
        }
        max_total_thickness = std::max(max_total_thickness, total_thickness);

        for (std::size_t i = 0; i < loop.node_ids.size(); ++i) {
            const std::size_t id_prev =
                loop.node_ids[(i + loop.node_ids.size() - 1) % loop.node_ids.size()];
            const std::size_t id_cur = loop.node_ids[i];
            const std::size_t id_next =
                loop.node_ids[(i + 1) % loop.node_ids.size()];

            const Point2D& p_prev = points[id_prev];
            const Point2D& p_cur = points[id_cur];
            const Point2D& p_next = points[id_next];

            const double ds_local = ComputeLocalLoopSpacing(loop, i, points);
            const double local_scale = Clamp(ds_local / reference_spacing, 0.75, 1.35);
            const double local_first_height = first_layer_height * local_scale;

            double t1x = p_cur.x - p_prev.x;
            double t1y = p_cur.y - p_prev.y;
            double t2x = p_next.x - p_cur.x;
            double t2y = p_next.y - p_cur.y;

            double n1x = -t1y;
            double n1y = t1x;
            double n2x = -t2y;
            double n2y = t2x;

            const double norm1 = std::sqrt(n1x * n1x + n1y * n1y);
            const double norm2 = std::sqrt(n2x * n2x + n2y * n2y);

            if (!(norm1 > 0.0) || !(norm2 > 0.0)) {
                continue;
            }

            n1x /= norm1;
            n1y /= norm1;
            n2x /= norm2;
            n2y /= norm2;

            double nx = n1x + n2x;
            double ny = n1y + n2y;

            const double norm = std::sqrt(nx * nx + ny * ny);
            if (!(norm > 0.0)) {
                continue;
            }

            nx /= norm;
            ny /= norm;

            if (is_ccw) {
                nx *= -1.0;
                ny *= -1.0;
            }

            double cumulative_offset = 0.0;

            for (int layer = 0; layer < n_layers; ++layer) {
                const double h_layer =
                    local_first_height * std::pow(growth, static_cast<double>(layer));
                cumulative_offset += h_layer;

                Point2D candidate;
                candidate.x = p_cur.x + cumulative_offset * nx;
                candidate.y = p_cur.y + cumulative_offset * ny;
                candidate.is_boundary = false;

                if (!IsPointInsideDomain(candidate, loops, points)) {
                    continue;
                }

                bool too_close = false;
                const double min_allowed = 0.55 * h_layer;

                for (const Point2D& p : points) {
                    if (EdgeLength(candidate, p) < min_allowed) {
                        too_close = true;
                        break;
                    }
                }

                if (too_close) {
                    continue;
                }

                candidate.id = points.size();
                points.push_back(candidate);
            }
        }
    }

    RenumberPoints(points);
    return max_total_thickness;
}

void PointCloudBuilder::AppendPointsBySizeFunction(
    std::vector<Point2D>& points,
    const std::vector<BoundaryLoop>& loops,
    const std::size_t outer_loop_id,
    const double exclusion_distance,
    const SizeFunctionSettings& settings,
    const double reference_spacing
) {
    if (points.empty() || loops.empty()) {
        return;
    }

    double min_x = points.front().x;
    double max_x = points.front().x;
    double min_y = points.front().y;
    double max_y = points.front().y;

    for (const Point2D& p : points) {
        min_x = std::min(min_x, p.x);
        max_x = std::max(max_x, p.x);
        min_y = std::min(min_y, p.y);
        max_y = std::max(max_y, p.y);
    }

    const double box_size = std::max(max_x - min_x, max_y - min_y);

    double h_min =
        (settings.h_min > 0.0) ? settings.h_min : std::max(1.05 * reference_spacing, 1e-3);

    if (exclusion_distance > 0.0) {
        h_min = std::max(h_min, 0.55 * exclusion_distance / std::max(1, 4));
    }

    const double h_max =
        (settings.h_max > 0.0) ? settings.h_max : std::max(8.0 * h_min, box_size / 18.0);

    const double transition_radius =
        (settings.transition_radius > 0.0) ? settings.transition_radius : std::max(25.0 * h_min, box_size / 6.0);

    const double alpha = (settings.alpha > 0.0) ? settings.alpha : 1.8;

    const double candidate_base_step = 0.85 * h_min;
    const double candidate_dy = std::sqrt(3.0) * 0.5 * candidate_base_step;

    const double grid_cell_size = h_max;

    auto make_grid_key = [min_x, min_y, grid_cell_size](const double x, const double y) -> GridKey {
        return GridKey{
            static_cast<int>(std::floor((x - min_x) / grid_cell_size)),
            static_cast<int>(std::floor((y - min_y) / grid_cell_size))
        };
    };

    std::map<GridKey, std::vector<std::size_t>> spatial_grid;

    for (std::size_t i = 0; i < points.size(); ++i) {
        spatial_grid[make_grid_key(points[i].x, points[i].y)].push_back(i);
    }

    for (int row = 0;; ++row) {
        const double y = min_y + static_cast<double>(row) * candidate_dy;
        if (y > max_y) {
            break;
        }

        const double x_shift = (row % 2 == 0) ? 0.0 : 0.5 * candidate_base_step;

        for (double x = min_x + x_shift; x <= max_x; x += candidate_base_step) {
            Point2D candidate;
            candidate.x = x;
            candidate.y = y;
            candidate.is_boundary = false;

            if (!IsPointInsideDomain(candidate, loops, points)) {
                continue;
            }

            double distance_to_body = transition_radius;
            if (loops.size() > 1) {
                distance_to_body = DistanceToInnerBoundaries(candidate, points, loops, outer_loop_id);
            }

            if (exclusion_distance > 0.0 && distance_to_body < exclusion_distance) {
                continue;
            }

            const double h_local = ComputeSizeByDistanceLaw(
                distance_to_body,
                h_min,
                h_max,
                transition_radius,
                alpha
            );

            const double min_allowed = 0.78 * h_local;
            const int reach = static_cast<int>(
                std::ceil(std::max(min_allowed, grid_cell_size) / grid_cell_size)
            );

            const GridKey center_key = make_grid_key(candidate.x, candidate.y);

            bool ok = true;

            for (int dix = -reach; dix <= reach && ok; ++dix) {
                for (int diy = -reach; diy <= reach && ok; ++diy) {
                    const GridKey key{center_key.ix + dix, center_key.iy + diy};
                    const auto it = spatial_grid.find(key);
                    if (it == spatial_grid.end()) {
                        continue;
                    }

                    for (const std::size_t point_id : it->second) {
                        double existing_distance = transition_radius;
                        if (loops.size() > 1) {
                            existing_distance = DistanceToInnerBoundaries(
                                points[point_id],
                                points,
                                loops,
                                outer_loop_id
                            );
                        }

                        const double h_existing = ComputeSizeByDistanceLaw(
                            existing_distance,
                            h_min,
                            h_max,
                            transition_radius,
                            alpha
                        );

                        const double allowed = 0.72 * std::min(h_local, h_existing);

                        if (EdgeLength(candidate, points[point_id]) < allowed) {
                            ok = false;
                            break;
                        }
                    }
                }
            }

            if (!ok) {
                continue;
            }

            candidate.id = points.size();
            points.push_back(candidate);
            spatial_grid[make_grid_key(candidate.x, candidate.y)].push_back(candidate.id);
        }
    }

    RenumberPoints(points);
}

void PointCloudBuilder::BuildPointCloud(
    std::vector<Point2D>& points,
    const std::vector<BoundaryLoop>& loops,
    const BoundaryLayerSettings& boundary_layer_settings,
    const SizeFunctionSettings& size_function_settings
) {
    if (points.empty() || loops.empty()) {
        return;
    }

    const std::size_t outer_loop_id = FindOuterLoopId(loops, points);

    double reference_spacing = 1.0;
    if (loops.size() > 1) {
        reference_spacing = ComputeMedianInnerBoundarySpacing(points, loops, outer_loop_id);
    }
    else {
        double min_x = points.front().x;
        double max_x = points.front().x;
        double min_y = points.front().y;
        double max_y = points.front().y;

        for (const Point2D& p : points) {
            min_x = std::min(min_x, p.x);
            max_x = std::max(max_x, p.x);
            min_y = std::min(min_y, p.y);
            max_y = std::max(max_y, p.y);
        }

        reference_spacing = std::max(std::max(max_x - min_x, max_y - min_y) / 100.0, 1e-3);
    }

    const double exclusion_distance = AppendInnerBoundaryLayerPoints(
        points,
        loops,
        outer_loop_id,
        boundary_layer_settings,
        reference_spacing
    );

    AppendPointsBySizeFunction(
        points,
        loops,
        outer_loop_id,
        exclusion_distance,
        size_function_settings,
        reference_spacing
    );
}
