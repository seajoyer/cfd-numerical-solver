#include "geometry/DelaunayMeshBuilder.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

#include <gmsh.h>

#include "geometry/Cell.hpp"
#include "geometry/DelaunayTriangulator.hpp"
#include "geometry/Face.hpp"
#include "geometry/Node.hpp"

namespace {
    [[nodiscard]] std::size_t FindOrAddPoint(
        const double x,
        const double y,
        const bool is_boundary,
        std::vector<Point2D>& points,
        const double tol = 1e-12
    ) {
        for (std::size_t i = 0; i < points.size(); ++i) {
            const double dx = points[i].x - x;
            const double dy = points[i].y - y;
            if (dx * dx + dy * dy <= tol * tol) {
                if (is_boundary) {
                    points[i].is_boundary = true;
                }
                return i;
            }
        }

        Point2D point;
        point.id = points.size();
        point.x = x;
        point.y = y;
        point.is_boundary = is_boundary;
        points.push_back(point);
        return point.id;
    }
}

DelaunayMeshBuilder::EdgeKey::EdgeKey(const std::size_t node_a, const std::size_t node_b) {
    if (node_a < node_b) {
        a = node_a;
        b = node_b;
    }
    else {
        a = node_b;
        b = node_a;
    }
}

bool DelaunayMeshBuilder::EdgeKey::operator<(const EdgeKey& other) const {
    if (a != other.a) {
        return a < other.a;
    }
    return b < other.b;
}

bool DelaunayMeshBuilder::EdgeKey::operator==(const EdgeKey& other) const {
    return a == other.a && b == other.b;
}

DelaunayMeshBuilder::GmshSessionGuard::GmshSessionGuard(const bool finalize_on_destroy)
    : finalize_on_destroy_(finalize_on_destroy) {}

DelaunayMeshBuilder::GmshSessionGuard::~GmshSessionGuard() {
    if (finalize_on_destroy_) {
        gmsh::finalize();
    }
}

void DelaunayMeshBuilder::ValidateInputDimension(const int dim) {
    if (dim != 2) {
        throw std::invalid_argument("DelaunayMeshBuilder: only dim=2 is supported");
    }
}

double DelaunayMeshBuilder::Orient2D(
    const Point2D& a,
    const Point2D& b,
    const Point2D& c
) {
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

double DelaunayMeshBuilder::EdgeLength(const Point2D& a, const Point2D& b) {
    const double dx = b.x - a.x;
    const double dy = b.y - a.y;
    return std::sqrt(dx * dx + dy * dy);
}

void DelaunayMeshBuilder::ReadPhysicalCurves(std::vector<GeometricCurve>& curves) {
    curves.clear();

    std::vector<std::pair<int, int>> physical_groups;
    gmsh::model::getPhysicalGroups(physical_groups, 1);

    for (const auto& [dim, physical_tag] : physical_groups) {
        if (dim != 1) {
            continue;
        }

        std::vector<int> entity_tags;
        gmsh::model::getEntitiesForPhysicalGroup(1, physical_tag, entity_tags);

        for (const int curve_tag : entity_tags) {
            curves.push_back(GeometricCurve{
                .curve_tag = curve_tag,
                .physical_tag = physical_tag
            });
        }
    }

    if (curves.empty()) {
        throw std::runtime_error("DelaunayMeshBuilder: no physical curves found");
    }
}

void DelaunayMeshBuilder::SamplePhysicalCurve(
    const int curve_tag,
    const int physical_tag,
    std::vector<Point2D>& points,
    std::vector<BoundarySegment>& boundary_segments,
    std::map<EdgeKey, int>& existing_boundary_edges
) {
    std::vector<double> param_min;
    std::vector<double> param_max;
    gmsh::model::getParametrizationBounds(1, curve_tag, param_min, param_max);

    if (param_min.empty() || param_max.empty()) {
        throw std::runtime_error("DelaunayMeshBuilder: empty curve parametrization bounds");
    }

    const double u0 = param_min[0];
    const double u1 = param_max[0];

    const int integration_samples = 200;
    double length = 0.0;

    {
        std::vector<double> prev_coord;
        gmsh::model::getValue(1, curve_tag, {u0}, prev_coord);

        for (int i = 1; i <= integration_samples; ++i) {
            const double ui = u0 + (u1 - u0) * static_cast<double>(i) /
                static_cast<double>(integration_samples);

            std::vector<double> coord;
            gmsh::model::getValue(1, curve_tag, {ui}, coord);

            const double dx = coord[0] - prev_coord[0];
            const double dy = coord[1] - prev_coord[1];
            length += std::sqrt(dx * dx + dy * dy);

            prev_coord = coord;
        }
    }

    if (!(length > 0.0)) {
        throw std::runtime_error("DelaunayMeshBuilder: curve length must be positive");
    }

    const double target_step = std::max(length / 40.0, 2e-3);
    const int n_segments = std::max(2, static_cast<int>(std::ceil(length / target_step)));

    std::vector<std::size_t> ids;
    ids.reserve(static_cast<std::size_t>(n_segments + 1));

    for (int i = 0; i <= n_segments; ++i) {
        const double s = static_cast<double>(i) / static_cast<double>(n_segments);
        const double u = u0 + s * (u1 - u0);

        std::vector<double> coord;
        gmsh::model::getValue(1, curve_tag, {u}, coord);

        ids.push_back(FindOrAddPoint(coord[0], coord[1], true, points));
    }

    for (std::size_t i = 1; i < ids.size(); ++i) {
        if (ids[i - 1] == ids[i]) {
            continue;
        }

        const EdgeKey key(ids[i - 1], ids[i]);
        if (existing_boundary_edges.find(key) != existing_boundary_edges.end()) {
            continue;
        }

        boundary_segments.push_back(BoundarySegment{
            .a = ids[i - 1],
            .b = ids[i],
            .boundary_tag = physical_tag
        });

        existing_boundary_edges[key] = physical_tag;
    }
}

void DelaunayMeshBuilder::ReadGeometryFromGeoFile(
    const std::string& file_path,
    std::vector<Point2D>& points,
    std::vector<BoundarySegment>& boundary_segments
) {
    const bool was_initialized = gmsh::isInitialized();
    if (!was_initialized) {
        gmsh::initialize();
    }

    GmshSessionGuard session_guard(!was_initialized);

    gmsh::clear();
    gmsh::open(file_path);
    gmsh::model::geo::synchronize();

    points.clear();
    boundary_segments.clear();

    std::vector<GeometricCurve> curves;
    ReadPhysicalCurves(curves);

    std::map<EdgeKey, int> existing_boundary_edges;

    for (const GeometricCurve& curve : curves) {
        SamplePhysicalCurve(
            curve.curve_tag,
            curve.physical_tag,
            points,
            boundary_segments,
            existing_boundary_edges
        );
    }

    if (points.size() < 3) {
        throw std::runtime_error("DelaunayMeshBuilder: too few boundary points");
    }

    if (boundary_segments.empty()) {
        throw std::runtime_error("DelaunayMeshBuilder: no boundary segments sampled");
    }
}

std::vector<BoundaryLoop> DelaunayMeshBuilder::BuildBoundaryLoops(
    const std::vector<BoundarySegment>& boundary_segments
) {
    std::map<std::size_t, std::vector<std::size_t>> adjacency;

    for (const BoundarySegment& segment : boundary_segments) {
        adjacency[segment.a].push_back(segment.b);
        adjacency[segment.b].push_back(segment.a);
    }

    for (const auto& [node_id, neighbors] : adjacency) {
        if (neighbors.size() != 2) {
            throw std::runtime_error(
                "DelaunayMeshBuilder: boundary graph is not a collection of simple loops"
            );
        }
    }

    std::set<EdgeKey> visited_edges;
    std::vector<BoundaryLoop> loops;

    for (const BoundarySegment& segment : boundary_segments) {
        const EdgeKey start_edge(segment.a, segment.b);
        if (visited_edges.find(start_edge) != visited_edges.end()) {
            continue;
        }

        BoundaryLoop loop;
        loop.node_ids.push_back(segment.a);

        std::size_t prev = segment.a;
        std::size_t cur = segment.b;

        visited_edges.insert(EdgeKey(prev, cur));

        while (cur != loop.node_ids.front()) {
            loop.node_ids.push_back(cur);

            const auto& neighbors = adjacency[cur];
            const std::size_t next = (neighbors[0] == prev) ? neighbors[1] : neighbors[0];

            prev = cur;
            cur = next;

            visited_edges.insert(EdgeKey(prev, cur));
        }

        if (loop.node_ids.size() < 3) {
            throw std::runtime_error("DelaunayMeshBuilder: invalid boundary loop");
        }

        loops.push_back(loop);
    }

    if (loops.empty()) {
        throw std::runtime_error("DelaunayMeshBuilder: no boundary loops reconstructed");
    }

    return loops;
}

double DelaunayMeshBuilder::ComputeLoopSignedArea(
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

bool DelaunayMeshBuilder::IsPointOnSegment(
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

bool DelaunayMeshBuilder::IsPointInsidePolygon(
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

bool DelaunayMeshBuilder::IsPointInsideDomain(
    const Point2D& p,
    const std::vector<BoundaryLoop>& loops,
    const std::vector<Point2D>& points
) {
    if (loops.empty()) {
        return false;
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

std::vector<DelaunayMeshBuilder::Triangle> DelaunayMeshBuilder::FilterTrianglesInsideDomain(
    const std::vector<Point2D>& points,
    const std::vector<Triangle>& triangles,
    const std::vector<BoundaryLoop>& loops
) {
    std::vector<Triangle> filtered;
    filtered.reserve(triangles.size());

    for (const Triangle& triangle : triangles) {
        Point2D center;
        center.x =
        (points[triangle.node_ids[0]].x +
            points[triangle.node_ids[1]].x +
            points[triangle.node_ids[2]].x) / 3.0;
        center.y =
        (points[triangle.node_ids[0]].y +
            points[triangle.node_ids[1]].y +
            points[triangle.node_ids[2]].y) / 3.0;

        if (IsPointInsideDomain(center, loops, points)) {
            filtered.push_back(triangle);
        }
    }

    if (filtered.empty()) {
        throw std::runtime_error("DelaunayMeshBuilder: no triangles remained inside domain");
    }

    return filtered;
}

double DelaunayMeshBuilder::DistancePointToSegment(
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

    const double t = std::max(0.0, std::min(1.0, (apx * abx + apy * aby) / ab2));

    const double qx = a.x + t * abx;
    const double qy = a.y + t * aby;

    const double dx = p.x - qx;
    const double dy = p.y - qy;
    return std::sqrt(dx * dx + dy * dy);
}

int DelaunayMeshBuilder::FindBoundaryTagForFace(
    const Point2D& a,
    const Point2D& b,
    const std::vector<Point2D>& points,
    const std::vector<BoundarySegment>& boundary_segments
) {
    const Point2D mid{
        .id = 0,
        .x = 0.5 * (a.x + b.x),
        .y = 0.5 * (a.y + b.y),
        .is_boundary = false
    };

    const double face_len = EdgeLength(a, b);
    int best_tag = -1;
    double best_dist = std::numeric_limits<double>::max();

    for (const BoundarySegment& segment : boundary_segments) {
        const Point2D& s0 = points[segment.a];
        const Point2D& s1 = points[segment.b];

        const double d = DistancePointToSegment(mid, s0, s1);
        const double tol = std::max(0.35 * face_len, 1e-8);

        if (d < tol && d < best_dist) {
            best_dist = d;
            best_tag = segment.boundary_tag;
        }
    }

    return best_tag;
}

Mesh DelaunayMeshBuilder::BuildMeshFromTriangles(
    const std::vector<Point2D>& points,
    const std::vector<Triangle>& triangles,
    const std::vector<BoundarySegment>& boundary_segments
) {
    Mesh mesh(2);

    auto& nodes = mesh.Nodes();
    auto& cells = mesh.Cells();
    auto& faces = mesh.Faces();

    nodes.reserve(points.size());
    for (const Point2D& point : points) {
        Node node;
        node.id = nodes.size();
        node.x = point.x;
        node.y = point.y;
        node.z = 0.0;
        nodes.push_back(node);
    }

    std::map<EdgeKey, std::size_t> edge_to_face_id;

    for (const Triangle& triangle : triangles) {
        Cell cell;
        cell.id = cells.size();
        cell.node_ids = {
            triangle.node_ids[0],
            triangle.node_ids[1],
            triangle.node_ids[2]
        };

        const Point2D& a = points[triangle.node_ids[0]];
        const Point2D& b = points[triangle.node_ids[1]];
        const Point2D& c = points[triangle.node_ids[2]];

        cell.center_x = (a.x + b.x + c.x) / 3.0;
        cell.center_y = (a.y + b.y + c.y) / 3.0;
        cell.center_z = 0.0;
        cell.volume = 0.5 * std::abs(Orient2D(a, b, c));

        cells.push_back(cell);

        const std::array<std::pair<std::size_t, std::size_t>, 3> edges = {
            {
                {triangle.node_ids[0], triangle.node_ids[1]},
                {triangle.node_ids[1], triangle.node_ids[2]},
                {triangle.node_ids[2], triangle.node_ids[0]}
            }
        };

        for (const auto& [na, nb] : edges) {
            const EdgeKey key(na, nb);
            const auto face_it = edge_to_face_id.find(key);

            if (face_it == edge_to_face_id.end()) {
                Face face;
                face.id = faces.size();
                face.node_ids = {na, nb};
                face.owner_cell_id = cell.id;
                face.neighbor_cell_id = Face::k_invalid_cell_id;
                face.center_x = 0.5 * (points[na].x + points[nb].x);
                face.center_y = 0.5 * (points[na].y + points[nb].y);
                face.center_z = 0.0;
                face.measure = EdgeLength(points[na], points[nb]);

                faces.push_back(face);
                edge_to_face_id[key] = face.id;
                cells[cell.id].face_ids.push_back(face.id);
            }
            else {
                Face& face = faces[face_it->second];

                if (face.neighbor_cell_id != Face::k_invalid_cell_id) {
                    throw std::runtime_error("DelaunayMeshBuilder: non-manifold edge detected");
                }

                face.neighbor_cell_id = cell.id;
                cells[cell.id].face_ids.push_back(face.id);
            }
        }
    }

    for (Face& face : faces) {
        if (face.IsBoundary()) {
            const Point2D a{
                .id = 0,
                .x = nodes[face.node_ids[0]].x,
                .y = nodes[face.node_ids[0]].y,
                .is_boundary = false
            };
            const Point2D b{
                .id = 0,
                .x = nodes[face.node_ids[1]].x,
                .y = nodes[face.node_ids[1]].y,
                .is_boundary = false
            };

            face.boundary_tag = FindBoundaryTagForFace(a, b, points, boundary_segments);
        }
        else {
            face.boundary_tag = -1;
        }
    }

    FinalizeFaceNormals(mesh);
    mesh.Validate();
    return mesh;
}

void DelaunayMeshBuilder::FinalizeFaceNormals(Mesh& mesh) {
    auto& faces = mesh.Faces();
    const auto& cells = mesh.Cells();
    const auto& nodes = mesh.Nodes();

    for (Face& face : faces) {
        const Node& n0 = nodes[face.node_ids[0]];
        const Node& n1 = nodes[face.node_ids[1]];

        double nx = -(n1.y - n0.y);
        double ny = (n1.x - n0.x);

        const double owner_dx = face.center_x - cells[face.owner_cell_id].center_x;
        const double owner_dy = face.center_y - cells[face.owner_cell_id].center_y;

        if (nx * owner_dx + ny * owner_dy < 0.0) {
            nx *= -1.0;
            ny *= -1.0;
        }

        if (face.IsInternal()) {
            const double dx =
                cells[face.neighbor_cell_id].center_x - cells[face.owner_cell_id].center_x;
            const double dy =
                cells[face.neighbor_cell_id].center_y - cells[face.owner_cell_id].center_y;

            if (nx * dx + ny * dy < 0.0) {
                nx *= -1.0;
                ny *= -1.0;
            }
        }

        const double norm = std::sqrt(nx * nx + ny * ny);
        if (!(norm > 0.0)) {
            throw std::runtime_error("DelaunayMeshBuilder: zero face normal");
        }

        face.normal_x = nx / norm;
        face.normal_y = ny / norm;
        face.normal_z = 0.0;
    }
}

Mesh DelaunayMeshBuilder::BuildFromGeoFile(const std::string& file_path, const int dim) {
    ValidateInputDimension(dim);

    std::vector<Point2D> points;
    std::vector<BoundarySegment> boundary_segments;

    ReadGeometryFromGeoFile(file_path, points, boundary_segments);

    const std::vector<BoundaryLoop> loops = BuildBoundaryLoops(boundary_segments);

    PointCloudBuilder::BoundaryLayerSettings bl_settings;
    bl_settings.first_layer_height = 0.0; // auto
    bl_settings.growth = 1.12;
    bl_settings.n_layers = 7;

    PointCloudBuilder::SizeFunctionSettings size_settings;
    size_settings.h_min = 0.0; // auto
    size_settings.h_max = 0.0; // auto
    size_settings.transition_radius = 0.0; // auto
    size_settings.alpha = 1.8;

    PointCloudBuilder::BuildPointCloud(points,
                                       loops,
                                       bl_settings,
                                       size_settings);

    std::vector<DelaunayTriangulator::Point2D> tri_points;
    tri_points.reserve(points.size());

    for (const Point2D& p : points) {
        tri_points.push_back(DelaunayTriangulator::Point2D{
            .id = p.id,
            .x = p.x,
            .y = p.y
        });
    }

    const std::vector<DelaunayTriangulator::Triangle> raw_triangles =
        DelaunayTriangulator::Triangulate(tri_points);

    std::vector<Triangle> triangles;
    triangles.reserve(raw_triangles.size());

    for (const auto& t : raw_triangles) {
        triangles.push_back(Triangle{{t.node_ids[0], t.node_ids[1], t.node_ids[2]}});
    }

    triangles = FilterTrianglesInsideDomain(points, triangles, loops);

    return BuildMeshFromTriangles(points, triangles, boundary_segments);
}
