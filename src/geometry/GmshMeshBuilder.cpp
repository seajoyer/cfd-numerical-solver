#include "geometry/GmshMeshBuilder.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

#include <gmsh.h>

#include "geometry/Cell.hpp"
#include "geometry/Face.hpp"
#include "geometry/Node.hpp"

bool GmshMeshBuilder::FaceKey::operator<(const FaceKey& other) const {
    return node_ids_sorted < other.node_ids_sorted;
}

GmshMeshBuilder::GmshSessionGuard::GmshSessionGuard(const bool finalize_on_destroy)
    : finalize_on_destroy_(finalize_on_destroy) {}

GmshMeshBuilder::GmshSessionGuard::~GmshSessionGuard() {
    if (finalize_on_destroy_) {
        gmsh::finalize();
    }
}

void GmshMeshBuilder::ValidateInputDimension(const int dim) {
    if (dim != 2 && dim != 3) {
        throw std::invalid_argument("GmshMeshBuilder: dim must be 2 or 3");
    }
}

void GmshMeshBuilder::EnsureGmshInitialized() {
    if (!gmsh::isInitialized()) {
        throw std::runtime_error("GmshMeshBuilder: Gmsh is not initialized");
    }
}

bool GmshMeshBuilder::IsSupportedCellElementType(const int dim, const int element_type) {
    if (dim == 2) {
        return element_type == 2 || element_type == 3;
    }

    if (dim == 3) {
        return element_type == 4 || element_type == 5 || element_type == 6 || element_type == 7;
    }

    return false;
}

bool GmshMeshBuilder::IsSupportedBoundaryElementType(const int dim, const int element_type) {
    if (dim == 2) {
        return element_type == 1;
    }

    if (dim == 3) {
        return element_type == 2 || element_type == 3;
    }

    return false;
}

int GmshMeshBuilder::GetExpectedNodeCountForElementType(const int element_type) {
    switch (element_type) {
    case 1:
        return 2;
    case 2:
        return 3;
    case 3:
        return 4;
    case 4:
        return 4;
    case 5:
        return 8;
    case 6:
        return 6;
    case 7:
        return 5;
    default:
        throw std::runtime_error("GmshMeshBuilder: unsupported element type");
    }
}

std::vector<std::vector<int>> GmshMeshBuilder::GetLocalFacesForElementType(const int element_type) {
    switch (element_type) {
    case 2:
        return {
            {0, 1},
            {1, 2},
            {2, 0}
        };

    case 3:
        return {
            {0, 1},
            {1, 2},
            {2, 3},
            {3, 0}
        };

    case 4:
        return {
            {0, 2, 1},
            {0, 1, 3},
            {1, 2, 3},
            {2, 0, 3}
        };

    case 5:
        return {
            {0, 3, 2, 1},
            {4, 5, 6, 7},
            {0, 1, 5, 4},
            {1, 2, 6, 5},
            {2, 3, 7, 6},
            {3, 0, 4, 7}
        };

    case 6:
        return {
            {0, 2, 1},
            {3, 4, 5},
            {0, 1, 4, 3},
            {1, 2, 5, 4},
            {2, 0, 3, 5}
        };

    case 7:
        return {
            {0, 3, 2, 1},
            {0, 1, 4},
            {1, 2, 4},
            {2, 3, 4},
            {3, 0, 4}
        };

    default:
        throw std::runtime_error("GmshMeshBuilder: unsupported cell element type");
    }
}

std::vector<std::size_t> GmshMeshBuilder::ConvertElementConnectivity(
    const std::vector<std::size_t>& element_nodes,
    const std::vector<int>& local_face_pattern
) {
    std::vector<std::size_t> face_nodes;
    face_nodes.reserve(local_face_pattern.size());

    for (const int local_id : local_face_pattern) {
        if (local_id < 0 || static_cast<std::size_t>(local_id) >= element_nodes.size()) {
            throw std::runtime_error("GmshMeshBuilder: invalid local face node index");
        }

        face_nodes.push_back(element_nodes[static_cast<std::size_t>(local_id)]);
    }

    return face_nodes;
}

GmshMeshBuilder::FaceKey GmshMeshBuilder::MakeFaceKey(const std::vector<std::size_t>& node_ids) {
    FaceKey key;
    key.node_ids_sorted = node_ids;
    std::sort(key.node_ids_sorted.begin(), key.node_ids_sorted.end());
    return key;
}

GmshMeshBuilder::Vec3 GmshMeshBuilder::Add(const Vec3& a, const Vec3& b) {
    return Vec3{a.x + b.x, a.y + b.y, a.z + b.z};
}

GmshMeshBuilder::Vec3 GmshMeshBuilder::Subtract(const Vec3& a, const Vec3& b) {
    return Vec3{a.x - b.x, a.y - b.y, a.z - b.z};
}

GmshMeshBuilder::Vec3 GmshMeshBuilder::Multiply(const double scalar, const Vec3& v) {
    return Vec3{scalar * v.x, scalar * v.y, scalar * v.z};
}

GmshMeshBuilder::Vec3 GmshMeshBuilder::Divide(const Vec3& v, const double scalar) {
    return Vec3{v.x / scalar, v.y / scalar, v.z / scalar};
}

double GmshMeshBuilder::Dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

GmshMeshBuilder::Vec3 GmshMeshBuilder::Cross(const Vec3& a, const Vec3& b) {
    return Vec3{
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
}

double GmshMeshBuilder::Norm(const Vec3& v) {
    return std::sqrt(Dot(v, v));
}

GmshMeshBuilder::Vec3 GmshMeshBuilder::Normalize(const Vec3& v) {
    const double norm = Norm(v);
    if (!(norm > 0.0)) {
        throw std::runtime_error("GmshMeshBuilder: cannot normalize zero vector");
    }

    return Divide(v, norm);
}

GmshMeshBuilder::Vec3 GmshMeshBuilder::ToVec3(const Node& node) {
    return Vec3{node.x, node.y, node.z};
}

GmshMeshBuilder::Vec3 GmshMeshBuilder::ComputePolygonCenter(
    const std::vector<std::size_t>& node_ids,
    const std::vector<Node>& nodes
) {
    Vec3 center{0.0, 0.0, 0.0};

    for (const std::size_t node_id : node_ids) {
        center = Add(center, ToVec3(nodes[node_id]));
    }

    return Divide(center, static_cast<double>(node_ids.size()));
}

double GmshMeshBuilder::ComputeEdgeLength(
    const std::vector<std::size_t>& node_ids,
    const std::vector<Node>& nodes
) {
    if (node_ids.size() != 2) {
        throw std::runtime_error("GmshMeshBuilder: edge must have exactly 2 nodes");
    }

    return Norm(Subtract(ToVec3(nodes[node_ids[1]]), ToVec3(nodes[node_ids[0]])));
}

double GmshMeshBuilder::ComputePolygonArea3D(
    const std::vector<std::size_t>& node_ids,
    const std::vector<Node>& nodes
) {
    if (node_ids.size() < 3) {
        throw std::runtime_error("GmshMeshBuilder: polygon must have at least 3 nodes");
    }

    const Vec3 face_center = ComputePolygonCenter(node_ids, nodes);
    double area = 0.0;

    for (std::size_t i = 0; i < node_ids.size(); ++i) {
        const Vec3 a = ToVec3(nodes[node_ids[i]]);
        const Vec3 b = ToVec3(nodes[node_ids[(i + 1) % node_ids.size()]]);
        const Vec3 ac = Subtract(a, face_center);
        const Vec3 bc = Subtract(b, face_center);

        area += 0.5 * Norm(Cross(ac, bc));
    }

    return area;
}

double GmshMeshBuilder::ComputeCellArea2D(
    const std::vector<std::size_t>& node_ids,
    const std::vector<Node>& nodes
) {
    double twice_area = 0.0;

    for (std::size_t i = 0; i < node_ids.size(); ++i) {
        const Node& a = nodes[node_ids[i]];
        const Node& b = nodes[node_ids[(i + 1) % node_ids.size()]];
        twice_area += a.x * b.y - b.x * a.y;
    }

    return 0.5 * std::abs(twice_area);
}

double GmshMeshBuilder::ComputeTetraVolume(
    const Vec3& a,
    const Vec3& b,
    const Vec3& c,
    const Vec3& d
) {
    return std::abs(Dot(Subtract(a, d), Cross(Subtract(b, d), Subtract(c, d)))) / 6.0;
}

double GmshMeshBuilder::ComputeCellVolume3D(
    const int element_type,
    const std::vector<std::size_t>& node_ids,
    const std::vector<Node>& nodes
) {
    if (node_ids.size() < 4) {
        throw std::runtime_error("GmshMeshBuilder: 3D cell must have at least 4 nodes");
    }

    if (element_type == 4) {
        return ComputeTetraVolume(
            ToVec3(nodes[node_ids[0]]),
            ToVec3(nodes[node_ids[1]]),
            ToVec3(nodes[node_ids[2]]),
            ToVec3(nodes[node_ids[3]])
        );
    }

    const Vec3 cell_center = ComputePolygonCenter(node_ids, nodes);
    const std::vector<std::vector<int>> local_faces = GetLocalFacesForElementType(element_type);

    double volume = 0.0;

    for (const std::vector<int>& local_face_pattern : local_faces) {
        const std::vector<std::size_t> face_nodes = ConvertElementConnectivity(node_ids, local_face_pattern);
        const Vec3 face_center = ComputePolygonCenter(face_nodes, nodes);

        for (std::size_t i = 0; i < face_nodes.size(); ++i) {
            const Vec3 a = ToVec3(nodes[face_nodes[i]]);
            const Vec3 b = ToVec3(nodes[face_nodes[(i + 1) % face_nodes.size()]]);
            volume += ComputeTetraVolume(cell_center, a, b, face_center);
        }
    }

    return volume;
}

std::map<GmshMeshBuilder::FaceKey, GmshMeshBuilder::BoundaryFaceRecord>
GmshMeshBuilder::BuildBoundaryFaceMap(
    const int dim,
    const std::unordered_map<std::size_t, std::size_t>& gmsh_to_internal_node
) {
    std::map<FaceKey, BoundaryFaceRecord> boundary_faces;

    std::vector<std::pair<int, int>> physical_groups;
    gmsh::model::getPhysicalGroups(physical_groups, dim - 1);

    for (const auto& [physical_dim, physical_tag] : physical_groups) {
        if (physical_dim != dim - 1) {
            continue;
        }

        std::vector<int> entity_tags;
        gmsh::model::getEntitiesForPhysicalGroup(physical_dim, physical_tag, entity_tags);

        for (const int entity_tag : entity_tags) {
            std::vector<int> element_types;
            std::vector<std::vector<std::size_t>> element_tags;
            std::vector<std::vector<std::size_t>> element_node_tags;

            gmsh::model::mesh::getElements(
                element_types,
                element_tags,
                element_node_tags,
                physical_dim,
                entity_tag
            );

            for (std::size_t block_id = 0; block_id < element_types.size(); ++block_id) {
                const int element_type = element_types[block_id];

                if (!IsSupportedBoundaryElementType(dim, element_type)) {
                    continue;
                }

                const int node_count = GetExpectedNodeCountForElementType(element_type);
                const std::vector<std::size_t>& flat_nodes = element_node_tags[block_id];

                if (flat_nodes.size() % static_cast<std::size_t>(node_count) != 0) {
                    throw std::runtime_error(
                        "GmshMeshBuilder: invalid boundary element connectivity block"
                    );
                }

                const std::size_t element_count =
                    flat_nodes.size() / static_cast<std::size_t>(node_count);

                for (std::size_t element_id = 0; element_id < element_count; ++element_id) {
                    std::vector<std::size_t> face_nodes;
                    face_nodes.reserve(static_cast<std::size_t>(node_count));

                    for (int local_id = 0; local_id < node_count; ++local_id) {
                        const std::size_t gmsh_node_tag =
                            flat_nodes[element_id * static_cast<std::size_t>(node_count) +
                                static_cast<std::size_t>(local_id)];

                        const auto it = gmsh_to_internal_node.find(gmsh_node_tag);
                        if (it == gmsh_to_internal_node.end()) {
                            throw std::runtime_error(
                                "GmshMeshBuilder: boundary node not found in node map"
                            );
                        }

                        face_nodes.push_back(it->second);
                    }

                    boundary_faces[MakeFaceKey(face_nodes)] = BoundaryFaceRecord{physical_tag};
                }
            }
        }
    }

    return boundary_faces;
}

void GmshMeshBuilder::BuildNodes(
    Mesh& mesh,
    std::unordered_map<std::size_t, std::size_t>& gmsh_to_internal_node
) {
    std::vector<std::size_t> node_tags;
    std::vector<double> coords;
    std::vector<double> params;

    gmsh::model::mesh::getNodes(node_tags, coords, params, -1, -1, false, false);

    if (coords.size() != 3 * node_tags.size()) {
        throw std::runtime_error("GmshMeshBuilder: invalid node coordinate array size");
    }

    auto& nodes = mesh.Nodes();
    nodes.reserve(node_tags.size());

    for (std::size_t i = 0; i < node_tags.size(); ++i) {
        Node node;
        node.id = nodes.size();
        node.x = coords[3 * i + 0];
        node.y = coords[3 * i + 1];
        node.z = coords[3 * i + 2];

        gmsh_to_internal_node[node_tags[i]] = node.id;
        nodes.push_back(node);
    }

    if (nodes.empty()) {
        throw std::runtime_error("GmshMeshBuilder: no nodes found in Gmsh model");
    }
}

void GmshMeshBuilder::BuildCellsAndFaces(
    Mesh& mesh,
    const int dim,
    const std::unordered_map<std::size_t, std::size_t>& gmsh_to_internal_node,
    const std::map<FaceKey, BoundaryFaceRecord>& boundary_faces
) {
    std::vector<int> element_types;
    std::vector<std::vector<std::size_t>> element_tags;
    std::vector<std::vector<std::size_t>> element_node_tags;

    gmsh::model::mesh::getElements(element_types, element_tags, element_node_tags, dim, -1);

    auto& cells = mesh.Cells();
    auto& faces = mesh.Faces();
    const auto& nodes = mesh.Nodes();

    std::map<FaceKey, std::size_t> face_key_to_face_id;

    for (std::size_t block_id = 0; block_id < element_types.size(); ++block_id) {
        const int element_type = element_types[block_id];

        if (!IsSupportedCellElementType(dim, element_type)) {
            continue;
        }

        const int node_count = GetExpectedNodeCountForElementType(element_type);
        const std::vector<std::size_t>& flat_nodes = element_node_tags[block_id];

        if (flat_nodes.size() % static_cast<std::size_t>(node_count) != 0) {
            throw std::runtime_error("GmshMeshBuilder: invalid cell element connectivity block");
        }

        const std::size_t element_count =
            flat_nodes.size() / static_cast<std::size_t>(node_count);

        const std::vector<std::vector<int>> local_faces =
            GetLocalFacesForElementType(element_type);

        for (std::size_t element_id = 0; element_id < element_count; ++element_id) {
            Cell cell;
            cell.id = cells.size();
            cell.node_ids.reserve(static_cast<std::size_t>(node_count));

            for (int local_id = 0; local_id < node_count; ++local_id) {
                const std::size_t gmsh_node_tag =
                    flat_nodes[element_id * static_cast<std::size_t>(node_count) +
                        static_cast<std::size_t>(local_id)];

                const auto it = gmsh_to_internal_node.find(gmsh_node_tag);
                if (it == gmsh_to_internal_node.end()) {
                    throw std::runtime_error("GmshMeshBuilder: cell node not found in node map");
                }

                cell.node_ids.push_back(it->second);
            }

            const Vec3 cell_center = ComputePolygonCenter(cell.node_ids, nodes);
            cell.center_x = cell_center.x;
            cell.center_y = cell_center.y;
            cell.center_z = cell_center.z;
            cell.volume = (dim == 2)
                              ? ComputeCellArea2D(cell.node_ids, nodes)
                              : ComputeCellVolume3D(element_type, cell.node_ids, nodes);

            cells.push_back(cell);

            for (const std::vector<int>& local_face_pattern : local_faces) {
                const std::vector<std::size_t> face_nodes =
                    ConvertElementConnectivity(cell.node_ids, local_face_pattern);

                const FaceKey key = MakeFaceKey(face_nodes);
                const auto face_it = face_key_to_face_id.find(key);

                if (face_it == face_key_to_face_id.end()) {
                    Face face;
                    face.id = faces.size();
                    face.node_ids = face_nodes;
                    face.owner_cell_id = cell.id;
                    face.neighbor_cell_id = Face::k_invalid_cell_id;

                    const Vec3 face_center = ComputePolygonCenter(face.node_ids, nodes);
                    face.center_x = face_center.x;
                    face.center_y = face_center.y;
                    face.center_z = face_center.z;

                    face.measure = (dim == 2)
                                       ? ComputeEdgeLength(face.node_ids, nodes)
                                       : ComputePolygonArea3D(face.node_ids, nodes);

                    const auto boundary_it = boundary_faces.find(key);
                    face.boundary_tag = (boundary_it != boundary_faces.end())
                                            ? boundary_it->second.boundary_tag
                                            : -1;

                    faces.push_back(face);
                    face_key_to_face_id[key] = face.id;
                    cells[cell.id].face_ids.push_back(face.id);
                }
                else {
                    Face& face = faces[face_it->second];

                    if (face.neighbor_cell_id != Face::k_invalid_cell_id) {
                        throw std::runtime_error("GmshMeshBuilder: non-manifold face detected");
                    }

                    face.neighbor_cell_id = cell.id;
                    face.boundary_tag = -1;
                    cells[cell.id].face_ids.push_back(face.id);
                }
            }
        }
    }

    if (cells.empty()) {
        throw std::runtime_error("GmshMeshBuilder: no supported cell elements found");
    }
}

void GmshMeshBuilder::FinalizeFaceNormals(Mesh& mesh) {
    auto& faces = mesh.Faces();
    const auto& cells = mesh.Cells();
    const auto& nodes = mesh.Nodes();

    for (Face& face : faces) {
        const Vec3 owner_center{
            cells[face.owner_cell_id].center_x,
            cells[face.owner_cell_id].center_y,
            cells[face.owner_cell_id].center_z
        };

        Vec3 direction;
        if (face.IsInternal()) {
            const Vec3 neighbor_center{
                cells[face.neighbor_cell_id].center_x,
                cells[face.neighbor_cell_id].center_y,
                cells[face.neighbor_cell_id].center_z
            };
            direction = Subtract(neighbor_center, owner_center);
        }
        else {
            const Vec3 face_center{.x = face.center_x, .y = face.center_y, .z = face.center_z};
            direction = Subtract(face_center, owner_center);
        }

        if (face.node_ids.size() == 2) {
            const Vec3 a = ToVec3(nodes[face.node_ids[0]]);
            const Vec3 b = ToVec3(nodes[face.node_ids[1]]);
            const Vec3 tangent = Subtract(b, a);

            Vec3 normal{.x = -tangent.y, .y = tangent.x, .z = 0.0};
            if (Dot(normal, direction) < 0.0) {
                normal = Multiply(-1.0, normal);
            }

            normal = Normalize(normal);

            face.normal_x = normal.x;
            face.normal_y = normal.y;
            face.normal_z = normal.z;
        }
        else {
            const Vec3 face_center = ComputePolygonCenter(face.node_ids, nodes);
            Vec3 area_vector{0.0, 0.0, 0.0};

            for (std::size_t i = 0; i < face.node_ids.size(); ++i) {
                const Vec3 a = Subtract(ToVec3(nodes[face.node_ids[i]]), face_center);
                const Vec3 b = Subtract(
                    ToVec3(nodes[face.node_ids[(i + 1) % face.node_ids.size()]]),
                    face_center
                );

                area_vector = Add(area_vector, Cross(a, b));
            }

            if (Dot(area_vector, direction) < 0.0) {
                area_vector = Multiply(-1.0, area_vector);
            }

            const Vec3 normal = Normalize(area_vector);

            face.normal_x = normal.x;
            face.normal_y = normal.y;
            face.normal_z = normal.z;
        }

        if (face.IsBoundary() && face.boundary_tag < 0) {
            throw std::runtime_error(
                "GmshMeshBuilder: boundary face has no physical boundary tag"
            );
        }
    }
}

Mesh GmshMeshBuilder::BuildFromCurrentModel(const int dim) {
    ValidateInputDimension(dim);
    EnsureGmshInitialized();

    Mesh mesh(dim);

    std::unordered_map<std::size_t, std::size_t> gmsh_to_internal_node;
    BuildNodes(mesh, gmsh_to_internal_node);

    const std::map<FaceKey, BoundaryFaceRecord> boundary_faces =
        BuildBoundaryFaceMap(dim, gmsh_to_internal_node);

    BuildCellsAndFaces(mesh, dim, gmsh_to_internal_node, boundary_faces);
    FinalizeFaceNormals(mesh);

    mesh.Validate();
    return mesh;
}

Mesh GmshMeshBuilder::BuildFromFile(const std::string& file_path, const int dim) {
    ValidateInputDimension(dim);

    const bool was_initialized = gmsh::isInitialized();
    if (!was_initialized) {
        gmsh::initialize();
    }

    GmshSessionGuard session_guard(!was_initialized);

    gmsh::clear();
    gmsh::open(file_path);

    return BuildFromCurrentModel(dim);
}

Mesh GmshMeshBuilder::BuildFromGeoFile(const std::string& file_path, const int dim) {
    ValidateInputDimension(dim);

    const bool was_initialized = gmsh::isInitialized();
    if (!was_initialized) {
        gmsh::initialize();
    }

    GmshSessionGuard session_guard(!was_initialized);

    gmsh::clear();
    gmsh::open(file_path);
    gmsh::model::mesh::generate(dim);

    return BuildFromCurrentModel(dim);
}
