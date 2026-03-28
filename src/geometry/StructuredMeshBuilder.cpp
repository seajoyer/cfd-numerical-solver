#include "geometry/StructuredMeshBuilder.hpp"

#include <array>
#include <stdexcept>
#include <vector>

namespace {
    std::size_t NodeIndex2D(const int i, const int j, const int nx, const int ny) {
        (void)ny;
        return static_cast<std::size_t>(j * (nx + 1) + i);
    }

    std::size_t CellIndex2D(const int i, const int j, const int nx, const int ny) {
        (void)ny;
        return static_cast<std::size_t>(j * nx + i);
    }

    std::size_t NodeIndex3D(const int i, const int j, const int k,
                            const int nx, const int ny, const int nz) {
        (void)nz;
        return static_cast<std::size_t>(k * (ny + 1) * (nx + 1) + j * (nx + 1) + i);
    }

    std::size_t CellIndex3D(const int i, const int j, const int k,
                            const int nx, const int ny, const int nz) {
        (void)nz;
        return static_cast<std::size_t>(k * ny * nx + j * nx + i);
    }
} // namespace

void StructuredMeshBuilder::Validate2DInput(const int nx, const int ny,
                                            const double x_min, const double x_max,
                                            const double y_min, const double y_max) {
    if (nx <= 0) {
        throw std::invalid_argument("StructuredMeshBuilder::BuildUniformCartesian2D: nx must be > 0");
    }
    if (ny <= 0) {
        throw std::invalid_argument("StructuredMeshBuilder::BuildUniformCartesian2D: ny must be > 0");
    }
    if (!(x_max > x_min)) {
        throw std::invalid_argument("StructuredMeshBuilder::BuildUniformCartesian2D: x_max must be > x_min");
    }
    if (!(y_max > y_min)) {
        throw std::invalid_argument("StructuredMeshBuilder::BuildUniformCartesian2D: y_max must be > y_min");
    }
}

void StructuredMeshBuilder::Validate3DInput(const int nx, const int ny, const int nz,
                                            const double x_min, const double x_max,
                                            const double y_min, const double y_max,
                                            const double z_min, const double z_max) {
    if (nx <= 0) {
        throw std::invalid_argument("StructuredMeshBuilder::BuildUniformCartesian3D: nx must be > 0");
    }
    if (ny <= 0) {
        throw std::invalid_argument("StructuredMeshBuilder::BuildUniformCartesian3D: ny must be > 0");
    }
    if (nz <= 0) {
        throw std::invalid_argument("StructuredMeshBuilder::BuildUniformCartesian3D: nz must be > 0");
    }
    if (!(x_max > x_min)) {
        throw std::invalid_argument("StructuredMeshBuilder::BuildUniformCartesian3D: x_max must be > x_min");
    }
    if (!(y_max > y_min)) {
        throw std::invalid_argument("StructuredMeshBuilder::BuildUniformCartesian3D: y_max must be > y_min");
    }
    if (!(z_max > z_min)) {
        throw std::invalid_argument("StructuredMeshBuilder::BuildUniformCartesian3D: z_max must be > z_min");
    }
}

Mesh StructuredMeshBuilder::BuildUniformCartesian2D(const int nx, const int ny,
                                                    const double x_min, const double x_max,
                                                    const double y_min, const double y_max) {
    Validate2DInput(nx, ny, x_min, x_max, y_min, y_max);

    Mesh mesh(2);

    const double dx = (x_max - x_min) / static_cast<double>(nx);
    const double dy = (y_max - y_min) / static_cast<double>(ny);
    const double cell_area = dx * dy;

    auto& nodes = mesh.Nodes();
    auto& faces = mesh.Faces();
    auto& cells = mesh.Cells();

    nodes.reserve(static_cast<std::size_t>((nx + 1) * (ny + 1)));
    cells.reserve(static_cast<std::size_t>(nx * ny));
    faces.reserve(static_cast<std::size_t>((nx + 1) * ny + nx * (ny + 1)));

    for (int j = 0; j <= ny; ++j) {
        const double y = y_min + static_cast<double>(j) * dy;
        for (int i = 0; i <= nx; ++i) {
            const double x = x_min + static_cast<double>(i) * dx;

            Node node;
            node.id = nodes.size();
            node.x = x;
            node.y = y;
            node.z = 0.0;

            nodes.push_back(node);
        }
    }

    for (int j = 0; j < ny; ++j) {
        const double yc = y_min + (static_cast<double>(j) + 0.5) * dy;
        for (int i = 0; i < nx; ++i) {
            const double xc = x_min + (static_cast<double>(i) + 0.5) * dx;

            Cell cell;
            cell.id = cells.size();
            cell.node_ids = {
                NodeIndex2D(i, j, nx, ny),
                NodeIndex2D(i + 1, j, nx, ny),
                NodeIndex2D(i + 1, j + 1, nx, ny),
                NodeIndex2D(i, j + 1, nx, ny)
            };
            cell.center_x = xc;
            cell.center_y = yc;
            cell.center_z = 0.0;
            cell.volume = cell_area;

            cells.push_back(cell);
        }
    }

    // Vertical faces: normal in +/- x direction
    for (int j = 0; j < ny; ++j) {
        const double yc = y_min + (static_cast<double>(j) + 0.5) * dy;
        for (int i = 0; i <= nx; ++i) {
            Face face;
            face.id = faces.size();
            face.node_ids = {
                NodeIndex2D(i, j, nx, ny),
                NodeIndex2D(i, j + 1, nx, ny)
            };
            face.center_x = x_min + static_cast<double>(i) * dx;
            face.center_y = yc;
            face.center_z = 0.0;
            face.measure = dy;

            if (i == 0) {
                face.owner_cell_id = CellIndex2D(0, j, nx, ny);
                face.neighbor_cell_id = Face::k_invalid_cell_id;
                face.normal_x = -1.0;
                face.normal_y = 0.0;
                face.normal_z = 0.0;
                face.boundary_tag = k_xmin_tag;
            } else if (i == nx) {
                face.owner_cell_id = CellIndex2D(nx - 1, j, nx, ny);
                face.neighbor_cell_id = Face::k_invalid_cell_id;
                face.normal_x = 1.0;
                face.normal_y = 0.0;
                face.normal_z = 0.0;
                face.boundary_tag = k_xmax_tag;
            } else {
                face.owner_cell_id = CellIndex2D(i - 1, j, nx, ny);
                face.neighbor_cell_id = CellIndex2D(i, j, nx, ny);
                face.normal_x = 1.0;
                face.normal_y = 0.0;
                face.normal_z = 0.0;
                face.boundary_tag = -1;
            }

            faces.push_back(face);
            cells[face.owner_cell_id].face_ids.push_back(face.id);
            if (face.IsInternal()) {
                cells[face.neighbor_cell_id].face_ids.push_back(face.id);
            }
        }
    }

    // Horizontal faces: normal in +/- y direction
    for (int j = 0; j <= ny; ++j) {
        const double y = y_min + static_cast<double>(j) * dy;
        for (int i = 0; i < nx; ++i) {
            Face face;
            face.id = faces.size();
            face.node_ids = {
                NodeIndex2D(i, j, nx, ny),
                NodeIndex2D(i + 1, j, nx, ny)
            };
            face.center_x = x_min + (static_cast<double>(i) + 0.5) * dx;
            face.center_y = y;
            face.center_z = 0.0;
            face.measure = dx;

            if (j == 0) {
                face.owner_cell_id = CellIndex2D(i, 0, nx, ny);
                face.neighbor_cell_id = Face::k_invalid_cell_id;
                face.normal_x = 0.0;
                face.normal_y = -1.0;
                face.normal_z = 0.0;
                face.boundary_tag = k_ymin_tag;
            } else if (j == ny) {
                face.owner_cell_id = CellIndex2D(i, ny - 1, nx, ny);
                face.neighbor_cell_id = Face::k_invalid_cell_id;
                face.normal_x = 0.0;
                face.normal_y = 1.0;
                face.normal_z = 0.0;
                face.boundary_tag = k_ymax_tag;
            } else {
                face.owner_cell_id = CellIndex2D(i, j - 1, nx, ny);
                face.neighbor_cell_id = CellIndex2D(i, j, nx, ny);
                face.normal_x = 0.0;
                face.normal_y = 1.0;
                face.normal_z = 0.0;
                face.boundary_tag = -1;
            }

            faces.push_back(face);
            cells[face.owner_cell_id].face_ids.push_back(face.id);
            if (face.IsInternal()) {
                cells[face.neighbor_cell_id].face_ids.push_back(face.id);
            }
        }
    }

    mesh.Validate();
    return mesh;
}

Mesh StructuredMeshBuilder::BuildUniformCartesian3D(const int nx, const int ny, const int nz,
                                                    const double x_min, const double x_max,
                                                    const double y_min, const double y_max,
                                                    const double z_min, const double z_max) {
    Validate3DInput(nx, ny, nz, x_min, x_max, y_min, y_max, z_min, z_max);

    Mesh mesh(3);

    const double dx = (x_max - x_min) / static_cast<double>(nx);
    const double dy = (y_max - y_min) / static_cast<double>(ny);
    const double dz = (z_max - z_min) / static_cast<double>(nz);

    const double cell_volume = dx * dy * dz;
    const double yz_area = dy * dz;
    const double xz_area = dx * dz;
    const double xy_area = dx * dy;

    auto& nodes = mesh.Nodes();
    auto& faces = mesh.Faces();
    auto& cells = mesh.Cells();

    nodes.reserve(static_cast<std::size_t>((nx + 1) * (ny + 1) * (nz + 1)));
    cells.reserve(static_cast<std::size_t>(nx * ny * nz));
    faces.reserve(static_cast<std::size_t>(
                      (nx + 1) * ny * nz +
                      nx * (ny + 1) * nz +
                      nx * ny * (nz + 1)
                  ));

    for (int k = 0; k <= nz; ++k) {
        const double z = z_min + static_cast<double>(k) * dz;
        for (int j = 0; j <= ny; ++j) {
            const double y = y_min + static_cast<double>(j) * dy;
            for (int i = 0; i <= nx; ++i) {
                const double x = x_min + static_cast<double>(i) * dx;

                Node node;
                node.id = nodes.size();
                node.x = x;
                node.y = y;
                node.z = z;

                nodes.push_back(node);
            }
        }
    }

    for (int k = 0; k < nz; ++k) {
        const double zc = z_min + (static_cast<double>(k) + 0.5) * dz;
        for (int j = 0; j < ny; ++j) {
            const double yc = y_min + (static_cast<double>(j) + 0.5) * dy;
            for (int i = 0; i < nx; ++i) {
                const double xc = x_min + (static_cast<double>(i) + 0.5) * dx;

                Cell cell;
                cell.id = cells.size();
                cell.node_ids = {
                    NodeIndex3D(i, j, k, nx, ny, nz),
                    NodeIndex3D(i + 1, j, k, nx, ny, nz),
                    NodeIndex3D(i + 1, j + 1, k, nx, ny, nz),
                    NodeIndex3D(i, j + 1, k, nx, ny, nz),
                    NodeIndex3D(i, j, k + 1, nx, ny, nz),
                    NodeIndex3D(i + 1, j, k + 1, nx, ny, nz),
                    NodeIndex3D(i + 1, j + 1, k + 1, nx, ny, nz),
                    NodeIndex3D(i, j + 1, k + 1, nx, ny, nz)
                };
                cell.center_x = xc;
                cell.center_y = yc;
                cell.center_z = zc;
                cell.volume = cell_volume;

                cells.push_back(cell);
            }
        }
    }

    // X-normal faces
    for (int k = 0; k < nz; ++k) {
        const double zc = z_min + (static_cast<double>(k) + 0.5) * dz;
        for (int j = 0; j < ny; ++j) {
            const double yc = y_min + (static_cast<double>(j) + 0.5) * dy;
            for (int i = 0; i <= nx; ++i) {
                Face face;
                face.id = faces.size();
                face.node_ids = {
                    NodeIndex3D(i, j, k, nx, ny, nz),
                    NodeIndex3D(i, j + 1, k, nx, ny, nz),
                    NodeIndex3D(i, j + 1, k + 1, nx, ny, nz),
                    NodeIndex3D(i, j, k + 1, nx, ny, nz)
                };
                face.center_x = x_min + static_cast<double>(i) * dx;
                face.center_y = yc;
                face.center_z = zc;
                face.measure = yz_area;

                if (i == 0) {
                    face.owner_cell_id = CellIndex3D(0, j, k, nx, ny, nz);
                    face.neighbor_cell_id = Face::k_invalid_cell_id;
                    face.normal_x = -1.0;
                    face.normal_y = 0.0;
                    face.normal_z = 0.0;
                    face.boundary_tag = k_xmin_tag;
                } else if (i == nx) {
                    face.owner_cell_id = CellIndex3D(nx - 1, j, k, nx, ny, nz);
                    face.neighbor_cell_id = Face::k_invalid_cell_id;
                    face.normal_x = 1.0;
                    face.normal_y = 0.0;
                    face.normal_z = 0.0;
                    face.boundary_tag = k_xmax_tag;
                } else {
                    face.owner_cell_id = CellIndex3D(i - 1, j, k, nx, ny, nz);
                    face.neighbor_cell_id = CellIndex3D(i, j, k, nx, ny, nz);
                    face.normal_x = 1.0;
                    face.normal_y = 0.0;
                    face.normal_z = 0.0;
                    face.boundary_tag = -1;
                }

                faces.push_back(face);
                cells[face.owner_cell_id].face_ids.push_back(face.id);
                if (face.IsInternal()) {
                    cells[face.neighbor_cell_id].face_ids.push_back(face.id);
                }
            }
        }
    }

    // Y-normal faces
    for (int k = 0; k < nz; ++k) {
        const double zc = z_min + (static_cast<double>(k) + 0.5) * dz;
        for (int j = 0; j <= ny; ++j) {
            const double y = y_min + static_cast<double>(j) * dy;
            for (int i = 0; i < nx; ++i) {
                Face face;
                face.id = faces.size();
                face.node_ids = {
                    NodeIndex3D(i, j, k, nx, ny, nz),
                    NodeIndex3D(i + 1, j, k, nx, ny, nz),
                    NodeIndex3D(i + 1, j, k + 1, nx, ny, nz),
                    NodeIndex3D(i, j, k + 1, nx, ny, nz)
                };
                face.center_x = x_min + (static_cast<double>(i) + 0.5) * dx;
                face.center_y = y;
                face.center_z = zc;
                face.measure = xz_area;

                if (j == 0) {
                    face.owner_cell_id = CellIndex3D(i, 0, k, nx, ny, nz);
                    face.neighbor_cell_id = Face::k_invalid_cell_id;
                    face.normal_x = 0.0;
                    face.normal_y = -1.0;
                    face.normal_z = 0.0;
                    face.boundary_tag = k_ymin_tag;
                } else if (j == ny) {
                    face.owner_cell_id = CellIndex3D(i, ny - 1, k, nx, ny, nz);
                    face.neighbor_cell_id = Face::k_invalid_cell_id;
                    face.normal_x = 0.0;
                    face.normal_y = 1.0;
                    face.normal_z = 0.0;
                    face.boundary_tag = k_ymax_tag;
                } else {
                    face.owner_cell_id = CellIndex3D(i, j - 1, k, nx, ny, nz);
                    face.neighbor_cell_id = CellIndex3D(i, j, k, nx, ny, nz);
                    face.normal_x = 0.0;
                    face.normal_y = 1.0;
                    face.normal_z = 0.0;
                    face.boundary_tag = -1;
                }

                faces.push_back(face);
                cells[face.owner_cell_id].face_ids.push_back(face.id);
                if (face.IsInternal()) {
                    cells[face.neighbor_cell_id].face_ids.push_back(face.id);
                }
            }
        }
    }

    // Z-normal faces
    for (int k = 0; k <= nz; ++k) {
        const double z = z_min + static_cast<double>(k) * dz;
        for (int j = 0; j < ny; ++j) {
            const double yc = y_min + (static_cast<double>(j) + 0.5) * dy;
            for (int i = 0; i < nx; ++i) {
                Face face;
                face.id = faces.size();
                face.node_ids = {
                    NodeIndex3D(i, j, k, nx, ny, nz),
                    NodeIndex3D(i + 1, j, k, nx, ny, nz),
                    NodeIndex3D(i + 1, j + 1, k, nx, ny, nz),
                    NodeIndex3D(i, j + 1, k, nx, ny, nz)
                };
                face.center_x = x_min + (static_cast<double>(i) + 0.5) * dx;
                face.center_y = yc;
                face.center_z = z;
                face.measure = xy_area;

                if (k == 0) {
                    face.owner_cell_id = CellIndex3D(i, j, 0, nx, ny, nz);
                    face.neighbor_cell_id = Face::k_invalid_cell_id;
                    face.normal_x = 0.0;
                    face.normal_y = 0.0;
                    face.normal_z = -1.0;
                    face.boundary_tag = k_zmin_tag;
                } else if (k == nz) {
                    face.owner_cell_id = CellIndex3D(i, j, nz - 1, nx, ny, nz);
                    face.neighbor_cell_id = Face::k_invalid_cell_id;
                    face.normal_x = 0.0;
                    face.normal_y = 0.0;
                    face.normal_z = 1.0;
                    face.boundary_tag = k_zmax_tag;
                } else {
                    face.owner_cell_id = CellIndex3D(i, j, k - 1, nx, ny, nz);
                    face.neighbor_cell_id = CellIndex3D(i, j, k, nx, ny, nz);
                    face.normal_x = 0.0;
                    face.normal_y = 0.0;
                    face.normal_z = 1.0;
                    face.boundary_tag = -1;
                }

                faces.push_back(face);
                cells[face.owner_cell_id].face_ids.push_back(face.id);
                if (face.IsInternal()) {
                    cells[face.neighbor_cell_id].face_ids.push_back(face.id);
                }
            }
        }
    }

    mesh.Validate();
    return mesh;
}
