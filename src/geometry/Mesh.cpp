#include "geometry/Mesh.hpp"

#include <cmath>
#include <string>

Mesh::Mesh(const int dim) : dim_(dim) {
    ValidateDimension();
}

int Mesh::GetDim() const {
    return dim_;
}

std::size_t Mesh::GetNodeCount() const {
    return nodes_.size();
}

std::size_t Mesh::GetFaceCount() const {
    return faces_.size();
}

std::size_t Mesh::GetCellCount() const {
    return cells_.size();
}

std::vector<Node>& Mesh::Nodes() {
    return nodes_;
}

const std::vector<Node>& Mesh::Nodes() const {
    return nodes_;
}

std::vector<Face>& Mesh::Faces() {
    return faces_;
}

const std::vector<Face>& Mesh::Faces() const {
    return faces_;
}

std::vector<Cell>& Mesh::Cells() {
    return cells_;
}

const std::vector<Cell>& Mesh::Cells() const {
    return cells_;
}

Node& Mesh::GetNode(const std::size_t node_id) {
    if (node_id >= nodes_.size()) {
        throw std::out_of_range("Mesh::GetNode: node_id is out of range");
    }
    return nodes_[node_id];
}

const Node& Mesh::GetNode(const std::size_t node_id) const {
    if (node_id >= nodes_.size()) {
        throw std::out_of_range("Mesh::GetNode: node_id is out of range");
    }
    return nodes_[node_id];
}

Face& Mesh::GetFace(const std::size_t face_id) {
    if (face_id >= faces_.size()) {
        throw std::out_of_range("Mesh::GetFace: face_id is out of range");
    }
    return faces_[face_id];
}

const Face& Mesh::GetFace(const std::size_t face_id) const {
    if (face_id >= faces_.size()) {
        throw std::out_of_range("Mesh::GetFace: face_id is out of range");
    }
    return faces_[face_id];
}

Cell& Mesh::GetCell(const std::size_t cell_id) {
    if (cell_id >= cells_.size()) {
        throw std::out_of_range("Mesh::GetCell: cell_id is out of range");
    }
    return cells_[cell_id];
}

const Cell& Mesh::GetCell(const std::size_t cell_id) const {
    if (cell_id >= cells_.size()) {
        throw std::out_of_range("Mesh::GetCell: cell_id is out of range");
    }
    return cells_[cell_id];
}

bool Mesh::IsBoundaryFace(const std::size_t face_id) const {
    return GetFace(face_id).IsBoundary();
}

bool Mesh::IsInternalFace(const std::size_t face_id) const {
    return GetFace(face_id).IsInternal();
}

void Mesh::Clear() {
    nodes_.clear();
    faces_.clear();
    cells_.clear();
}

void Mesh::Validate() const {
    ValidateDimension();
    ValidateNodeIds();
    ValidateFaceIds();
    ValidateCellIds();
    ValidateFaceConnectivity();
    ValidateCellConnectivity();
    ValidateGeometry();
}

void Mesh::ValidateDimension() const {
    if (dim_ < 1 || dim_ > 3) {
        throw std::runtime_error("Mesh::Validate: dim must be 1..3");
    }
}

void Mesh::ValidateNodeIds() const {
    for (std::size_t node_id = 0; node_id < nodes_.size(); ++node_id) {
        if (nodes_[node_id].id != node_id) {
            throw std::runtime_error("Mesh::Validate: node id does not match storage index");
        }
    }
}

void Mesh::ValidateFaceIds() const {
    for (std::size_t face_id = 0; face_id < faces_.size(); ++face_id) {
        if (faces_[face_id].id != face_id) {
            throw std::runtime_error("Mesh::Validate: face id does not match storage index");
        }
    }
}

void Mesh::ValidateCellIds() const {
    for (std::size_t cell_id = 0; cell_id < cells_.size(); ++cell_id) {
        if (cells_[cell_id].id != cell_id) {
            throw std::runtime_error("Mesh::Validate: cell id does not match storage index");
        }
    }
}

void Mesh::ValidateFaceConnectivity() const {
    for (const Face& face : faces_) {
        if (face.node_ids.empty()) {
            throw std::runtime_error("Mesh::Validate: face has empty node_ids");
        }

        if (face.owner_cell_id == Face::k_invalid_cell_id) {
            throw std::runtime_error("Mesh::Validate: face has invalid owner_cell_id");
        }

        if (face.owner_cell_id >= cells_.size()) {
            throw std::runtime_error("Mesh::Validate: face owner_cell_id is out of range");
        }

        if (face.IsInternal() && face.neighbor_cell_id >= cells_.size()) {
            throw std::runtime_error("Mesh::Validate: internal face neighbor_cell_id is out of range");
        }

        if (face.IsInternal() && face.neighbor_cell_id == face.owner_cell_id) {
            throw std::runtime_error("Mesh::Validate: internal face owner and neighbor are identical");
        }

        if (face.IsBoundary() && face.boundary_tag < 0) {
            throw std::runtime_error("Mesh::Validate: boundary face must have non-negative boundary_tag");
        }

        if (face.IsInternal() && face.boundary_tag >= 0) {
            throw std::runtime_error("Mesh::Validate: internal face must not have boundary_tag");
        }

        for (const std::size_t node_id : face.node_ids) {
            if (node_id >= nodes_.size()) {
                throw std::runtime_error("Mesh::Validate: face node_id is out of range");
            }
        }
    }
}

void Mesh::ValidateCellConnectivity() const {
    for (const Cell& cell : cells_) {
        if (cell.node_ids.empty()) {
            throw std::runtime_error("Mesh::Validate: cell has empty node_ids");
        }

        if (cell.face_ids.empty()) {
            throw std::runtime_error("Mesh::Validate: cell has empty face_ids");
        }

        for (const std::size_t node_id : cell.node_ids) {
            if (node_id >= nodes_.size()) {
                throw std::runtime_error("Mesh::Validate: cell node_id is out of range");
            }
        }

        for (const std::size_t face_id : cell.face_ids) {
            if (face_id >= faces_.size()) {
                throw std::runtime_error("Mesh::Validate: cell face_id is out of range");
            }

            const Face& face = faces_[face_id];
            if (face.owner_cell_id != cell.id && face.neighbor_cell_id != cell.id) {
                throw std::runtime_error("Mesh::Validate: cell does not belong to one of its listed faces");
            }
        }
    }
}

void Mesh::ValidateGeometry() const {
    constexpr double eps = 1e-14;

    for (const Face& face : faces_) {
        if (!(face.measure > 0.0)) {
            throw std::runtime_error("Mesh::Validate: face measure must be > 0");
        }

        const double normal_norm =
            std::sqrt(face.normal_x * face.normal_x +
                      face.normal_y * face.normal_y +
                      face.normal_z * face.normal_z);

        if (std::abs(normal_norm - 1.0) > 1e-10) {
            throw std::runtime_error("Mesh::Validate: face normal must be unit-length");
        }
    }

    for (const Cell& cell : cells_) {
        if (!(cell.volume > 0.0)) {
            throw std::runtime_error("Mesh::Validate: cell volume must be > 0");
        }

        if (cell.face_ids.size() < 2) {
            throw std::runtime_error("Mesh::Validate: cell must have at least 2 faces");
        }

        if (dim_ == 2 && cell.face_ids.size() < 3) {
            throw std::runtime_error("Mesh::Validate: 2D cell must have at least 3 faces");
        }

        if (dim_ == 3 && cell.face_ids.size() < 4) {
            throw std::runtime_error("Mesh::Validate: 3D cell must have at least 4 faces");
        }

        const bool center_is_finite =
            std::isfinite(cell.center_x) &&
            std::isfinite(cell.center_y) &&
            std::isfinite(cell.center_z);

        if (!center_is_finite) {
            throw std::runtime_error("Mesh::Validate: cell center contains non-finite value");
        }
    }

    for (const Node& node : nodes_) {
        const bool coords_are_finite =
            std::isfinite(node.x) &&
            std::isfinite(node.y) &&
            std::isfinite(node.z);

        if (!coords_are_finite) {
            throw std::runtime_error("Mesh::Validate: node coordinates contain non-finite value");
        }
    }

    for (const Face& face : faces_) {
        const bool center_is_finite =
            std::isfinite(face.center_x) &&
            std::isfinite(face.center_y) &&
            std::isfinite(face.center_z);

        if (!center_is_finite) {
            throw std::runtime_error("Mesh::Validate: face center contains non-finite value");
        }

        const bool normal_is_finite =
            std::isfinite(face.normal_x) &&
            std::isfinite(face.normal_y) &&
            std::isfinite(face.normal_z);

        if (!normal_is_finite) {
            throw std::runtime_error("Mesh::Validate: face normal contains non-finite value");
        }

        if (std::abs(face.normal_x) < eps &&
            std::abs(face.normal_y) < eps &&
            std::abs(face.normal_z) < eps) {
            throw std::runtime_error("Mesh::Validate: face normal must be non-zero");
        }
    }
}
