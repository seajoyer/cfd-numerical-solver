#include "parallel/MeshSerialization.hpp"

#include <cstring>
#include <stdexcept>
#include <type_traits>

template <typename T>
void MeshSerialization::PackPod(std::vector<char>& buffer, const T& value) {
    static_assert(std::is_trivially_copyable_v<T>);
    const char* ptr = reinterpret_cast<const char*>(&value);
    buffer.insert(buffer.end(), ptr, ptr + sizeof(T));
}

template <typename T>
T MeshSerialization::UnpackPod(const std::vector<char>& buffer, std::size_t& offset) {
    static_assert(std::is_trivially_copyable_v<T>);

    if (offset + sizeof(T) > buffer.size()) {
        throw std::runtime_error("MeshSerialization: buffer underflow while reading POD");
    }

    T value{};
    std::memcpy(&value, buffer.data() + offset, sizeof(T));
    offset += sizeof(T);
    return value;
}

void MeshSerialization::PackSizeTVector(std::vector<char>& buffer,
                                        const std::vector<std::size_t>& values) {
    const std::size_t n = values.size();
    PackPod(buffer, n);

    if (!values.empty()) {
        const char* ptr = reinterpret_cast<const char*>(values.data());
        buffer.insert(buffer.end(), ptr, ptr + sizeof(std::size_t) * values.size());
    }
}

std::vector<std::size_t> MeshSerialization::UnpackSizeTVector(const std::vector<char>& buffer,
                                                              std::size_t& offset) {
    const std::size_t n = UnpackPod<std::size_t>(buffer, offset);
    std::vector<std::size_t> values(n);

    if (n > 0) {
        const std::size_t bytes = sizeof(std::size_t) * n;
        if (offset + bytes > buffer.size()) {
            throw std::runtime_error("MeshSerialization: buffer underflow while reading vector payload");
        }

        std::memcpy(values.data(), buffer.data() + offset, bytes);
        offset += bytes;
    }

    return values;
}

void MeshSerialization::PackNode(std::vector<char>& buffer, const Node& node) {
    PackPod(buffer, node.id);
    PackPod(buffer, node.x);
    PackPod(buffer, node.y);
    PackPod(buffer, node.z);
}

Node MeshSerialization::UnpackNode(const std::vector<char>& buffer, std::size_t& offset) {
    Node node;
    node.id = UnpackPod<std::size_t>(buffer, offset);
    node.x = UnpackPod<double>(buffer, offset);
    node.y = UnpackPod<double>(buffer, offset);
    node.z = UnpackPod<double>(buffer, offset);
    return node;
}

void MeshSerialization::PackFace(std::vector<char>& buffer, const Face& face) {
    PackPod(buffer, face.id);
    PackSizeTVector(buffer, face.node_ids);
    PackPod(buffer, face.owner_cell_id);
    PackPod(buffer, face.neighbor_cell_id);
    PackPod(buffer, face.kind);
    PackPod(buffer, face.remote_rank);
    PackPod(buffer, face.remote_cell_id);
    PackPod(buffer, face.center_x);
    PackPod(buffer, face.center_y);
    PackPod(buffer, face.center_z);
    PackPod(buffer, face.measure);
    PackPod(buffer, face.normal_x);
    PackPod(buffer, face.normal_y);
    PackPod(buffer, face.normal_z);
    PackPod(buffer, face.boundary_tag);
}

Face MeshSerialization::UnpackFace(const std::vector<char>& buffer, std::size_t& offset) {
    Face face;
    face.id = UnpackPod<std::size_t>(buffer, offset);
    face.node_ids = UnpackSizeTVector(buffer, offset);
    face.owner_cell_id = UnpackPod<std::size_t>(buffer, offset);
    face.neighbor_cell_id = UnpackPod<std::size_t>(buffer, offset);
    face.kind = UnpackPod<FaceKind>(buffer, offset);
    face.remote_rank = UnpackPod<int>(buffer, offset);
    face.remote_cell_id = UnpackPod<std::size_t>(buffer, offset);
    face.center_x = UnpackPod<double>(buffer, offset);
    face.center_y = UnpackPod<double>(buffer, offset);
    face.center_z = UnpackPod<double>(buffer, offset);
    face.measure = UnpackPod<double>(buffer, offset);
    face.normal_x = UnpackPod<double>(buffer, offset);
    face.normal_y = UnpackPod<double>(buffer, offset);
    face.normal_z = UnpackPod<double>(buffer, offset);
    face.boundary_tag = UnpackPod<int>(buffer, offset);
    return face;
}

void MeshSerialization::PackCell(std::vector<char>& buffer, const Cell& cell) {
    PackPod(buffer, cell.id);
    PackPod(buffer, cell.local_id);
    PackSizeTVector(buffer, cell.node_ids);
    PackSizeTVector(buffer, cell.face_ids);
    PackPod(buffer, cell.center_x);
    PackPod(buffer, cell.center_y);
    PackPod(buffer, cell.center_z);
    PackPod(buffer, cell.volume);
}

Cell MeshSerialization::UnpackCell(const std::vector<char>& buffer, std::size_t& offset) {
    Cell cell;
    cell.id = UnpackPod<std::size_t>(buffer, offset);
    cell.local_id = UnpackPod<std::size_t>(buffer, offset);
    cell.node_ids = UnpackSizeTVector(buffer, offset);
    cell.face_ids = UnpackSizeTVector(buffer, offset);
    cell.center_x = UnpackPod<double>(buffer, offset);
    cell.center_y = UnpackPod<double>(buffer, offset);
    cell.center_z = UnpackPod<double>(buffer, offset);
    cell.volume = UnpackPod<double>(buffer, offset);
    return cell;
}

std::vector<char> MeshSerialization::Pack(const Mesh& mesh) {
    std::vector<char> buffer;
    buffer.reserve(4096);

    PackPod(buffer, mesh.GetDim());
    PackPod(buffer, mesh.GetOwnedCellCount());
    PackPod(buffer, mesh.GetGhostCellCount());

    const std::size_t n_nodes = mesh.GetNodeCount();
    const std::size_t n_faces = mesh.GetFaceCount();
    const std::size_t n_cells = mesh.GetCellCount();

    PackPod(buffer, n_nodes);
    PackPod(buffer, n_faces);
    PackPod(buffer, n_cells);

    for (const Node& node : mesh.Nodes()) {
        PackNode(buffer, node);
    }

    for (const Face& face : mesh.Faces()) {
        PackFace(buffer, face);
    }

    for (const Cell& cell : mesh.Cells()) {
        PackCell(buffer, cell);
    }

    return buffer;
}

Mesh MeshSerialization::Unpack(const std::vector<char>& buffer) {
    std::size_t offset = 0;

    const int dim = UnpackPod<int>(buffer, offset);
    const std::size_t owned_cell_count = UnpackPod<std::size_t>(buffer, offset);
    const std::size_t ghost_cell_count = UnpackPod<std::size_t>(buffer, offset);

    const std::size_t n_nodes = UnpackPod<std::size_t>(buffer, offset);
    const std::size_t n_faces = UnpackPod<std::size_t>(buffer, offset);
    const std::size_t n_cells = UnpackPod<std::size_t>(buffer, offset);

    Mesh mesh(dim);
    mesh.SetOwnedCellCount(owned_cell_count);
    mesh.SetGhostCellCount(ghost_cell_count);

    auto& nodes = mesh.Nodes();
    auto& faces = mesh.Faces();
    auto& cells = mesh.Cells();

    nodes.reserve(n_nodes);
    faces.reserve(n_faces);
    cells.reserve(n_cells);

    for (std::size_t i = 0; i < n_nodes; ++i) {
        nodes.push_back(UnpackNode(buffer, offset));
    }

    for (std::size_t i = 0; i < n_faces; ++i) {
        faces.push_back(UnpackFace(buffer, offset));
    }

    for (std::size_t i = 0; i < n_cells; ++i) {
        cells.push_back(UnpackCell(buffer, offset));
    }

    if (offset != buffer.size()) {
        throw std::runtime_error("MeshSerialization: trailing bytes after unpack");
    }

    return mesh;
}

// explicit instantiations
template void MeshSerialization::PackPod<int>(std::vector<char>&, const int&);
template void MeshSerialization::PackPod<std::size_t>(std::vector<char>&, const std::size_t&);
template void MeshSerialization::PackPod<double>(std::vector<char>&, const double&);
template void MeshSerialization::PackPod<FaceKind>(std::vector<char>&, const FaceKind&);

template int MeshSerialization::UnpackPod<int>(const std::vector<char>&, std::size_t&);
template std::size_t MeshSerialization::UnpackPod<std::size_t>(const std::vector<char>&, std::size_t&);
template double MeshSerialization::UnpackPod<double>(const std::vector<char>&, std::size_t&);
template FaceKind MeshSerialization::UnpackPod<FaceKind>(const std::vector<char>&, std::size_t&);
