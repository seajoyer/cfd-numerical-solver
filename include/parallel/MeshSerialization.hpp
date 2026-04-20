#ifndef MESHSERIALIZATION_HPP
#define MESHSERIALIZATION_HPP

#include <cstddef>
#include <vector>

#include "geometry/Mesh.hpp"

class MeshSerialization final {
public:
    [[nodiscard]] static std::vector<char> Pack(const Mesh& mesh);
    [[nodiscard]] static Mesh Unpack(const std::vector<char>& buffer);

private:
    template <typename T>
    static void PackPod(std::vector<char>& buffer, const T& value);

    template <typename T>
    static T UnpackPod(const std::vector<char>& buffer, std::size_t& offset);

    static void PackNode(std::vector<char>& buffer, const Node& node);
    static Node UnpackNode(const std::vector<char>& buffer, std::size_t& offset);

    static void PackFace(std::vector<char>& buffer, const Face& face);
    static Face UnpackFace(const std::vector<char>& buffer, std::size_t& offset);

    static void PackCell(std::vector<char>& buffer, const Cell& cell);
    static Cell UnpackCell(const std::vector<char>& buffer, std::size_t& offset);

    static void PackSizeTVector(std::vector<char>& buffer,
                                const std::vector<std::size_t>& values);

    static std::vector<std::size_t> UnpackSizeTVector(const std::vector<char>& buffer,
                                                      std::size_t& offset);
};

#endif  // MESHSERIALIZATION_HPP
