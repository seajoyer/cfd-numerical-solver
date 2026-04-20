#ifndef HALOSERIALIZATION_HPP
#define HALOSERIALIZATION_HPP

#include <cstddef>
#include <vector>

#include "parallel/DomainDecomposition.hpp"

class HaloSerialization final {
public:
    [[nodiscard]] static std::vector<char> Pack(
        const std::vector<DomainDecomposition::NeighborHalo>& halos
    );

    [[nodiscard]] static std::vector<DomainDecomposition::NeighborHalo> Unpack(
        const std::vector<char>& buffer
    );

private:
    template <typename T>
    static void PackPod(std::vector<char>& buffer, const T& value);

    template <typename T>
    static T UnpackPod(const std::vector<char>& buffer, std::size_t& offset);

    static void PackSizeTVector(std::vector<char>& buffer,
                                const std::vector<std::size_t>& values);

    static std::vector<std::size_t> UnpackSizeTVector(const std::vector<char>& buffer,
                                                      std::size_t& offset);
};

#endif  // HALOSERIALIZATION_HPP
