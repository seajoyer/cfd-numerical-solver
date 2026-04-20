#include "parallel/HaloSerialization.hpp"

#include <cstring>
#include <stdexcept>
#include <type_traits>

template <typename T>
void HaloSerialization::PackPod(std::vector<char>& buffer, const T& value) {
    static_assert(std::is_trivially_copyable_v<T>);
    const char* ptr = reinterpret_cast<const char*>(&value);
    buffer.insert(buffer.end(), ptr, ptr + sizeof(T));
}

template <typename T>
T HaloSerialization::UnpackPod(const std::vector<char>& buffer, std::size_t& offset) {
    static_assert(std::is_trivially_copyable_v<T>);

    if (offset + sizeof(T) > buffer.size()) {
        throw std::runtime_error("HaloSerialization: buffer underflow while reading POD");
    }

    T value{};
    std::memcpy(&value, buffer.data() + offset, sizeof(T));
    offset += sizeof(T);
    return value;
}

void HaloSerialization::PackSizeTVector(std::vector<char>& buffer,
                                        const std::vector<std::size_t>& values) {
    const std::size_t n = values.size();
    PackPod(buffer, n);

    if (!values.empty()) {
        const char* ptr = reinterpret_cast<const char*>(values.data());
        buffer.insert(buffer.end(), ptr, ptr + sizeof(std::size_t) * values.size());
    }
}

std::vector<std::size_t> HaloSerialization::UnpackSizeTVector(const std::vector<char>& buffer,
                                                              std::size_t& offset) {
    const std::size_t n = UnpackPod<std::size_t>(buffer, offset);
    std::vector<std::size_t> values(n);

    if (n > 0) {
        const std::size_t bytes = sizeof(std::size_t) * n;
        if (offset + bytes > buffer.size()) {
            throw std::runtime_error("HaloSerialization: buffer underflow while reading vector payload");
        }

        std::memcpy(values.data(), buffer.data() + offset, bytes);
        offset += bytes;
    }

    return values;
}

std::vector<char> HaloSerialization::Pack(
    const std::vector<DomainDecomposition::NeighborHalo>& halos
) {
    std::vector<char> buffer;
    buffer.reserve(512);

    const std::size_t n_halos = halos.size();
    PackPod(buffer, n_halos);

    for (const DomainDecomposition::NeighborHalo& halo : halos) {
        PackPod(buffer, halo.remote_rank);
        PackSizeTVector(buffer, halo.send_local_ids);
        PackSizeTVector(buffer, halo.recv_local_ids);
    }

    return buffer;
}

std::vector<DomainDecomposition::NeighborHalo> HaloSerialization::Unpack(
    const std::vector<char>& buffer
) {
    std::size_t offset = 0;

    const std::size_t n_halos = UnpackPod<std::size_t>(buffer, offset);
    std::vector<DomainDecomposition::NeighborHalo> halos;
    halos.reserve(n_halos);

    for (std::size_t i = 0; i < n_halos; ++i) {
        DomainDecomposition::NeighborHalo halo;
        halo.remote_rank = UnpackPod<int>(buffer, offset);
        halo.send_local_ids = UnpackSizeTVector(buffer, offset);
        halo.recv_local_ids = UnpackSizeTVector(buffer, offset);
        halos.push_back(std::move(halo));
    }

    if (offset != buffer.size()) {
        throw std::runtime_error("HaloSerialization: trailing bytes after unpack");
    }

    return halos;
}

// explicit instantiations
template void HaloSerialization::PackPod<int>(std::vector<char>&, const int&);
template void HaloSerialization::PackPod<std::size_t>(std::vector<char>&, const std::size_t&);

template int HaloSerialization::UnpackPod<int>(const std::vector<char>&, std::size_t&);
template std::size_t HaloSerialization::UnpackPod<std::size_t>(const std::vector<char>&, std::size_t&);
