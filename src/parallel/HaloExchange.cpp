#include "parallel/HaloExchange.hpp"

#include <stdexcept>

HaloExchange::HaloExchange(MPI_Comm comm, int size, bool exchange_reactant_mass_fraction)
    : comm_(comm),
      size_(size),
      exchange_reactant_mass_fraction_(exchange_reactant_mass_fraction) {}

void HaloExchange::Exchange(DataLayer& layer, const Mesh& mesh) const {
    const int ng = mesh.GetPadding();
    if (ng <= 0 || size_ <= 1) {
        return;
    }

    ExchangeX(layer, mesh);

    if (mesh.GetDim() >= 2) {
        ExchangeY(layer, mesh);
    }

    if (mesh.GetDim() >= 3) {
        ExchangeZ(layer, mesh);
    }
}

void HaloExchange::ExchangeX(DataLayer& layer, const Mesh& mesh) const {
    const int left_rank = mesh.GetNeighborRank(Axis::X, Side::Left);
    const int right_rank = mesh.GetNeighborRank(Axis::X, Side::Right);

    const int i0 = mesh.GetCoreStartX();
    const int i1 = mesh.GetCoreEndExclusiveX();
    const int ng = mesh.GetPadding();
    (void)ng;

    if (left_rank >= 0) {
        const std::vector<double> send_left = PackX(layer, mesh, i0);
        std::vector<double> recv_left(send_left.size(), 0.0);

        MPI_Sendrecv(send_left.data(),
                     static_cast<int>(send_left.size()),
                     MPI_DOUBLE,
                     left_rank,
                     100,
                     recv_left.data(),
                     static_cast<int>(recv_left.size()),
                     MPI_DOUBLE,
                     left_rank,
                     101,
                     comm_,
                     MPI_STATUS_IGNORE);

        UnpackX(layer, mesh, i0 - mesh.GetPadding(), recv_left);
    }

    if (right_rank >= 0) {
        const std::vector<double> send_right = PackX(layer, mesh, i1 - mesh.GetPadding());
        std::vector<double> recv_right(send_right.size(), 0.0);

        MPI_Sendrecv(send_right.data(),
                     static_cast<int>(send_right.size()),
                     MPI_DOUBLE,
                     right_rank,
                     101,
                     recv_right.data(),
                     static_cast<int>(recv_right.size()),
                     MPI_DOUBLE,
                     right_rank,
                     100,
                     comm_,
                     MPI_STATUS_IGNORE);

        UnpackX(layer, mesh, i1, recv_right);
    }
}

void HaloExchange::ExchangeY(DataLayer& layer, const Mesh& mesh) const {
    const int left_rank = mesh.GetNeighborRank(Axis::Y, Side::Left);
    const int right_rank = mesh.GetNeighborRank(Axis::Y, Side::Right);

    const int j0 = mesh.GetCoreStartY();
    const int j1 = mesh.GetCoreEndExclusiveY();

    if (left_rank >= 0) {
        const std::vector<double> send_left = PackY(layer, mesh, j0);
        std::vector<double> recv_left(send_left.size(), 0.0);

        MPI_Sendrecv(send_left.data(),
                     static_cast<int>(send_left.size()),
                     MPI_DOUBLE,
                     left_rank,
                     200,
                     recv_left.data(),
                     static_cast<int>(recv_left.size()),
                     MPI_DOUBLE,
                     left_rank,
                     201,
                     comm_,
                     MPI_STATUS_IGNORE);

        UnpackY(layer, mesh, j0 - mesh.GetPadding(), recv_left);
    }

    if (right_rank >= 0) {
        const std::vector<double> send_right = PackY(layer, mesh, j1 - mesh.GetPadding());
        std::vector<double> recv_right(send_right.size(), 0.0);

        MPI_Sendrecv(send_right.data(),
                     static_cast<int>(send_right.size()),
                     MPI_DOUBLE,
                     right_rank,
                     201,
                     recv_right.data(),
                     static_cast<int>(recv_right.size()),
                     MPI_DOUBLE,
                     right_rank,
                     200,
                     comm_,
                     MPI_STATUS_IGNORE);

        UnpackY(layer, mesh, j1, recv_right);
    }
}

void HaloExchange::ExchangeZ(DataLayer& layer, const Mesh& mesh) const {
    const int left_rank = mesh.GetNeighborRank(Axis::Z, Side::Left);
    const int right_rank = mesh.GetNeighborRank(Axis::Z, Side::Right);

    const int k0 = mesh.GetCoreStartZ();
    const int k1 = mesh.GetCoreEndExclusiveZ();

    if (left_rank >= 0) {
        const std::vector<double> send_left = PackZ(layer, mesh, k0);
        std::vector<double> recv_left(send_left.size(), 0.0);

        MPI_Sendrecv(send_left.data(),
                     static_cast<int>(send_left.size()),
                     MPI_DOUBLE,
                     left_rank,
                     300,
                     recv_left.data(),
                     static_cast<int>(recv_left.size()),
                     MPI_DOUBLE,
                     left_rank,
                     301,
                     comm_,
                     MPI_STATUS_IGNORE);

        UnpackZ(layer, mesh, k0 - mesh.GetPadding(), recv_left);
    }

    if (right_rank >= 0) {
        const std::vector<double> send_right = PackZ(layer, mesh, k1 - mesh.GetPadding());
        std::vector<double> recv_right(send_right.size(), 0.0);

        MPI_Sendrecv(send_right.data(),
                     static_cast<int>(send_right.size()),
                     MPI_DOUBLE,
                     right_rank,
                     301,
                     recv_right.data(),
                     static_cast<int>(recv_right.size()),
                     MPI_DOUBLE,
                     right_rank,
                     300,
                     comm_,
                     MPI_STATUS_IGNORE);

        UnpackZ(layer, mesh, k1, recv_right);
    }
}

std::vector<double> HaloExchange::PackX(const DataLayer& layer, const Mesh& mesh, const int i_begin) const {
    const int ng = mesh.GetPadding();
    const int sy = mesh.GetSy();
    const int sz = mesh.GetSz();

    const std::size_t n_fields =
        DataLayer::k_nvar + (exchange_reactant_mass_fraction_ ? 1u : 0u);

    const std::size_t count =
        n_fields *
        static_cast<std::size_t>(ng) *
        static_cast<std::size_t>(sy) *
        static_cast<std::size_t>(sz);

    std::vector<double> buffer(count, 0.0);

    const auto& U = layer.U();
    std::size_t p = 0;

    for (std::size_t v = 0; v < DataLayer::k_nvar; ++v) {
        for (int i = i_begin; i < i_begin + ng; ++i) {
            for (int j = 0; j < sy; ++j) {
                for (int k = 0; k < sz; ++k) {
                    buffer[p++] = U(v, i, j, k);
                }
            }
        }
    }

    if (exchange_reactant_mass_fraction_) {
        const auto& lambda = layer.ReactantMassFraction();

        for (int i = i_begin; i < i_begin + ng; ++i) {
            for (int j = 0; j < sy; ++j) {
                for (int k = 0; k < sz; ++k) {
                    buffer[p++] = lambda(i, j, k);
                }
            }
        }
    }

    return buffer;
}

std::vector<double> HaloExchange::PackY(const DataLayer& layer, const Mesh& mesh, const int j_begin) const {
    const int ng = mesh.GetPadding();
    const int sx = mesh.GetSx();
    const int sz = mesh.GetSz();

    const std::size_t n_fields =
        DataLayer::k_nvar + (exchange_reactant_mass_fraction_ ? 1u : 0u);

    const std::size_t count =
        n_fields *
        static_cast<std::size_t>(sx) *
        static_cast<std::size_t>(ng) *
        static_cast<std::size_t>(sz);

    std::vector<double> buffer(count, 0.0);

    const auto& U = layer.U();
    std::size_t p = 0;

    for (std::size_t v = 0; v < DataLayer::k_nvar; ++v) {
        for (int i = 0; i < sx; ++i) {
            for (int j = j_begin; j < j_begin + ng; ++j) {
                for (int k = 0; k < sz; ++k) {
                    buffer[p++] = U(v, i, j, k);
                }
            }
        }
    }

    if (exchange_reactant_mass_fraction_) {
        const auto& lambda = layer.ReactantMassFraction();

        for (int i = 0; i < sx; ++i) {
            for (int j = j_begin; j < j_begin + ng; ++j) {
                for (int k = 0; k < sz; ++k) {
                    buffer[p++] = lambda(i, j, k);
                }
            }
        }
    }

    return buffer;
}

std::vector<double> HaloExchange::PackZ(const DataLayer& layer, const Mesh& mesh, const int k_begin) const {
    const int ng = mesh.GetPadding();
    const int sx = mesh.GetSx();
    const int sy = mesh.GetSy();

    const std::size_t n_fields =
        DataLayer::k_nvar + (exchange_reactant_mass_fraction_ ? 1u : 0u);

    const std::size_t count =
        n_fields *
        static_cast<std::size_t>(sx) *
        static_cast<std::size_t>(sy) *
        static_cast<std::size_t>(ng);

    std::vector<double> buffer(count, 0.0);

    const auto& U = layer.U();
    std::size_t p = 0;

    for (std::size_t v = 0; v < DataLayer::k_nvar; ++v) {
        for (int i = 0; i < sx; ++i) {
            for (int j = 0; j < sy; ++j) {
                for (int k = k_begin; k < k_begin + ng; ++k) {
                    buffer[p++] = U(v, i, j, k);
                }
            }
        }
    }

    if (exchange_reactant_mass_fraction_) {
        const auto& lambda = layer.ReactantMassFraction();

        for (int i = 0; i < sx; ++i) {
            for (int j = 0; j < sy; ++j) {
                for (int k = k_begin; k < k_begin + ng; ++k) {
                    buffer[p++] = lambda(i, j, k);
                }
            }
        }
    }

    return buffer;
}

void HaloExchange::UnpackX(DataLayer& layer, const Mesh& mesh, const int i_begin,
                           const std::vector<double>& buffer) const {
    const int ng = mesh.GetPadding();
    const int sy = mesh.GetSy();
    const int sz = mesh.GetSz();

    const std::size_t n_fields =
        DataLayer::k_nvar + (exchange_reactant_mass_fraction_ ? 1u : 0u);

    const std::size_t expected =
        n_fields *
        static_cast<std::size_t>(ng) *
        static_cast<std::size_t>(sy) *
        static_cast<std::size_t>(sz);

    if (buffer.size() != expected) {
        throw std::runtime_error("HaloExchange::UnpackX: invalid buffer size");
    }

    auto& U = layer.U();
    std::size_t p = 0;

    for (std::size_t v = 0; v < DataLayer::k_nvar; ++v) {
        for (int i = i_begin; i < i_begin + ng; ++i) {
            for (int j = 0; j < sy; ++j) {
                for (int k = 0; k < sz; ++k) {
                    U(v, i, j, k) = buffer[p++];
                }
            }
        }
    }

    if (exchange_reactant_mass_fraction_) {
        auto& lambda = layer.ReactantMassFraction();

        for (int i = i_begin; i < i_begin + ng; ++i) {
            for (int j = 0; j < sy; ++j) {
                for (int k = 0; k < sz; ++k) {
                    lambda(i, j, k) = buffer[p++];
                }
            }
        }
    }
}

void HaloExchange::UnpackY(DataLayer& layer, const Mesh& mesh, const int j_begin,
                           const std::vector<double>& buffer) const {
    const int ng = mesh.GetPadding();
    const int sx = mesh.GetSx();
    const int sz = mesh.GetSz();

    const std::size_t n_fields =
        DataLayer::k_nvar + (exchange_reactant_mass_fraction_ ? 1u : 0u);

    const std::size_t expected =
        n_fields *
        static_cast<std::size_t>(sx) *
        static_cast<std::size_t>(ng) *
        static_cast<std::size_t>(sz);

    if (buffer.size() != expected) {
        throw std::runtime_error("HaloExchange::UnpackY: invalid buffer size");
    }

    auto& U = layer.U();
    std::size_t p = 0;

    for (std::size_t v = 0; v < DataLayer::k_nvar; ++v) {
        for (int i = 0; i < sx; ++i) {
            for (int j = j_begin; j < j_begin + ng; ++j) {
                for (int k = 0; k < sz; ++k) {
                    U(v, i, j, k) = buffer[p++];
                }
            }
        }
    }

    if (exchange_reactant_mass_fraction_) {
        auto& lambda = layer.ReactantMassFraction();

        for (int i = 0; i < sx; ++i) {
            for (int j = j_begin; j < j_begin + ng; ++j) {
                for (int k = 0; k < sz; ++k) {
                    lambda(i, j, k) = buffer[p++];
                }
            }
        }
    }
}

void HaloExchange::UnpackZ(DataLayer& layer, const Mesh& mesh, const int k_begin,
                           const std::vector<double>& buffer) const {
    const int ng = mesh.GetPadding();
    const int sx = mesh.GetSx();
    const int sy = mesh.GetSy();

    const std::size_t n_fields =
        DataLayer::k_nvar + (exchange_reactant_mass_fraction_ ? 1u : 0u);

    const std::size_t expected =
        n_fields *
        static_cast<std::size_t>(sx) *
        static_cast<std::size_t>(sy) *
        static_cast<std::size_t>(ng);

    if (buffer.size() != expected) {
        throw std::runtime_error("HaloExchange::UnpackZ: invalid buffer size");
    }

    auto& U = layer.U();
    std::size_t p = 0;

    for (std::size_t v = 0; v < DataLayer::k_nvar; ++v) {
        for (int i = 0; i < sx; ++i) {
            for (int j = 0; j < sy; ++j) {
                for (int k = k_begin; k < k_begin + ng; ++k) {
                    U(v, i, j, k) = buffer[p++];
                }
            }
        }
    }

    if (exchange_reactant_mass_fraction_) {
        auto& lambda = layer.ReactantMassFraction();

        for (int i = 0; i < sx; ++i) {
            for (int j = 0; j < sy; ++j) {
                for (int k = k_begin; k < k_begin + ng; ++k) {
                    lambda(i, j, k) = buffer[p++];
                }
            }
        }
    }
}