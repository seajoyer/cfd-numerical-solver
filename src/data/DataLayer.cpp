#include "data/DataLayer.hpp"

DataLayer::DataLayer(const int N, const int padding)
    : n_(N), n_ghost_cells_(padding) {
    RecomputeSizes();
    Allocate1D();
}

DataLayer::DataLayer(const int N, const int padding, const int dim)
    : n_(N), n_ghost_cells_(padding), dimention_(dim) {
    RecomputeSizes();
    Allocate1D();
}

void DataLayer::SetN(const int new_N) {
    n_ = new_N;
    RecomputeSizes();
    Allocate1D();
}

void DataLayer::SetPadding(const int new_padding) {
    n_ghost_cells_ = new_padding;
    RecomputeSizes();
    Allocate1D();
}

void DataLayer::SetDim(const int new_dim) {
    dimention_ = new_dim;
    RecomputeSizes();
    Allocate1D();
}

void DataLayer::RecomputeSizes() {
    total_size_ = n_ + 2 * n_ghost_cells_;
    if (total_size_ < 0) total_size_ = 0;
}

void DataLayer::Allocate1D() {
    std::vector<std::size_t> shape = {static_cast<std::size_t>(total_size_)};
    rho.resize(shape);
    u.resize(shape);
    P.resize(shape);
    p.resize(shape);
    e.resize(shape);
    U.resize(shape);
    V.resize(shape);
    m.resize(shape);
    xb.resize(shape);
    xc.resize(shape);
}

auto DataLayer::GetCoreEndExclusive(const int axis) const -> int {
    (void)axis;
    return n_ghost_cells_ + n_;
}

auto DataLayer::GetCoreStart(const int axis) const -> int {
    (void)axis;
    return n_ghost_cells_;
}
