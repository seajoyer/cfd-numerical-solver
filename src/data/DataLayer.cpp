#include "data/DataLayer.hpp"
#include <stdexcept>
#include <vector>

DataLayer::DataLayer(const int N, const int padding)
    : n_(N), n_ghost_cells_(padding)  {
    if (N <= 0) {
        throw std::invalid_argument("Number of cells N must be positive");
    }
    if (padding < 0) {
        throw std::invalid_argument("Padding must be non-negative");
    }
    RecomputeSizes();
    Allocate1D();
}

DataLayer::DataLayer(const int N, const int padding, const int dim)
    : n_(N), n_ghost_cells_(padding), dimension_(dim) {
    if (N <= 0) {
        throw std::invalid_argument("Number of cells N must be positive");
    }
    if (padding < 0) {
        throw std::invalid_argument("Padding must be non-negative");
    }
    if (dim < 1 || dim > 3) {
        throw std::invalid_argument("Dimension must be 1, 2, or 3");
    }
    RecomputeSizes();
    Allocate1D();
}

void DataLayer::SetN(const int new_N) {
    if (new_N <= 0) {
        throw std::invalid_argument("Number of cells N must be positive");
    }
    n_ = new_N;
    RecomputeSizes();
    Allocate1D();
}

void DataLayer::SetPadding(const int new_padding) {
    if (new_padding < 0) {
        throw std::invalid_argument("Padding must be non-negative");
    }
    n_ghost_cells_ = new_padding;
    RecomputeSizes();
    Allocate1D();
}

void DataLayer::SetDim(const int new_dim) {
    if (new_dim < 1 || new_dim > 3) {
        throw std::invalid_argument("Dimension must be 1, 2, or 3");
    }
    dimension_ = new_dim;
    RecomputeSizes();
    Allocate1D();
}

void DataLayer::RecomputeSizes() {
    total_size_ = n_ + 2 * n_ghost_cells_;
    if (total_size_ < 0) {
        total_size_ = 0;
    }
}

void DataLayer::Allocate1D() {
    const std::vector<std::size_t> shape = {static_cast<std::size_t>(total_size_)};
    
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
    
    // Initialize arrays to zero
    rho.fill(0.0);
    u.fill(0.0);
    P.fill(0.0);
    p.fill(0.0);
    e.fill(0.0);
    U.fill(0.0);
    V.fill(0.0);
    m.fill(0.0);
    xb.fill(0.0);
    xc.fill(0.0);
}

auto DataLayer::GetCoreStart(const int axis) const -> int {
    (void)axis;  // TODO: Currently unused for 1D
    return n_ghost_cells_;
}

auto DataLayer::GetCoreEndExclusive(const int axis) const -> int {
    (void)axis;  // TODO: Currently unused for 1D
    return n_ghost_cells_ + n_;
}
