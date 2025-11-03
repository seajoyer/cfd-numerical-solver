#include "data/DataLayer.hpp"

DataLayer::DataLayer(const int N, const int padding)
    : n(N), nGhostCells(padding) {
    RecomputeSizes();
    Allocate1D();
}

DataLayer::DataLayer(const int N, const int padding, const int dim)
    : n(N), nGhostCells(padding), dimention(dim) {
    RecomputeSizes();
    Allocate1D();
}

void DataLayer::SetN(const int new_N) {
    n = new_N;
    RecomputeSizes();
    Allocate1D();
}

void DataLayer::SetPadding(const int new_padding) {
    nGhostCells = new_padding;
    RecomputeSizes();
    Allocate1D();
}

void DataLayer::SetDim(const int new_dim) {
    dimention = new_dim;
    RecomputeSizes();
    Allocate1D();
}

void DataLayer::RecomputeSizes() {
    totalSize = n + 2 * nGhostCells;
    if (totalSize < 0) totalSize = 0;
}

void DataLayer::Allocate1D() {
    std::vector<std::size_t> shape = {static_cast<std::size_t>(totalSize)};
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
