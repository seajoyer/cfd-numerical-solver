#include "data/DataLayer.hpp"

#include "geometry/Mesh.hpp"

DataLayer::DataLayer(const std::size_t n_cells) {
    Resize(n_cells);
}

void DataLayer::Resize(const std::size_t n_cells) {
    if (n_cells == n_cells_ && IsAllocated()) {
        return;
    }

    Allocate(n_cells);
}

void DataLayer::ResizeFrom(const Mesh& mesh) {
    Resize(mesh.GetCellCount());
}

xt::xtensor<double, 2>& DataLayer::U() {
    return U_;
}

const xt::xtensor<double, 2>& DataLayer::U() const {
    return U_;
}

xt::xtensor<double, 1>& DataLayer::ReactantMassFraction() {
    return reactant_mass_fraction_;
}

const xt::xtensor<double, 1>& DataLayer::ReactantMassFraction() const {
    return reactant_mass_fraction_;
}

std::size_t DataLayer::GetCellCount() const {
    return n_cells_;
}

bool DataLayer::IsAllocated() const {
    return
        U_.dimension() == 2 &&
        reactant_mass_fraction_.dimension() == 1 &&
        U_.shape()[0] == n_cells_ &&
        U_.shape()[1] == k_nvar &&
        reactant_mass_fraction_.shape()[0] == n_cells_;
}

void DataLayer::Allocate(const std::size_t n_cells) {
    n_cells_ = n_cells;

    U_ = xt::zeros<double>({n_cells_, k_nvar});
    reactant_mass_fraction_ = xt::zeros<double>({n_cells_});
}
