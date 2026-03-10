#include "data/Workspace.hpp"

#include "data/Mesh.hpp"

void Workspace::ResizeFrom(const Mesh& mesh) {
    const std::size_t sx = static_cast<std::size_t>(mesh.GetSx());
    const std::size_t sy = static_cast<std::size_t>(mesh.GetSy());
    const std::size_t sz = static_cast<std::size_t>(mesh.GetSz());

    if (sx == sx_ && sy == sy_ && sz == sz_ && IsAllocated()) {
        return;
    }

    Allocate(sx, sy, sz);
}

xt::xtensor<double, 4>& Workspace::W() {
    return W_;
}

const xt::xtensor<double, 4>& Workspace::W() const {
    return W_;
}

xt::xtensor<double, 4>& Workspace::Rhs() {
    return rhs_;
}

const xt::xtensor<double, 4>& Workspace::Rhs() const {
    return rhs_;
}

void Workspace::ZeroRhs() {
    rhs_.fill(0.0);
}

void Workspace::ZeroW() {
    W_.fill(0.0);
}

bool Workspace::IsAllocated() const {
    return W_.dimension() == 4 && rhs_.dimension() == 4 && sx_ > 0 && sy_ > 0 && sz_ > 0;
}

void Workspace::Allocate(const std::size_t sx, const std::size_t sy, const std::size_t sz) {
    sx_ = sx;
    sy_ = sy;
    sz_ = sz;

    W_ = xt::zeros<double>({k_nvar, sx_, sy_, sz_});
    rhs_ = xt::zeros<double>({k_nvar, sx_, sy_, sz_});
}