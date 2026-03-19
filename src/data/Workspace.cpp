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

// -------------------- cell-centered getters --------------------

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

xt::xtensor<double, 3>& Workspace::Temperature() {
    return temperature_;
}

const xt::xtensor<double, 3>& Workspace::Temperature() const {
    return temperature_;
}

xt::xtensor<double, 3>& Workspace::InternalEnergy() {
    return internal_energy_;
}

const xt::xtensor<double, 3>& Workspace::InternalEnergy() const {
    return internal_energy_;
}

xt::xtensor<double, 3>& Workspace::Q() {
    return q_;
}

const xt::xtensor<double, 3>& Workspace::Q() const {
    return q_;
}

xt::xtensor<double, 4>& Workspace::D() {
    return D_;
}

const xt::xtensor<double, 4>& Workspace::D() const {
    return D_;
}

// -------------------- face-centered getters --------------------

xt::xtensor<double, 3>& Workspace::Ux() {
    return ux_;
}

const xt::xtensor<double, 3>& Workspace::Ux() const {
    return ux_;
}

xt::xtensor<double, 3>& Workspace::Vy() {
    return vy_;
}

const xt::xtensor<double, 3>& Workspace::Vy() const {
    return vy_;
}

xt::xtensor<double, 3>& Workspace::Wz() {
    return wz_;
}

const xt::xtensor<double, 3>& Workspace::Wz() const {
    return wz_;
}

xt::xtensor<double, 3>& Workspace::UxOld() {
    return ux_old_;
}

const xt::xtensor<double, 3>& Workspace::UxOld() const {
    return ux_old_;
}

xt::xtensor<double, 3>& Workspace::VyOld() {
    return vy_old_;
}

const xt::xtensor<double, 3>& Workspace::VyOld() const {
    return vy_old_;
}

xt::xtensor<double, 3>& Workspace::WzOld() {
    return wz_old_;
}

const xt::xtensor<double, 3>& Workspace::WzOld() const {
    return wz_old_;
}

// -------------------- zero helpers --------------------

void Workspace::ZeroWc() {
    W_.fill(0.0);
}

void Workspace::ZeroRhs() {
    rhs_.fill(0.0);
}

void Workspace::ZeroTemperature() {
    temperature_.fill(0.0);
}

void Workspace::ZeroInternalEnergy() {
    internal_energy_.fill(0.0);
}

void Workspace::ZeroQ() {
    q_.fill(0.0);
}

void Workspace::ZeroD() {
    D_.fill(0.0);
}

void Workspace::ZeroUx() {
    ux_.fill(0.0);
}

void Workspace::ZeroVy() {
    vy_.fill(0.0);
}

void Workspace::ZeroWz() {
    wz_.fill(0.0);
}

void Workspace::ZeroUxOld() {
    ux_old_.fill(0.0);
}

void Workspace::ZeroVyOld() {
    vy_old_.fill(0.0);
}

void Workspace::ZeroWzOld() {
    wz_old_.fill(0.0);
}

void Workspace::ZeroAll() {
    ZeroWc();
    ZeroRhs();
    ZeroTemperature();
    ZeroInternalEnergy();
    ZeroQ();
    ZeroD();
    ZeroUx();
    ZeroVy();
    ZeroWz();
    ZeroUxOld();
    ZeroVyOld();
    ZeroWzOld();
}

bool Workspace::IsAllocated() const {
    return
        W_.dimension() == 4 &&
        rhs_.dimension() == 4 &&
        temperature_.dimension() == 3 &&
        internal_energy_.dimension() == 3 &&
        q_.dimension() == 3 &&
        D_.dimension() == 4 &&
        ux_.dimension() == 3 &&
        vy_.dimension() == 3 &&
        wz_.dimension() == 3 &&
        ux_old_.dimension() == 3 &&
        vy_old_.dimension() == 3 &&
        wz_old_.dimension() == 3 &&
        sx_ > 0 && sy_ > 0 && sz_ > 0;
}

void Workspace::Allocate(const std::size_t sx, const std::size_t sy, const std::size_t sz) {
    sx_ = sx;
    sy_ = sy;
    sz_ = sz;

    // cell-centered
    W_ = xt::zeros<double>({k_nvar, sx_, sy_, sz_});
    rhs_ = xt::zeros<double>({k_nvar, sx_, sy_, sz_});

    temperature_ = xt::zeros<double>({sx_, sy_, sz_});
    internal_energy_ = xt::zeros<double>({sx_, sy_, sz_});
    q_ = xt::zeros<double>({sx_, sy_, sz_});
    D_ = xt::zeros<double>({k_ndelta, sx_, sy_, sz_});

    // face-centered velocities
    ux_ = xt::zeros<double>({sx_ + 1, sy_, sz_});
    vy_ = xt::zeros<double>({sx_, sy_ + 1, sz_});
    wz_ = xt::zeros<double>({sx_, sy_, sz_ + 1});

    ux_old_ = xt::zeros<double>({sx_ + 1, sy_, sz_});
    vy_old_ = xt::zeros<double>({sx_, sy_ + 1, sz_});
    wz_old_ = xt::zeros<double>({sx_, sy_, sz_ + 1});
}
