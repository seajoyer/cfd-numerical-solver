#include "data/Workspace.hpp"

#include "geometry/Mesh.hpp"

void Workspace::ResizeFrom(const Mesh& mesh) {
    const std::size_t n_cells = mesh.GetCellCount();
    const std::size_t n_faces = mesh.GetFaceCount();

    if (n_cells == 0) {
        throw std::invalid_argument("Workspace::ResizeFrom: mesh has zero cells");
    }
    if (n_faces == 0) {
        throw std::invalid_argument("Workspace::ResizeFrom: mesh has zero faces");
    }

    if (n_cells == n_cells_ && n_faces == n_faces_ && IsAllocated()) {
        return;
    }

    Allocate(n_cells, n_faces);
}

// -------------------- cell-centered getters --------------------

xt::xtensor<double, 2>& Workspace::W() {
    return W_;
}

const xt::xtensor<double, 2>& Workspace::W() const {
    return W_;
}

xt::xtensor<double, 2>& Workspace::Rhs() {
    return rhs_;
}

const xt::xtensor<double, 2>& Workspace::Rhs() const {
    return rhs_;
}

xt::xtensor<double, 1>& Workspace::Temperature() {
    return temperature_;
}

const xt::xtensor<double, 1>& Workspace::Temperature() const {
    return temperature_;
}

xt::xtensor<double, 1>& Workspace::InternalEnergy() {
    return internal_energy_;
}

const xt::xtensor<double, 1>& Workspace::InternalEnergy() const {
    return internal_energy_;
}

xt::xtensor<double, 1>& Workspace::Q() {
    return q_;
}

const xt::xtensor<double, 1>& Workspace::Q() const {
    return q_;
}

xt::xtensor<double, 2>& Workspace::D() {
    return D_;
}

const xt::xtensor<double, 2>& Workspace::D() const {
    return D_;
}

// -------------------- face-centered getters --------------------

xt::xtensor<double, 1>& Workspace::Ux() {
    return ux_;
}

const xt::xtensor<double, 1>& Workspace::Ux() const {
    return ux_;
}

xt::xtensor<double, 1>& Workspace::Vy() {
    return vy_;
}

const xt::xtensor<double, 1>& Workspace::Vy() const {
    return vy_;
}

xt::xtensor<double, 1>& Workspace::Wz() {
    return wz_;
}

const xt::xtensor<double, 1>& Workspace::Wz() const {
    return wz_;
}

xt::xtensor<double, 1>& Workspace::UxOld() {
    return ux_old_;
}

const xt::xtensor<double, 1>& Workspace::UxOld() const {
    return ux_old_;
}

xt::xtensor<double, 1>& Workspace::VyOld() {
    return vy_old_;
}

const xt::xtensor<double, 1>& Workspace::VyOld() const {
    return vy_old_;
}

xt::xtensor<double, 1>& Workspace::WzOld() {
    return wz_old_;
}

const xt::xtensor<double, 1>& Workspace::WzOld() const {
    return wz_old_;
}

// -------------------- zero helpers --------------------

void Workspace::ZeroW() {
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
    ZeroW();
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
        W_.dimension() == 2 &&
        rhs_.dimension() == 2 &&
        temperature_.dimension() == 1 &&
        internal_energy_.dimension() == 1 &&
        q_.dimension() == 1 &&
        D_.dimension() == 2 &&
        ux_.dimension() == 1 &&
        vy_.dimension() == 1 &&
        wz_.dimension() == 1 &&
        ux_old_.dimension() == 1 &&
        vy_old_.dimension() == 1 &&
        wz_old_.dimension() == 1 &&
        n_cells_ > 0 &&
        n_faces_ > 0;
}

std::size_t Workspace::GetCellCount() const {
    return n_cells_;
}

std::size_t Workspace::GetFaceCount() const {
    return n_faces_;
}

void Workspace::Allocate(const std::size_t n_cells, const std::size_t n_faces) {
    n_cells_ = n_cells;
    n_faces_ = n_faces;

    // cell-centered
    W_ = xt::zeros<double>({n_cells_, k_nvar});
    rhs_ = xt::zeros<double>({n_cells_, k_nvar});

    temperature_ = xt::zeros<double>({n_cells_});
    internal_energy_ = xt::zeros<double>({n_cells_});
    q_ = xt::zeros<double>({n_cells_});
    D_ = xt::zeros<double>({n_cells_, k_ndelta});

    // face-centered velocities
    ux_ = xt::zeros<double>({n_faces_});
    vy_ = xt::zeros<double>({n_faces_});
    wz_ = xt::zeros<double>({n_faces_});

    ux_old_ = xt::zeros<double>({n_faces_});
    vy_old_ = xt::zeros<double>({n_faces_});
    wz_old_ = xt::zeros<double>({n_faces_});
}
