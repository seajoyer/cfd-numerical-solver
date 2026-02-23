#include "data/DataLayer.hpp"

#include <algorithm>
#include <cmath>

DataLayer::DataLayer(const int nx, const int ny, const int nz, const int padding, const int dim)
    : nx_(nx), ny_(ny), nz_(nz), ng_(padding), dim_(dim) {
    Validate();

    // enforce dimensional contraction (so metadata stays consistent everywhere)
    if (dim_ < 2) {
        ny_ = 1;
    }
    if (dim_ < 3) {
        nz_ = 1;
    }

    RecomputeSizes();
    Allocate();
}

int DataLayer::GetDim() const {
    return dim_;
}

int DataLayer::GetPadding() const {
    return ng_;
}

int DataLayer::GetNx() const {
    return nx_;
}

int DataLayer::GetNy() const {
    return ny_;
}

int DataLayer::GetNz() const {
    return nz_;
}

int DataLayer::GetSx() const {
    return sx_;
}

int DataLayer::GetSy() const {
    return sy_;
}

int DataLayer::GetSz() const {
    return sz_;
}

int DataLayer::GetCoreStartX() const {
    return ng_;
}

int DataLayer::GetCoreStartY() const {
    return dim_ >= 2 ? ng_ : 0;
}

int DataLayer::GetCoreStartZ() const {
    return dim_ >= 3 ? ng_ : 0;
}

int DataLayer::GetCoreEndExclusiveX() const {
    return ng_ + nx_;
}

int DataLayer::GetCoreEndExclusiveY() const {
    return dim_ >= 2 ? ng_ + ny_ : 1;
}

int DataLayer::GetCoreEndExclusiveZ() const {
    return dim_ >= 3 ? ng_ + nz_ : 1;
}

int DataLayer::GetCoreNx() const {
    return nx_;
}

int DataLayer::GetCoreNy() const {
    return ny_;
}

int DataLayer::GetCoreNz() const {
    return nz_;
}

xt::xtensor<double, 4>& DataLayer::U() {
    return U_;
}

const xt::xtensor<double, 4>& DataLayer::U() const {
    return U_;
}

xt::xtensor<double, 1>& DataLayer::Xb() {
    return xb_;
}

const xt::xtensor<double, 1>& DataLayer::Xb() const {
    return xb_;
}

xt::xtensor<double, 1>& DataLayer::Xc() {
    return xc_;
}

const xt::xtensor<double, 1>& DataLayer::Xc() const {
    return xc_;
}

xt::xtensor<double, 1>& DataLayer::Yb() {
    return yb_;
}

const xt::xtensor<double, 1>& DataLayer::Yb() const {
    return yb_;
}

xt::xtensor<double, 1>& DataLayer::Yc() {
    return yc_;
}

const xt::xtensor<double, 1>& DataLayer::Yc() const {
    return yc_;
}

xt::xtensor<double, 1>& DataLayer::Zb() {
    return zb_;
}

const xt::xtensor<double, 1>& DataLayer::Zb() const {
    return zb_;
}

xt::xtensor<double, 1>& DataLayer::Zc() {
    return zc_;
}

const xt::xtensor<double, 1>& DataLayer::Zc() const {
    return zc_;
}

const xt::xtensor<double, 1>& DataLayer::Dx() const {
    return dx_;
}

const xt::xtensor<double, 1>& DataLayer::Dy() const {
    return dy_;
}

const xt::xtensor<double, 1>& DataLayer::Dz() const {
    return dz_;
}

const xt::xtensor<double, 1>& DataLayer::InvDx() const {
    return inv_dx_;
}

const xt::xtensor<double, 1>& DataLayer::InvDy() const {
    return inv_dy_;
}

const xt::xtensor<double, 1>& DataLayer::InvDz() const {
    return inv_dz_;
}

void DataLayer::Validate() const {
    if (nx_ <= 0) throw std::invalid_argument("nx must be > 0");
    if (dim_ < 1 || dim_ > 3) throw std::invalid_argument("dim must be 1..3");
    if (ng_ < 0) throw std::invalid_argument("padding must be >= 0");

    // ny/nz are allowed to be anything on input; we will force them to 1 for lower dims,
    // but for dim>=2/3 they must be valid:
    if (dim_ >= 2 && ny_ <= 0) throw std::invalid_argument("ny must be > 0 for dim>=2");
    if (dim_ >= 3 && nz_ <= 0) throw std::invalid_argument("nz must be > 0 for dim>=3");
}

void DataLayer::RecomputeSizes() {
    sx_ = nx_ + 2 * ng_;
    sy_ = dim_ >= 2 ? ny_ + 2 * ng_ : 1;
    sz_ = dim_ >= 3 ? nz_ + 2 * ng_ : 1;
}

void DataLayer::Allocate() {
    U_ = xt::zeros<double>({
        k_nvar,
        static_cast<std::size_t>(sx_),
        static_cast<std::size_t>(sy_),
        static_cast<std::size_t>(sz_)
    });

    // boundaries are size (s + 1), centers and metrics are size (s)
    xb_ = xt::zeros<double>({static_cast<std::size_t>(sx_ + 1)});
    xc_ = xt::zeros<double>({static_cast<std::size_t>(sx_)});
    dx_ = xt::zeros<double>({static_cast<std::size_t>(sx_)});
    inv_dx_ = xt::zeros<double>({static_cast<std::size_t>(sx_)});

    yb_ = xt::zeros<double>({static_cast<std::size_t>(sy_ + 1)});
    yc_ = xt::zeros<double>({static_cast<std::size_t>(sy_)});
    dy_ = xt::zeros<double>({static_cast<std::size_t>(sy_)});
    inv_dy_ = xt::zeros<double>({static_cast<std::size_t>(sy_)});

    zb_ = xt::zeros<double>({static_cast<std::size_t>(sz_ + 1)});
    zc_ = xt::zeros<double>({static_cast<std::size_t>(sz_)});
    dz_ = xt::zeros<double>({static_cast<std::size_t>(sz_)});
    inv_dz_ = xt::zeros<double>({static_cast<std::size_t>(sz_)});
}

void DataLayer::UpdateMetricsFromCoordinates() {
    // X
    for (int i = 0; i < sx_; ++i) {
        const double xl = xb_(static_cast<std::size_t>(i));
        const double xr = xb_(static_cast<std::size_t>(i + 1));
        const double d = xr - xl;

        xc_(static_cast<std::size_t>(i)) = 0.5 * (xl + xr);
        dx_(static_cast<std::size_t>(i)) = d;
        inv_dx_(static_cast<std::size_t>(i)) = (d != 0.0) ? (1.0 / d) : 0.0;
    }

    // Y
    for (int j = 0; j < sy_; ++j) {
        const double yl = yb_(static_cast<std::size_t>(j));
        const double yr = yb_(static_cast<std::size_t>(j + 1));
        const double d = yr - yl;

        yc_(static_cast<std::size_t>(j)) = 0.5 * (yl + yr);
        dy_(static_cast<std::size_t>(j)) = d;
        inv_dy_(static_cast<std::size_t>(j)) = (d != 0.0) ? (1.0 / d) : 0.0;
    }

    // Z
    for (int k = 0; k < sz_; ++k) {
        const double zl = zb_(static_cast<std::size_t>(k));
        const double zr = zb_(static_cast<std::size_t>(k + 1));
        const double d = zr - zl;

        zc_(static_cast<std::size_t>(k)) = 0.5 * (zl + zr);
        dz_(static_cast<std::size_t>(k)) = d;
        inv_dz_(static_cast<std::size_t>(k)) = (d != 0.0) ? (1.0 / d) : 0.0;
    }
}

bool DataLayer::IsGlobalBoundary(const Axis axis, const Side side) const {
    const auto a = static_cast<std::size_t>(axis);
    const auto s = static_cast<std::size_t>(side);

    // Axis absent for lower dims => treat as non-global (no physical BC there)
    if (a == 1 && dim_ < 2) return false;
    if (a == 2 && dim_ < 3) return false;

    return is_global_boundary_[a][s];
}

void DataLayer::SetGlobalBoundary(const Axis axis, const Side side, const bool is_global) {
    const auto a = static_cast<std::size_t>(axis);
    const auto s = static_cast<std::size_t>(side);

    if (a == 1 && dim_ < 2) {
        throw std::invalid_argument("SetGlobalBoundary: Axis Y is not active for dim<2");
    }
    if (a == 2 && dim_ < 3) {
        throw std::invalid_argument("SetGlobalBoundary: Axis Z is not active for dim<3");
    }

    is_global_boundary_[a][s] = is_global;
}

void DataLayer::SetAllGlobalBoundaries(const bool is_global) {
    for (int a = 0; a < 3; ++a) {
        is_global_boundary_[a][0] = is_global;
        is_global_boundary_[a][1] = is_global;
    }
}
