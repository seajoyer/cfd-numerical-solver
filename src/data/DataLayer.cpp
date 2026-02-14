#include "data/DataLayer.hpp"

// ---------- 1D constructor (backward compatible) ----------

DataLayer::DataLayer(const int N, const int padding)
    : nx_(N), ny_(0), n_ghost_cells_(padding), dimension_(1) {
    if (N <= 0) throw std::invalid_argument("Number of cells N must be positive");
    if (padding < 0) throw std::invalid_argument("Padding must be non-negative");
    RecomputeSizes();
    Allocate1D();
}

// ---------- N x N constructor with dim ----------

DataLayer::DataLayer(const int N, const int padding, const int dim)
    : nx_(N), n_ghost_cells_(padding), dimension_(dim) {
    if (N <= 0) throw std::invalid_argument("Number of cells N must be positive");
    if (padding < 0) throw std::invalid_argument("Padding must be non-negative");
    if (dim < 1 || dim > 3) throw std::invalid_argument("Dimension must be 1, 2, or 3");

    ny_ = (dim >= 2) ? N : 0;

    RecomputeSizes();
    if (dim == 1) {
        Allocate1D();
    } else {
        Allocate2D();
    }
}

// ---------- Nx x Ny constructor for 2D ----------

DataLayer::DataLayer(const int Nx, const int Ny, const int padding, const int dim)
    : nx_(Nx), ny_(Ny), n_ghost_cells_(padding), dimension_(dim) {
    if (Nx <= 0) throw std::invalid_argument("Nx must be positive");
    if (dim >= 2 && Ny <= 0) throw std::invalid_argument("Ny must be positive for 2D");
    if (padding < 0) throw std::invalid_argument("Padding must be non-negative");
    if (dim < 1 || dim > 3) throw std::invalid_argument("Dimension must be 1, 2, or 3");

    if (dim == 1) ny_ = 0;

    RecomputeSizes();
    if (dim == 1) {
        Allocate1D();
    } else {
        Allocate2D();
    }
}

// ---------- Recompute sizes ----------

void DataLayer::RecomputeSizes() {
    total_size_x_ = nx_ + 2 * n_ghost_cells_;
    if (total_size_x_ < 0) total_size_x_ = 0;

    if (dimension_ >= 2) {
        total_size_y_ = ny_ + 2 * n_ghost_cells_;
        if (total_size_y_ < 0) total_size_y_ = 0;
    } else {
        total_size_y_ = 0;
    }
}

// ---------- 1D allocation ----------

void DataLayer::Allocate1D() {
    const auto size = static_cast<std::size_t>(total_size_x_);

    rho = xt::zeros<double>({size});
    u   = xt::zeros<double>({size});
    P   = xt::zeros<double>({size});
    p   = xt::zeros<double>({size});
    e   = xt::zeros<double>({size});
    U   = xt::zeros<double>({size});
    V   = xt::zeros<double>({size});
    m   = xt::zeros<double>({size});
    xb  = xt::zeros<double>({size});
    xc  = xt::zeros<double>({size});

    // v, q, yb, yc not allocated for 1D (remain empty)
    v  = xt::xarray<double>();
    q  = xt::xarray<double>();
    yb = xt::xarray<double>();
    yc = xt::xarray<double>();
}

// ---------- 2D allocation ----------

void DataLayer::Allocate2D() {
    const auto sx = static_cast<std::size_t>(total_size_x_);
    const auto sy = static_cast<std::size_t>(total_size_y_);

    rho = xt::zeros<double>({sx, sy});
    u   = xt::zeros<double>({sx, sy});
    v   = xt::zeros<double>({sx, sy});
    P   = xt::zeros<double>({sx, sy});
    p   = xt::zeros<double>({sx, sy});  // rho*u
    q   = xt::zeros<double>({sx, sy});  // rho*v
    e   = xt::zeros<double>({sx, sy});
    U   = xt::zeros<double>({sx, sy});
    V   = xt::zeros<double>({sx, sy});
    m   = xt::zeros<double>({sx, sy});

    // 1D coordinate arrays along each axis
    xb = xt::zeros<double>({sx});
    xc = xt::zeros<double>({sx});
    yb = xt::zeros<double>({sy});
    yc = xt::zeros<double>({sy});
}

// ---------- Setters ----------

void DataLayer::SetN(const int new_N) {
    if (new_N <= 0) throw std::invalid_argument("Number of cells N must be positive");
    nx_ = new_N;
    if (dimension_ >= 2) ny_ = new_N;
    RecomputeSizes();
    if (dimension_ == 1) Allocate1D(); else Allocate2D();
}

void DataLayer::SetPadding(const int new_padding) {
    if (new_padding < 0) throw std::invalid_argument("Padding must be non-negative");
    n_ghost_cells_ = new_padding;
    RecomputeSizes();
    if (dimension_ == 1) Allocate1D(); else Allocate2D();
}

void DataLayer::SetDim(const int new_dim) {
    if (new_dim < 1 || new_dim > 3) throw std::invalid_argument("Dimension must be 1, 2, or 3");
    dimension_ = new_dim;
    if (new_dim >= 2 && ny_ == 0) ny_ = nx_;
    if (new_dim == 1) ny_ = 0;
    RecomputeSizes();
    if (new_dim == 1) Allocate1D(); else Allocate2D();
}

// ---------- Core start/end ----------

auto DataLayer::GetCoreStart(const int axis) const -> int {
    (void)axis;
    return n_ghost_cells_;
}

auto DataLayer::GetCoreEndExclusive(const int axis) const -> int {
    if (axis == 0) return n_ghost_cells_ + nx_;
    if (axis == 1) return n_ghost_cells_ + ny_;
    return n_ghost_cells_ + nx_;
}

auto DataLayer::GetTotalSize(int axis) const -> int {
    if (axis == 0) return total_size_x_;
    if (axis == 1) return total_size_y_;
    return total_size_x_;
}

// ---------- 1D state access ----------

auto DataLayer::GetPrimitive(const int i) const -> Primitive {
    Primitive w;
    w.rho = rho(i);
    w.u = u(i);
    w.v = 0.0;
    w.P = P(i);
    return w;
}

void DataLayer::SetPrimitive(const int i, const Primitive& w) {
    rho(i) = w.rho;
    u(i) = w.u;
    P(i) = w.P;
}

// ---------- 2D state access ----------

auto DataLayer::GetPrimitive2D(int i, int j) const -> Primitive {
    return {rho(i, j), u(i, j), v(i, j), P(i, j)};
}

void DataLayer::SetPrimitive2D(int i, int j, const Primitive& w) {
    rho(i, j) = w.rho;
    u(i, j) = w.u;
    v(i, j) = w.v;
    P(i, j) = w.P;
}

auto DataLayer::GetConservative2D(int i, int j, double gamma) const -> Conservative {
    const double rho_val = rho(i, j);
    const double u_val = u(i, j);
    const double v_val = v(i, j);
    const double P_val = P(i, j);

    const double rhoU = rho_val * u_val;
    const double rhoV = rho_val * v_val;
    const double E = P_val / (gamma - 1.0) + 0.5 * rho_val * (u_val * u_val + v_val * v_val);

    return {rho_val, rhoU, rhoV, E};
}

void DataLayer::SetConservative2D(int i, int j, const Conservative& uc,
                                  double gamma, double dx, double dy) {
    const double rho_val = uc.rho;
    const double u_val = (rho_val > 0.0) ? uc.rhoU / rho_val : 0.0;
    const double v_val = (rho_val > 0.0) ? uc.rhoV / rho_val : 0.0;
    const double kinetic = 0.5 * rho_val * (u_val * u_val + v_val * v_val);
    const double P_val = (gamma - 1.0) * (uc.E - kinetic);

    rho(i, j) = rho_val;
    u(i, j) = u_val;
    v(i, j) = v_val;
    P(i, j) = P_val;

    p(i, j) = uc.rhoU;
    q(i, j) = uc.rhoV;
    V(i, j) = (rho_val > 0.0) ? 1.0 / rho_val : 0.0;

    const double Eint = uc.E - kinetic;
    const double eint = (rho_val > 0.0) ? Eint / rho_val : 0.0;

    U(i, j) = eint;
    e(i, j) = (rho_val > 0.0) ? uc.E : 0.0;
    m(i, j) = rho_val * dx * dy;
}
