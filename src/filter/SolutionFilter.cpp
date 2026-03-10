#include "filter/SolutionFilter.hpp"

#include <algorithm>
#include <cmath>

#include "data/DataLayer.hpp"
#include "data/Mesh.hpp"

SolutionFilter::SolutionFilter(const Settings& settings) {
    (void)settings;

    eps_diff_ = 0.05;
    eps_anti_ = 0.03;
    rho_floor_ = 1e-10;
    p_floor_ = 1e-10;
    enable_antidiffusion_ = true;
}

void SolutionFilter::ResizeFrom(const Mesh& mesh) const {
    const int sx = mesh.GetSx();
    const int sy = mesh.GetSy();
    const int sz = mesh.GetSz();

    if (sx == sx_ && sy == sy_ && sz == sz_ && rho0_.size() > 0) {
        return;
    }

    sx_ = sx;
    sy_ = sy;
    sz_ = sz;

    rho0_ = xt::zeros<double>({
        static_cast<std::size_t>(sx_),
        static_cast<std::size_t>(sy_),
        static_cast<std::size_t>(sz_)
    });
    p0_ = xt::zeros<double>(rho0_.shape());
    rho1_ = xt::zeros<double>(rho0_.shape());
    p1_ = xt::zeros<double>(rho0_.shape());
}

void SolutionFilter::ExtractRhoP(const DataLayer& layer, const Mesh& mesh, const double gamma) const {
    const auto& U = layer.U();

    for (int k = 0; k < sz_; ++k) {
        for (int j = 0; j < sy_; ++j) {
            for (int i = 0; i < sx_; ++i) {
                if (!mesh.IsFluidCell(i, j, k)) {
                    rho0_(i, j, k) = 0.0;
                    p0_(i, j, k) = 0.0;
                    rho1_(i, j, k) = 0.0;
                    p1_(i, j, k) = 0.0;
                    continue;
                }

                const double rho_in = U(DataLayer::k_rho, i, j, k);
                const double rho = rho_in > rho_floor_ ? rho_in : rho_floor_;
                const double inv_rho = 1.0 / rho;

                const double u = U(DataLayer::k_rhoU, i, j, k) * inv_rho;
                const double v = U(DataLayer::k_rhoV, i, j, k) * inv_rho;
                const double w = U(DataLayer::k_rhoW, i, j, k) * inv_rho;

                const double kinetic = 0.5 * rho * (u * u + v * v + w * w);
                const double eint = U(DataLayer::k_E, i, j, k) - kinetic;
                const double P_raw = (gamma - 1.0) * eint;
                const double P = P_raw > p_floor_ ? P_raw : p_floor_;

                rho0_(i, j, k) = rho;
                p0_(i, j, k) = P;

                rho1_(i, j, k) = rho;
                p1_(i, j, k) = P;
            }
        }
    }
}

void SolutionFilter::ApplyAxis(const Mesh& mesh, const Axis axis) const {
    const AxisStride st = AxisStride::FromAxis(axis);

    for (int k = 0; k < sz_; ++k) {
        for (int j = 0; j < sy_; ++j) {
            for (int i = 0; i < sx_; ++i) {
                if (!mesh.IsFluidCell(i, j, k)) {
                    rho1_(i, j, k) = rho0_(i, j, k);
                    p1_(i, j, k) = p0_(i, j, k);
                    continue;
                }

                const int ip = i + st.di;
                const int jp = j + st.dj;
                const int kp = k + st.dk;

                const int im = i - st.di;
                const int jm = j - st.dj;
                const int km = k - st.dk;

                const bool ok =
                    (ip >= 0 && ip < sx_ && jp >= 0 && jp < sy_ && kp >= 0 && kp < sz_) &&
                    (im >= 0 && im < sx_ && jm >= 0 && jm < sy_ && km >= 0 && km < sz_) &&
                    mesh.IsFluidCell(ip, jp, kp) &&
                    mesh.IsFluidCell(im, jm, km);

                if (!ok) {
                    rho1_(i, j, k) = rho0_(i, j, k);
                    p1_(i, j, k) = p0_(i, j, k);
                    continue;
                }

                rho1_(i, j, k) =
                    rho0_(i, j, k) +
                    eps_diff_ * (rho0_(ip, jp, kp) - 2.0 * rho0_(i, j, k) + rho0_(im, jm, km));

                p1_(i, j, k) =
                    p0_(i, j, k) +
                    eps_diff_ * (p0_(ip, jp, kp) - 2.0 * p0_(i, j, k) + p0_(im, jm, km));
            }
        }
    }

    if (!enable_antidiffusion_ || !(eps_anti_ > 0.0)) {
        rho0_ = rho1_;
        p0_ = p1_;
        return;
    }

    for (int k = 0; k < sz_; ++k) {
        for (int j = 0; j < sy_; ++j) {
            for (int i = 0; i < sx_; ++i) {
                if (!mesh.IsFluidCell(i, j, k)) {
                    rho0_(i, j, k) = rho1_(i, j, k);
                    p0_(i, j, k) = p1_(i, j, k);
                    continue;
                }

                const int ip = i + st.di;
                const int jp = j + st.dj;
                const int kp = k + st.dk;

                const int im = i - st.di;
                const int jm = j - st.dj;
                const int km = k - st.dk;

                const bool ok =
                    (ip >= 0 && ip < sx_ && jp >= 0 && jp < sy_ && kp >= 0 && kp < sz_) &&
                    (im >= 0 && im < sx_ && jm >= 0 && jm < sy_ && km >= 0 && km < sz_) &&
                    mesh.IsFluidCell(ip, jp, kp) &&
                    mesh.IsFluidCell(im, jm, km);

                if (!ok) {
                    rho0_(i, j, k) = rho1_(i, j, k);
                    p0_(i, j, k) = p1_(i, j, k);
                    continue;
                }

                rho0_(i, j, k) =
                    rho1_(i, j, k) -
                    eps_anti_ * (rho1_(ip, jp, kp) - 2.0 * rho1_(i, j, k) + rho1_(im, jm, km));

                p0_(i, j, k) =
                    p1_(i, j, k) -
                    eps_anti_ * (p1_(ip, jp, kp) - 2.0 * p1_(i, j, k) + p1_(im, jm, km));
            }
        }
    }
}

void SolutionFilter::WriteBackConservative(DataLayer& layer, const Mesh& mesh, const double gamma) const {
    auto& U = layer.U();

    const int i0 = mesh.GetCoreStartX();
    const int i1 = mesh.GetCoreEndExclusiveX();
    const int j0 = mesh.GetCoreStartY();
    const int j1 = mesh.GetCoreEndExclusiveY();
    const int k0 = mesh.GetCoreStartZ();
    const int k1 = mesh.GetCoreEndExclusiveZ();

    for (int k = k0; k < k1; ++k) {
        for (int j = j0; j < j1; ++j) {
            for (int i = i0; i < i1; ++i) {
                if (!mesh.IsFluidCell(i, j, k)) {
                    continue;
                }

                const double rho_old = U(DataLayer::k_rho, i, j, k);
                const double rho = std::max(rho0_(i, j, k), rho_floor_);
                const double P = std::max(p0_(i, j, k), p_floor_);

                const double rho_safe = rho_old > rho_floor_ ? rho_old : rho_floor_;
                const double inv_rho_old = 1.0 / rho_safe;

                const double u = U(DataLayer::k_rhoU, i, j, k) * inv_rho_old;
                const double v = U(DataLayer::k_rhoV, i, j, k) * inv_rho_old;
                const double w = U(DataLayer::k_rhoW, i, j, k) * inv_rho_old;

                U(DataLayer::k_rho, i, j, k) = rho;
                U(DataLayer::k_rhoU, i, j, k) = rho * u;
                U(DataLayer::k_rhoV, i, j, k) = rho * v;
                U(DataLayer::k_rhoW, i, j, k) = rho * w;

                const double kinetic = 0.5 * rho * (u * u + v * v + w * w);
                const double E = P / (gamma - 1.0) + kinetic;
                U(DataLayer::k_E, i, j, k) = E;
            }
        }
    }
}

void SolutionFilter::Apply(DataLayer& layer, const Mesh& mesh, const double gamma) const {
    const int ng = mesh.GetPadding();
    if (ng < 1) {
        return;
    }

    ResizeFrom(mesh);
    ExtractRhoP(layer, mesh, gamma);

    ApplyAxis(mesh, Axis::X);
    if (mesh.GetDim() >= 2) {
        ApplyAxis(mesh, Axis::Y);
    }
    if (mesh.GetDim() >= 3) {
        ApplyAxis(mesh, Axis::Z);
    }

    for (int k = 0; k < sz_; ++k) {
        for (int j = 0; j < sy_; ++j) {
            for (int i = 0; i < sx_; ++i) {
                if (!mesh.IsFluidCell(i, j, k)) {
                    continue;
                }
                rho0_(i, j, k) = std::max(rho0_(i, j, k), rho_floor_);
                p0_(i, j, k) = std::max(p0_(i, j, k), p_floor_);
            }
        }
    }

    WriteBackConservative(layer, mesh, gamma);
}