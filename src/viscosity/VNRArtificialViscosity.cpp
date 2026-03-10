#include "viscosity/VNRArtificialViscosity.hpp"

#include "data/DataLayer.hpp"
#include "data/Mesh.hpp"

VNRArtificialViscosity::VNRArtificialViscosity(const Settings& settings,
                                               const double C1,
                                               const double C2)
    : C1_(C1), C2_(C2), settings_(settings) {}

PrimitiveCell VNRArtificialViscosity::LoadPrimitive(const xt::xtensor<double, 4>& W,
                                                    const int i,
                                                    const int j,
                                                    const int k) {
    PrimitiveCell w;
    w.rho = W(var::u_rho, i, j, k);
    w.u = W(var::u_u, i, j, k);
    w.v = W(var::u_v, i, j, k);
    w.w = W(var::u_w, i, j, k);
    w.P = W(var::u_P, i, j, k);
    return w;
}

void VNRArtificialViscosity::ResizeFrom(const Mesh& mesh) const {
    const int sx = mesh.GetSx();
    const int sy = mesh.GetSy();
    const int sz = mesh.GetSz();

    const bool ok =
        (sx == sx_) && (sy == sy_) && (sz == sz_) &&
        (qx_.size() > 0 || sx_ <= 1) &&
        (qy_.size() > 0 || sy_ <= 1) &&
        (qz_.size() > 0 || sz_ <= 1);

    if (ok) {
        return;
    }

    sx_ = sx;
    sy_ = sy;
    sz_ = sz;

    const std::size_t Sx = static_cast<std::size_t>(sx_);
    const std::size_t Sy = static_cast<std::size_t>(sy_);
    const std::size_t Sz = static_cast<std::size_t>(sz_);

    qx_ = xt::zeros<double>({(Sx > 0 ? Sx - 1 : 0), Sy, Sz});
    qy_ = xt::zeros<double>({Sx, (Sy > 0 ? Sy - 1 : 0), Sz});
    qz_ = xt::zeros<double>({Sx, Sy, (Sz > 0 ? Sz - 1 : 0)});
}

xt::xtensor<double, 3>& VNRArtificialViscosity::QFace(const Axis axis) const {
    if (axis == Axis::X) {
        return qx_;
    }
    if (axis == Axis::Y) {
        return qy_;
    }
    return qz_;
}

const xt::xtensor<double, 1>& VNRArtificialViscosity::InvMetric(const Mesh& mesh, const Axis axis) const {
    if (axis == Axis::X) {
        return mesh.InvDx();
    }
    if (axis == Axis::Y) {
        return mesh.InvDy();
    }
    return mesh.InvDz();
}

std::size_t VNRArtificialViscosity::MomVarIndex(const Axis axis) {
    if (axis == Axis::X) {
        return DataLayer::k_rhoU;
    }
    if (axis == Axis::Y) {
        return DataLayer::k_rhoV;
    }
    return DataLayer::k_rhoW;
}

std::size_t VNRArtificialViscosity::VelVarIndex(const Axis axis) {
    if (axis == Axis::X) {
        return var::u_u;
    }
    if (axis == Axis::Y) {
        return var::u_v;
    }
    return var::u_w;
}

void VNRArtificialViscosity::ComputeQFaces(const Mesh& mesh,
                                           const xt::xtensor<double, 4>& W,
                                           const double gamma,
                                           const Axis axis) const {
    auto& qface = QFace(axis);
    qface.fill(0.0);

    const AxisStride st = AxisStride::FromAxis(axis);

    const int sx = mesh.GetSx();
    const int sy = mesh.GetSy();
    const int sz = mesh.GetSz();

    const int fx = axis == Axis::X ? (sx - 1) : sx;
    const int fy = axis == Axis::Y ? (sy - 1) : sy;
    const int fz = axis == Axis::Z ? (sz - 1) : sz;

    for (int k = 0; k < fz; ++k) {
        for (int j = 0; j < fy; ++j) {
            for (int i = 0; i < fx; ++i) {
                const int ir = i + st.di;
                const int jr = j + st.dj;
                const int kr = k + st.dk;

                const PrimitiveCell wl = LoadPrimitive(W, i, j, k);
                const PrimitiveCell wr = LoadPrimitive(W, ir, jr, kr);

                const double un_l = axis == Axis::X ? wl.u : (axis == Axis::Y ? wl.v : wl.w);
                const double un_r = axis == Axis::X ? wr.u : (axis == Axis::Y ? wr.v : wr.w);
                const double du = un_r - un_l;

                if (du >= 0.0) {
                    continue;
                }

                const double rho_bar = 0.5 * (wl.rho + wr.rho);
                if (!(rho_bar > 0.0) || !std::isfinite(rho_bar)) {
                    continue;
                }

                const double c_l = SoundSpeed(wl, gamma);
                const double c_r = SoundSpeed(wr, gamma);
                const double c_bar = 0.5 * (c_l + c_r);

                const double du_abs = std::fabs(du);
                const double q = rho_bar * (C2_ * du * du + C1_ * c_bar * du_abs);

                if (q > 0.0 && std::isfinite(q)) {
                    qface(static_cast<std::size_t>(i),
                          static_cast<std::size_t>(j),
                          static_cast<std::size_t>(k)) = q;
                }
            }
        }
    }
}

void VNRArtificialViscosity::AddAxisContribution(const Mesh& mesh,
                                                 const xt::xtensor<double, 4>& W,
                                                 xt::xtensor<double, 4>& rhs,
                                                 const Axis axis) const {
    const AxisStride st = AxisStride::FromAxis(axis);
    const auto& inv_h = InvMetric(mesh, axis);
    const auto& qface = QFace(axis);

    const std::size_t mom = MomVarIndex(axis);
    const std::size_t vvel = VelVarIndex(axis);

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

                const double q_plus = qface(static_cast<std::size_t>(i),
                                            static_cast<std::size_t>(j),
                                            static_cast<std::size_t>(k));

                const int im = i - st.di;
                const int jm = j - st.dj;
                const int km = k - st.dk;

                const double q_minus = qface(static_cast<std::size_t>(im),
                                             static_cast<std::size_t>(jm),
                                             static_cast<std::size_t>(km));

                const double v0 = W(vvel, i, j, k);
                const double vp = W(vvel, i + st.di, j + st.dj, k + st.dk);
                const double vm = W(vvel, i - st.di, j - st.dj, k - st.dk);

                const double vel_plus = 0.5 * (v0 + vp);
                const double vel_minus = 0.5 * (vm + v0);

                const int idx = axis == Axis::X ? i : (axis == Axis::Y ? j : k);
                const double inv = inv_h(static_cast<std::size_t>(idx));

                rhs(mom, i, j, k) += -(q_plus - q_minus) * inv;
                rhs(DataLayer::k_E, i, j, k) += -((vel_plus * q_plus) - (vel_minus * q_minus)) * inv;
            }
        }
    }
}

void VNRArtificialViscosity::AddToRhs(const DataLayer& layer,
                                      const Mesh& mesh,
                                      const xt::xtensor<double, 4>& W,
                                      const double gamma,
                                      const double dt,
                                      xt::xtensor<double, 4>& rhs) const {
    (void)layer;
    (void)dt;

    ResizeFrom(mesh);

    ComputeQFaces(mesh, W, gamma, Axis::X);
    AddAxisContribution(mesh, W, rhs, Axis::X);

    if (mesh.GetDim() >= 2) {
        ComputeQFaces(mesh, W, gamma, Axis::Y);
        AddAxisContribution(mesh, W, rhs, Axis::Y);
    }

    if (mesh.GetDim() >= 3) {
        ComputeQFaces(mesh, W, gamma, Axis::Z);
        AddAxisContribution(mesh, W, rhs, Axis::Z);
    }
}