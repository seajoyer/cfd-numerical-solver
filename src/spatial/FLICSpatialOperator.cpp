#include "spatial/FLICSpatialOperator.hpp"

#include <cmath>
#include <stdexcept>

#include "bc/BoundaryManager.hpp"
#include "data/Variables.hpp"
#include "viscosity/ArtificialViscosity.hpp"
#include "viscosity/VNRArtificialViscosity.hpp"

FLICSpatialOperator::FLICSpatialOperator(const Settings& settings,
                                         std::shared_ptr<BoundaryManager> boundary_manager)
    : SpatialOperator(std::move(boundary_manager)) {
    viscosity_ = nullptr;
    if (settings.viscosity) {
        viscosity_ = std::make_shared<VNRArtificialViscosity>(settings);
    }
}

const xt::xtensor<double, 1>& FLICSpatialOperator::InvMetric(const Mesh& mesh, const Axis axis) const {
    if (axis == Axis::X) {
        return mesh.InvDx();
    }
    if (axis == Axis::Y) {
        return mesh.InvDy();
    }
    return mesh.InvDz();
}

double FLICSpatialOperator::CellVelocityComponent(const xt::xtensor<double, 4>& U_star,
                                                  const int i,
                                                  const int j,
                                                  const int k,
                                                  const Axis axis,
                                                  const double rho_floor) {
    const double rho_in = U_star(DataLayer::k_rho, i, j, k);
    const double rho = rho_in > rho_floor ? rho_in : rho_floor;

    if (axis == Axis::X) {
        return U_star(DataLayer::k_rhoU, i, j, k) / rho;
    }
    if (axis == Axis::Y) {
        return U_star(DataLayer::k_rhoV, i, j, k) / rho;
    }
    return U_star(DataLayer::k_rhoW, i, j, k) / rho;
}

void FLICSpatialOperator::ComputeLagrangianRhs(const Mesh& mesh,
                                               const xt::xtensor<double, 4>& W,
                                               xt::xtensor<double, 4>& rhs) const {
    const auto& inv_dx = mesh.InvDx();
    const auto& inv_dy = mesh.InvDy();
    const auto& inv_dz = mesh.InvDz();

    rhs.fill(0.0);

    const int dim = mesh.GetDim();
    const int i0 = 1;
    const int i1 = mesh.GetSx() - 1;
    const int j0 = dim >= 2 ? 1 : 0;
    const int j1 = dim >= 2 ? mesh.GetSy() - 1 : 1;
    const int k0 = dim >= 3 ? 1 : 0;
    const int k1 = dim >= 3 ? mesh.GetSz() - 1 : 1;

    for (int i = i0; i < i1; ++i) {
        const double idx = inv_dx(std::size_t(i));

        for (int j = j0; j < j1; ++j) {
            const double idy = dim >= 2 ? inv_dy(std::size_t(j)) : 0.0;

            for (int k = k0; k < k1; ++k) {
                if (!mesh.IsFluidCell(i, j, k)) {
                    continue;
                }

                const double idz = dim >= 3 ? inv_dz(std::size_t(k)) : 0.0;

                const double dPdx = 0.5 * (W(var::u_P, i + 1, j, k) - W(var::u_P, i - 1, j, k)) * idx;
                rhs(DataLayer::k_rhoU, i, j, k) -= dPdx;

                double div_Pv = 0.0;

                const double Pu_ip1 = W(var::u_P, i + 1, j, k) * W(var::u_u, i + 1, j, k);
                const double Pu_im1 = W(var::u_P, i - 1, j, k) * W(var::u_u, i - 1, j, k);
                div_Pv += 0.5 * (Pu_ip1 - Pu_im1) * idx;

                if (dim >= 2) {
                    const double dPdy = 0.5 * (W(var::u_P, i, j + 1, k) - W(var::u_P, i, j - 1, k)) * idy;
                    rhs(DataLayer::k_rhoV, i, j, k) -= dPdy;

                    const double Pv_jp1 = W(var::u_P, i, j + 1, k) * W(var::u_v, i, j + 1, k);
                    const double Pv_jm1 = W(var::u_P, i, j - 1, k) * W(var::u_v, i, j - 1, k);
                    div_Pv += 0.5 * (Pv_jp1 - Pv_jm1) * idy;
                }

                if (dim >= 3) {
                    const double dPdz = 0.5 * (W(var::u_P, i, j, k + 1) - W(var::u_P, i, j, k - 1)) * idz;
                    rhs(DataLayer::k_rhoW, i, j, k) -= dPdz;

                    const double Pw_kp1 = W(var::u_P, i, j, k + 1) * W(var::u_w, i, j, k + 1);
                    const double Pw_km1 = W(var::u_P, i, j, k - 1) * W(var::u_w, i, j, k - 1);
                    div_Pv += 0.5 * (Pw_kp1 - Pw_km1) * idz;
                }

                rhs(DataLayer::k_E, i, j, k) -= div_Pv;
            }
        }
    }
}

void FLICSpatialOperator::BuildUStar(const DataLayer& layer,
                                     const Mesh& mesh,
                                     const xt::xtensor<double, 4>& U,
                                     const xt::xtensor<double, 4>& rhs_lag,
                                     const double dt,
                                     xt::xtensor<double, 4>& U_star) const {
    const int sx = mesh.GetSx();
    const int sy = mesh.GetSy();
    const int sz = mesh.GetSz();

    for (std::size_t v = 0; v < DataLayer::k_nvar; ++v) {
        for (int i = 0; i < sx; ++i) {
            for (int j = 0; j < sy; ++j) {
                for (int k = 0; k < sz; ++k) {
                    U_star(v, i, j, k) = U(v, i, j, k) + dt * rhs_lag(v, i, j, k);
                }
            }
        }
    }
}

void FLICSpatialOperator::AccumulateAxisAdvection(const Mesh& mesh,
                                                  const xt::xtensor<double, 4>& U_star,
                                                  xt::xtensor<double, 4>& rhs,
                                                  const Axis axis) const {
    const AxisStride st = AxisStride::FromAxis(axis);
    const auto& inv_h = InvMetric(mesh, axis);

    const int i0 = mesh.GetCoreStartX();
    const int i1 = mesh.GetCoreEndExclusiveX();
    const int j0 = mesh.GetCoreStartY();
    const int j1 = mesh.GetCoreEndExclusiveY();
    const int k0 = mesh.GetCoreStartZ();
    const int k1 = mesh.GetCoreEndExclusiveZ();

    const int i_start = i0 - st.di;
    const int j_start = j0 - st.dj;
    const int k_start = k0 - st.dk;

    for (int i = i_start; i < i1; ++i) {
        for (int j = j_start; j < j1; ++j) {
            for (int k = k_start; k < k1; ++k) {
                const int ip = i + st.di;
                const int jp = j + st.dj;
                const int kp = k + st.dk;

                const double uL = CellVelocityComponent(U_star, i, j, k, axis);
                const double uR = CellVelocityComponent(U_star, ip, jp, kp, axis);
                const double u_face = 0.5 * (uL + uR);

                const int up_i = u_face >= 0.0 ? i : ip;
                const int up_j = u_face >= 0.0 ? j : jp;
                const int up_k = u_face >= 0.0 ? k : kp;

                for (std::size_t v = 0; v < DataLayer::k_nvar; ++v) {
                    const double flux = u_face * U_star(v, up_i, up_j, up_k);

                    if (i >= i0 && j >= j0 && k >= k0 && mesh.IsFluidCell(i, j, k)) {
                        const double invL = inv_h((axis == Axis::X) ? i : (axis == Axis::Y) ? j : k);
                        rhs(v, i, j, k) -= flux * invL;
                    }

                    if (ip < i1 && jp < j1 && kp < k1 && mesh.IsFluidCell(ip, jp, kp)) {
                        const double invR = inv_h((axis == Axis::X) ? ip : (axis == Axis::Y) ? jp : kp);
                        rhs(v, ip, jp, kp) += flux * invR;
                    }
                }
            }
        }
    }
}

void FLICSpatialOperator::ComputeEulerianAdvectionRhs(const Mesh& mesh,
                                                      const xt::xtensor<double, 4>& U_star,
                                                      xt::xtensor<double, 4>& rhs) const {
    AccumulateAxisAdvection(mesh, U_star, rhs, Axis::X);
    if (mesh.GetDim() >= 2) {
        AccumulateAxisAdvection(mesh, U_star, rhs, Axis::Y);
    }
    if (mesh.GetDim() >= 3) {
        AccumulateAxisAdvection(mesh, U_star, rhs, Axis::Z);
    }
}

void FLICSpatialOperator::ComputeRHS(DataLayer& layer,
                                     const Mesh& mesh,
                                     Workspace& workspace,
                                     const double gamma,
                                     const double dt) const {
    if (!boundary_manager_) {
        throw std::runtime_error("FLICSpatialOperator: boundary_manager_ is null");
    }
    if (dt <= 0.0) {
        throw std::runtime_error("FLICSpatialOperator: dt must be > 0");
    }

    const int ng = mesh.GetPadding();
    if (ng < 1) {
        throw std::runtime_error("FLICSpatialOperator: requires at least 1 ghost cell (ng>=1)");
    }

    workspace.ResizeFrom(mesh);

    boundary_manager_->UpdateHalo(layer, mesh);
    boundary_manager_->ApplyPhysicalBc(layer, mesh);

    ConvertUtoW(layer.U(), workspace.W(), gamma,
                0, mesh.GetSx(),
                0, mesh.GetSy(),
                0, mesh.GetSz());

    workspace.ZeroRhs();
    auto& rhs = workspace.Rhs();

    ComputeLagrangianRhs(mesh, workspace.W(), rhs);

    xt::xtensor<double, 4> U_star = xt::empty_like(layer.U());
    BuildUStar(layer, mesh, layer.U(), rhs, dt, U_star);

    xt::xtensor<double, 4> U_saved = xt::eval(layer.U());

    layer.U() = U_star;
    boundary_manager_->UpdateHalo(layer, mesh);
    boundary_manager_->ApplyPhysicalBc(layer, mesh);

    const auto& U_star_bc = layer.U();

    ComputeEulerianAdvectionRhs(mesh, U_star_bc, rhs);

    layer.U() = U_saved;

    if (viscosity_) {
        viscosity_->AddToRhs(layer, mesh, workspace.W(), gamma, dt, rhs);
    }
}
