#include "spatial/MaderSpatialOperator.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "bc/BoundaryManager.hpp"
#include "config/Settings.hpp"
#include "data/DataLayer.hpp"
#include "data/Mesh.hpp"
#include "data/Variables.hpp"
#include "data/Workspace.hpp"

namespace {
    constexpr double k_eps = 1e-14;

    inline double GetRhoFloor(const EOS& eos) {
        if (eos.GetType() == EosType::IdealGas) {
            return eos.GetIdealGasParameters().rho_floor;
        }
        if (eos.GetType() == EosType::HugoniotGruneisen) {
            return eos.GetHugoniotGruneisenParameters().rho_floor;
        }
        return 1e-14;
    }

    inline bool IsCoreCell(const Mesh& mesh, const int i, const int j, const int k) {
        return i >= mesh.GetCoreStartX() && i < mesh.GetCoreEndExclusiveX() &&
            j >= mesh.GetCoreStartY() && j < mesh.GetCoreEndExclusiveY() &&
            k >= mesh.GetCoreStartZ() && k < mesh.GetCoreEndExclusiveZ();
    }

    inline void AddTransportContributionIfCore(xt::xtensor<double, 4>& D,
                                               const Mesh& mesh,
                                               const int i,
                                               const int j,
                                               const int k,
                                               const double dm,
                                               const double de,
                                               const double dw,
                                               const double dpu,
                                               const double dpv) {
        if (!IsCoreCell(mesh, i, j, k)) {
            return;
        }

        D(Workspace::k_dm, i, j, k) += dm;
        D(Workspace::k_de, i, j, k) += de;
        D(Workspace::k_dw, i, j, k) += dw;
        D(Workspace::k_dpu, i, j, k) += dpu;
        D(Workspace::k_dpv, i, j, k) += dpv;
    }

    inline void RebuildFaceVelocitiesFromCellCentered(const Mesh& mesh,
                                                      Workspace& workspace) {
        auto& W = workspace.W();
        auto& Ux = workspace.Ux();
        auto& Vy = workspace.Vy();
        auto& Wz = workspace.Wz();

        const int sx = mesh.GetSx();
        const int sy = mesh.GetSy();
        const int sz = mesh.GetSz();

        for (int k = 0; k < sz; ++k) {
            for (int j = 0; j < sy; ++j) {
                Ux(0, j, k) = W(Workspace::k_u, 0, j, k);
                for (int i = 1; i < sx; ++i) {
                    Ux(i, j, k) = 0.5 * (W(Workspace::k_u, i - 1, j, k) +
                        W(Workspace::k_u, i, j, k));
                }
                Ux(sx, j, k) = W(Workspace::k_u, sx - 1, j, k);
            }
        }

        for (int k = 0; k < sz; ++k) {
            for (int i = 0; i < sx; ++i) {
                Vy(i, 0, k) = W(Workspace::k_v, i, 0, k);
                for (int j = 1; j < sy; ++j) {
                    Vy(i, j, k) = 0.5 * (W(Workspace::k_v, i, j - 1, k) +
                        W(Workspace::k_v, i, j, k));
                }
                Vy(i, sy, k) = W(Workspace::k_v, i, sy - 1, k);
            }
        }

        for (int j = 0; j < sy; ++j) {
            for (int i = 0; i < sx; ++i) {
                Wz(i, j, 0) = W(Workspace::k_w, i, j, 0);
                for (int k = 1; k < sz; ++k) {
                    Wz(i, j, k) = 0.5 * (W(Workspace::k_w, i, j, k - 1) +
                        W(Workspace::k_w, i, j, k));
                }
                Wz(i, j, sz) = W(Workspace::k_w, i, j, sz - 1);
            }
        }
    }

    inline double ComputeTransferredReactantFractionShargatovImpl(
        const Mesh& mesh,
        const DataLayer& layer,
        const int donor_i,
        const int donor_j,
        const int donor_k,
        const int accept_i,
        const int accept_j,
        const int accept_k) {
        const auto& reactant = layer.ReactantMassFraction();

        const auto is_pure_a = [](const double w) {
            return std::abs(w) <= 1e-12;
        };
        const auto is_pure_b = [](const double w) {
            return std::abs(w - 1.0) <= 1e-12;
        };
        const auto is_mixed = [](const double w) {
            return w > 1e-12 && w < 1.0 - 1e-12;
        };

        const double donor_w = std::clamp(reactant(donor_i, donor_j, donor_k), 0.0, 1.0);
        if (!is_mixed(donor_w)) {
            return donor_w;
        }

        bool has_pure_a_neighbor = false;
        bool has_pure_b_neighbor = false;

        const int sx = mesh.GetSx();
        const int sy = mesh.GetSy();
        const int sz = mesh.GetSz();

        const int ni[4] = {donor_i - 1, donor_i + 1, donor_i, donor_i};
        const int nj[4] = {donor_j, donor_j, donor_j - 1, donor_j + 1};

        for (int n = 0; n < 4; ++n) {
            const int ii = ni[n];
            const int jj = nj[n];
            const int kk = donor_k;

            if (ii < 0 || ii >= sx || jj < 0 || jj >= sy || kk < 0 || kk >= sz) {
                continue;
            }
            if (ii == accept_i && jj == accept_j && kk == accept_k) {
                continue;
            }
            if (!mesh.IsFluidCell(ii, jj, kk)) {
                continue;
            }

            const double wn = std::clamp(reactant(ii, jj, kk), 0.0, 1.0);
            if (is_pure_a(wn)) {
                has_pure_a_neighbor = true;
            }
            if (is_pure_b(wn)) {
                has_pure_b_neighbor = true;
            }
        }

        if (has_pure_a_neighbor && !has_pure_b_neighbor) {
            return 1.0;
        }
        if (has_pure_b_neighbor && !has_pure_a_neighbor) {
            return 0.0;
        }

        return donor_w;
    }
} // namespace

MaderSpatialOperator::MaderSpatialOperator(const Settings& settings, std::shared_ptr<BoundaryManager> boundary_manager)
    : SpatialOperator(std::move(boundary_manager)) {
    MaderChemistryParameters chem;
    chem.enabled = settings.chemistry_enabled;
    chem.z_freq = settings.chemistry_z_freq;
    chem.activation_energy = settings.chemistry_activation_energy;
    chem.gas_constant = settings.chemistry_gas_constant;
    chem.heat_release = settings.chemistry_heat_release;
    settings_ = settings;

    SetChemistryParameters(chem);
}

void MaderSpatialOperator::SetEos(const EOS& eos) {
    eos_ = eos;
}

const EOS& MaderSpatialOperator::GetEos() const {
    return eos_;
}

void MaderSpatialOperator::SetChemistryParameters(const MaderChemistryParameters& params) {
    chemistry_params_ = params;
}

const MaderChemistryParameters& MaderSpatialOperator::GetChemistryParameters() const {
    return chemistry_params_;
}

void MaderSpatialOperator::SetViscosityParameters(const MaderViscosityParameters& params) {
    viscosity_params_ = params;
}

const MaderViscosityParameters& MaderSpatialOperator::GetViscosityParameters() const {
    return viscosity_params_;
}

void MaderSpatialOperator::SetTransportParameters(const MaderTransportParameters& params) {
    transport_params_ = params;
}

const MaderTransportParameters& MaderSpatialOperator::GetTransportParameters() const {
    return transport_params_;
}

void MaderSpatialOperator::ComputeRHS(DataLayer&,
                                      const Mesh&,
                                      Workspace&,
                                      const double,
                                      const double) const {
    // Mader method is phase-split and is not implemented as a single RHS operator.
}

void MaderSpatialOperator::Phase1_EosAndChemistry(DataLayer& layer,
                                                  const Mesh& mesh,
                                                  Workspace& workspace,
                                                  const double gamma,
                                                  const double dt) const {
    (void)gamma;
    ApplyBoundaryConditions(layer, mesh);
    ComputeCellCenteredThermodynamics(layer, mesh, workspace, gamma, dt);
    InitializeFaceVelocitiesIfNeeded(layer, mesh, workspace);
    ++cycle_counter_;
}

void MaderSpatialOperator::Phase2_PressureForces(DataLayer& layer,
                                                 const Mesh& mesh,
                                                 Workspace& workspace,
                                                 const double gamma,
                                                 const double dt) const {
    (void)gamma;

    ApplyBoundaryConditions(layer, mesh);
    SaveOldFaceVelocities(workspace);
    ComputeArtificialViscosity(mesh, workspace);
    UpdateUxFaces(mesh, workspace, dt);
    UpdateVyFaces(mesh, workspace, dt);

    auto& W = workspace.W();
    auto& Ux = workspace.Ux();
    auto& Vy = workspace.Vy();

    const int sx = mesh.GetSx();
    const int sy = mesh.GetSy();
    const int sz = mesh.GetSz();

    for (int k = 0; k < sz; ++k) {
        for (int j = 0; j < sy; ++j) {
            for (int i = 0; i < sx; ++i) {
                W(Workspace::k_u, i, j, k) = 0.5 * (Ux(i, j, k) + Ux(i + 1, j, k));
                W(Workspace::k_v, i, j, k) = 0.5 * (Vy(i, j, k) + Vy(i, j + 1, k));
            }
        }
    }

    ApplyBoundaryConditions(layer, mesh);
}

void MaderSpatialOperator::Phase3_ZipEnergy(DataLayer& layer,
                                            const Mesh& mesh,
                                            Workspace& workspace,
                                            const double gamma,
                                            const double dt) const {
    (void)gamma;
    ApplyBoundaryConditions(layer, mesh);
    UpdateInternalEnergyZip(mesh, workspace, dt);
    RebuildConservativeEnergy(layer, mesh, workspace);
    ApplyBoundaryConditions(layer, mesh);
}

void MaderSpatialOperator::Phase4_Transport(DataLayer& layer,
                                            const Mesh& mesh,
                                            Workspace& workspace,
                                            const double gamma,
                                            const double dt) const {
    (void)gamma;

    ApplyBoundaryConditions(layer, mesh);
    ZeroTransportAccumulators(workspace);

    if (transport_params_.method != MaderTransportMethod::DonorAcceptor) {
        throw std::logic_error("MaderSpatialOperator: only DonorAcceptor transport is implemented");
    }

    TransportDonorAcceptor(mesh, workspace, layer, dt);
}

void MaderSpatialOperator::Phase5_Finalize(DataLayer& layer,
                                           const Mesh& mesh,
                                           Workspace& workspace,
                                           const double gamma,
                                           const double dt) const {
    (void)dt;
    ApplyTransportAccumulators(layer, mesh, workspace, gamma);
    ReconstructCellFieldsFromConservative(layer, mesh, workspace, gamma);
    ApplyBoundaryConditions(layer, mesh);
}

void MaderSpatialOperator::ApplyBoundaryConditions(DataLayer& layer, const Mesh& mesh) const {
    if (!boundary_manager_) {
        return;
    }

    boundary_manager_->UpdateHalo(layer, mesh);
    boundary_manager_->ApplyPhysicalBc(layer, mesh);
}

void MaderSpatialOperator::ComputeCellCenteredThermodynamics(DataLayer& layer,
                                                             const Mesh& mesh,
                                                             Workspace& workspace,
                                                             const double gamma,
                                                             const double dt) const {
    (void)gamma;

    auto& U = layer.U();
    auto& W = workspace.W();
    auto& T = workspace.Temperature();
    auto& I = workspace.InternalEnergy();
    auto& reactant = layer.ReactantMassFraction();

    const int sx = mesh.GetSx();
    const int sy = mesh.GetSy();
    const int sz = mesh.GetSz();

    const double rho_floor = GetRhoFloor(eos_);

    for (int k = 0; k < sz; ++k) {
        for (int j = 0; j < sy; ++j) {
            for (int i = 0; i < sx; ++i) {
                if (!mesh.IsFluidCell(i, j, k)) {
                    W(Workspace::k_rho, i, j, k) = 0.0;
                    W(Workspace::k_u, i, j, k) = 0.0;
                    W(Workspace::k_v, i, j, k) = 0.0;
                    W(Workspace::k_w, i, j, k) = 0.0;
                    W(Workspace::k_p, i, j, k) = 0.0;
                    T(i, j, k) = 0.0;
                    I(i, j, k) = 0.0;
                    reactant(i, j, k) = chemistry_params_.reactant_floor;
                    continue;
                }

                const double rho_in = U(var::rho, i, j, k);
                const double rhoU_in = U(var::rhoU, i, j, k);
                const double rhoV_in = U(var::rhoV, i, j, k);
                const double rhoW_in = U(var::rhoW, i, j, k);
                const double E_in = U(var::E, i, j, k);

                const double rho = std::max(rho_in, rho_floor);
                const double inv_rho = 1.0 / rho;

                const double u = rhoU_in * inv_rho;
                const double v = rhoV_in * inv_rho;
                const double w = rhoW_in * inv_rho;

                const double kinetic = 0.5 * (u * u + v * v + w * w);
                double I_cell = std::max(E_in * inv_rho - kinetic, 0.0);

                const double lambda_old = std::clamp(
                    reactant(i, j, k),
                    chemistry_params_.reactant_floor,
                    1.0
                );

                const EosCellInput eos_in_1{
                    .rho = rho,
                    .I = I_cell,
                    .lambda = lambda_old
                };
                const EosCellOutput eos_out_1 = eos_.Evaluate(eos_in_1);

                double lambda_new = lambda_old;

                if (chemistry_params_.enabled) {
                    const bool hot_enough = eos_out_1.T > chemistry_params_.min_temperature;
                    const bool enough_reactant = lambda_old > chemistry_params_.min_reactant;
                    const bool delay_passed = cycle_counter_ >= chemistry_params_.delay_cycles;

                    if (hot_enough && enough_reactant && delay_passed) {
                        const double exponent =
                            -chemistry_params_.activation_energy /
                            (chemistry_params_.gas_constant * std::max(eos_out_1.T, k_eps));

                        const double reaction_rate =
                            chemistry_params_.z_freq * lambda_old * std::exp(exponent);

                        lambda_new = lambda_old - dt * reaction_rate;
                        lambda_new = std::max(lambda_new, chemistry_params_.reactant_floor);

                        if (lambda_new < chemistry_params_.min_reactant) {
                            lambda_new = chemistry_params_.reactant_floor;
                        }

                        lambda_new = std::clamp(lambda_new, chemistry_params_.reactant_floor, 1.0);
                    }
                }

                const double d_lambda = lambda_old - lambda_new;
                I_cell += chemistry_params_.heat_release * d_lambda;
                I_cell = std::max(I_cell, 0.0);

                const EosCellInput eos_in_2{
                    .rho = rho,
                    .I = I_cell,
                    .lambda = lambda_new
                };
                const EosCellOutput eos_out_2 = eos_.Evaluate(eos_in_2);

                W(Workspace::k_rho, i, j, k) = rho;
                W(Workspace::k_u, i, j, k) = u;
                W(Workspace::k_v, i, j, k) = v;
                W(Workspace::k_w, i, j, k) = w;
                W(Workspace::k_p, i, j, k) = eos_out_2.P;

                I(i, j, k) = I_cell;
                T(i, j, k) = eos_out_2.T;
                reactant(i, j, k) = lambda_new;

                U(var::rho, i, j, k) = rho;
                U(var::rhoU, i, j, k) = rho * u;
                U(var::rhoV, i, j, k) = rho * v;
                U(var::rhoW, i, j, k) = rho * w;
                U(var::E, i, j, k) = rho * (I_cell + kinetic);
            }
        }
    }
}

void MaderSpatialOperator::InitializeFaceVelocitiesIfNeeded(DataLayer&,
                                                            const Mesh& mesh,
                                                            Workspace& workspace) const {
    RebuildFaceVelocitiesFromCellCentered(mesh, workspace);
}

void MaderSpatialOperator::SaveOldFaceVelocities(Workspace& workspace) const {
    workspace.UxOld() = workspace.Ux();
    workspace.VyOld() = workspace.Vy();
    workspace.WzOld() = workspace.Wz();
}

void MaderSpatialOperator::ComputeArtificialViscosity(const Mesh& mesh,
                                                      Workspace& workspace) const {
    auto& W = workspace.W();
    auto& Ux = workspace.Ux();
    auto& Vy = workspace.Vy();
    auto& Wz = workspace.Wz();
    auto& Q = workspace.Q();

    const int sx = mesh.GetSx();
    const int sy = mesh.GetSy();
    const int sz = mesh.GetSz();

    Q.fill(0.0);

    if (!viscosity_params_.enabled || viscosity_params_.coefficient <= 0.0) {
        return;
    }

    for (int k = 0; k < sz; ++k) {
        for (int j = 0; j < sy; ++j) {
            for (int i = 0; i < sx; ++i) {
                if (!mesh.IsFluidCell(i, j, k)) {
                    continue;
                }

                const double rho = std::max(W(Workspace::k_rho, i, j, k), k_eps);

                const double comp_x = Ux(i, j, k) - Ux(i + 1, j, k);
                const double cx = std::max(comp_x, 0.0);

                const double comp_y = Vy(i, j, k) - Vy(i, j + 1, k);
                const double cy = std::max(comp_y, 0.0);

                const double comp_z = Wz(i, j, k) - Wz(i, j, k + 1);
                const double cz = std::max(comp_z, 0.0);

                if (cx > 0.0 || cy > 0.0 || cz > 0.0) {
                    Q(i, j, k) = viscosity_params_.coefficient * rho * (cx * cx + cy * cy + cz * cz);
                }
            }
        }
    }
}

void MaderSpatialOperator::UpdateUxFaces(const Mesh& mesh,
                                         Workspace& workspace,
                                         const double dt) const {
    auto& W = workspace.W();
    auto& Q = workspace.Q();
    auto& Ux = workspace.Ux();

    const auto& dx = mesh.Dx();

    const int sy = mesh.GetSy();
    const int sz = mesh.GetSz();

    for (int k = 0; k < sz; ++k) {
        for (int j = 0; j < sy; ++j) {
            for (int iface = 1; iface < mesh.GetSx(); ++iface) {
                const int il = iface - 1;
                const int ir = iface;

                if (!mesh.IsFluidCell(il, j, k) || !mesh.IsFluidCell(ir, j, k)) {
                    continue;
                }

                const double rho_l = std::max(W(Workspace::k_rho, il, j, k), k_eps);
                const double rho_r = std::max(W(Workspace::k_rho, ir, j, k), k_eps);

                const double denom = 0.5 * (rho_l * dx(il) + rho_r * dx(ir));
                if (denom <= 0.0) {
                    continue;
                }

                const double p_l = W(Workspace::k_p, il, j, k);
                const double p_r = W(Workspace::k_p, ir, j, k);

                const double q_l = Q(il, j, k);
                const double q_r = Q(ir, j, k);

                Ux(iface, j, k) += dt * ((p_l - p_r) + (q_l - q_r)) / denom;
            }
        }
    }
}

void MaderSpatialOperator::UpdateVyFaces(const Mesh& mesh,
                                         Workspace& workspace,
                                         const double dt) const {
    auto& W = workspace.W();
    auto& Q = workspace.Q();
    auto& Vy = workspace.Vy();

    const auto& dy = mesh.Dy();

    const int sx = mesh.GetSx();
    const int sz = mesh.GetSz();

    for (int k = 0; k < sz; ++k) {
        for (int jface = 1; jface < mesh.GetSy(); ++jface) {
            const int jb = jface - 1;
            const int jt = jface;

            for (int i = 0; i < sx; ++i) {
                if (!mesh.IsFluidCell(i, jb, k) || !mesh.IsFluidCell(i, jt, k)) {
                    continue;
                }

                const double rho_b = std::max(W(Workspace::k_rho, i, jb, k), k_eps);
                const double rho_t = std::max(W(Workspace::k_rho, i, jt, k), k_eps);

                const double denom = 0.5 * (rho_b * dy(jb) + rho_t * dy(jt));
                if (denom <= 0.0) {
                    continue;
                }

                const double p_b = W(Workspace::k_p, i, jb, k);
                const double p_t = W(Workspace::k_p, i, jt, k);

                const double q_b = Q(i, jb, k);
                const double q_t = Q(i, jt, k);

                Vy(i, jface, k) += dt * ((p_b - p_t) + (q_b - q_t)) / denom;
            }
        }
    }
}

void MaderSpatialOperator::UpdateInternalEnergyZip(const Mesh& mesh,
                                                   Workspace& workspace,
                                                   const double dt) const {
    auto& W = workspace.W();
    auto& I = workspace.InternalEnergy();
    auto& Q = workspace.Q();
    auto& Ux = workspace.Ux();
    auto& Vy = workspace.Vy();
    auto& UxOld = workspace.UxOld();
    auto& VyOld = workspace.VyOld();

    const auto& dx = mesh.Dx();
    const auto& dy = mesh.Dy();

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

                const double rho = std::max(W(Workspace::k_rho, i, j, k), k_eps);

                const double ux_l_half = 0.5 * (UxOld(i, j, k) + Ux(i, j, k));
                const double ux_r_half = 0.5 * (UxOld(i + 1, j, k) + Ux(i + 1, j, k));

                const double vy_b_half = 0.5 * (VyOld(i, j, k) + Vy(i, j, k));
                const double vy_t_half = 0.5 * (VyOld(i, j + 1, k) + Vy(i, j + 1, k));

                const double div_half =
                    (ux_r_half - ux_l_half) / std::max(dx(i), k_eps) +
                    (vy_t_half - vy_b_half) / std::max(dy(j), k_eps);

                const double p = W(Workspace::k_p, i, j, k);
                const double q = Q(i, j, k);

                I(i, j, k) = std::max(I(i, j, k) - dt * (p + q) * div_half / rho, 0.0);
            }
        }
    }
}

void MaderSpatialOperator::RebuildConservativeEnergy(DataLayer& layer,
                                                     const Mesh& mesh,
                                                     Workspace& workspace) const {
    auto& U = layer.U();
    auto& W = workspace.W();
    auto& I = workspace.InternalEnergy();
    auto& Ux = workspace.Ux();
    auto& Vy = workspace.Vy();

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

                const double rho = std::max(W(Workspace::k_rho, i, j, k), k_eps);

                const double u = 0.5 * (Ux(i, j, k) + Ux(i + 1, j, k));
                const double v = 0.5 * (Vy(i, j, k) + Vy(i, j + 1, k));
                const double w = W(Workspace::k_w, i, j, k);

                W(Workspace::k_u, i, j, k) = u;
                W(Workspace::k_v, i, j, k) = v;

                U(var::rho, i, j, k) = rho;
                U(var::rhoU, i, j, k) = rho * u;
                U(var::rhoV, i, j, k) = rho * v;
                U(var::rhoW, i, j, k) = rho * w;
                U(var::E, i, j, k) = rho * (I(i, j, k) + 0.5 * (u * u + v * v + w * w));
            }
        }
    }
}

void MaderSpatialOperator::ZeroTransportAccumulators(Workspace& workspace) const {
    workspace.ZeroD();
}

void MaderSpatialOperator::TransportDonorAcceptor(const Mesh& mesh,
                                                  Workspace& workspace,
                                                  DataLayer& layer,
                                                  const double dt) const {
    TransportDonorAcceptorR(mesh, workspace, layer, dt);
    TransportDonorAcceptorZ(mesh, workspace, layer, dt);
}

void MaderSpatialOperator::TransportDonorAcceptorR(const Mesh& mesh,
                                                   Workspace& workspace,
                                                   DataLayer& layer,
                                                   const double dt) const {
    auto& W = workspace.W();
    auto& I = workspace.InternalEnergy();
    auto& D = workspace.D();
    auto& Ux = workspace.Ux();
    auto& reactant = layer.ReactantMassFraction();

    const auto& dx = mesh.Dx();

    const int i0 = mesh.GetCoreStartX();
    const int i1 = mesh.GetCoreEndExclusiveX();
    const int j0 = mesh.GetCoreStartY();
    const int j1 = mesh.GetCoreEndExclusiveY();
    const int k0 = mesh.GetCoreStartZ();
    const int k1 = mesh.GetCoreEndExclusiveZ();

    for (int k = k0; k < k1; ++k) {
        for (int j = j0; j < j1; ++j) {
            for (int iface = i0; iface <= i1; ++iface) {
                const int il = iface - 1;
                const int ir = iface;

                if (il < 0 || ir >= mesh.GetSx()) {
                    continue;
                }
                if (!mesh.IsFluidCell(il, j, k) || !mesh.IsFluidCell(ir, j, k)) {
                    continue;
                }

                const double u_face = Ux(iface, j, k);
                if (std::abs(u_face) <= k_eps) {
                    continue;
                }

                const int donor_i = (u_face >= 0.0) ? il : ir;
                const int accept_i = (u_face >= 0.0) ? ir : il;

                const double rho_d = std::max(W(Workspace::k_rho, donor_i, j, k), k_eps);
                const double alpha = u_face * dt / std::max(dx(donor_i), k_eps);
                const double dmass = rho_d * std::abs(alpha);

                if (dmass <= 0.0) {
                    continue;
                }

                const double u_d = W(Workspace::k_u, donor_i, j, k);
                const double v_d = W(Workspace::k_v, donor_i, j, k);
                const double w_d = W(Workspace::k_w, donor_i, j, k);
                const double I_d = I(donor_i, j, k);
                const double E_d = I_d + 0.5 * (u_d * u_d + v_d * v_d + w_d * w_d);

                double w_transfer = 0.0;
                if (transport_params_.composition_mode == MaderCompositionTransportMode::Standard) {
                    w_transfer = ComputeTransferredReactantFractionStandard(
                        reactant(donor_i, j, k));
                }
                else {
                    w_transfer = ComputeTransferredReactantFractionShargatovR(
                        mesh, layer,
                        donor_i, j, k,
                        accept_i, j, k);
                }

                const double de = E_d * dmass;
                const double dw = w_transfer * dmass;
                const double dpu = u_d * dmass;
                const double dpv = v_d * dmass;

                AddTransportContributionIfCore(D, mesh, accept_i, j, k, +dmass, +de, +dw, +dpu, +dpv);
                AddTransportContributionIfCore(D, mesh, donor_i, j, k, -dmass, -de, -dw, -dpu, -dpv);
            }
        }
    }
}

void MaderSpatialOperator::TransportDonorAcceptorZ(const Mesh& mesh,
                                                   Workspace& workspace,
                                                   DataLayer& layer,
                                                   const double dt) const {
    auto& W = workspace.W();
    auto& I = workspace.InternalEnergy();
    auto& D = workspace.D();
    auto& Vy = workspace.Vy();
    auto& reactant = layer.ReactantMassFraction();

    const auto& dy = mesh.Dy();

    const int i0 = mesh.GetCoreStartX();
    const int i1 = mesh.GetCoreEndExclusiveX();
    const int j0 = mesh.GetCoreStartY();
    const int j1 = mesh.GetCoreEndExclusiveY();
    const int k0 = mesh.GetCoreStartZ();
    const int k1 = mesh.GetCoreEndExclusiveZ();

    for (int k = k0; k < k1; ++k) {
        for (int jface = j0; jface <= j1; ++jface) {
            const int jb = jface - 1;
            const int jt = jface;

            if (jb < 0 || jt >= mesh.GetSy()) {
                continue;
            }

            for (int i = i0; i < i1; ++i) {
                if (!mesh.IsFluidCell(i, jb, k) || !mesh.IsFluidCell(i, jt, k)) {
                    continue;
                }

                const double v_face = Vy(i, jface, k);
                if (std::abs(v_face) <= k_eps) {
                    continue;
                }

                const int donor_j = (v_face >= 0.0) ? jb : jt;
                const int accept_j = (v_face >= 0.0) ? jt : jb;

                const double rho_d = std::max(W(Workspace::k_rho, i, donor_j, k), k_eps);
                const double beta = v_face * dt / std::max(dy(donor_j), k_eps);
                const double dmass = rho_d * std::abs(beta);

                if (dmass <= 0.0) {
                    continue;
                }

                const double u_d = W(Workspace::k_u, i, donor_j, k);
                const double v_d = W(Workspace::k_v, i, donor_j, k);
                const double w_d = W(Workspace::k_w, i, donor_j, k);
                const double I_d = I(i, donor_j, k);
                const double E_d = I_d + 0.5 * (u_d * u_d + v_d * v_d + w_d * w_d);

                double w_transfer = 0.0;
                if (transport_params_.composition_mode == MaderCompositionTransportMode::Standard) {
                    w_transfer = ComputeTransferredReactantFractionStandard(
                        reactant(i, donor_j, k));
                }
                else {
                    w_transfer = ComputeTransferredReactantFractionShargatovZ(
                        mesh, layer,
                        i, donor_j, k,
                        i, accept_j, k);
                }

                const double de = E_d * dmass;
                const double dw = w_transfer * dmass;
                const double dpu = u_d * dmass;
                const double dpv = v_d * dmass;

                AddTransportContributionIfCore(D, mesh, i, accept_j, k, +dmass, +de, +dw, +dpu, +dpv);
                AddTransportContributionIfCore(D, mesh, i, donor_j, k, -dmass, -de, -dw, -dpu, -dpv);
            }
        }
    }
}

void MaderSpatialOperator::ApplyTransportAccumulators(DataLayer& layer,
                                                      const Mesh& mesh,
                                                      Workspace& workspace,
                                                      const double gamma) const {
    auto& U = layer.U();
    auto& W = workspace.W();
    auto& I = workspace.InternalEnergy();
    auto& T = workspace.Temperature();
    auto& reactant = layer.ReactantMassFraction();
    auto& D = workspace.D();

    const double rho_floor = GetRhoFloor(eos_);

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

                const double rho_old = std::max(W(Workspace::k_rho, i, j, k), rho_floor);
                const double u_old = W(Workspace::k_u, i, j, k);
                const double v_old = W(Workspace::k_v, i, j, k);
                const double w_old = W(Workspace::k_w, i, j, k);
                const double I_old = I(i, j, k);
                const double lambda_old = reactant(i, j, k);

                const double E_old = I_old + 0.5 * (u_old * u_old + v_old * v_old + w_old * w_old);

                const double dm = D(Workspace::k_dm, i, j, k);
                const double de = D(Workspace::k_de, i, j, k);
                const double dw = D(Workspace::k_dw, i, j, k);
                const double dpu = D(Workspace::k_dpu, i, j, k);
                const double dpv = D(Workspace::k_dpv, i, j, k);

                const double rho_new = std::max(rho_old + dm, rho_floor);

                const double rho_u_new = rho_old * u_old + dpu;
                const double rho_v_new = rho_old * v_old + dpv;
                const double rho_E_new = rho_old * E_old + de;
                const double rho_W_new = rho_old * lambda_old + dw;

                const double u_new = rho_u_new / rho_new;
                const double v_new = rho_v_new / rho_new;
                const double w_new = w_old;

                double lambda_new = rho_W_new / rho_new;
                lambda_new = std::clamp(lambda_new, chemistry_params_.reactant_floor, 1.0);

                double I_new =
                    rho_E_new / rho_new -
                    0.5 * (u_new * u_new + v_new * v_new + w_new * w_new);
                I_new = std::max(I_new, 0.0);

                const EosCellInput eos_in{
                    .rho = rho_new,
                    .I = I_new,
                    .lambda = lambda_new
                };
                const EosCellOutput eos_out = eos_.Evaluate(eos_in);

                W(Workspace::k_rho, i, j, k) = rho_new;
                W(Workspace::k_u, i, j, k) = u_new;
                W(Workspace::k_v, i, j, k) = v_new;
                W(Workspace::k_w, i, j, k) = w_new;
                W(Workspace::k_p, i, j, k) = eos_out.P;

                I(i, j, k) = I_new;
                T(i, j, k) = eos_out.T;
                reactant(i, j, k) = lambda_new;

                U(var::rho, i, j, k) = rho_new;
                U(var::rhoU, i, j, k) = rho_new * u_new;
                U(var::rhoV, i, j, k) = rho_new * v_new;
                U(var::rhoW, i, j, k) = rho_new * w_new;
                U(var::E, i, j, k) = rho_new * (I_new + 0.5 * (u_new * u_new + v_new * v_new + w_new * w_new));
            }
        }
    }
}

void MaderSpatialOperator::ReconstructCellFieldsFromConservative(DataLayer& layer,
                                                                 const Mesh& mesh,
                                                                 Workspace& workspace,
                                                                 const double gamma) const {
    (void)gamma;

    auto& U = layer.U();
    auto& W = workspace.W();
    auto& T = workspace.Temperature();
    auto& I = workspace.InternalEnergy();
    auto& reactant = layer.ReactantMassFraction();

    const int sx = mesh.GetSx();
    const int sy = mesh.GetSy();
    const int sz = mesh.GetSz();

    const double rho_floor = GetRhoFloor(eos_);

    for (int k = 0; k < sz; ++k) {
        for (int j = 0; j < sy; ++j) {
            for (int i = 0; i < sx; ++i) {
                if (!mesh.IsFluidCell(i, j, k)) {
                    W(Workspace::k_rho, i, j, k) = 0.0;
                    W(Workspace::k_u, i, j, k) = 0.0;
                    W(Workspace::k_v, i, j, k) = 0.0;
                    W(Workspace::k_w, i, j, k) = 0.0;
                    W(Workspace::k_p, i, j, k) = 0.0;
                    I(i, j, k) = 0.0;
                    T(i, j, k) = 0.0;
                    reactant(i, j, k) = chemistry_params_.reactant_floor;
                    continue;
                }

                const double rho = std::max(U(var::rho, i, j, k), rho_floor);
                const double inv_rho = 1.0 / rho;

                const double u = U(var::rhoU, i, j, k) * inv_rho;
                const double v = U(var::rhoV, i, j, k) * inv_rho;
                const double w = U(var::rhoW, i, j, k) * inv_rho;

                const double kinetic = 0.5 * (u * u + v * v + w * w);
                const double I_cell = std::max(U(var::E, i, j, k) * inv_rho - kinetic, 0.0);
                const double lambda = std::clamp(
                    reactant(i, j, k),
                    chemistry_params_.reactant_floor,
                    1.0
                );

                const EosCellInput eos_in{
                    .rho = rho,
                    .I = I_cell,
                    .lambda = lambda
                };
                const EosCellOutput eos_out = eos_.Evaluate(eos_in);

                W(Workspace::k_rho, i, j, k) = rho;
                W(Workspace::k_u, i, j, k) = u;
                W(Workspace::k_v, i, j, k) = v;
                W(Workspace::k_w, i, j, k) = w;
                W(Workspace::k_p, i, j, k) = eos_out.P;

                I(i, j, k) = I_cell;
                T(i, j, k) = eos_out.T;
                reactant(i, j, k) = lambda;
            }
        }
    }

    RebuildFaceVelocitiesFromCellCentered(mesh, workspace);
}

bool MaderSpatialOperator::IsPureA(const double w) const {
    return std::abs(w) <= 1e-12;
}

bool MaderSpatialOperator::IsPureB(const double w) const {
    return std::abs(w - 1.0) <= 1e-12;
}

bool MaderSpatialOperator::IsMixed(const double w) const {
    return w > 1e-12 && w < 1.0 - 1e-12;
}

double MaderSpatialOperator::ComputeTransferredReactantFractionStandard(
    const double donor_w) const {
    return std::clamp(donor_w, 0.0, 1.0);
}

double MaderSpatialOperator::ComputeTransferredReactantFractionShargatovR(
    const Mesh& mesh,
    const DataLayer& layer,
    const int donor_i, const int donor_j, const int donor_k,
    const int accept_i, const int accept_j, const int accept_k) const {
    return ComputeTransferredReactantFractionShargatovImpl(
        mesh, layer,
        donor_i, donor_j, donor_k,
        accept_i, accept_j, accept_k);
}

double MaderSpatialOperator::ComputeTransferredReactantFractionShargatovZ(
    const Mesh& mesh,
    const DataLayer& layer,
    const int donor_i, const int donor_j, const int donor_k,
    const int accept_i, const int accept_j, const int accept_k) const {
    return ComputeTransferredReactantFractionShargatovImpl(
        mesh, layer,
        donor_i, donor_j, donor_k,
        accept_i, accept_j, accept_k);
}
