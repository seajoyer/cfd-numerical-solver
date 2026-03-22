#include "bc/NonReflectiveBoundary.hpp"

#include <cmath>
#include <stdexcept>

#include "data/DataLayer.hpp"
#include "data/Variables.hpp"
#include "data/Mesh.hpp"

namespace {
    auto ConservativeToPrimitive(const FarfieldConservative& U, double gamma) -> PrimitiveCell {
        PrimitiveCell W;
        W.rho = U.rho;

        if (!(W.rho > 0.0)) {
            throw std::runtime_error("NonReflectiveBoundary: far-field rho must be positive");
        }

        W.u = U.rhoU / W.rho;
        W.v = U.rhoV / W.rho;
        W.w = U.rhoW / W.rho;

        const double kinetic = 0.5 * W.rho * (W.u * W.u + W.v * W.v + W.w * W.w);
        W.P = (gamma - 1.0) * (U.E - kinetic);

        if (!(W.P > 0.0)) {
            throw std::runtime_error("NonReflectiveBoundary: far-field pressure must be positive");
        }

        return W;
    }

    auto ConservativeToPrimitive(const xt::xtensor<double, 4>& U,
                                 int i, int j, int k,
                                 double gamma,
                                 double rho_floor = 1e-14,
                                 double p_floor = 1e-14) -> PrimitiveCell {
        PrimitiveCell W;

        W.rho = std::max(U(DataLayer::k_rho, i, j, k), rho_floor);
        W.u = U(DataLayer::k_rhoU, i, j, k) / W.rho;
        W.v = U(DataLayer::k_rhoV, i, j, k) / W.rho;
        W.w = U(DataLayer::k_rhoW, i, j, k) / W.rho;

        const double kinetic = 0.5 * W.rho * (W.u * W.u + W.v * W.v + W.w * W.w);
        W.P = std::max((gamma - 1.0) * (U(DataLayer::k_E, i, j, k) - kinetic), p_floor);

        return W;
    }

    auto PrimitiveToConservative(const PrimitiveCell& W, double gamma) -> FarfieldConservative {
        FarfieldConservative U;
        U.rho = W.rho;
        U.rhoU = W.rho * W.u;
        U.rhoV = W.rho * W.v;
        U.rhoW = W.rho * W.w;

        const double kinetic = 0.5 * W.rho * (W.u * W.u + W.v * W.v + W.w * W.w);
        U.E = W.P / (gamma - 1.0) + kinetic;

        return U;
    }

    void SetPlaneToStateX(xt::xtensor<double, 4>& U, int dst_i, const FarfieldConservative& s) {
        xt::view(U, DataLayer::k_rho, dst_i, xt::all(), xt::all()) = s.rho;
        xt::view(U, DataLayer::k_rhoU, dst_i, xt::all(), xt::all()) = s.rhoU;
        xt::view(U, DataLayer::k_rhoV, dst_i, xt::all(), xt::all()) = s.rhoV;
        xt::view(U, DataLayer::k_rhoW, dst_i, xt::all(), xt::all()) = s.rhoW;
        xt::view(U, DataLayer::k_E, dst_i, xt::all(), xt::all()) = s.E;
    }

    void SetPlaneToStateY(xt::xtensor<double, 4>& U, int dst_j, const FarfieldConservative& s) {
        xt::view(U, DataLayer::k_rho, xt::all(), dst_j, xt::all()) = s.rho;
        xt::view(U, DataLayer::k_rhoU, xt::all(), dst_j, xt::all()) = s.rhoU;
        xt::view(U, DataLayer::k_rhoV, xt::all(), dst_j, xt::all()) = s.rhoV;
        xt::view(U, DataLayer::k_rhoW, xt::all(), dst_j, xt::all()) = s.rhoW;
        xt::view(U, DataLayer::k_E, xt::all(), dst_j, xt::all()) = s.E;
    }

    void SetPlaneToStateZ(xt::xtensor<double, 4>& U, int dst_k, const FarfieldConservative& s) {
        xt::view(U, DataLayer::k_rho, xt::all(), xt::all(), dst_k) = s.rho;
        xt::view(U, DataLayer::k_rhoU, xt::all(), xt::all(), dst_k) = s.rhoU;
        xt::view(U, DataLayer::k_rhoV, xt::all(), xt::all(), dst_k) = s.rhoV;
        xt::view(U, DataLayer::k_rhoW, xt::all(), xt::all(), dst_k) = s.rhoW;
        xt::view(U, DataLayer::k_E, xt::all(), xt::all(), dst_k) = s.E;
    }
} // namespace

NonReflectiveBoundary::NonReflectiveBoundary(const FarfieldConservative& farfield_U, double gamma)
    : farfield_U_(farfield_U), gamma_(gamma) {}

void NonReflectiveBoundary::Apply(DataLayer& layer, const Mesh& mesh, const Axis axis, const Side side) const {
    const int ng = mesh.GetPadding();
    if (ng == 0) {
        return;
    }

    const int dim = mesh.GetDim();
    if (axis == Axis::Y && dim < 2) {
        return;
    }
    if (axis == Axis::Z && dim < 3) {
        return;
    }

    auto& U = layer.U();

    const int i0 = mesh.GetCoreStartX();
    const int i1 = mesh.GetCoreEndExclusiveX();
    const int j0 = mesh.GetCoreStartY();
    const int j1 = mesh.GetCoreEndExclusiveY();
    const int k0 = mesh.GetCoreStartZ();
    const int k1 = mesh.GetCoreEndExclusiveZ();

    const int outward_sign = side == Side::Left ? -1 : 1;

    int ii = i0;
    int jj = j0;
    int kk = k0;

    if (axis == Axis::X) {
        ii = side == Side::Left ? i0 : (i1 - 1);
    }
    else if (axis == Axis::Y) {
        jj = side == Side::Left ? j0 : (j1 - 1);
    }
    else {
        kk = side == Side::Left ? k0 : (k1 - 1);
    }

    const PrimitiveCell Wi = ConservativeToPrimitive(U, ii, jj, kk, gamma_);
    const PrimitiveCell Wf = ConservativeToPrimitive(farfield_U_, gamma_);

    double ui_n = 0.0;
    double uf_n = 0.0;
    double ui_t1 = 0.0;
    double ui_t2 = 0.0;
    double uf_t1 = 0.0;
    double uf_t2 = 0.0;

    if (axis == Axis::X) {
        ui_n = outward_sign * Wi.u;
        uf_n = outward_sign * Wf.u;
        ui_t1 = Wi.v;
        ui_t2 = Wi.w;
        uf_t1 = Wf.v;
        uf_t2 = Wf.w;
    }
    else if (axis == Axis::Y) {
        ui_n = outward_sign * Wi.v;
        uf_n = outward_sign * Wf.v;
        ui_t1 = Wi.u;
        ui_t2 = Wi.w;
        uf_t1 = Wf.u;
        uf_t2 = Wf.w;
    }
    else {
        ui_n = outward_sign * Wi.w;
        uf_n = outward_sign * Wf.w;
        ui_t1 = Wi.u;
        ui_t2 = Wi.v;
        uf_t1 = Wf.u;
        uf_t2 = Wf.v;
    }

    const double ci = SoundSpeed(Wi, gamma_);
    const double cf = SoundSpeed(Wf, gamma_);

    const double mach_out = ui_n / ci;

    FarfieldConservative ghost_state{};

    if (mach_out >= 1.0) {
        ghost_state = PrimitiveToConservative(Wi, gamma_);
    }
    else if (mach_out <= -1.0) {
        ghost_state = farfield_U_;
    }
    else {
        const double J_plus_i = ui_n + 2.0 * ci / (gamma_ - 1.0);
        const double J_minus_f = uf_n - 2.0 * cf / (gamma_ - 1.0);

        const double ub_n = 0.5 * (J_plus_i + J_minus_f);
        const double cb = 0.25 * (gamma_ - 1.0) * (J_plus_i - J_minus_f);

        if (!(cb > 0.0)) {
            ghost_state = ui_n >= 0.0 ? PrimitiveToConservative(Wi, gamma_) : farfield_U_;
        }
        else {
            const bool outflow = ui_n >= 0.0;
            const PrimitiveCell& Wref = outflow ? Wi : Wf;

            const double entropy = Wref.P / std::pow(Wref.rho, gamma_);
            const double rho_b = std::pow((cb * cb) / (gamma_ * entropy), 1.0 / (gamma_ - 1.0));
            const double p_b = entropy * std::pow(rho_b, gamma_);

            PrimitiveCell Wb;
            Wb.rho = rho_b;
            Wb.P = p_b;

            if (axis == Axis::X) {
                Wb.u = outward_sign * ub_n;
                Wb.v = outflow ? ui_t1 : uf_t1;
                Wb.w = outflow ? ui_t2 : uf_t2;
            }
            else if (axis == Axis::Y) {
                Wb.v = outward_sign * ub_n;
                Wb.u = outflow ? ui_t1 : uf_t1;
                Wb.w = outflow ? ui_t2 : uf_t2;
            }
            else {
                Wb.w = outward_sign * ub_n;
                Wb.u = outflow ? ui_t1 : uf_t1;
                Wb.v = outflow ? ui_t2 : uf_t2;
            }

            ghost_state = PrimitiveToConservative(Wb, gamma_);
        }
    }

    for (int g = 0; g < ng; ++g) {
        if (axis == Axis::X) {
            const int dst_i = side == Side::Left ? (i0 - 1 - g) : (i1 + g);
            SetPlaneToStateX(U, dst_i, ghost_state);
        }
        else if (axis == Axis::Y) {
            const int dst_j = side == Side::Left ? (j0 - 1 - g) : (j1 + g);
            SetPlaneToStateY(U, dst_j, ghost_state);
        }
        else {
            const int dst_k = side == Side::Left ? (k0 - 1 - g) : (k1 + g);
            SetPlaneToStateZ(U, dst_k, ghost_state);
        }
    }
}
