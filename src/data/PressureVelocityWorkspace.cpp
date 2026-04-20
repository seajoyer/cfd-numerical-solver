#include "data/PressureVelocityWorkspace.hpp"

#include "data/Mesh.hpp"

void PressureVelocityWorkspace::ResizeFrom(const Mesh& mesh) {
    const std::size_t sx = static_cast<std::size_t>(mesh.GetSx());
    const std::size_t sy = static_cast<std::size_t>(mesh.GetSy());
    const std::size_t sz = static_cast<std::size_t>(mesh.GetSz());

    if (sx == sx_ && sy == sy_ && sz == sz_ && IsAllocated()) {
        return;
    }

    Allocate(sx, sy, sz);
}

xt::xtensor<double, 3>& PressureVelocityWorkspace::PressureCorrection() {
    return pressure_correction_;
}

const xt::xtensor<double, 3>& PressureVelocityWorkspace::PressureCorrection() const {
    return pressure_correction_;
}

xt::xtensor<double, 3>& PressureVelocityWorkspace::PressureRhs() {
    return pressure_rhs_;
}

const xt::xtensor<double, 3>& PressureVelocityWorkspace::PressureRhs() const {
    return pressure_rhs_;
}

xt::xtensor<double, 3>& PressureVelocityWorkspace::ApE() {
    return ap_e_;
}

const xt::xtensor<double, 3>& PressureVelocityWorkspace::ApE() const {
    return ap_e_;
}

xt::xtensor<double, 3>& PressureVelocityWorkspace::ApW() {
    return ap_w_;
}

const xt::xtensor<double, 3>& PressureVelocityWorkspace::ApW() const {
    return ap_w_;
}

xt::xtensor<double, 3>& PressureVelocityWorkspace::ApN() {
    return ap_n_;
}

const xt::xtensor<double, 3>& PressureVelocityWorkspace::ApN() const {
    return ap_n_;
}

xt::xtensor<double, 3>& PressureVelocityWorkspace::ApS() {
    return ap_s_;
}

const xt::xtensor<double, 3>& PressureVelocityWorkspace::ApS() const {
    return ap_s_;
}

xt::xtensor<double, 3>& PressureVelocityWorkspace::ApT() {
    return ap_t_;
}

const xt::xtensor<double, 3>& PressureVelocityWorkspace::ApT() const {
    return ap_t_;
}

xt::xtensor<double, 3>& PressureVelocityWorkspace::ApB() {
    return ap_b_;
}

const xt::xtensor<double, 3>& PressureVelocityWorkspace::ApB() const {
    return ap_b_;
}

xt::xtensor<double, 3>& PressureVelocityWorkspace::ApP() {
    return ap_p_;
}

const xt::xtensor<double, 3>& PressureVelocityWorkspace::ApP() const {
    return ap_p_;
}

xt::xtensor<double, 3>& PressureVelocityWorkspace::UxStar() {
    return ux_star_;
}

const xt::xtensor<double, 3>& PressureVelocityWorkspace::UxStar() const {
    return ux_star_;
}

xt::xtensor<double, 3>& PressureVelocityWorkspace::VyStar() {
    return vy_star_;
}

const xt::xtensor<double, 3>& PressureVelocityWorkspace::VyStar() const {
    return vy_star_;
}

xt::xtensor<double, 3>& PressureVelocityWorkspace::WzStar() {
    return wz_star_;
}

const xt::xtensor<double, 3>& PressureVelocityWorkspace::WzStar() const {
    return wz_star_;
}

xt::xtensor<double, 3>& PressureVelocityWorkspace::APu() {
    return a_pu_;
}

const xt::xtensor<double, 3>& PressureVelocityWorkspace::APu() const {
    return a_pu_;
}

xt::xtensor<double, 3>& PressureVelocityWorkspace::APv() {
    return a_pv_;
}

const xt::xtensor<double, 3>& PressureVelocityWorkspace::APv() const {
    return a_pv_;
}

xt::xtensor<double, 3>& PressureVelocityWorkspace::APw() {
    return a_pw_;
}

const xt::xtensor<double, 3>& PressureVelocityWorkspace::APw() const {
    return a_pw_;
}

void PressureVelocityWorkspace::ZeroPressureCorrection() {
    pressure_correction_.fill(0.0);
}

void PressureVelocityWorkspace::ZeroPressureRhs() {
    pressure_rhs_.fill(0.0);
}

void PressureVelocityWorkspace::ZeroPressureCoefficients() {
    ap_e_.fill(0.0);
    ap_w_.fill(0.0);
    ap_n_.fill(0.0);
    ap_s_.fill(0.0);
    ap_t_.fill(0.0);
    ap_b_.fill(0.0);
    ap_p_.fill(0.0);
}

void PressureVelocityWorkspace::ZeroUxStar() {
    ux_star_.fill(0.0);
}

void PressureVelocityWorkspace::ZeroVyStar() {
    vy_star_.fill(0.0);
}

void PressureVelocityWorkspace::ZeroWzStar() {
    wz_star_.fill(0.0);
}

void PressureVelocityWorkspace::ZeroAPu() {
    a_pu_.fill(0.0);
}

void PressureVelocityWorkspace::ZeroAPv() {
    a_pv_.fill(0.0);
}

void PressureVelocityWorkspace::ZeroAPw() {
    a_pw_.fill(0.0);
}

void PressureVelocityWorkspace::ZeroAll() {
    ZeroPressureCorrection();
    ZeroPressureRhs();
    ZeroPressureCoefficients();
    ZeroMomentumCoefficients();
    ZeroUxStar();
    ZeroVyStar();
    ZeroWzStar();
    ZeroAPu();
    ZeroAPv();
    ZeroAPw();
}

bool PressureVelocityWorkspace::IsAllocated() const {
    return
        pressure_correction_.dimension() == 3 &&
        pressure_rhs_.dimension() == 3 &&
        ap_e_.dimension() == 3 &&
        ap_w_.dimension() == 3 &&
        ap_n_.dimension() == 3 &&
        ap_s_.dimension() == 3 &&
        ap_t_.dimension() == 3 &&
        ap_b_.dimension() == 3 &&
        ap_p_.dimension() == 3 &&
        ux_star_.dimension() == 3 &&
        vy_star_.dimension() == 3 &&
        wz_star_.dimension() == 3 &&
        a_pu_.dimension() == 3 &&
        a_pv_.dimension() == 3 &&
        a_pw_.dimension() == 3 &&
        au_e_.dimension() == 3 &&
        au_w_.dimension() == 3 &&
        au_n_.dimension() == 3 &&
        au_s_.dimension() == 3 &&
        av_e_.dimension() == 3 &&
        av_w_.dimension() == 3 &&
        av_n_.dimension() == 3 &&
        av_s_.dimension() == 3 &&
        sx_ > 0 && sy_ > 0 && sz_ > 0;
}

std::size_t PressureVelocityWorkspace::GetSx() const {
    return sx_;
}

std::size_t PressureVelocityWorkspace::GetSy() const {
    return sy_;
}

std::size_t PressureVelocityWorkspace::GetSz() const {
    return sz_;
}

xt::xtensor<double, 3>& PressureVelocityWorkspace::AuE() {
    return au_e_;
}

const xt::xtensor<double, 3>& PressureVelocityWorkspace::AuE() const {
    return au_e_;
}

xt::xtensor<double, 3>& PressureVelocityWorkspace::AuW() {
    return au_w_;
}

const xt::xtensor<double, 3>& PressureVelocityWorkspace::AuW() const {
    return au_w_;
}

xt::xtensor<double, 3>& PressureVelocityWorkspace::AuN() {
    return au_n_;
}

const xt::xtensor<double, 3>& PressureVelocityWorkspace::AuN() const {
    return au_n_;
}

xt::xtensor<double, 3>& PressureVelocityWorkspace::AuS() {
    return au_s_;
}

const xt::xtensor<double, 3>& PressureVelocityWorkspace::AuS() const {
    return au_s_;
}

xt::xtensor<double, 3>& PressureVelocityWorkspace::AvE() {
    return av_e_;
}

const xt::xtensor<double, 3>& PressureVelocityWorkspace::AvE() const {
    return av_e_;
}

xt::xtensor<double, 3>& PressureVelocityWorkspace::AvW() {
    return av_w_;
}

const xt::xtensor<double, 3>& PressureVelocityWorkspace::AvW() const {
    return av_w_;
}

xt::xtensor<double, 3>& PressureVelocityWorkspace::AvN() {
    return av_n_;
}

const xt::xtensor<double, 3>& PressureVelocityWorkspace::AvN() const {
    return av_n_;
}

xt::xtensor<double, 3>& PressureVelocityWorkspace::AvS() {
    return av_s_;
}

const xt::xtensor<double, 3>& PressureVelocityWorkspace::AvS() const {
    return av_s_;
}

void PressureVelocityWorkspace::ZeroMomentumCoefficients() {
    au_e_.fill(0.0);
    au_w_.fill(0.0);
    au_n_.fill(0.0);
    au_s_.fill(0.0);

    av_e_.fill(0.0);
    av_w_.fill(0.0);
    av_n_.fill(0.0);
    av_s_.fill(0.0);
}

void PressureVelocityWorkspace::Allocate(const std::size_t sx,
                                         const std::size_t sy,
                                         const std::size_t sz) {
    sx_ = sx;
    sy_ = sy;
    sz_ = sz;

    pressure_correction_ = xt::zeros<double>({sx_, sy_, sz_});
    pressure_rhs_ = xt::zeros<double>({sx_, sy_, sz_});

    ap_e_ = xt::zeros<double>({sx_, sy_, sz_});
    ap_w_ = xt::zeros<double>({sx_, sy_, sz_});
    ap_n_ = xt::zeros<double>({sx_, sy_, sz_});
    ap_s_ = xt::zeros<double>({sx_, sy_, sz_});
    ap_t_ = xt::zeros<double>({sx_, sy_, sz_});
    ap_b_ = xt::zeros<double>({sx_, sy_, sz_});
    ap_p_ = xt::zeros<double>({sx_, sy_, sz_});

    ux_star_ = xt::zeros<double>({sx_ + 1, sy_, sz_});
    vy_star_ = xt::zeros<double>({sx_, sy_ + 1, sz_});
    wz_star_ = xt::zeros<double>({sx_, sy_, sz_ + 1});

    a_pu_ = xt::zeros<double>({sx_ + 1, sy_, sz_});
    a_pv_ = xt::zeros<double>({sx_, sy_ + 1, sz_});
    a_pw_ = xt::zeros<double>({sx_, sy_, sz_ + 1});

    au_e_ = xt::zeros<double>({sx_ + 1, sy_, sz_});
    au_w_ = xt::zeros<double>({sx_ + 1, sy_, sz_});
    au_n_ = xt::zeros<double>({sx_ + 1, sy_, sz_});
    au_s_ = xt::zeros<double>({sx_ + 1, sy_, sz_});

    av_e_ = xt::zeros<double>({sx_, sy_ + 1, sz_});
    av_w_ = xt::zeros<double>({sx_, sy_ + 1, sz_});
    av_n_ = xt::zeros<double>({sx_, sy_ + 1, sz_});
    av_s_ = xt::zeros<double>({sx_, sy_ + 1, sz_});
}
