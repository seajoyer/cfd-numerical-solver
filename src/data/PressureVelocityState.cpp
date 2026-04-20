#include "data/PressureVelocityState.hpp"

#include "data/Mesh.hpp"

void PressureVelocityState::ResizeFrom(const Mesh& mesh) {
    const std::size_t sx = static_cast<std::size_t>(mesh.GetSx());
    const std::size_t sy = static_cast<std::size_t>(mesh.GetSy());
    const std::size_t sz = static_cast<std::size_t>(mesh.GetSz());

    if (sx == sx_ && sy == sy_ && sz == sz_ && IsAllocated()) {
        return;
    }

    Allocate(sx, sy, sz);
}

xt::xtensor<double, 3>& PressureVelocityState::Pressure() {
    return pressure_;
}

const xt::xtensor<double, 3>& PressureVelocityState::Pressure() const {
    return pressure_;
}

xt::xtensor<double, 3>& PressureVelocityState::PressureOld() {
    return pressure_old_;
}

const xt::xtensor<double, 3>& PressureVelocityState::PressureOld() const {
    return pressure_old_;
}

xt::xtensor<double, 3>& PressureVelocityState::Ux() {
    return ux_;
}

const xt::xtensor<double, 3>& PressureVelocityState::Ux() const {
    return ux_;
}

xt::xtensor<double, 3>& PressureVelocityState::Vy() {
    return vy_;
}

const xt::xtensor<double, 3>& PressureVelocityState::Vy() const {
    return vy_;
}

xt::xtensor<double, 3>& PressureVelocityState::Wz() {
    return wz_;
}

const xt::xtensor<double, 3>& PressureVelocityState::Wz() const {
    return wz_;
}

xt::xtensor<double, 3>& PressureVelocityState::UxOld() {
    return ux_old_;
}

const xt::xtensor<double, 3>& PressureVelocityState::UxOld() const {
    return ux_old_;
}

xt::xtensor<double, 3>& PressureVelocityState::VyOld() {
    return vy_old_;
}

const xt::xtensor<double, 3>& PressureVelocityState::VyOld() const {
    return vy_old_;
}

xt::xtensor<double, 3>& PressureVelocityState::WzOld() {
    return wz_old_;
}

const xt::xtensor<double, 3>& PressureVelocityState::WzOld() const {
    return wz_old_;
}

void PressureVelocityState::CopyCurrentToOld() {
    pressure_old_ = pressure_;
    ux_old_ = ux_;
    vy_old_ = vy_;
    wz_old_ = wz_;
}

void PressureVelocityState::ZeroPressure() {
    pressure_.fill(0.0);
}

void PressureVelocityState::ZeroPressureOld() {
    pressure_old_.fill(0.0);
}

void PressureVelocityState::ZeroUx() {
    ux_.fill(0.0);
}

void PressureVelocityState::ZeroVy() {
    vy_.fill(0.0);
}

void PressureVelocityState::ZeroWz() {
    wz_.fill(0.0);
}

void PressureVelocityState::ZeroUxOld() {
    ux_old_.fill(0.0);
}

void PressureVelocityState::ZeroVyOld() {
    vy_old_.fill(0.0);
}

void PressureVelocityState::ZeroWzOld() {
    wz_old_.fill(0.0);
}

void PressureVelocityState::ZeroAll() {
    ZeroPressure();
    ZeroPressureOld();
    ZeroUx();
    ZeroVy();
    ZeroWz();
    ZeroUxOld();
    ZeroVyOld();
    ZeroWzOld();
}

bool PressureVelocityState::IsAllocated() const {
    return
        pressure_.dimension() == 3 &&
        pressure_old_.dimension() == 3 &&
        ux_.dimension() == 3 &&
        vy_.dimension() == 3 &&
        wz_.dimension() == 3 &&
        ux_old_.dimension() == 3 &&
        vy_old_.dimension() == 3 &&
        wz_old_.dimension() == 3 &&
        sx_ > 0 && sy_ > 0 && sz_ > 0;
}

std::size_t PressureVelocityState::GetSx() const {
    return sx_;
}

std::size_t PressureVelocityState::GetSy() const {
    return sy_;
}

std::size_t PressureVelocityState::GetSz() const {
    return sz_;
}

void PressureVelocityState::Allocate(const std::size_t sx,
                                     const std::size_t sy,
                                     const std::size_t sz) {
    sx_ = sx;
    sy_ = sy;
    sz_ = sz;

    pressure_ = xt::zeros<double>({sx_, sy_, sz_});
    pressure_old_ = xt::zeros<double>({sx_, sy_, sz_});

    ux_ = xt::zeros<double>({sx_ + 1, sy_, sz_});
    vy_ = xt::zeros<double>({sx_, sy_ + 1, sz_});
    wz_ = xt::zeros<double>({sx_, sy_, sz_ + 1});

    ux_old_ = xt::zeros<double>({sx_ + 1, sy_, sz_});
    vy_old_ = xt::zeros<double>({sx_, sy_ + 1, sz_});
    wz_old_ = xt::zeros<double>({sx_, sy_, sz_ + 1});
}
