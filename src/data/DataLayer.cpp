#include "data/DataLayer.hpp"

DataLayer::DataLayer(const int sx, const int sy, const int sz) {
    Resize(sx, sy, sz);
}

void DataLayer::Resize(const int sx, const int sy, const int sz) {
    ValidateSizes(sx, sy, sz);

    if (sx == sx_ && sy == sy_ && sz == sz_ && IsAllocated()) {
        return;
    }

    Allocate(sx, sy, sz);
}

xt::xtensor<double, 4>& DataLayer::U() {
    return U_;
}

const xt::xtensor<double, 4>& DataLayer::U() const {
    return U_;
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

bool DataLayer::IsAllocated() const {
    return U_.dimension() == 4 && sx_ > 0 && sy_ > 0 && sz_ > 0;
}

void DataLayer::ValidateSizes(const int sx, const int sy, const int sz) const {
    if (sx <= 0) {
        throw std::invalid_argument("DataLayer: sx must be > 0");
    }
    if (sy <= 0) {
        throw std::invalid_argument("DataLayer: sy must be > 0");
    }
    if (sz <= 0) {
        throw std::invalid_argument("DataLayer: sz must be > 0");
    }
}

void DataLayer::Allocate(const int sx, const int sy, const int sz) {
    sx_ = sx;
    sy_ = sy;
    sz_ = sz;

    U_ = xt::zeros<double>({
        k_nvar,
        static_cast<std::size_t>(sx_),
        static_cast<std::size_t>(sy_),
        static_cast<std::size_t>(sz_)
    });
}
