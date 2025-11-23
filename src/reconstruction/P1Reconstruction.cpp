#include "reconstruction/P1Reconstruction.hpp"

#include <algorithm>
#include <cmath>

auto P1Reconstruction::Minmod(const double a, const double b) -> double {
    if (a * b <= 0.0) {
        return 0.0;
    }
    return std::abs(a) < std::abs(b) ? a : b;
}

auto P1Reconstruction::Mc(const double a, const double b) -> double {
    if (a * b <= 0.0) {
        return 0.0;
    }
    const double sign_a = a > 0.0 ? 1.0 : -1.0;
    const double abs_a = std::abs(a);
    const double abs_b = std::abs(b);

    return sign_a * std::min({2.0 * abs_a, 2.0 * abs_b, 0.5 * (abs_a + abs_b)});
}

auto P1Reconstruction::Superbee(const double a, const double b) -> double {
    if (a * b <= 0.0) {
        return 0.0;
    }
    const double sign_a = a > 0.0 ? 1.0 : -1.0;
    const double abs_a = std::abs(a);
    const double abs_b = std::abs(b);

    const double term1 = std::min(2.0 * abs_a, abs_b);
    const double term2 = std::min(abs_a, 2.0 * abs_b);

    return sign_a * std::max(term1, term2);
}

auto P1Reconstruction::ApplyLimiter(const double a, const double b) const -> double {
    switch (limiter_type_) {
        case LimiterType::kMc:
            return Mc(a, b);
        case LimiterType::kSuperbee:
            return Superbee(a, b);
        case LimiterType::kMinmod:
        default:
            return Minmod(a, b);
    }
}

void P1Reconstruction::ReconstructStates(const DataLayer& layer,
                                         xt::xarray<Primitive>& left_states,
                                         xt::xarray<Primitive>& right_states) const {
    const int total_size = layer.GetTotalSize();
    const int n_interfaces = total_size > 0 ? total_size - 1 : 0;

    left_states = xt::xarray<Primitive>::from_shape(
        {static_cast<std::size_t>(n_interfaces)});
    right_states = xt::xarray<Primitive>::from_shape(
        {static_cast<std::size_t>(n_interfaces)});

    if (total_size <= 1) {
        return;
    }

    xt::xarray<Primitive> slopes =
        xt::xarray<Primitive>::from_shape({static_cast<std::size_t>(total_size)});

    for (int k = 0; k < total_size; ++k) {
        const int k_l = std::max(0, k - 1);
        const int k_r = std::min(total_size - 1, k + 1);

        const Primitive w_l = layer.GetPrimitive(k_l);
        const Primitive w_c = layer.GetPrimitive(k);
        const Primitive w_r = layer.GetPrimitive(k_r);

        slopes(k).rho = ApplyLimiter(w_c.rho - w_l.rho, w_r.rho - w_c.rho);
        slopes(k).u = ApplyLimiter(w_c.u - w_l.u, w_r.u - w_c.u);
        slopes(k).P = ApplyLimiter(w_c.P - w_l.P, w_r.P - w_c.P);
    }

    for (int i = 0; i < n_interfaces; ++i) {
        const Primitive w_i = layer.GetPrimitive(i);
        const Primitive w_ip1 = layer.GetPrimitive(i + 1);

        left_states(i) = w_i + 0.5 * slopes(i);
        right_states(i) = w_ip1 - 0.5 * slopes(i + 1);
    }
}