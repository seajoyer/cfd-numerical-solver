#include "reconstruction/P1Reconstruction.hpp"

#include <algorithm>
#include <cmath>

auto P1Reconstruction::minmod(const double a, const double b) -> double {
    if (a * b <= 0.0) {
        return 0.0;
    }
    return std::abs(a) < std::abs(b) ? a : b;
}

auto P1Reconstruction::mc(const double a, const double b) -> double {
    if (a * b <= 0.0) {
        return 0.0;
    }
    const double sign_a = a > 0.0 ? 1.0 : -1.0;
    const double abs_a = std::abs(a);
    const double abs_b = std::abs(b);

    return sign_a * std::min({2.0 * abs_a,
                              2.0 * abs_b,
                              0.5 * (abs_a + abs_b)});
}

auto P1Reconstruction::superbee(const double a, const double b) -> double {
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

double P1Reconstruction::apply_limiter(const double a, const double b) const {
    switch (limiter_type_) {
        case LimiterType::MC:
            return mc(a, b);
        case LimiterType::SUPERBEE:
            return superbee(a, b);
        case LimiterType::MINMOD:
        default:
            return minmod(a, b);
    }
}

void P1Reconstruction::ComputeInterfaceStates(const DataLayer& layer,
                                              const int interface_index,
                                              Primitive& left_state,
                                              Primitive& right_state) const {
    const int total_size = layer.GetTotalSize();
    const int i = interface_index;
    const int j = interface_index + 1;

    if (i < 0 || j >= total_size) {
        left_state = layer.GetPrimitive(i);
        right_state = layer.GetPrimitive(j);
        return;
    }

    const int i_l = std::max(0, i - 1);
    const int i_r = std::min(total_size - 1, i + 1);

    const int j_l = std::max(0, j - 1);
    const int j_r = std::min(total_size - 1, j + 1);

    const Primitive w_i_l = layer.GetPrimitive(i_l);
    const Primitive w_i = layer.GetPrimitive(i);
    const Primitive w_i_r = layer.GetPrimitive(i_r);

    const Primitive w_j_l = layer.GetPrimitive(j_l);
    const Primitive w_j = layer.GetPrimitive(j);
    const Primitive w_j_r = layer.GetPrimitive(j_r);

    Primitive slope_i, slope_j;

    slope_i.rho = apply_limiter(w_i.rho - w_i_l.rho, w_i_r.rho - w_i.rho);
    slope_i.u = apply_limiter(w_i.u - w_i_l.u, w_i_r.u - w_i.u);
    slope_i.P = apply_limiter(w_i.P - w_i_l.P, w_i_r.P - w_i.P);

    slope_j.rho = apply_limiter(w_j.rho - w_j_l.rho, w_j_r.rho - w_j.rho);
    slope_j.u = apply_limiter(w_j.u - w_j_l.u, w_j_r.u - w_j.u);
    slope_j.P = apply_limiter(w_j.P - w_j_l.P, w_j_r.P - w_j.P);

    left_state = w_i + 0.5 * slope_i;
    right_state = w_j - 0.5 * slope_j;
}