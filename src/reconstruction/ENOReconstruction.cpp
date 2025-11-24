#include "reconstruction/ENOReconstruction.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace {
double DividedDifference(const xt::xarray<double>& f,
                         const std::vector<int>& idx) {
    const std::size_t m = idx.size();
    if (m == 0) {
        return 0.0;
    }
    if (m == 1) {
        return f(idx[0]);
    }

    std::vector<double> g(m);
    for (std::size_t i = 0; i < m; ++i) {
        g[i] = f(idx[i]);
    }

    for (std::size_t k = 1; k < m; ++k) {
        for (std::size_t i = 0; i < m - k; ++i) {
            const double xj = idx[i];
            const double xjk = idx[i + k];
            const double denom = xjk - xj;
            if (denom != 0.0) {
                g[i] = (g[i + 1] - g[i]) / denom;
            } else {
                g[i] = 0.0;
            }
        }
    }
    return g[0];
}

void BuildENOStencil(const xt::xarray<double>& f,
                     int base_index,
                     int order,
                     int N,
                     std::vector<int>& stencil_out) {
    stencil_out.clear();
    if (N <= 0 || base_index < 0 || base_index >= N || order <= 0) {
        if (base_index >= 0 && base_index < N) {
            stencil_out.push_back(base_index);
        }
        return;
    }

    stencil_out.push_back(base_index);
    int left = base_index;
    int right = base_index;

    const int max_order = std::min(order, N);

    for (int level = 1; level < max_order; ++level) {
        const int cand_left = left - 1;
        const int cand_right = right + 1;

        const bool has_left = cand_left >= 0;
        const bool has_right = cand_right < N;

        if (!has_left && !has_right) {
            break;
        }

        int chosen_side = 0;

        if (has_left && !has_right) {
            chosen_side = -1;
        } else if (!has_left && has_right) {
            chosen_side = +1;
        } else {
            std::vector<int> stencil_left;
            stencil_left.reserve(level + 1);
            stencil_left.push_back(cand_left);
            stencil_left.insert(stencil_left.end(), stencil_out.begin(),
                                stencil_out.end());
            const double dd_left = std::abs(DividedDifference(f, stencil_left));

            std::vector<int> stencil_right = stencil_out;
            stencil_right.push_back(cand_right);
            const double dd_right = std::abs(DividedDifference(f, stencil_right));

            if (dd_left < dd_right) {
                chosen_side = -1;
            } else if (dd_right < dd_left) {
                chosen_side = +1;
            } else {
                const double center = base_index;
                const double dl = std::abs(static_cast<double>(cand_left) - center);
                const double dr = std::abs(static_cast<double>(cand_right) - center);
                chosen_side = dl <= dr ? -1 : +1;
            }
        }

        if (chosen_side < 0) {
            left = cand_left;
            stencil_out.insert(stencil_out.begin(), left);
        } else {
            right = cand_right;
            stencil_out.push_back(right);
        }
    }
}

double LagrangeInterpolate(const xt::xarray<double>& f,
                           const std::vector<int>& stencil,
                           double x_eval) {
    const std::size_t m = stencil.size();
    if (m == 0) {
        return 0.0;
    }
    if (m == 1) {
        return f(stencil[0]);
    }

    double result = 0.0;
    for (std::size_t k = 0; k < m; ++k) {
        const double xk = stencil[k];
        double Lk = 1.0;
        for (std::size_t j = 0; j < m; ++j) {
            if (j == k) {
                continue;
            }
            const double xj = stencil[j];
            const double denom = xk - xj;
            if (denom != 0.0) {
                Lk *= (x_eval - xj) / denom;
            }
        }
        result += f(stencil[k]) * Lk;
    }
    return result;
}

double ENOReconstructScalar(const xt::xarray<double>& f,
                            int base_index,
                            double x_eval,
                            int order,
                            int N) {
    std::vector<int> stencil;
    BuildENOStencil(f, base_index, order, N, stencil);
    return LagrangeInterpolate(f, stencil, x_eval);
}

Primitive ENOReconstructPrimitive(const DataLayer& layer,
                                  int base_index,
                                  double x_eval,
                                  int order) {
    const int N = layer.GetTotalSize();

    Primitive w;
    w.rho = ENOReconstructScalar(layer.rho, base_index, x_eval, order, N);
    w.u = ENOReconstructScalar(layer.u, base_index, x_eval, order, N);
    w.P = ENOReconstructScalar(layer.P, base_index, x_eval, order, N);
    return w;
}
} // namespace

ENOReconstruction::ENOReconstruction(const int order) : order_(order) {
}

void ENOReconstruction::SetOrder(const int order) {
    order_ = order;
}

void ENOReconstruction::ReconstructStates(const DataLayer& layer,
                                          xt::xarray<Primitive>& left_states,
                                          xt::xarray<Primitive>& right_states) const {
    const int total_size = layer.GetTotalSize();
    const int n_interfaces = std::max(0, total_size - 1);
    const std::size_t n_int = static_cast<std::size_t>(n_interfaces);

    left_states = xt::xarray<Primitive>::from_shape({n_int});
    right_states = xt::xarray<Primitive>::from_shape({n_int});

    if (n_interfaces <= 0) {
        return;
    }

    for (int i = 0; i < n_interfaces; ++i) {
        const double x_half = static_cast<double>(i) + 0.5;

        left_states(i) = ENOReconstructPrimitive(layer, i, x_half, order_);
        right_states(i) = ENOReconstructPrimitive(layer, i + 1, x_half, order_);
    }
}