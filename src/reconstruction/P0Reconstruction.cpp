#include "reconstruction/P0Reconstruction.hpp"

void P0Reconstruction::ReconstructStates(const DataLayer& layer,
                                         xt::xarray<Primitive>& left_states,
                                         xt::xarray<Primitive>& right_states) const {
    const int total_size = layer.GetTotalSize();
    const int n_interfaces = total_size > 0 ? total_size - 1 : 0;

    left_states = xt::xarray<Primitive>::from_shape(
        {static_cast<std::size_t>(n_interfaces)});
    right_states = xt::xarray<Primitive>::from_shape(
        {static_cast<std::size_t>(n_interfaces)});

    for (int i = 0; i < n_interfaces; ++i) {
        left_states(i) = layer.GetPrimitive(i);
        right_states(i) = layer.GetPrimitive(i + 1);
    }
}