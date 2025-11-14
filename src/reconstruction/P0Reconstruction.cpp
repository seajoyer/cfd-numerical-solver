#include "reconstruction/P0Reconstruction.hpp"

void P0Reconstruction::Reconstruct(const std::vector<Primitive> &cell_states,
                                   std::vector<Primitive> &left_states,
                                   std::vector<Primitive> &right_states) const {
    const std::size_t n_interfaces = left_states.size();

    for (std::size_t i = 0; i < n_interfaces; ++i) {
        left_states[i] = cell_states[i];
        right_states[i] = cell_states[i + 1];
    }
}
