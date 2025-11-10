#include "reconstruction/P0Reconstruction.hpp"

void P0Reconstruction::Reconstruct(const std::vector<Primitive> &cellStates,
                                   std::vector<Primitive> &leftStates,
                                   std::vector<Primitive> &rightStates) const {
    const std::size_t nInterfaces = leftStates.size();

    for (std::size_t i = 0; i < nInterfaces; ++i) {
        leftStates[i] = cellStates[i];
        rightStates[i] = cellStates[i + 1];
    }
}
