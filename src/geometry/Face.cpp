#include "geometry/Face.hpp"

bool Face::IsBoundary() const {
    return neighbor_cell_id == k_invalid_cell_id;
}

bool Face::IsInternal() const {
    return !IsBoundary();
}
