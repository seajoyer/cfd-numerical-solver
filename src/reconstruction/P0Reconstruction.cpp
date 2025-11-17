#include "reconstruction/P0Reconstruction.hpp"

void P0Reconstruction::ComputeInterfaceStates(const DataLayer& layer,
                                              const int interface_index,
                                              Primitive& left_state,
                                              Primitive& right_state) const {
    left_state = layer.GetPrimitive(interface_index);
    right_state = layer.GetPrimitive(interface_index + 1);
}