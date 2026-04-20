#include "geometry/Face.hpp"

bool Face::IsPhysicalBoundary() const {
    return kind == FaceKind::PhysicalBoundary;
}

bool Face::IsInternal() const {
    return kind == FaceKind::Interior;
}

bool Face::IsMPIBoundary() const {
    return kind == FaceKind::MPIBoundary;
}
