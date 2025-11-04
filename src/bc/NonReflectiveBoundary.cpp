#include "bc/NonReflectiveBoundary.hpp"
#include "data/DataLayer.hpp"

void NonReflectiveBoundary::Apply(DataLayer &layer, int axis, Side side) const {
    (void) layer;
    (void) axis;
    (void) side;
    // TODO: непроницаемая стенка: нормальная скорость = 0, скалярные зеркально/копирование
}
