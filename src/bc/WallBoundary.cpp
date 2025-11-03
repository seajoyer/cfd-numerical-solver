#include "bc/WallBoundary.hpp"
#include "data/DataLayer.hpp"

void WallBoundary::Apply(DataLayer &layer, int axis, Side side) const {
    (void) layer;
    (void) axis;
    (void) side;
    // TODO: непроницаемая стенка: нормальная скорость = 0, скалярные зеркально/копирование
}