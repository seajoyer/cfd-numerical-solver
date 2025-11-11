#include "solver/TimeStepCalculator.hpp"
#include "solver/EOS.hpp"
#include <cmath>

auto TimeStepCalculator::ComputeDt(const std::vector<Primitive> &cell_states,
                                     double dx,
                                     double cfl,
                                     double gamma) -> double {
    double max_speed = 0.0;

    for (const Primitive &w: cell_states) {
        if (w.rho <= 0.0 || w.P <= 0.0) {
            continue;
        }

        const double a = EOS::SoundSpeed(w, gamma);
        const double speed = std::fabs(w.u) + a;

        if (speed > max_speed) {
            max_speed = speed;
        }
    }

    if (max_speed <= 0.0) {
        return 0.0;
    }

    return cfl * dx / max_speed;
}
