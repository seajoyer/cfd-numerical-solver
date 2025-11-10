#include "solver/TimeStepCalculator.hpp"
#include <cmath>

double TimeStepCalculator::ComputeDt(const std::vector<Primitive> &cellStates,
                                     double dx,
                                     double cfl,
                                     double gamma) {
    double maxSpeed = 0.0;

    for (const Primitive &w: cellStates) {
        if (w.rho <= 0.0 || w.P <= 0.0) {
            continue;
        }

        const double a = EOS::SoundSpeed(w, gamma);
        const double speed = std::fabs(w.u) + a;

        if (speed > maxSpeed) {
            maxSpeed = speed;
        }
    }

    if (maxSpeed <= 0.0) {
        return 0.0;
    }

    return cfl * dx / maxSpeed;
}
