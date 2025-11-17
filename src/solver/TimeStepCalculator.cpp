#include "solver/TimeStepCalculator.hpp"


double TimeStepCalculator::ComputeDt(const DataLayer &layer,
                                     const double dx,
                                     const double cfl,
                                     const double gamma) {
    const int coreStart = layer.GetCoreStart(0);
    const int coreEnd   = layer.GetCoreEndExclusive(0);

    double maxSpeed = 0.0;

    for (int i = coreStart; i < coreEnd; ++i) {
        const double rho = layer.rho(i);
        const double u   = layer.u(i);
        const double P   = layer.P(i);

        if (rho <= 0.0 || P <= 0.0) {
            continue;
        }

        const double c  = std::sqrt(gamma * P / rho);
        const double s  = std::abs(u) + c;
        if (s > maxSpeed) {
            maxSpeed = s;
        }
    }

    if (maxSpeed <= 0.0) {
        return 0.0;
    }

    return cfl * dx / maxSpeed;
}
