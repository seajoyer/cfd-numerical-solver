#ifndef GLOBALLIMITER_HPP
#define GLOBALLIMITER_HPP

#include "config/Settings.hpp"
#include "data/DataLayer.hpp"

/**
 * @class GlobalLimiter
 * @brief Simple post-step global limiter for suppressing non-physical oscillations.
 *
 * This limiter applies a very simple local-extremum smoothing to the total energy
 * density field e. For each core cell i (excluding immediate neighbors), if e(i)
 * is a local extremum, it is relaxed toward the average of its two neighbors:
 *
 *   e_i <- e_i + strength * (0.5*(e_{i-1}+e_{i+1}) - e_i).
 *
 * After updating e, the limiter recomputes pressure P and specific internal
 * energy U consistently via EOS and enforces positivity via PositivityLimiter.
 *
 * @note This limiter is intentionally primitive and is intended as a "global crutch".
 */
class GlobalLimiter {
public:
    /**
     * @brief Constructs limiter from settings.
     */
    explicit GlobalLimiter() = default;

    /**
     * @brief Applies limiter in-place to DataLayer.
     *
     * @param layer    Data layer to be modified.
     * @param dx       Cell size.
     * @param settings Global simulation settings (gamma etc.).
     */
    void Apply(DataLayer& layer, double dx, const Settings& settings) const;

private:
    double strength_{1.0};
};

#endif  // GLOBALLIMITER_HPP
