#ifndef SOLUTIONFILTER_HPP
#define SOLUTIONFILTER_HPP

#include "data/DataLayer.hpp"
#include "config/Settings.hpp"

/**
 * @class SolutionFilter
 * @brief Explicit diffusion / anti-diffusion filter for 1D solutions.
 *
 * Implements a simple second-difference based filter:
 *  - diffusion (smoothing)
 *  - optional anti-diffusion (partial sharpening)
 *
 * Operates directly on DataLayer primitive variables.
 */
class SolutionFilter {
public:
    SolutionFilter() = default;

    /**
     * @brief Applies diffusion + anti-diffusion filter.
     *
     * @param layer    DataLayer to be filtered.
     */
    void Apply(DataLayer& layer) const;
};

#endif  // SOLUTIONFILTER_HPP
