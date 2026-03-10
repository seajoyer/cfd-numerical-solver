#ifndef SOLUTIONFILTER_HPP
#define SOLUTIONFILTER_HPP

#include <cstddef>

#include "config/Settings.hpp"
#include "data/Variables.hpp"

/**
 * @class SolutionFilter
 * @brief Explicit diffusion / anti-diffusion filter for conservative Euler state.
 *
 * Applies a second-difference smoothing and optional anti-diffusion (sharpening)
 * to selected primitive scalars (rho and P), then maps back to conservative U.
 *
 * Notes:
 *  - Operates on fluid core cells only.
 *  - Uses axis-split filtering (X/Y/Z) depending on mesh dimension.
 *  - Keeps velocities (u,v,w) unchanged; updates E consistently with filtered P and rho.
 */
class DataLayer;
class Mesh;

class SolutionFilter final {
public:
    explicit SolutionFilter(const Settings& settings);

    /** @brief Apply the filter in-place to DataLayer::U() on fluid core cells only. */
    void Apply(DataLayer& layer, const Mesh& mesh, double gamma) const;

private:
    double eps_diff_ = 0.05;
    double eps_anti_ = 0.03;
    double rho_floor_ = 1e-10;
    double p_floor_ = 1e-10;
    bool enable_antidiffusion_ = true;

    mutable xt::xtensor<double, 3> rho0_;
    mutable xt::xtensor<double, 3> p0_;
    mutable xt::xtensor<double, 3> rho1_;
    mutable xt::xtensor<double, 3> p1_;
    mutable int sx_ = 0;
    mutable int sy_ = 0;
    mutable int sz_ = 0;

    void ResizeFrom(const Mesh& mesh) const;

    void ExtractRhoP(const DataLayer& layer, const Mesh& mesh, double gamma) const;
    void WriteBackConservative(DataLayer& layer, const Mesh& mesh, double gamma) const;

    void ApplyAxis(const Mesh& mesh, Axis axis) const;
};

#endif  // SOLUTIONFILTER_HPP