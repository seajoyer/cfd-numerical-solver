#ifndef P1RECONSTRUCTION_HPP
#define P1RECONSTRUCTION_HPP

#include "reconstruction/Reconstruction.hpp"

/**
 * @enum LimiterType
 * @brief Types of slope limiters for piecewise-linear reconstruction.
 */
enum class LimiterType : std::uint8_t { kMinmod, kMc, kSuperbee };

/**
 * @class P1Reconstruction
 * @brief Piecewise-linear reconstruction with slope limiters.
 *
 * This scheme computes limited linear slopes for each primitive component
 * in each cell using a chosen limiter (minmod, MC, or Superbee) applied
 * to backward and forward differences. Interface states are then built as
 *
 *  w_L(i+1/2) = w_i   + 0.5 * slope_i
 *  w_R(i+1/2) = w_{i+1} - 0.5 * slope_{i+1}
 *
 * where ghost cells provide the necessary neighbors near boundaries.
 */
class P1Reconstruction : public Reconstruction {
public:
    P1Reconstruction() = default;
    ~P1Reconstruction() override = default;

    /**
     * @brief Computes left/right primitive states at all interfaces.
     *
     * For a 1D layout with total_size cells (including ghosts), there are
     * total_size - 1 interfaces. Slopes are first computed for every cell
     * using a chosen limiter, then used to reconstruct interface states.
     *
     * @param layer        DataLayer with primitive fields.
     * @param left_states  Output: left primitive states at all interfaces.
     * @param right_states Output: right primitive states at all interfaces.
     */
    void ReconstructStates(const DataLayer& layer,
                           xt::xarray<Primitive>& left_states,
                           xt::xarray<Primitive>& right_states) const override;

    /**
     * @brief Sets the slope limiter type.
     *
     * @param type LimiterType to use (kMinmod, kMc, kSuperbee).
     */
    void SetLimiter(LimiterType type) { limiter_type_ = type; }

private:
    LimiterType limiter_type_{LimiterType::kMinmod};

    /**
     * @brief Minmod limiter for two arguments.
     *
     * Returns 0 if a and b have opposite signs or either is zero,
     * otherwise returns the argument with smaller magnitude.
     *
     * @param a First slope candidate.
     * @param b Second slope candidate.
     * @return Limited slope.
     */
    static auto Minmod(double a, double b) -> double;

    /**
     * @brief Monotonized central (MC) limiter for two arguments.
     *
     * @param a First slope candidate.
     * @param b Second slope candidate.
     * @return Limited slope.
     */
    static auto Mc(double a, double b) -> double;

    /**
     * @brief Superbee limiter for two arguments.
     *
     * @param a First slope candidate.
     * @param b Second slope candidate.
     * @return Limited slope.
     */
    static auto Superbee(double a, double b) -> double;

    /**
     * @brief Applies the currently selected limiter to two slope candidates.
     *
     * @param a First slope candidate.
     * @param b Second slope candidate.
     * @return Limited slope.
     */
    [[nodiscard]] auto ApplyLimiter(double a, double b) const -> double;
};

#endif  // P1RECONSTRUCTION_HPP