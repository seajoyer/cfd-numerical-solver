#ifndef P1RECONSTRUCTION_HPP
#define P1RECONSTRUCTION_HPP

#include <cstdint>

#include "reconstruction/Reconstruction.hpp"

/**
 * @enum LimiterType
 * @brief Types of slope limiters for piecewise-linear reconstruction.
 */
enum class LimiterType : std::uint8_t { kMinmod = 0, kMc = 1, kSuperbee = 2 };

/**
 * @class P1Reconstruction
 * @brief Piecewise-linear reconstruction with slope limiters (MUSCL-type).
 *
 * For a face between left cell (i,j,k) and right cell (+1 along axis):
 *   WL = W_L + 0.5 * slope(L)
 *   WR = W_R - 0.5 * slope(R)
 *
 * Slopes are computed component-wise using a limiter applied to backward and forward differences:
 *   slope(i) = limiter( W(i) - W(i-1), W(i+1) - W(i) )
 *
 * Contract:
 *  - No allocations in hot path.
 *  - Uses only local stencil and ghost cells for boundary support.
 */
class P1Reconstruction final : public Reconstruction {
public:
    P1Reconstruction() = default;
    ~P1Reconstruction() override = default;

    void ReconstructFace(const xt::xtensor<double, 4>& W,
                         Axis axis,
                         int i, int j, int k,
                         PrimitiveCell& WL,
                         PrimitiveCell& WR) const override;

    void SetLimiter(LimiterType type) { limiter_type_ = type; }

private:
    LimiterType limiter_type_{LimiterType::kMinmod};

    static auto Minmod(double a, double b) -> double;
    static auto Mc(double a, double b) -> double;
    static auto Superbee(double a, double b) -> double;

    [[nodiscard]] auto ApplyLimiter(double a, double b) const -> double;
};

#endif  // P1RECONSTRUCTION_HPP