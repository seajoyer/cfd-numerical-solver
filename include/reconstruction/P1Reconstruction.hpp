#ifndef P1RECONSTRUCTION_HPP
#define P1RECONSTRUCTION_HPP

#include "reconstruction/Reconstruction.hpp"

enum class LimiterType : std::uint8_t { kMinmod, kMc, kSuperbee };

/**
 * @class P1Reconstruction
 * @brief Piecewise-linear reconstruction with a minmod slope limiter.
 *
 * This scheme computes limited linear slopes for each primitive component
 * in each cell using the minmod function applied to backward and forward
 * differences. Interface states are then built as
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
     * @brief Computes left/right primitive states at interface i.
     *
     * For interface i between cell i and cell i+1, this method:
     *  - computes limited slopes for both cell i and cell i+1 using
     *    minmod(backward difference, forward difference),
     *  - reconstructs left/right states at the interface.
     *
     * Cell indices are clamped to [0, totalSize-1] to safely work near
     * the domain boundaries using ghost cells.
     *
     * @param layer          DataLayer with primitive fields.
     * @param interface_index Interface index (0 ≤ i < totalSize - 1).
     * @param left_state      Output: left primitive state at the interface.
     * @param right_state     Output: right primitive state at the interface.
     */
    void Extracted() const;
    void ComputeInterfaceStates(const DataLayer& layer, int interface_index,
                                Primitive& left_state,
                                Primitive& right_state) const override;

    void SetLimiter(LimiterType type) { limiter_type_ = type; }

   private:
    LimiterType limiter_type_{LimiterType::kMinmod};

    /**
     * @brief Minmod limiter for two arguments.
     *
     * Returns:
     *  - 0 if a and b have opposite signs or either is zero,
     *  - the argument with smaller magnitude otherwise.
     *
     * @param a First slope candidate.
     * @param b Second slope candidate.
     * @return Limited slope.
     */
    static auto Minmod(double a, double b) -> double;
    static auto Mc(double a, double b) -> double;
    static auto Superbee(double a, double b) -> double;

    [[nodiscard]] auto ApplyLimiter(double a, double b) const -> double;
};

#endif  // P1RECONSTRUCTION_HPP
