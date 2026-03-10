#ifndef DATALAYER_HPP
#define DATALAYER_HPP

#include <cstddef>
#include <stdexcept>
#include <xtensor.hpp>

/**
 * @brief Storage owner for conservative state only.
 *
 * DataLayer owns only conservative variables:
 *   U(var, i, j, k), var in {rho, rhoU, rhoV, rhoW, E}.
 *
 * Geometry, grid metrics, padding, dimension, core ranges, and boundary flags
 * are intentionally not stored here and must live in Mesh.
 *
 * Storage layout:
 *   U(var, i, j, k), shape = (5, sx, sy, sz).
 */
class DataLayer final {
public:
    static constexpr std::size_t k_nvar = 5;
    static constexpr std::size_t k_rho = 0;
    static constexpr std::size_t k_rhoU = 1;
    static constexpr std::size_t k_rhoV = 2;
    static constexpr std::size_t k_rhoW = 3;
    static constexpr std::size_t k_E = 4;

    DataLayer() = default;

    /**
     * @brief Construct conservative storage with given padded sizes.
     * @param sx Total number of cells in x, including ghosts.
     * @param sy Total number of cells in y, including ghosts.
     * @param sz Total number of cells in z, including ghosts.
     */
    DataLayer(int sx, int sy, int sz);

    /**
     * @brief Resize conservative storage to the given padded sizes.
     * @details Reallocates only if shape changed.
     * @param sx Total number of cells in x, including ghosts.
     * @param sy Total number of cells in y, including ghosts.
     * @param sz Total number of cells in z, including ghosts.
     */
    void Resize(int sx, int sy, int sz);

    /** @brief Conservative state array U(var,i,j,k). Shape (5,sx,sy,sz). */
    [[nodiscard]] xt::xtensor<double, 4>& U();
    [[nodiscard]] const xt::xtensor<double, 4>& U() const;

    /** @brief Total padded size in x. */
    [[nodiscard]] int GetSx() const;

    /** @brief Total padded size in y. */
    [[nodiscard]] int GetSy() const;

    /** @brief Total padded size in z. */
    [[nodiscard]] int GetSz() const;

    /** @brief Check whether conservative storage is allocated. */
    [[nodiscard]] bool IsAllocated() const;

private:
    int sx_ = 0;
    int sy_ = 0;
    int sz_ = 0;

    xt::xtensor<double, 4> U_;

    void ValidateSizes(int sx, int sy, int sz) const;
    void Allocate(int sx, int sy, int sz);
};

#endif  // DATALAYER_HPP
