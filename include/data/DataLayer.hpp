#ifndef DATALAYER_HPP
#define DATALAYER_HPP

#include <xtensor/containers/xarray.hpp>
#include <xtensor/io/xio.hpp>

#include "data/Variables.hpp"

/**
 * @struct DataLayer
 * @brief Storage container for all physical and computational fields using xtensor.
 *
 * This class contains all primary and auxiliary fields used in numerical
 * schemes (density, velocity, pressure, energy, etc.) along with boundary
 * and cell center coordinates. It manages memory allocation and provides
 * convenient accessors for grid size and ghost cell information.
 *
 * Fields stored:
 * - rho: density
 * - u: velocity (x-component)
 * - P: pressure (capital P for thermodynamic pressure)
 * - p: momentum
 * - e: total energy density
 * - U: internal energy
 * - V: specific volume
 * - m: mass
 * - xb: cell boundary coordinates
 * - xc: cell center coordinates
 */
struct DataLayer {
    xt::xarray<double> rho;  ///< Density
    xt::xarray<double> u;    ///< Velocity
    xt::xarray<double> P;    ///< Thermodynamic pressure
    xt::xarray<double> p;    ///< Momentum
    xt::xarray<double> e;    ///< Total energy density
    xt::xarray<double> U;    ///< Specific internal energy
    xt::xarray<double> V;    ///< Specific volume
    xt::xarray<double> m;    ///< Cell mass
    xt::xarray<double> xb;   ///< Cell boundary coordinates
    xt::xarray<double> xc;   ///< Cell center coordinates

    DataLayer() = default;

    /**
     * @brief Constructs containers for all variables with ghost cells (1D).
     *
     * @param N Number of core (physical) cells in the computational domain.
     * @param padding Number of ghost cells on each boundary side.
     */
    DataLayer(int N, int padding);

    /**
     * @brief Constructs containers for a given dimension.
     *
     * @param N Number of core cells per direction.
     * @param padding Number of ghost cells per direction.
     * @param dim Dimension of the problem (1, 2, or 3).
     */
    DataLayer(int N, int padding, int dim);

    ~DataLayer() = default;

    /**
     * @brief Returns starting index of the core (physical) domain.
     * @param axis Spatial axis (0=x, 1=y, 2=z; currently ignored for 1D)
     * @return Index of first physical cell
     */
    [[nodiscard]] auto GetCoreStart(int axis = 0) const -> int;

    /**
     * @brief Returns one-past-last index of the core domain.
     * @param axis Spatial axis (0=x, 1=y, 2=z; currently ignored for 1D)
     * @return End index (exclusive)
     */
    [[nodiscard]] auto GetCoreEndExclusive(int axis = 0) const -> int;

    /**
     * @brief Returns the number of physical cells (without ghost zones).
     * @return Number of cells (N)
     */
    [[nodiscard]] auto GetN() const -> int { return n_; }

    /**
     * @brief Returns the number of ghost cells on each side.
     * @return Number of ghost cells (padding)
     */
    [[nodiscard]] auto GetNGhostCells() const -> int { return n_ghost_cells_; }

    /**
     * @brief Returns the number of ghost cells on each side (alias).
     * @return Number of ghost cells (padding)
     */
    [[nodiscard]] auto GetPadding() const -> int { return n_ghost_cells_; }

    /**
     * @brief Returns the current problem dimension.
     * @return Dimension of the grid (1, 2, or 3)
     */
    [[nodiscard]] auto GetDimension() const -> int { return dimension_; }

    /**
     * @brief Returns the current problem dimension (alias).
     * @return Dimension of the grid (1, 2, or 3)
     */
    [[nodiscard]] auto GetDim() const -> int { return dimension_; }

    /**
     * @brief Returns total size including ghost cells.
     * @return Total array size
     */
    [[nodiscard]] auto GetTotalSize() const -> int { return total_size_; }

    /**
     * @brief Changes number of cells and reallocates all internal arrays.
     * @param new_N New number of physical cells
     */
    void SetN(int new_N);

    /**
     * @brief Changes number of ghost cells and reallocates arrays.
     * @param new_padding New ghost cell count
     */
    void SetPadding(int new_padding);

    /**
     * @brief Sets spatial dimension and reallocates arrays.
     * @param new_dim Dimension (1, 2, or 3)
     */
    void SetDim(int new_dim);

    /**
     * @brief Returns primitive state (rho, u, P) in a given 1D cell.
     *
     * This is a convenience accessor that bundles the SoA storage
     * (rho(i), u(i), P(i)) into a single Primitive structure.
     * Only 1D index is currently supported; 2D/3D overloads can be
     * added later without changing existing code.
     *
     * @param i Cell index in the 1D layout (including ghost cells).
     * @return Primitive state at cell i.
     */
    [[nodiscard]] auto GetPrimitive(int i) const -> Primitive;

    /**
     * @brief Sets primitive state (rho, u, P) in a given 1D cell.
     *
     * Writes the components of the provided Primitive structure into
     * the SoA arrays (rho, u, P) at index i. Auxiliary quantities
     * (e, U, V, m) are left unchanged by this method.
     *
     * @param i Cell index in the 1D layout (including ghost cells).
     * @param w Primitive state to store.
     */
    void SetPrimitive(int i, const Primitive& w);

   private:
    int n_ = 0;              ///< Number of core cells
    int n_ghost_cells_ = 0;  ///< Number of ghost cells on each side
    int dimension_ = 1;      ///< Spatial dimension (1, 2, or 3)
    int total_size_ = 0;     ///< Total array size (core + ghost)

    /**
     * @brief Recomputes total sizes after n or padding are changed.
     */
    void RecomputeSizes();

    /**
     * @brief Allocates 1D arrays for all stored fields using xtensor.
     *
     * Called internally by constructors and setters.
     */
    void Allocate1D();
};

#endif  // DATALAYER_HPP
