#ifndef DATALAYER_HPP
#define DATALAYER_HPP

#include <xtensor/io/xio.hpp>
#include <xtensor/views/xview.hpp>
#include <xtensor/containers/xarray.hpp>

/**
 * @class DataLayer
 * @brief Stores physical variables and mesh-related data for a single time layer.
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
 * - e: specific internal energy
 * - U: conserved energy
 * - V: volume
 * - m: mass
 * - xb: cell boundary coordinates
 * - xc: cell center coordinates
 */
struct DataLayer {
    xt::xarray<double> rho;  // Density
    xt::xarray<double> u;    // Velocity
    xt::xarray<double> P;    // Pressure
    xt::xarray<double> p;    // Momentum
    xt::xarray<double> e;    // Specific internal energy
    xt::xarray<double> U;    // Conserved energy
    xt::xarray<double> V;    // Volume
    xt::xarray<double> m;    // Mass
    xt::xarray<double> xb;   // Cell boundary coordinates
    xt::xarray<double> xc;   // Cell center coordinates

    DataLayer() = default;

    /**
     * @brief Constructs containers for all variables with ghost cells.
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
     * @param axis Spatial axis (0=x, 1=y, 2=z; ignored for 1D)
     * @return Index of first physical cell
     */
    [[nodiscard]] auto GetCoreStart(int axis = 0) const -> int;
    
    /**
     * @brief Returns one-past-last index of the core domain.
     * @param axis Spatial axis (0=x, 1=y, 2=z; ignored for 1D)
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
    [[nodiscard]] auto GetPadding() const -> int { return n_ghost_cells_; }

    /**
     * @brief Returns the current problem dimension.
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

private:
    int n_ = 0;
    int n_ghost_cells_ = 0;
    int dimension_ = 1;
    int total_size_ = 0;

    /**
     * @brief Recomputes total sizes after n or padding are changed.
     */
    void RecomputeSizes();

    /**
     * @brief Allocates 1D arrays for all stored fields.
     *
     * Called internally by constructors and setters.
     */
    void Allocate1D();
};

#endif // DATALAYER_HPP
