#ifndef DATALAYER_HPP
#define DATALAYER_HPP

#include <xtensor/io/xio.hpp>
#include <xtensor/views/xview.hpp>
#include <xtensor/containers/xarray.hpp>

/**
 * @class DataLayer
 * @brief Stores physical variables and mesh-related data for a single time layer.
 *
 * This class contains all primary and auxiliary fields used in 1D numerical
 * schemes (density, velocity, pressure, energy, etc.) along with boundary
 * and cell center coordinates. It also manages memory allocation and provides
 * convenient accessors for grid size and ghost cell information.
 */
struct DataLayer {
    xt::xarray<double> rho;
    xt::xarray<double> u;
    xt::xarray<double> P;
    xt::xarray<double> p;
    xt::xarray<double> e;
    xt::xarray<double> U;
    xt::xarray<double> V;
    xt::xarray<double> m;
    xt::xarray<double> xb;
    xt::xarray<double> xc;

    DataLayer() = default;

    /**
     * @brief Constructs containers for all inner variables in 1D.
     *
     * @param n Number of core (physical) cells in the computational domain.
     * @param padding Number of ghost cells on each boundary side.
     */
    DataLayer(int N, int padding);

    /**
     * @brief Constructs containers for a given dimension (only 1D currently supported).
     *
     * @param n Number of core cells per direction.
     * @param padding Number of ghost cells per direction.
     * @param dim Dimension of the problem (1 for now).
     */
    DataLayer(int N, int padding, int dim);
    ~DataLayer() = default;

    /**
     * @brief Returns starting index of the core (physical) domain.
     * @param axis Spatial axis (ignored for 1D)
     * @return Index of first physical cell
     */
    int GetCoreStart(const int axis = 0) const { (void)axis; return nGhostCells; }
    /**
     * @brief Returns one-past-last index of the core domain.
     * @param axis Spatial axis (ignored for 1D)
     * @return End index (exclusive)
     */
    int GetCoreEndExclusive(const int axis = 0) const { (void)axis; return nGhostCells + n; }

    /**
     * @brief Returns the number of physical cells (without ghost zones).
     * @return Number of cells (N)
     */
    int GetN() const { return n; }

    /**
     * @brief Returns the number of ghost cells on each side.
     * @return Number of ghost cells (padding)
     */
    int GetPadding() const { return nGhostCells; }

    /**
     * @brief Returns the current problem dimension (1D by default).
     * @return Dimension of the grid
     */
    int GetDim() const { return dimention; }

    /**
     * @brief Changes number of cells and reallocates all internal arrays.
     * @param newN New number of physical cells
     */
    void SetN(int newN);

    /**
     * @brief Changes number of ghost cells and reallocates arrays.
     * @param newPadding New ghost cell count
     */
    void SetPadding(int newPadding);

    /**
     * @brief Sets spatial dimension (currently only dim=1 supported).
     * @param newDim Dimension (1)
     */
    void SetDim(int newDim);

private:
    int n;
    int nGhostCells;
    int dimention = 1;
    int totalSize = 0;

    /**
     * @brief Recomputes total sizes after n or padding are changed.
     */
    void RecomputeSizes();

    /**
     * @brief Allocates 1D arrays for all stored fields.
     *
     * Called internally by constructors and ResizeAll().
     */
    void Allocate1D();
};

#endif //DATALAYER_HPP
