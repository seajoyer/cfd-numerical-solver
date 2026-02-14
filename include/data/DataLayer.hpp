#ifndef DATALAYER_HPP
#define DATALAYER_HPP

#include <xtensor.hpp>
#include <xtensor/io/xio.hpp>

#include "data/Variables.hpp"

/**
 * @file DataLayer.hpp
 * @brief Data storage container for all simulation fields (1D and 2D)
 *
 * ## Storage Layout (1D)
 * For N physical cells with padding P:
 *   [ghost: P][core: N][ghost: P]  =>  total = N + 2*P
 *
 * ## Storage Layout (2D)
 * For Nx x Ny physical cells with padding P:
 *   Shape: (Nx + 2*P) x (Ny + 2*P)
 *   All field arrays are 2D xtensor arrays.
 *   Access: field(i, j) where i is x-index, j is y-index.
 *
 *   Ghost regions surround the core on all four sides.
 */
struct DataLayer {
    // ==================== Field Arrays ====================
    
    xt::xarray<double> rho;  ///< Density field
    xt::xarray<double> u;    ///< Velocity field in x-direction
    xt::xarray<double> v;    ///< Velocity field in y-direction (2D only; size 0 in 1D)
    xt::xarray<double> P;    ///< Thermodynamic pressure field
    xt::xarray<double> p;    ///< x-Momentum density (rho*u)
    xt::xarray<double> q;    ///< y-Momentum density (rho*v, 2D only)
    xt::xarray<double> e;    ///< Total energy density
    xt::xarray<double> U;    ///< Specific internal energy
    xt::xarray<double> V;    ///< Specific volume (1/rho)
    xt::xarray<double> m;    ///< Cell mass
    
    xt::xarray<double> xb;   ///< Cell boundary x-coordinates
    xt::xarray<double> xc;   ///< Cell center x-coordinates
    xt::xarray<double> yb;   ///< Cell boundary y-coordinates (2D only)
    xt::xarray<double> yc;   ///< Cell center y-coordinates (2D only)

    // ==================== Constructors ====================
    
    DataLayer() = default;

    /**
     * @brief Constructs 1D data layer with specified resolution
     */
    DataLayer(int N, int padding);

    /**
     * @brief Constructs data layer with specified dimension
     * For dim=1: allocates 1D arrays of size N+2*padding
     * For dim=2: allocates 2D arrays of size (N+2*padding) x (N+2*padding)
     */
    DataLayer(int N, int padding, int dim);

    /**
     * @brief Constructs 2D data layer with separate Nx, Ny
     */
    DataLayer(int Nx, int Ny, int padding, int dim);

    ~DataLayer() = default;

    // ==================== Grid Information ====================
    
    /**
     * @brief Returns starting index of core domain along given axis
     * @param axis 0=x, 1=y
     */
    [[nodiscard]] auto GetCoreStart(int axis = 0) const -> int;

    /**
     * @brief Returns one-past-last index of core domain (exclusive end)
     * @param axis 0=x, 1=y
     */
    [[nodiscard]] auto GetCoreEndExclusive(int axis = 0) const -> int;

    [[nodiscard]] auto GetNx() const -> int { return nx_; }
    [[nodiscard]] auto GetNy() const -> int { return ny_; }
    [[nodiscard]] auto GetN() const -> int { return nx_; }  // 1D compat
    [[nodiscard]] auto GetNGhostCells() const -> int { return n_ghost_cells_; }
    [[nodiscard]] auto GetPadding() const -> int { return n_ghost_cells_; }
    [[nodiscard]] auto GetDimension() const -> int { return dimension_; }
    [[nodiscard]] auto GetDim() const -> int { return dimension_; }

    /**
     * @brief Returns total array size (1D) or total size along x (2D)
     */
    [[nodiscard]] auto GetTotalSize() const -> int { return total_size_x_; }

    /**
     * @brief Returns total size along given axis
     */
    [[nodiscard]] auto GetTotalSize(int axis) const -> int;

    // ==================== Configuration ====================
    
    void SetN(int new_N);
    void SetPadding(int new_padding);
    void SetDim(int new_dim);

    // ==================== State Access (1D) ====================
    
    [[nodiscard]] auto GetPrimitive(int i) const -> Primitive;
    void SetPrimitive(int i, const Primitive& w);

    // ==================== State Access (2D) ====================
    
    [[nodiscard]] auto GetPrimitive2D(int i, int j) const -> Primitive;
    void SetPrimitive2D(int i, int j, const Primitive& w);

    /**
     * @brief Gets conservative state at (i,j) for 2D
     */
    [[nodiscard]] auto GetConservative2D(int i, int j, double gamma) const -> Conservative;

    /**
     * @brief Sets conservative state at (i,j) for 2D, updates all derived fields
     */
    void SetConservative2D(int i, int j, const Conservative& uc, double gamma, double dx, double dy);

private:
    int nx_ = 0;             ///< Physical cells in x
    int ny_ = 0;             ///< Physical cells in y (0 for 1D)
    int n_ghost_cells_ = 0;
    int dimension_ = 1;
    int total_size_x_ = 0;   ///< Total cells in x including ghosts
    int total_size_y_ = 0;   ///< Total cells in y including ghosts (0 for 1D)

    void RecomputeSizes();
    void Allocate1D();
    void Allocate2D();
};

#endif  // DATALAYER_HPP
