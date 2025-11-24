#ifndef DATALAYER_HPP
#define DATALAYER_HPP

#include <xtensor.hpp>
#include <xtensor/io/xio.hpp>

#include "data/Variables.hpp"

/**
 * @file DataLayer.hpp
 * @brief Data storage container for all simulation fields
 */

/**
 * @struct DataLayer
 * @brief Structure-of-arrays storage for computational fluid dynamics fields
 * 
 * DataLayer manages all physical and auxiliary variables required for finite-volume
 * simulations of the Euler equations. Fields are stored in xtensor arrays with
 * separate allocations for core cells and ghost cells (padding).
 * 
 * ## Storage Layout (1D)
 * For N physical cells with padding P:
 * ```
 * [ghost: P cells][core: N cells][ghost: P cells]
 * Total size: N + 2*P
 * ```
 * 
 * ## Field Categories
 * 
 * **Primitive variables** (directly physical):
 * - rho: density (ρ)
 * - u: velocity
 * - P: thermodynamic pressure
 * 
 * **Conservative variables** (integrated quantities):
 * - e: total energy density
 * - p: momentum density (ρu)
 * 
 * **Auxiliary variables** (derived):
 * - U: specific internal energy
 * - V: specific volume (1/ρ)
 * - m: cell mass (ρ·Δx)
 * 
 * **Geometric data**:
 * - xc: cell center coordinates
 * - xb: cell boundary coordinates
 * 
 * ## Usage Pattern
 * ```cpp
 * DataLayer layer(N, padding, dim);
 * 
 * // Access core cells
 * int start = layer.GetCoreStart();
 * int end = layer.GetCoreEndExclusive();
 * for (int i = start; i < end; ++i) {
 *     layer.rho(i) = 1.0;
 *     layer.u(i) = 0.0;
 *     layer.P(i) = 1.0;
 * }
 * 
 * // Convenient primitive state access
 * Primitive w = layer.GetPrimitive(i);
 * layer.SetPrimitive(i, w);
 * ```
 * 
 * @note Ghost cells must be filled via boundary conditions before solver steps
 * @note Currently optimized for 1D; 2D/3D support planned
 * @see BoundaryCondition, Primitive, Conservative
 */
struct DataLayer {
    // ==================== Field Arrays ====================
    
    /** @brief Density field (ρ) [kg/m³] */
    xt::xarray<double> rho;
    
    /** @brief Velocity field in x-direction (u) [m/s] */
    xt::xarray<double> u;
    
    /** @brief Thermodynamic pressure field (P) [Pa] */
    xt::xarray<double> P;
    
    /** @brief Momentum density field (p = ρu) [kg/(m²·s)] */
    xt::xarray<double> p;
    
    /** @brief Total energy density (e) [J/m³] */
    xt::xarray<double> e;
    
    /** @brief Specific internal energy (U = e_int/ρ) [J/kg] */
    xt::xarray<double> U;
    
    /** @brief Specific volume (V = 1/ρ) [m³/kg] */
    xt::xarray<double> V;
    
    /** @brief Cell mass (m = ρ·Δx) [kg] */
    xt::xarray<double> m;
    
    /** @brief Cell boundary coordinates [m]
     * @note Size: total_size
     */
    xt::xarray<double> xb;
    
    /** @brief Cell center coordinates [m]
     * @note Size: total_size
     */
    xt::xarray<double> xc;

    // ==================== Constructors ====================
    
    /**
     * @brief Default constructor creating empty layer
     * @note Must call SetN() and SetPadding() before use
     */
    DataLayer() = default;

    /**
     * @brief Constructs 1D data layer with specified resolution
     * 
     * Allocates arrays for N physical cells plus ghost cells.
     * All fields are zero-initialized.
     * 
     * @param N Number of physical (core) cells
     * @param padding Number of ghost cells on each boundary
     * @throw std::invalid_argument if N ≤ 0 or padding < 0
     */
    DataLayer(int N, int padding);

    /**
     * @brief Constructs data layer with specified dimension
     * 
     * @param N Number of cells per direction
     * @param padding Number of ghost cells per direction
     * @param dim Spatial dimension (1, 2, or 3)
     * @throw std::invalid_argument if parameters invalid or dim ∉ {1,2,3}
     * @note Only dim=1 is currently fully implemented
     */
    DataLayer(int N, int padding, int dim);

    /**
     * @brief Default destructor
     */
    ~DataLayer() = default;

    // ==================== Grid Information ====================
    
    /**
     * @brief Returns starting index of core (physical) domain
     * 
     * For 1D with padding P, returns P.
     * First physical cell is at index GetCoreStart().
     * 
     * @param axis Spatial axis (0=x, 1=y, 2=z)
     * @return Index of first physical cell
     * @note Currently axis parameter unused (1D only)
     */
    [[nodiscard]] auto GetCoreStart(int axis = 0) const -> int;

    /**
     * @brief Returns one-past-last index of core domain (exclusive end)
     * 
     * For 1D with N cells and padding P, returns P + N.
     * Last physical cell is at index GetCoreEndExclusive() - 1.
     * 
     * @param axis Spatial axis (0=x, 1=y, 2=z)
     * @return Exclusive end index of physical cells
     * @note Currently axis parameter unused (1D only)
     */
    [[nodiscard]] auto GetCoreEndExclusive(int axis = 0) const -> int;

    /**
     * @brief Returns number of physical cells (excluding ghosts)
     * @return Core grid size N
     */
    [[nodiscard]] auto GetN() const -> int { return n_; }

    /**
     * @brief Returns number of ghost cells on each boundary
     * @return Padding size
     */
    [[nodiscard]] auto GetNGhostCells() const -> int { return n_ghost_cells_; }

    /**
     * @brief Returns padding (alias for GetNGhostCells)
     * @return Number of ghost cells per side
     */
    [[nodiscard]] auto GetPadding() const -> int { return n_ghost_cells_; }

    /**
     * @brief Returns spatial dimension of the problem
     * @return Dimension (1, 2, or 3)
     */
    [[nodiscard]] auto GetDimension() const -> int { return dimension_; }

    /**
     * @brief Returns dimension (alias for GetDimension)
     * @return Spatial dimension
     */
    [[nodiscard]] auto GetDim() const -> int { return dimension_; }

    /**
     * @brief Returns total array size including ghost cells
     * 
     * For 1D: total_size = N + 2·padding
     * 
     * @return Total number of cells (physical + ghost)
     */
    [[nodiscard]] auto GetTotalSize() const -> int { return total_size_; }

    // ==================== Configuration ====================
    
    /**
     * @brief Changes grid resolution and reallocates all arrays
     * 
     * Discards existing data and creates new zero-initialized arrays
     * with the specified number of physical cells.
     * 
     * @param new_N New number of physical cells
     * @throw std::invalid_argument if new_N ≤ 0
     */
    void SetN(int new_N);

    /**
     * @brief Changes ghost cell count and reallocates arrays
     * 
     * @param new_padding New number of ghost cells per side
     * @throw std::invalid_argument if new_padding < 0
     */
    void SetPadding(int new_padding);

    /**
     * @brief Changes spatial dimension and reallocates arrays
     * 
     * @param new_dim New dimension (1, 2, or 3)
     * @throw std::invalid_argument if new_dim ∉ {1,2,3}
     * @warning 2D/3D support incomplete; use with caution
     */
    void SetDim(int new_dim);

    // ==================== State Access ====================
    
    /**
     * @brief Retrieves primitive state at specified cell
     * 
     * Convenience method that bundles (ρ, u, P) into a Primitive struct.
     * 
     * @param i Cell index (0-based, includes ghost cells)
     * @return Primitive state (rho, u, P) at cell i
     * @note No bounds checking; ensure i < GetTotalSize()
     */
    [[nodiscard]] auto GetPrimitive(int i) const -> Primitive;

    /**
     * @brief Sets primitive state at specified cell
     * 
     * Writes (ρ, u, P) components to arrays.
     * Auxiliary fields (e, U, V, m, p) are NOT updated.
     * 
     * @param i Cell index (0-based, includes ghost cells)
     * @param w Primitive state to store
     * @note No bounds checking; ensure i < GetTotalSize()
     * @note Caller responsible for maintaining consistency with auxiliary fields
     */
    void SetPrimitive(int i, const Primitive& w);

private:
    /** @brief Number of physical cells */
    int n_ = 0;
    
    /** @brief Number of ghost cells on each boundary */
    int n_ghost_cells_ = 0;
    
    /** @brief Spatial dimension (1, 2, or 3) */
    int dimension_ = 1;
    
    /** @brief Total cell count including ghosts */
    int total_size_ = 0;

    /**
     * @brief Recomputes total_size after changing n or n_ghost_cells
     */
    void RecomputeSizes();

    /**
     * @brief Allocates 1D arrays for all fields
     * 
     * Creates zero-initialized xtensor arrays of size total_size
     * for all physical, conservative, auxiliary, and geometric fields.
     */
    void Allocate1D();
};

#endif  // DATALAYER_HPP
