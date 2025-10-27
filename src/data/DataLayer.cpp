#include "../../include/data/DataLayer.hpp"


/**
 * @brief Constructs containers for initial values for all inner variables
 *
 * @param N Number of cells in scheme (1D case)
 * @param padding Number of ghost cells for scheme (1D case)
 */
DataLayer::DataLayer(const int N, const int padding) {
    rho = xt::zeros<double>({N + padding * 2});
    u = xt::zeros<double>({N + padding * 2});
    m = xt::zeros<double>({N + padding * 2});
    P = xt::zeros<double>({N + padding * 2});
    e = xt::zeros<double>({N + padding * 2});
    U = xt::zeros<double>({N + padding * 2});
    xb = xt::zeros<double>({N + padding * 2});
    xc = xt::zeros<double>({N + padding * 2});
}


/**
 * @brief Constructs containers for initial values for all inner variables(1D case)
 *
 * @param N Number of cells in scheme
 * @param padding Number of ghost cells for scheme
 * @param dim problem dimension
 */
DataLayer::DataLayer(const int N, const int padding, const int dim) {
    const xt::xarray<double> zero_dim = xt::ones<int>({dim}) * (N + padding * 2);
    rho = xt::zeros<double>(zero_dim);
    u = xt::zeros<double>(zero_dim);
    m = xt::zeros<double>(zero_dim);
    P = xt::zeros<double>(zero_dim);
    e = xt::zeros<double>(zero_dim);
    U = xt::zeros<double>(zero_dim);
    xb = xt::zeros<double>(zero_dim);
    xc = xt::zeros<double>(zero_dim);
}
