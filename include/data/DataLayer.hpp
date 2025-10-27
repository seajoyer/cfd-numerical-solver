#ifndef DATALAYER_HPP
#define DATALAYER_HPP

#include <xtensor/io/xio.hpp>
#include <xtensor/views/xview.hpp>
#include <xtensor/containers/xarray.hpp>

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
    DataLayer(int N, int padding);
    DataLayer(int N, int padding, int dim);
    ~DataLayer() = default;
};

#endif //DATALAYER_HPP
