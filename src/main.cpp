#include <iostream>
#include <xtensor/io/xio.hpp>
#include <xtensor/views/xview.hpp>
#include "../include/data/DataLayer.hpp"

int main() {
    const auto *init_cond = new DataLayer(10, 2);

    std::cout << xt::view(init_cond->rho, 0);
    return 0;
}
