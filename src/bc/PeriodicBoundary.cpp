#include "bc/PeriodicBoundary.hpp"
#include "data/DataLayer.hpp"


void PeriodicBoundary::Apply(DataLayer &layer, int axis, Side side) const {
    (void) axis;

    const int pad = layer.GetPadding();
    const int n = layer.GetN();
    const int coreStart = layer.GetPadding();
    const int coreEnd = layer.GetCoreEndExclusive();


    if (side == Side::Min) {
        auto LCopy = xt::view(layer.rho, xt::range(coreEnd - pad, coreEnd));
        auto UCopy = xt::view(layer.u, xt::range(coreEnd - pad, coreEnd));
        auto PCopy = xt::view(layer.P, xt::range(coreEnd - pad, coreEnd));
        auto pCopy = xt::view(layer.p, xt::range(coreEnd - pad, coreEnd));
        auto eCopy = xt::view(layer.e, xt::range(coreEnd - pad, coreEnd));
        auto UeCopy = xt::view(layer.U, xt::range(coreEnd - pad, coreEnd));
        auto VCopy = xt::view(layer.V, xt::range(coreEnd - pad, coreEnd));
        auto mCopy = xt::view(layer.m, xt::range(coreEnd - pad, coreEnd));
        auto xbCopy = xt::view(layer.xb, xt::range(coreEnd - pad, coreEnd));
        auto xcCopy = xt::view(layer.xc, xt::range(coreEnd - pad, coreEnd));

        xt::view(layer.rho, xt::range(0, pad)) = LCopy;
        xt::view(layer.u, xt::range(0, pad)) = UCopy;
        xt::view(layer.P, xt::range(0, pad)) = PCopy;
        xt::view(layer.p, xt::range(0, pad)) = pCopy;
        xt::view(layer.e, xt::range(0, pad)) = eCopy;
        xt::view(layer.U, xt::range(0, pad)) = UeCopy;
        xt::view(layer.V, xt::range(0, pad)) = VCopy;
        xt::view(layer.m, xt::range(0, pad)) = mCopy;
        xt::view(layer.xb, xt::range(0, pad)) = xbCopy;
        xt::view(layer.xc, xt::range(0, pad)) = xcCopy;
    } else {
        auto copyFromL = xt::view(layer.rho, xt::range(coreStart, coreStart + pad));
        auto copyFromU = xt::view(layer.u, xt::range(coreStart, coreStart + pad));
        auto copyFromP = xt::view(layer.P, xt::range(coreStart, coreStart + pad));
        auto copyFromp = xt::view(layer.p, xt::range(coreStart, coreStart + pad));
        auto copyFrome = xt::view(layer.e, xt::range(coreStart, coreStart + pad));
        auto copyFromUe = xt::view(layer.U, xt::range(coreStart, coreStart + pad));
        auto copyFromV = xt::view(layer.V, xt::range(coreStart, coreStart + pad));
        auto copyFromm = xt::view(layer.m, xt::range(coreStart, coreStart + pad));
        auto copyFromxb = xt::view(layer.xb, xt::range(coreStart, coreStart + pad));
        auto copyFromxc = xt::view(layer.xc, xt::range(coreStart, coreStart + pad));

        xt::view(layer.rho, xt::range(coreEnd, coreEnd + pad)) = copyFromL;
        xt::view(layer.u, xt::range(coreEnd, coreEnd + pad)) = copyFromU;
        xt::view(layer.P, xt::range(coreEnd, coreEnd + pad)) = copyFromP;
        xt::view(layer.p, xt::range(coreEnd, coreEnd + pad)) = copyFromp;
        xt::view(layer.e, xt::range(coreEnd, coreEnd + pad)) = copyFrome;
        xt::view(layer.U, xt::range(coreEnd, coreEnd + pad)) = copyFromUe;
        xt::view(layer.V, xt::range(coreEnd, coreEnd + pad)) = copyFromV;
        xt::view(layer.m, xt::range(coreEnd, coreEnd + pad)) = copyFromm;
        xt::view(layer.xb, xt::range(coreEnd, coreEnd + pad)) = copyFromxb;
        xt::view(layer.xc, xt::range(coreEnd, coreEnd + pad)) = copyFromxc;
    }

    (void) n;
}