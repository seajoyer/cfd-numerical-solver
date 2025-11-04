#include "bc/PeriodicBoundary.hpp"
#include "data/DataLayer.hpp"


void PeriodicBoundary::Apply(DataLayer &layer, int axis, Side side) const {
    (void) axis;

    const int pad = layer.GetPadding();
    const int n = layer.GetN();
    const int core_start = layer.GetPadding();
    const int core_end = layer.GetCoreEndExclusive();


    if (side == Side::kLeft) {
        auto l_copy = xt::view(layer.rho, xt::range(core_end - pad, core_end));
        auto u_copy = xt::view(layer.u, xt::range(core_end - pad, core_end));
        auto P_copy = xt::view(layer.P, xt::range(core_end - pad, core_end));
        auto p_copy = xt::view(layer.p, xt::range(core_end - pad, core_end));
        auto e_copy = xt::view(layer.e, xt::range(core_end - pad, core_end));
        auto ue_copy = xt::view(layer.U, xt::range(core_end - pad, core_end));
        auto v_copy = xt::view(layer.V, xt::range(core_end - pad, core_end));
        auto m_copy = xt::view(layer.m, xt::range(core_end - pad, core_end));
        auto xb_copy = xt::view(layer.xb, xt::range(core_end - pad, core_end));
        auto xc_copy = xt::view(layer.xc, xt::range(core_end - pad, core_end));

        xt::view(layer.rho, xt::range(0, pad)) = l_copy;
        xt::view(layer.u, xt::range(0, pad)) = u_copy;
        xt::view(layer.P, xt::range(0, pad)) = P_copy;
        xt::view(layer.p, xt::range(0, pad)) = p_copy;
        xt::view(layer.e, xt::range(0, pad)) = e_copy;
        xt::view(layer.U, xt::range(0, pad)) = ue_copy;
        xt::view(layer.V, xt::range(0, pad)) = v_copy;
        xt::view(layer.m, xt::range(0, pad)) = m_copy;
        xt::view(layer.xb, xt::range(0, pad)) = xb_copy;
        xt::view(layer.xc, xt::range(0, pad)) = xc_copy;
    } else {
        auto copy_from_l = xt::view(layer.rho, xt::range(core_start, core_start + pad));
        auto copy_from_u = xt::view(layer.u, xt::range(core_start, core_start + pad));
        auto copy_from_P = xt::view(layer.P, xt::range(core_start, core_start + pad));
        auto copy_from_p = xt::view(layer.p, xt::range(core_start, core_start + pad));
        auto copy_from_e = xt::view(layer.e, xt::range(core_start, core_start + pad));
        auto copy_from_ue = xt::view(layer.U, xt::range(core_start, core_start + pad));
        auto copy_from_v = xt::view(layer.V, xt::range(core_start, core_start + pad));
        auto copy_from_m = xt::view(layer.m, xt::range(core_start, core_start + pad));
        auto copy_from_xb = xt::view(layer.xb, xt::range(core_start, core_start + pad));
        auto copy_from_xc = xt::view(layer.xc, xt::range(core_start, core_start + pad));

        xt::view(layer.rho, xt::range(core_end, core_end + pad)) = copy_from_l;
        xt::view(layer.u, xt::range(core_end, core_end + pad)) = copy_from_u;
        xt::view(layer.P, xt::range(core_end, core_end + pad)) = copy_from_P;
        xt::view(layer.p, xt::range(core_end, core_end + pad)) = copy_from_p;
        xt::view(layer.e, xt::range(core_end, core_end + pad)) = copy_from_e;
        xt::view(layer.U, xt::range(core_end, core_end + pad)) = copy_from_ue;
        xt::view(layer.V, xt::range(core_end, core_end + pad)) = copy_from_v;
        xt::view(layer.m, xt::range(core_end, core_end + pad)) = copy_from_m;
        xt::view(layer.xb, xt::range(core_end, core_end + pad)) = copy_from_xb;
        xt::view(layer.xc, xt::range(core_end, core_end + pad)) = copy_from_xc;
    }

    (void) n;
}
