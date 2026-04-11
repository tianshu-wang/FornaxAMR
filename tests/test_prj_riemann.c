#include <math.h>

#include "prj.h"

static int prj_almost_equal(double a, double b)
{
    return fabs(a - b) < 1.0e-10;
}

int main(void)
{
    prj_eos eos;
    double W[PRJ_NVAR_PRIM] = {1.2, 0.7, -0.2, 0.1, 2.5, 0.05};
    double flux[PRJ_NVAR_CONS];
    double q[PRJ_EOS_NQUANT];
    double U[PRJ_NVAR_CONS];
    double expected[PRJ_NVAR_CONS];
    double WL[PRJ_NVAR_PRIM] = {1.0, 2.0, 0.0, 0.0, 1.5, 0.1};
    double WR[PRJ_NVAR_PRIM] = {1.0, -2.0, 0.0, 0.0, 1.5, 0.1};
    double flux_lr[PRJ_NVAR_CONS];
    double flux_rl[PRJ_NVAR_CONS];

    eos.filename[0] = '\0';

    prj_eos_rey(&eos, W[PRJ_PRIM_RHO], W[PRJ_PRIM_EINT], W[PRJ_PRIM_YE], q);
    prj_eos_prim2cons(&eos, W, U);
    expected[PRJ_CONS_RHO] = U[PRJ_CONS_RHO] * W[PRJ_PRIM_V1];
    expected[PRJ_CONS_MOM1] = U[PRJ_CONS_MOM1] * W[PRJ_PRIM_V1] + q[PRJ_EOS_PRESSURE];
    expected[PRJ_CONS_MOM2] = U[PRJ_CONS_MOM2] * W[PRJ_PRIM_V1];
    expected[PRJ_CONS_MOM3] = U[PRJ_CONS_MOM3] * W[PRJ_PRIM_V1];
    expected[PRJ_CONS_ETOT] = (U[PRJ_CONS_ETOT] + q[PRJ_EOS_PRESSURE]) * W[PRJ_PRIM_V1];
    expected[PRJ_CONS_YE] = U[PRJ_CONS_YE] * W[PRJ_PRIM_V1] / U[PRJ_CONS_RHO];

    prj_riemann_hllc(W, W, &eos, flux);
    if (!prj_almost_equal(flux[PRJ_CONS_RHO], expected[PRJ_CONS_RHO]) ||
        !prj_almost_equal(flux[PRJ_CONS_MOM1], expected[PRJ_CONS_MOM1]) ||
        !prj_almost_equal(flux[PRJ_CONS_ETOT], expected[PRJ_CONS_ETOT])) {
        return 1;
    }

    prj_riemann_hllc(WL, WR, &eos, flux_lr);
    WL[PRJ_PRIM_V1] = -WL[PRJ_PRIM_V1];
    WR[PRJ_PRIM_V1] = -WR[PRJ_PRIM_V1];
    prj_riemann_hllc(WR, WL, &eos, flux_rl);
    if (!prj_almost_equal(flux_lr[PRJ_CONS_RHO], -flux_rl[PRJ_CONS_RHO])) {
        return 2;
    }

    return 0;
}
