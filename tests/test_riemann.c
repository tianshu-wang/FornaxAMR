#include <math.h>

#include "prj.h"

static int prj_almost_equal(double a, double b)
{
    return fabs(a - b) < 1.0e-10;
}

int main(void)
{
    prj_eos eos;
    prj_coord coord = {-1.0, 1.0, -1.0, 1.0, -1.0, 1.0};
    prj_mesh mesh;
    double W[PRJ_NVAR_PRIM] = {1.2, 0.7, -0.2, 0.1, 2.5, 0.05};
    double flux[PRJ_NVAR_CONS];
    double q[PRJ_EOS_NQUANT];
    double U[PRJ_NVAR_CONS];
    double expected[PRJ_NVAR_CONS];
    double WL[PRJ_NVAR_PRIM] = {1.0, 2.0, 0.0, 0.0, 1.5, 0.1};
    double WR[PRJ_NVAR_PRIM] = {1.0, -2.0, 0.0, 0.0, 1.5, 0.1};
    double WS_left[PRJ_NVAR_PRIM] = {1.0, 2.0, 0.0, 0.0, 2.0, 0.1};
    double WS_right[PRJ_NVAR_PRIM] = {2.0, -2.0, 0.0, 0.0, 8.0, 0.1};
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
    expected[PRJ_CONS_YE] = U[PRJ_CONS_YE] * W[PRJ_PRIM_V1];

    prj_riemann_hllc(W, W, &eos, flux);
    if (!prj_almost_equal(flux[PRJ_CONS_RHO], expected[PRJ_CONS_RHO]) ||
        !prj_almost_equal(flux[PRJ_CONS_MOM1], expected[PRJ_CONS_MOM1]) ||
        !prj_almost_equal(flux[PRJ_CONS_MOM2], expected[PRJ_CONS_MOM2]) ||
        !prj_almost_equal(flux[PRJ_CONS_MOM3], expected[PRJ_CONS_MOM3]) ||
        !prj_almost_equal(flux[PRJ_CONS_ETOT], expected[PRJ_CONS_ETOT]) ||
        !prj_almost_equal(flux[PRJ_CONS_YE], expected[PRJ_CONS_YE])) {
        return 1;
    }

    prj_riemann_hllc(WL, WR, &eos, flux_lr);
    WL[PRJ_PRIM_V1] = -WL[PRJ_PRIM_V1];
    WR[PRJ_PRIM_V1] = -WR[PRJ_PRIM_V1];
    prj_riemann_hllc(WR, WL, &eos, flux_rl);
    if (!prj_almost_equal(flux_lr[PRJ_CONS_RHO], -flux_rl[PRJ_CONS_RHO])) {
        return 2;
    }

    if (prj_riemann_detect_shock(WS_left, WS_right) != X1DIR) {
        return 3;
    }

    if (prj_mesh_init(&mesh, 2, 2, 2, 1, &coord) != 0) {
        return 4;
    }
    prj_amr_refine_block(&mesh, 7);
    {
        prj_block *fine_parent = prj_mesh_get_block(&mesh, 7);
        prj_block *fine = prj_mesh_get_block(&mesh, fine_parent->children[0]);
        prj_block *coarse = prj_mesh_get_block(&mesh, 3);
        int v;

        prj_riemann_set_mesh(&mesh);
        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
            fine->flux[X1DIR][VIDX(v, 0, 0, 0)] = 100.0 + (double)v;
            fine->flux[X1DIR][VIDX(v, 0, 1, 0)] = 200.0 + (double)v;
            fine->flux[X1DIR][VIDX(v, 0, 0, 1)] = 300.0 + (double)v;
            fine->flux[X1DIR][VIDX(v, 0, 1, 1)] = 400.0 + (double)v;
            coarse->flux[X1DIR][VIDX(v, PRJ_BLOCK_SIZE, 0, 0)] = 0.0;
        }
        prj_riemann_flux_send(fine);
        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
            if (!prj_almost_equal(coarse->flux[X1DIR][VIDX(v, PRJ_BLOCK_SIZE, 0, 0)], 250.0 + (double)v)) {
                prj_mesh_destroy(&mesh);
                return 5;
            }
        }
    }
    prj_mesh_destroy(&mesh);

    return 0;
}
