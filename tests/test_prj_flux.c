#include <math.h>

#include "prj.h"

static int prj_almost_equal(double a, double b)
{
    return fabs(a - b) < 1.0e-9;
}

int main(void)
{
    prj_eos eos;
    prj_rad rad;
    prj_block block;
    double total[PRJ_NVAR_CONS];
    int v;
    int i;
    int j;
    int k;

    eos.filename[0] = '\0';
    prj_rad_init(&rad);
    block.id = -1;
    block.W = 0;
    block.W1 = 0;
    block.U = 0;
    block.dUdt = 0;
    block.flux[0] = 0;
    block.flux[1] = 0;
    block.flux[2] = 0;
    block.dx[0] = 1.0;
    block.dx[1] = 1.0;
    block.dx[2] = 1.0;
    block.area[0] = 1.0;
    block.area[1] = 1.0;
    block.area[2] = 1.0;
    block.vol = 1.0;

    if (prj_block_alloc_data(&block) != 0) {
        return 1;
    }

    for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
        for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
            for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                block.W[VIDX(PRJ_PRIM_RHO, i, j, k)] = 2.0;
                block.W[VIDX(PRJ_PRIM_V1, i, j, k)] = 0.0;
                block.W[VIDX(PRJ_PRIM_V2, i, j, k)] = 3.0;
                block.W[VIDX(PRJ_PRIM_V3, i, j, k)] = 0.0;
                block.W[VIDX(PRJ_PRIM_EINT, i, j, k)] = 1.5;
                block.W[VIDX(PRJ_PRIM_YE, i, j, k)] = 0.2;
            }
        }
    }
    prj_flux_update(&eos, &rad, block.W, block.flux);
    if (!prj_almost_equal(block.flux[X1DIR][VIDX(PRJ_CONS_RHO, 1, 1, 1)], 0.0)) {
        prj_block_free_data(&block);
        return 2;
    }
    if (!prj_almost_equal(block.flux[X2DIR][VIDX(PRJ_CONS_RHO, 1, 1, 1)], 6.0)) {
        prj_block_free_data(&block);
        return 3;
    }
    if (!prj_almost_equal(block.flux[X3DIR][VIDX(PRJ_CONS_RHO, 1, 1, 1)], 0.0)) {
        prj_block_free_data(&block);
        return 4;
    }

    for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
        for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
            for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                double x = ((double)i + 0.5) / (double)PRJ_BLOCK_SIZE;
                double pulse = exp(-50.0 * (x - 0.5) * (x - 0.5));
                double rho = 1.0 + 0.1 * pulse;

                block.W[VIDX(PRJ_PRIM_RHO, i, j, k)] = rho;
                block.W[VIDX(PRJ_PRIM_V1, i, j, k)] = 0.0;
                block.W[VIDX(PRJ_PRIM_V2, i, j, k)] = 0.0;
                block.W[VIDX(PRJ_PRIM_V3, i, j, k)] = 0.0;
                block.W[VIDX(PRJ_PRIM_EINT, i, j, k)] = 1.0 / rho;
                block.W[VIDX(PRJ_PRIM_YE, i, j, k)] = 0.1;
            }
        }
    }
    prj_flux_update(&eos, &rad, block.W, block.flux);
    prj_flux_div(block.flux, block.area, block.vol, block.dUdt);
    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        total[v] = 0.0;
    }
    for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
        for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
            for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                for (v = 0; v < PRJ_NVAR_CONS; ++v) {
                    total[v] += block.dUdt[VIDX(v, i, j, k)] * block.vol;
                }
            }
        }
    }
    if (!prj_almost_equal(total[PRJ_CONS_RHO], 0.0)) {
        prj_block_free_data(&block);
        return 5;
    }

    prj_block_free_data(&block);
    return 0;
}
