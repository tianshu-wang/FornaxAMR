#include <math.h>
#include <stdio.h>

#include "prj.h"

static double prj_total_cons(const prj_mesh *mesh, int var)
{
    double total = 0.0;
    int bidx;

    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        const prj_block *block = &mesh->blocks[bidx];
        int i;
        int j;
        int k;

        if (block->id < 0 || block->active != 1) {
            continue;
        }
        for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
            for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    total += block->U[VIDX(var, i, j, k)] * block->vol;
                }
            }
        }
    }
    return total;
}

static void prj_store_prim(double *dst, int i, int j, int k, const double *w)
{
    int v;

    for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
        dst[VIDX(v, i, j, k)] = w[v];
    }
}

static void prj_store_cons(double *dst, int i, int j, int k, const double *u)
{
    int v;

    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        dst[VIDX(v, i, j, k)] = u[v];
    }
}

static void prj_fill_uniform(prj_mesh *mesh, double rho, double eint)
{
    int bidx;

    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        int i;
        int j;
        int k;

        if (block->id < 0 || block->active != 1) {
            continue;
        }
        for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
            for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
                for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                    double W[PRJ_NVAR_PRIM];
                    double U[PRJ_NVAR_CONS];

                    W[PRJ_PRIM_RHO] = rho;
                    W[PRJ_PRIM_V1] = 0.0;
                    W[PRJ_PRIM_V2] = 0.0;
                    W[PRJ_PRIM_V3] = 0.0;
                    W[PRJ_PRIM_EINT] = eint;
                    W[PRJ_PRIM_YE] = 0.1;
                    prj_store_prim(block->W, i, j, k, W);
                    prj_eos_prim2cons((prj_eos *)0, W, U);
                    prj_store_cons(block->U, i, j, k, U);
                }
            }
        }
    }
}

int main(void)
{
    prj_coord coord = {-1.0, 1.0, -1.0, 1.0, -1.0, 1.0};
    prj_bc bc = {
        PRJ_BC_OUTFLOW, PRJ_BC_OUTFLOW,
        PRJ_BC_OUTFLOW, PRJ_BC_OUTFLOW,
        PRJ_BC_OUTFLOW, PRJ_BC_OUTFLOW
    };
    prj_eos eos;
    prj_rad rad;
    prj_mesh mesh;
    prj_mesh amr_mesh;
    double mass0;
    double energy0;
    double mass1;
    double energy1;
    double dt;
    int step;

    eos.filename[0] = '\0';
    prj_rad_init(&rad);

    if (prj_mesh_init(&mesh, 1, 1, 1, 0, &coord) != 0) {
        return 1;
    }
    prj_fill_uniform(&mesh, 1.0, 1.0);
    mass0 = prj_total_cons(&mesh, PRJ_CONS_RHO);
    energy0 = prj_total_cons(&mesh, PRJ_CONS_ETOT);
    for (step = 0; step < 10; ++step) {
        dt = prj_timeint_calc_dt(&mesh, &eos, 0.4);
        printf("uniform step %d boundary fluxes 0 0\n", step);
        prj_timeint_step(&mesh, &coord, &bc, &eos, &rad, dt);
    }
    mass1 = prj_total_cons(&mesh, PRJ_CONS_RHO);
    energy1 = prj_total_cons(&mesh, PRJ_CONS_ETOT);
    if (fabs((mass1 - mass0) / mass0) > 1.0e-12) {
        prj_mesh_destroy(&mesh);
        return 2;
    }
    if (fabs((energy1 - energy0) / energy0) > 1.0e-12) {
        prj_mesh_destroy(&mesh);
        return 3;
    }
    prj_mesh_destroy(&mesh);

    if (prj_mesh_init(&amr_mesh, 2, 2, 2, 1, &coord) != 0) {
        return 4;
    }
    prj_fill_uniform(&amr_mesh, 2.0, 0.5);
    prj_amr_refine_block(&amr_mesh, 7);
    prj_fill_uniform(&amr_mesh, 2.0, 0.5);
    mass0 = prj_total_cons(&amr_mesh, PRJ_CONS_RHO);
    energy0 = prj_total_cons(&amr_mesh, PRJ_CONS_ETOT);
    for (step = 0; step < 10; ++step) {
        dt = prj_timeint_calc_dt(&amr_mesh, &eos, 0.4);
        printf("amr step %d boundary fluxes 0 0\n", step);
        prj_timeint_step(&amr_mesh, &coord, &bc, &eos, &rad, dt);
    }
    mass1 = prj_total_cons(&amr_mesh, PRJ_CONS_RHO);
    energy1 = prj_total_cons(&amr_mesh, PRJ_CONS_ETOT);
    if (fabs((mass1 - mass0) / mass0) > 1.0e-12) {
        prj_mesh_destroy(&amr_mesh);
        return 5;
    }
    if (fabs((energy1 - energy0) / energy0) > 1.0e-12) {
        prj_mesh_destroy(&amr_mesh);
        return 6;
    }
    prj_mesh_destroy(&amr_mesh);
    return 0;
}
