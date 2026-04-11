#include <math.h>

#include "prj.h"

static double prj_src_radius_center(const prj_grav_mono *grav_mono, int idx)
{
    return 0.5 * (grav_mono->rf[idx] + grav_mono->rf[idx + 1]);
}

static double prj_src_interp_accel(const prj_grav_mono *grav_mono, double r)
{
    int idx;

    if (grav_mono == 0 || grav_mono->nbins <= 0 || grav_mono->accel == 0) {
        return 0.0;
    }
    if (r <= prj_src_radius_center(grav_mono, 0)) {
        return grav_mono->accel[0];
    }
    for (idx = 0; idx < grav_mono->nbins - 1; ++idx) {
        double r0 = prj_src_radius_center(grav_mono, idx);
        double r1 = prj_src_radius_center(grav_mono, idx + 1);

        if (r <= r1) {
            double weight = (r - r0) / (r1 - r0);

            return (1.0 - weight) * grav_mono->accel[idx] + weight * grav_mono->accel[idx + 1];
        }
    }
    return grav_mono->accel[grav_mono->nbins - 1];
}

void prj_src_geom(prj_eos *eos, double *W, double *dUdt)
{
    (void)eos;
    (void)W;
    (void)dUdt;
}

void prj_src_user(prj_eos *eos, double *W, double *dUdt)
{
    (void)eos;
    (void)W;
    (void)dUdt;
}

void prj_src_monopole_gravity(prj_mesh *mesh, const prj_grav_mono *grav_mono,
    double *restrict W, double *restrict U, double *restrict dUdt)
{
    int bidx;

    if (mesh == 0 || grav_mono == 0) {
        return;
    }

    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        int i;
        int j;
        int k;

        if (block->id < 0 || block->active != 1 || block->W != W || block->U != U || block->dUdt != dUdt) {
            continue;
        }
        for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
            for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    double x1 = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
                    double x2 = block->xmin[1] + ((double)j + 0.5) * block->dx[1];
                    double x3 = block->xmin[2] + ((double)k + 0.5) * block->dx[2];
                    double r = sqrt(x1 * x1 + x2 * x2 + x3 * x3);
                    double accel;
                    double g1;
                    double g2;
                    double g3;

                    if (r <= 0.0) {
                        continue;
                    }

                    accel = prj_src_interp_accel(grav_mono, r);
                    g1 = accel * x1 / r;
                    g2 = accel * x2 / r;
                    g3 = accel * x3 / r;

                    dUdt[VIDX(PRJ_CONS_MOM1, i, j, k)] += U[VIDX(PRJ_CONS_RHO, i, j, k)] * g1;
                    dUdt[VIDX(PRJ_CONS_MOM2, i, j, k)] += U[VIDX(PRJ_CONS_RHO, i, j, k)] * g2;
                    dUdt[VIDX(PRJ_CONS_MOM3, i, j, k)] += U[VIDX(PRJ_CONS_RHO, i, j, k)] * g3;
                    dUdt[VIDX(PRJ_CONS_ETOT, i, j, k)] +=
                        U[VIDX(PRJ_CONS_MOM1, i, j, k)] * g1 +
                        U[VIDX(PRJ_CONS_MOM2, i, j, k)] * g2 +
                        U[VIDX(PRJ_CONS_MOM3, i, j, k)] * g3;
                }
            }
        }
    }
}

void prj_src_update(prj_mesh *mesh, prj_eos *eos, double *restrict W,
    double *restrict U, double *restrict dUdt)
{
    prj_src_geom(eos, W, dUdt);
    prj_src_user(eos, W, dUdt);
    prj_src_monopole_gravity(mesh, prj_gravity_active_monopole(), W, U, dUdt);
}
