#include <math.h>

#include "prj.h"

#if PRJ_NRAD > 0
static double prj_src_interp_lapse(const prj_grav_mono *grav_mono, double r)
{
    int idx;

    if (grav_mono == 0 || grav_mono->nbins <= 0 || grav_mono->lapse == 0) {
        return 1.0;
    }
    if (r <= grav_mono->rf[0]) {
        return grav_mono->lapse[0];
    }
    for (idx = 0; idx < grav_mono->nbins; ++idx) {
        double r0 = grav_mono->rf[idx];
        double r1 = grav_mono->rf[idx + 1];

        if (r <= r1) {
            double weight = (r - r0) / (r1 - r0);

            return (1.0 - weight) * grav_mono->lapse[idx] + weight * grav_mono->lapse[idx + 1];
        }
    }
    return grav_mono->lapse[grav_mono->nbins];
}
#endif

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

                    accel = prj_gravity_interp_accel(grav_mono, r);
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
#if PRJ_NRAD > 0
                    {
                        double lapse = prj_src_interp_lapse(grav_mono, r);
                        int field;
                        int group;

                        for (field = 0; field < PRJ_NRAD; ++field) {
                            for (group = 0; group < PRJ_NEGROUP; ++group) {
                                dUdt[VIDX(PRJ_CONS_RAD_F1(field, group), i, j, k)] +=
                                    lapse * U[VIDX(PRJ_CONS_RAD_E(field, group), i, j, k)] * g1;
                                dUdt[VIDX(PRJ_CONS_RAD_F2(field, group), i, j, k)] +=
                                    lapse * U[VIDX(PRJ_CONS_RAD_E(field, group), i, j, k)] * g2;
                                dUdt[VIDX(PRJ_CONS_RAD_F3(field, group), i, j, k)] +=
                                    lapse * U[VIDX(PRJ_CONS_RAD_E(field, group), i, j, k)] * g3;
                            }
                        }
                    }
#endif
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
