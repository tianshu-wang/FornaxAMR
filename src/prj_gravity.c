#include <math.h>
#include <stddef.h>
#include <stdlib.h>

#include "prj.h"

#if defined(PRJ_ENABLE_MPI)
#include <mpi.h>
#endif

#define PRJ_GRAVITY_DEFAULT_NBINS 256

static prj_grav_mono *prj_gravity_active = 0;
static double prj_gravity_rmax = 0.0;

static double prj_gravity_min_double(double a, double b)
{
    return a < b ? a : b;
}

static int prj_gravity_block_is_local_active(const prj_block *block)
{
    prj_mpi *mpi = prj_mpi_current();

    return block != 0 && block->id >= 0 && block->active == 1 &&
        (mpi == 0 || block->rank == mpi->rank);
}

static int prj_gravity_bin_index(const prj_grav_mono *grav_mono, double r)
{
    double log_span;
    int idx;

    if (grav_mono == 0 || grav_mono->nbins <= 0 || r < grav_mono->dr_min || r >= prj_gravity_rmax) {
        return -1;
    }

    log_span = log(prj_gravity_rmax / grav_mono->dr_min);
    if (log_span <= 0.0) {
        return -1;
    }
    idx = (int)((double)grav_mono->nbins * log(r / grav_mono->dr_min) / log_span);
    if (idx < 0) {
        idx = 0;
    }
    if (idx >= grav_mono->nbins) {
        return -1;
    }
    return idx;
}

void prj_gravity_free(prj_grav_mono *grav_mono)
{
    if (grav_mono == 0) {
        return;
    }

    free(grav_mono->rf);
    free(grav_mono->ms);
    free(grav_mono->phi);
    free(grav_mono->accel);
    grav_mono->rf = 0;
    grav_mono->ms = 0;
    grav_mono->phi = 0;
    grav_mono->accel = 0;
    grav_mono->nbins = 0;
    grav_mono->dr_min = 0.0;
    if (prj_gravity_active == grav_mono) {
        prj_gravity_active = 0;
        prj_gravity_rmax = 0.0;
    }
}

void prj_gravity_init(prj_sim *sim)
{
    prj_grav_mono *grav_mono;
    double span1;
    double span2;
    double span3;
    double min_cell;
    double log_span;
    int i;

    if (sim == 0) {
        return;
    }

    grav_mono = &sim->monopole_grav;
    prj_gravity_free(grav_mono);

    span1 = fabs(sim->mesh.coord.x1max - sim->mesh.coord.x1min);
    span2 = fabs(sim->mesh.coord.x2max - sim->mesh.coord.x2min);
    span3 = fabs(sim->mesh.coord.x3max - sim->mesh.coord.x3min);
    prj_gravity_rmax = 0.5 * prj_gravity_min_double(span1, prj_gravity_min_double(span2, span3));
    min_cell = prj_gravity_min_double(sim->mesh.blocks[0].dx[0],
        prj_gravity_min_double(sim->mesh.blocks[0].dx[1], sim->mesh.blocks[0].dx[2]));

    grav_mono->nbins = PRJ_GRAVITY_DEFAULT_NBINS;
    grav_mono->dr_min = 0.5 * min_cell;
    if (grav_mono->dr_min <= 0.0 || prj_gravity_rmax <= grav_mono->dr_min) {
        grav_mono->dr_min = prj_gravity_rmax / (double)grav_mono->nbins;
    }

    grav_mono->rf = (double *)calloc((size_t)grav_mono->nbins + 1U, sizeof(*grav_mono->rf));
    grav_mono->ms = (double *)calloc((size_t)grav_mono->nbins, sizeof(*grav_mono->ms));
    grav_mono->phi = (double *)calloc((size_t)grav_mono->nbins + 1U, sizeof(*grav_mono->phi));
    grav_mono->accel = (double *)calloc((size_t)grav_mono->nbins, sizeof(*grav_mono->accel));
    if (grav_mono->rf == 0 || grav_mono->ms == 0 || grav_mono->phi == 0 || grav_mono->accel == 0) {
        prj_gravity_free(grav_mono);
        return;
    }

    log_span = log(prj_gravity_rmax / grav_mono->dr_min);
    for (i = 0; i <= grav_mono->nbins; ++i) {
        grav_mono->rf[i] = exp(log_span * (double)i / (double)grav_mono->nbins) * grav_mono->dr_min;
    }

    prj_gravity_active = grav_mono;
}

const prj_grav_mono *prj_gravity_active_monopole(void)
{
    return prj_gravity_active;
}

void prj_gravity_monopole_reduce(prj_mesh *mesh)
{
    prj_grav_mono *grav_mono = prj_gravity_active;
    int bidx;
    int idx;

    if (mesh == 0 || grav_mono == 0 || grav_mono->ms == 0) {
        return;
    }

    for (idx = 0; idx < grav_mono->nbins; ++idx) {
        grav_mono->ms[idx] = 0.0;
    }

    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        const prj_block *block = &mesh->blocks[bidx];
        int i;
        int j;
        int k;

        if (!prj_gravity_block_is_local_active(block)) {
            continue;
        }
        for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
            for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    double x1 = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
                    double x2 = block->xmin[1] + ((double)j + 0.5) * block->dx[1];
                    double x3 = block->xmin[2] + ((double)k + 0.5) * block->dx[2];
                    double r = sqrt(x1 * x1 + x2 * x2 + x3 * x3);

                    idx = prj_gravity_bin_index(grav_mono, r);
                    if (idx >= 0) {
                        grav_mono->ms[idx] += block->W[VIDX(PRJ_PRIM_RHO, i, j, k)] * block->vol;
                    }
                }
            }
        }
    }

#if defined(PRJ_ENABLE_MPI)
    {
        prj_mpi *mpi = prj_mpi_current();

        if (mpi != 0 && mpi->totrank > 1) {
            double *restrict global_ms;

            global_ms = (double *)calloc((size_t)grav_mono->nbins, sizeof(*global_ms));
            if (global_ms == 0) {
                return;
            }
            MPI_Allreduce(grav_mono->ms, global_ms, grav_mono->nbins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            for (idx = 0; idx < grav_mono->nbins; ++idx) {
                grav_mono->ms[idx] = global_ms[idx];
            }
            free(global_ms);
        }
    }
#endif
}

void prj_gravity_monopole_integrate(prj_mesh *mesh)
{
    prj_grav_mono *grav_mono = prj_gravity_active;
    double *restrict enclosed_face;
    int idx;

    (void)mesh;

    if (grav_mono == 0 || grav_mono->nbins <= 0) {
        return;
    }

    enclosed_face = (double *)calloc((size_t)grav_mono->nbins + 1U, sizeof(*enclosed_face));
    if (enclosed_face == 0) {
        return;
    }

    for (idx = 1; idx < grav_mono->nbins; ++idx) {
        grav_mono->ms[idx] += grav_mono->ms[idx - 1];
    }
    enclosed_face[0] = 0.0;
    for (idx = 0; idx < grav_mono->nbins; ++idx) {
        enclosed_face[idx + 1] = grav_mono->ms[idx];
    }

    grav_mono->phi[grav_mono->nbins] =
        -PRJ_GNEWT * enclosed_face[grav_mono->nbins] / prj_gravity_rmax;
    for (idx = grav_mono->nbins - 1; idx >= 0; --idx) {
        double r0 = grav_mono->rf[idx];
        double r1 = grav_mono->rf[idx + 1];
        double f0 = PRJ_GNEWT * enclosed_face[idx] / (r0 * r0);
        double f1 = PRJ_GNEWT * enclosed_face[idx + 1] / (r1 * r1);

        grav_mono->phi[idx] = grav_mono->phi[idx + 1] - 0.5 * (r1 - r0) * (f0 + f1);
    }
    for (idx = 0; idx < grav_mono->nbins; ++idx) {
        grav_mono->accel[idx] =
            -(grav_mono->phi[idx + 1] - grav_mono->phi[idx]) /
            (grav_mono->rf[idx + 1] - grav_mono->rf[idx]);
    }

    free(enclosed_face);
}

int prj_gravity_apply(void)
{
    return 0;
}
