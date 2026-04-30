#include <math.h>
#include <stddef.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "prj.h"

#if defined(PRJ_ENABLE_MPI)
#include <mpi.h>
#endif

#define PRJ_GRAVITY_DEFAULT_NBINS 1024
#define PRJ_GRAVITY_CACHE_INVALID (-1)
#define PRJ_GRAVITY_CACHE_LAST_VALUE 2.0
#define PRJ_GRAVITY_CACHE_SKIP_REDUCE 3.0

static prj_grav_mono *prj_gravity_active = 0;
static double prj_gravity_rmax = 0.0;

static double prj_gravity_min_double(double a, double b)
{
    return a < b ? a : b;
}

static double prj_gravity_radius_center(const prj_grav_mono *grav_mono, int idx)
{
    return 0.5 * (grav_mono->rf[idx] + grav_mono->rf[idx + 1]);
}

static double prj_gravity_abs_double(double a)
{
    return a < 0.0 ? -a : a;
}

static double prj_gravity_clamp_double(double x, double lo, double hi)
{
    if (x < lo) {
        return lo;
    }
    if (x > hi) {
        return hi;
    }
    return x;
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

    if (grav_mono == 0 || grav_mono->nbins <= 0 || r < 0.0 || r >= prj_gravity_rmax) {
        return -1;
    }
    if (r < grav_mono->dr_min) {
        return 0;
    }
    log_span = log(prj_gravity_rmax / grav_mono->dr_min);
    if (log_span <= 0.0 || grav_mono->nbins < 2) {
        return -1;
    }
    idx = 1 + (int)((double)(grav_mono->nbins - 1) * log(r / grav_mono->dr_min) / log_span);
    if (idx < 1) {
        idx = 1;
    }
    if (idx >= grav_mono->nbins) {
        return -1;
    }
    return idx;
}

static void prj_gravity_cache_entry(const prj_grav_mono *grav_mono, double r, int *ridx, double *fr)
{
    int idx;

    if (ridx == 0 || fr == 0) {
        return;
    }
    *ridx = PRJ_GRAVITY_CACHE_INVALID;
    *fr = 0.0;
    if (grav_mono == 0 || grav_mono->nbins <= 0 || grav_mono->rf == 0 || r < 0.0 ||
        prj_gravity_rmax <= 0.0) {
        return;
    }
    if (r >= prj_gravity_rmax) {
        *ridx = grav_mono->nbins - 1;
        *fr = PRJ_GRAVITY_CACHE_SKIP_REDUCE;
        return;
    }

    idx = prj_gravity_bin_index(grav_mono, r);
    if (idx < 0) {
        return;
    }

    *ridx = idx;
    if (grav_mono->nbins == 1) {
        double r0 = prj_gravity_radius_center(grav_mono, 0);

        if (r <= r0 && r0 > 0.0) {
            *fr = -prj_gravity_clamp_double(r / r0, 0.0, 1.0);
        } else {
            *fr = PRJ_GRAVITY_CACHE_LAST_VALUE;
        }
        return;
    }

    if (idx == 0) {
        double r0 = prj_gravity_radius_center(grav_mono, 0);
        double r1 = prj_gravity_radius_center(grav_mono, 1);

        if (r <= r0) {
            *fr = r0 > 0.0 ? -prj_gravity_clamp_double(r / r0, 0.0, 1.0) : 0.0;
        } else {
            *fr = 1.0 + prj_gravity_clamp_double((r - r0) / (r1 - r0), 0.0, 1.0);
        }
        return;
    }

    if (idx == grav_mono->nbins - 1) {
        double r0 = prj_gravity_radius_center(grav_mono, idx - 1);
        double r1 = prj_gravity_radius_center(grav_mono, idx);

        if (r <= r1) {
            *fr = prj_gravity_clamp_double((r - r0) / (r1 - r0), 0.0, 1.0);
        } else {
            *fr = PRJ_GRAVITY_CACHE_LAST_VALUE;
        }
        return;
    }

    {
        double rc = prj_gravity_radius_center(grav_mono, idx);

        if (r <= rc) {
            double rl = prj_gravity_radius_center(grav_mono, idx - 1);

            *fr = prj_gravity_clamp_double((r - rl) / (rc - rl), 0.0, 1.0);
        } else {
            double rr = prj_gravity_radius_center(grav_mono, idx + 1);

            *fr = 1.0 + prj_gravity_clamp_double((r - rc) / (rr - rc), 0.0, 1.0);
        }
    }
}

static void prj_gravity_cache_block_clear(prj_block *block)
{
    int i;

    if (block == 0 || block->ridx == 0 || block->fr == 0) {
        return;
    }
    for (i = 0; i < PRJ_BLOCK_NCELLS; ++i) {
        block->ridx[i] = PRJ_GRAVITY_CACHE_INVALID;
        block->fr[i] = 0.0;
    }
}

static double prj_gravity_min_cell_size(const prj_mesh *mesh)
{
    double min_cell = 0.0;
    int b;

    if (mesh != 0) {
        for (b = 0; b < mesh->nblocks; ++b) {
            const prj_block *block = &mesh->blocks[b];
            double dx;

            if (!prj_gravity_block_is_local_active(block)) {
                continue;
            }
            dx = prj_gravity_min_double(block->dx[0],
                prj_gravity_min_double(block->dx[1], block->dx[2]));
            if (min_cell == 0.0 || dx < min_cell) {
                min_cell = dx;
            }
        }
    }
#if defined(PRJ_ENABLE_MPI)
    {
        double sentinel = min_cell > 0.0 ? min_cell : 1.0e300;
        double global_min = sentinel;

        MPI_Allreduce(&sentinel, &global_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        if (global_min < 1.0e300) {
            min_cell = global_min;
        }
    }
#endif
    return min_cell;
}

static void prj_gravity_build_rf(prj_grav_mono *grav_mono, const prj_mesh *mesh)
{
    double min_cell;
    double log_span;
    int i;

    if (grav_mono == 0 || grav_mono->nbins <= 0 || grav_mono->rf == 0 ||
        prj_gravity_rmax <= 0.0) {
        return;
    }
    min_cell = prj_gravity_min_cell_size(mesh);
    grav_mono->dr_min = 1.5 * min_cell;
    if (grav_mono->dr_min <= 0.0 || prj_gravity_rmax <= grav_mono->dr_min) {
        grav_mono->dr_min = prj_gravity_rmax / (double)grav_mono->nbins;
    }
    log_span = log(prj_gravity_rmax / grav_mono->dr_min);
    grav_mono->rf[0] = 0.0;
    if (grav_mono->nbins == 1) {
        grav_mono->rf[1] = prj_gravity_rmax;
        return;
    }
    for (i = 1; i <= grav_mono->nbins; ++i) {
        grav_mono->rf[i] = grav_mono->dr_min *
            exp(log_span * (double)(i - 1) / (double)(grav_mono->nbins - 1));
    }
}

void prj_gravity_cache_block(prj_block *block)
{
    int i;
    int j;
    int k;

    if (block == 0 || block->ridx == 0 || block->fr == 0) {
        return;
    }
    if (block->id < 0 || block->dx[0] <= 0.0 || block->dx[1] <= 0.0 || block->dx[2] <= 0.0) {
        prj_gravity_cache_block_clear(block);
        return;
    }

    for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
        for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
            for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                double x1 = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
                double x2 = block->xmin[1] + ((double)j + 0.5) * block->dx[1];
                double x3 = block->xmin[2] + ((double)k + 0.5) * block->dx[2];
                double r = sqrt(x1 * x1 + x2 * x2 + x3 * x3);
                int cache_idx = prj_block_cache_index(i, j, k);

                prj_gravity_cache_entry(prj_gravity_active, r, &block->ridx[cache_idx], &block->fr[cache_idx]);
            }
        }
    }
}

void prj_gravity_cache_mesh(prj_mesh *mesh)
{
    int bidx;

    if (mesh == 0) {
        return;
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_gravity_cache_block(&mesh->blocks[bidx]);
    }
}

void prj_gravity_rebuild_grid(prj_sim *sim)
{
    int i;

    if (sim == 0) {
        return;
    }
    prj_gravity_build_rf(&sim->monopole_grav, &sim->mesh);
    for (i = 0; i <= sim->monopole_grav.nbins; ++i) {
        sim->monopole_grav.lapse[i] = 1.0;
    }
    prj_gravity_cache_mesh(&sim->mesh);
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
    free(grav_mono->lapse);
    free(grav_mono->vol);
    free(grav_mono->rho_avg);
    free(grav_mono->vr_avg);
    free(grav_mono->pgas_avg);
    free(grav_mono->uavg_int);
    free(grav_mono->erad_avg);
    free(grav_mono->prad_avg);
    free(grav_mono->vdotF_avg);
    grav_mono->rf = 0;
    grav_mono->ms = 0;
    grav_mono->phi = 0;
    grav_mono->accel = 0;
    grav_mono->lapse = 0;
    grav_mono->vol = 0;
    grav_mono->rho_avg = 0;
    grav_mono->vr_avg = 0;
    grav_mono->pgas_avg = 0;
    grav_mono->uavg_int = 0;
    grav_mono->erad_avg = 0;
    grav_mono->prad_avg = 0;
    grav_mono->vdotF_avg = 0;
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

    grav_mono->nbins = PRJ_GRAVITY_DEFAULT_NBINS;

    grav_mono->rf = (double *)calloc((size_t)grav_mono->nbins + 1U, sizeof(*grav_mono->rf));
    grav_mono->ms = (double *)calloc((size_t)grav_mono->nbins, sizeof(*grav_mono->ms));
    grav_mono->phi = (double *)calloc((size_t)grav_mono->nbins + 1U, sizeof(*grav_mono->phi));
    grav_mono->accel = (double *)calloc((size_t)grav_mono->nbins, sizeof(*grav_mono->accel));
    grav_mono->lapse = (double *)calloc((size_t)grav_mono->nbins + 1U, sizeof(*grav_mono->lapse));
    grav_mono->vol = (double *)calloc((size_t)grav_mono->nbins, sizeof(*grav_mono->vol));
    grav_mono->rho_avg = (double *)calloc((size_t)grav_mono->nbins, sizeof(*grav_mono->rho_avg));
    grav_mono->vr_avg = (double *)calloc((size_t)grav_mono->nbins, sizeof(*grav_mono->vr_avg));
    grav_mono->pgas_avg = (double *)calloc((size_t)grav_mono->nbins, sizeof(*grav_mono->pgas_avg));
    grav_mono->uavg_int = (double *)calloc((size_t)grav_mono->nbins, sizeof(*grav_mono->uavg_int));
    grav_mono->erad_avg = (double *)calloc((size_t)grav_mono->nbins, sizeof(*grav_mono->erad_avg));
    grav_mono->prad_avg = (double *)calloc((size_t)grav_mono->nbins, sizeof(*grav_mono->prad_avg));
    grav_mono->vdotF_avg = (double *)calloc((size_t)grav_mono->nbins, sizeof(*grav_mono->vdotF_avg));
    if (grav_mono->rf == 0 || grav_mono->ms == 0 || grav_mono->phi == 0 || grav_mono->accel == 0 ||
        grav_mono->lapse == 0 || grav_mono->vol == 0 || grav_mono->rho_avg == 0 ||
        grav_mono->vr_avg == 0 || grav_mono->pgas_avg == 0 || grav_mono->uavg_int == 0 ||
        grav_mono->erad_avg == 0 || grav_mono->prad_avg == 0 || grav_mono->vdotF_avg == 0) {
        prj_gravity_free(grav_mono);
        return;
    }

    prj_gravity_build_rf(grav_mono, &sim->mesh);
    for (i = 0; i <= grav_mono->nbins; ++i) {
        grav_mono->lapse[i] = 1.0;
    }

    prj_gravity_active = grav_mono;
    prj_gravity_cache_mesh(&sim->mesh);
}

const prj_grav_mono *prj_gravity_active_monopole(void)
{
    return prj_gravity_active;
}

double prj_gravity_interp_accel(const prj_grav_mono *grav_mono, double r)
{
    int idx;

    if (grav_mono == 0 || grav_mono->nbins <= 0 || grav_mono->accel == 0) {
        return 0.0;
    }
    if (r <= prj_gravity_radius_center(grav_mono, 0)) {
        double r0 = prj_gravity_radius_center(grav_mono, 0);

        if (r0 <= 0.0) {
            return grav_mono->accel[0];
        }
        return grav_mono->accel[0] * (r / r0);
    }
    for (idx = 0; idx < grav_mono->nbins - 1; ++idx) {
        double r0 = prj_gravity_radius_center(grav_mono, idx);
        double r1 = prj_gravity_radius_center(grav_mono, idx + 1);

        if (r <= r1) {
            double weight = (r - r0) / (r1 - r0);

            return (1.0 - weight) * grav_mono->accel[idx] + weight * grav_mono->accel[idx + 1];
        }
    }
    return grav_mono->accel[grav_mono->nbins - 1];
}

double prj_gravity_interp_lapse(const prj_grav_mono *grav_mono, double r)
{
    int idx;

    if (grav_mono == 0 || grav_mono->nbins <= 0 || grav_mono->lapse == 0) {
        return 1.0;
    }
    if (r <= prj_gravity_radius_center(grav_mono, 0)) {
        return grav_mono->lapse[0];
    }
    for (idx = 0; idx < grav_mono->nbins - 1; ++idx) {
        double r0 = prj_gravity_radius_center(grav_mono, idx);
        double r1 = prj_gravity_radius_center(grav_mono, idx + 1);

        if (r <= r1) {
            double weight = (r - r0) / (r1 - r0);

            return (1.0 - weight) * grav_mono->lapse[idx] + weight * grav_mono->lapse[idx + 1];
        }
    }
    return grav_mono->lapse[grav_mono->nbins - 1];
}

double prj_gravity_block_accel_at(const prj_block *block, int i, int j, int k)
{
    const prj_grav_mono *grav_mono = prj_gravity_active;
    int cache_idx;
    int idx;
    double fr;

    if (block == 0 || grav_mono == 0 || grav_mono->nbins <= 0 ||
        grav_mono->accel == 0 || block->ridx == 0 || block->fr == 0) {
        return 0.0;
    }
    cache_idx = prj_block_cache_index(i, j, k);
    idx = block->ridx[cache_idx];
    fr = block->fr[cache_idx];
    if (idx == PRJ_GRAVITY_CACHE_INVALID) {
        return 0.0;
    }
    if (grav_mono->nbins == 1) {
        return fr < 0.0 ? grav_mono->accel[0] * (-fr) : grav_mono->accel[0];
    }
    if (fr < 0.0) {
        return grav_mono->accel[0] * (-fr);
    }
    if (fr <= 1.0) {
        int left = idx > 0 ? idx - 1 : 0;

        return (1.0 - fr) * grav_mono->accel[left] + fr * grav_mono->accel[left + 1];
    }
    if (fr >= PRJ_GRAVITY_CACHE_LAST_VALUE || idx >= grav_mono->nbins - 1) {
        return grav_mono->accel[grav_mono->nbins - 1];
    }
    fr -= 1.0;
    return (1.0 - fr) * grav_mono->accel[idx] + fr * grav_mono->accel[idx + 1];
}

double prj_gravity_block_lapse_at(const prj_block *block, int i, int j, int k)
{
    const prj_grav_mono *grav_mono = prj_gravity_active;
    int cache_idx;
    int idx;
    double fr;

    if (block == 0 || grav_mono == 0 || grav_mono->nbins <= 0 ||
        grav_mono->lapse == 0 || block->ridx == 0 || block->fr == 0) {
        return 1.0;
    }
    cache_idx = prj_block_cache_index(i, j, k);
    idx = block->ridx[cache_idx];
    fr = block->fr[cache_idx];
    if (idx == PRJ_GRAVITY_CACHE_INVALID) {
        return 1.0;
    }
    if (grav_mono->nbins == 1) {
        return grav_mono->lapse[0];
    }
    if (fr < 0.0) {
        return grav_mono->lapse[0];
    }
    if (fr <= 1.0) {
        int left = idx > 0 ? idx - 1 : 0;

        return (1.0 - fr) * grav_mono->lapse[left] + fr * grav_mono->lapse[left + 1];
    }
    if (fr >= PRJ_GRAVITY_CACHE_LAST_VALUE || idx >= grav_mono->nbins - 1) {
        return grav_mono->lapse[grav_mono->nbins - 1];
    }
    fr -= 1.0;
    return (1.0 - fr) * grav_mono->lapse[idx] + fr * grav_mono->lapse[idx + 1];
}

void prj_gravity_monopole_reduce(prj_mesh *mesh, int stage)
{
    prj_grav_mono *grav_mono = prj_gravity_active;
    int bidx;
    int idx;

    if (mesh == 0 || grav_mono == 0 || grav_mono->ms == 0 || grav_mono->vol == 0 ||
        grav_mono->rho_avg == 0 || grav_mono->vr_avg == 0 || grav_mono->pgas_avg == 0 ||
        grav_mono->uavg_int == 0 || grav_mono->erad_avg == 0 || grav_mono->prad_avg == 0 ||
        grav_mono->vdotF_avg == 0) {
        return;
    }

    for (idx = 0; idx < grav_mono->nbins; ++idx) {
        grav_mono->ms[idx] = 0.0;
        grav_mono->vol[idx] = 0.0;
        grav_mono->rho_avg[idx] = 0.0;
        grav_mono->vr_avg[idx] = 0.0;
        grav_mono->pgas_avg[idx] = 0.0;
        grav_mono->uavg_int[idx] = 0.0;
        grav_mono->erad_avg[idx] = 0.0;
        grav_mono->prad_avg[idx] = 0.0;
        grav_mono->vdotF_avg[idx] = 0.0;
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
                    int cache_idx = prj_block_cache_index(i, j, k);
                    double x1 = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
                    double x2 = block->xmin[1] + ((double)j + 0.5) * block->dx[1];
                    double x3 = block->xmin[2] + ((double)k + 0.5) * block->dx[2];
                    double r = sqrt(x1 * x1 + x2 * x2 + x3 * x3);
                    double fr = block->fr != 0 ? block->fr[cache_idx] : 0.0;
                    double rho;
                    double v1;
                    double v2;
                    double v3;
                    double eint;
                    double pgas;
                    double vr;
                    double erad;
                    double prad;
                    double vdotF;
                    double *W = stage == 2 ? block->W1 : block->W;

                    idx = block->ridx != 0 ? block->ridx[cache_idx] : PRJ_GRAVITY_CACHE_INVALID;
                    if (idx >= 0 && fr < PRJ_GRAVITY_CACHE_SKIP_REDUCE) {
                        rho = W[VIDX(PRJ_PRIM_RHO, i, j, k)];
                        v1 = W[VIDX(PRJ_PRIM_V1, i, j, k)];
                        v2 = W[VIDX(PRJ_PRIM_V2, i, j, k)];
                        v3 = W[VIDX(PRJ_PRIM_V3, i, j, k)];
                        eint = W[VIDX(PRJ_PRIM_EINT, i, j, k)];
                        vr = r > 0.0 ? (v1 * x1 + v2 * x2 + v3 * x3) / r : 0.0;
                        pgas = block->eosvar[EIDX(PRJ_EOSVAR_PRESSURE, i, j, k)];
                        erad = 0.0;
                        prad = 0.0;
                        vdotF = 0.0;
#if PRJ_NRAD > 0
                        {
                            int field;
                            int group;

                            for (field = 0; field < PRJ_NRAD; ++field) {
                                for (group = 0; group < PRJ_NEGROUP; ++group) {
                                    double e_rad = W[VIDX(PRJ_PRIM_RAD_E(field, group), i, j, k)];
                                    double f1 = W[VIDX(PRJ_PRIM_RAD_F1(field, group), i, j, k)];
                                    double f2 = W[VIDX(PRJ_PRIM_RAD_F2(field, group), i, j, k)];
                                    double f3 = W[VIDX(PRJ_PRIM_RAD_F3(field, group), i, j, k)];
                                    double fr = r > 0.0 ? (f1 * x1 + f2 * x2 + f3 * x3) / r : 0.0;

                                    erad += e_rad / (PRJ_CLIGHT * PRJ_CLIGHT);
                                    prad += (e_rad / 3.0) / (PRJ_CLIGHT * PRJ_CLIGHT * PRJ_CLIGHT);
                                    vdotF += (v1 * f1 + v2 * f2 + v3 * f3) /
                                        (PRJ_CLIGHT * PRJ_CLIGHT * PRJ_CLIGHT);
                                    (void)fr;
                                }
                            }
                        }
#endif
                        grav_mono->vol[idx] += block->vol;
                        grav_mono->rho_avg[idx] += block->vol * rho;
                        grav_mono->vr_avg[idx] += block->vol * vr;
                        grav_mono->pgas_avg[idx] += block->vol * pgas;
                        grav_mono->uavg_int[idx] += block->vol * rho * eint;
                        grav_mono->erad_avg[idx] += block->vol * erad;
                        grav_mono->prad_avg[idx] += block->vol * prad;
                        grav_mono->vdotF_avg[idx] += block->vol * vdotF;
                        grav_mono->ms[idx] += rho * block->vol;
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
            double *restrict global_vol;
            double *restrict global_rho_avg;
            double *restrict global_vr_avg;
            double *restrict global_pgas_avg;
            double *restrict global_uavg_int;
            double *restrict global_erad_avg;
            double *restrict global_prad_avg;
            double *restrict global_vdotF_avg;

            global_ms = (double *)calloc((size_t)grav_mono->nbins, sizeof(*global_ms));
            global_vol = (double *)calloc((size_t)grav_mono->nbins, sizeof(*global_vol));
            global_rho_avg = (double *)calloc((size_t)grav_mono->nbins, sizeof(*global_rho_avg));
            global_vr_avg = (double *)calloc((size_t)grav_mono->nbins, sizeof(*global_vr_avg));
            global_pgas_avg = (double *)calloc((size_t)grav_mono->nbins, sizeof(*global_pgas_avg));
            global_uavg_int = (double *)calloc((size_t)grav_mono->nbins, sizeof(*global_uavg_int));
            global_erad_avg = (double *)calloc((size_t)grav_mono->nbins, sizeof(*global_erad_avg));
            global_prad_avg = (double *)calloc((size_t)grav_mono->nbins, sizeof(*global_prad_avg));
            global_vdotF_avg = (double *)calloc((size_t)grav_mono->nbins, sizeof(*global_vdotF_avg));
            if (global_ms == 0 || global_vol == 0 || global_rho_avg == 0 || global_vr_avg == 0 ||
                global_pgas_avg == 0 || global_uavg_int == 0 || global_erad_avg == 0 ||
                global_prad_avg == 0 || global_vdotF_avg == 0) {
                free(global_vdotF_avg);
                free(global_prad_avg);
                free(global_erad_avg);
                free(global_uavg_int);
                free(global_pgas_avg);
                free(global_vr_avg);
                free(global_rho_avg);
                free(global_vol);
                free(global_ms);
                return;
            }
            MPI_Allreduce(grav_mono->ms, global_ms, grav_mono->nbins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(grav_mono->vol, global_vol, grav_mono->nbins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(grav_mono->rho_avg, global_rho_avg, grav_mono->nbins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(grav_mono->vr_avg, global_vr_avg, grav_mono->nbins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(grav_mono->pgas_avg, global_pgas_avg, grav_mono->nbins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(grav_mono->uavg_int, global_uavg_int, grav_mono->nbins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(grav_mono->erad_avg, global_erad_avg, grav_mono->nbins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(grav_mono->prad_avg, global_prad_avg, grav_mono->nbins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(grav_mono->vdotF_avg, global_vdotF_avg, grav_mono->nbins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            for (idx = 0; idx < grav_mono->nbins; ++idx) {
                grav_mono->ms[idx] = global_ms[idx];
                grav_mono->vol[idx] = global_vol[idx];
                grav_mono->rho_avg[idx] = global_rho_avg[idx];
                grav_mono->vr_avg[idx] = global_vr_avg[idx];
                grav_mono->pgas_avg[idx] = global_pgas_avg[idx];
                grav_mono->uavg_int[idx] = global_uavg_int[idx];
                grav_mono->erad_avg[idx] = global_erad_avg[idx];
                grav_mono->prad_avg[idx] = global_prad_avg[idx];
                grav_mono->vdotF_avg[idx] = global_vdotF_avg[idx];
            }
            free(global_vdotF_avg);
            free(global_prad_avg);
            free(global_erad_avg);
            free(global_uavg_int);
            free(global_pgas_avg);
            free(global_vr_avg);
            free(global_rho_avg);
            free(global_vol);
            free(global_ms);
        }
    }
#endif

    for (idx = 0; idx < grav_mono->nbins; ++idx) {
        double r0 = grav_mono->rf[idx];
        double r1 = grav_mono->rf[idx + 1];
        double shell_vol = (4.0 / 3.0) * M_PI * (r1 * r1 * r1 - r0 * r0 * r0);

        if (shell_vol > 0.0) {
            grav_mono->rho_avg[idx] /= shell_vol;
            grav_mono->vr_avg[idx] /= shell_vol;
            grav_mono->pgas_avg[idx] /= shell_vol;
            grav_mono->uavg_int[idx] /= shell_vol;
            grav_mono->erad_avg[idx] /= shell_vol;
            grav_mono->prad_avg[idx] /= shell_vol;
            grav_mono->vdotF_avg[idx] /= shell_vol;
        }
    }
}

void prj_gravity_monopole_integrate(prj_mesh *mesh)
{
    prj_grav_mono *grav_mono = prj_gravity_active;
    double *restrict enclosed_face;
    double *restrict baryon_mass_face;
    double *restrict gamma_face;
    int idx;

    (void)mesh;

    if (grav_mono == 0 || grav_mono->nbins <= 0) {
        return;
    }

    enclosed_face = (double *)calloc((size_t)grav_mono->nbins + 1U, sizeof(*enclosed_face));
    baryon_mass_face = (double *)calloc((size_t)grav_mono->nbins + 1U, sizeof(*baryon_mass_face));
    gamma_face = (double *)calloc((size_t)grav_mono->nbins + 1U, sizeof(*gamma_face));
    if (enclosed_face == 0 || baryon_mass_face == 0 || gamma_face == 0) {
        free(gamma_face);
        free(baryon_mass_face);
        free(enclosed_face);
        return;
    }

    for (idx = 1; idx < grav_mono->nbins; ++idx) {
        grav_mono->ms[idx] += grav_mono->ms[idx - 1];
    }
    enclosed_face[0] = 0.0;
    for (idx = 0; idx < grav_mono->nbins; ++idx) {
        enclosed_face[idx + 1] = grav_mono->ms[idx];
    }

#if PRJ_GRAVITY_USE_GR
    {
        double c2 = PRJ_CLIGHT * PRJ_CLIGHT;
        int iter;
        double mlast = 0.0;

        baryon_mass_face[0] = 0.0;
        for (idx = 0; idx < grav_mono->nbins; ++idx) {
            baryon_mass_face[idx + 1] = enclosed_face[idx + 1];
        }

        /* Warm-start gamma from baryon mass only. */
        gamma_face[0] = 1.0;
        for (idx = 0; idx < grav_mono->nbins; ++idx) {
            double r1 = grav_mono->rf[idx + 1];
            double gsq = 1.0 - 2.0 * PRJ_GNEWT * baryon_mass_face[idx + 1] / (r1 * c2);

            if (gsq < 1.0e-12) {
                gsq = 1.0e-12;
            }
            gamma_face[idx + 1] = sqrt(gsq);
        }

        /* Iterate m_TOV and gamma to self-consistency. */
        enclosed_face[0] = 0.0;
        for (iter = 0; iter < 20; ++iter) {
            double m_tov = 0.0;

            for (idx = 0; idx < grav_mono->nbins; ++idx) {
                double r0 = grav_mono->rf[idx];
                double r1 = grav_mono->rf[idx + 1];
                double shell_vol = (4.0 / 3.0) * M_PI * (r1 * r1 * r1 - r0 * r0 * r0);
                double vedge = idx < grav_mono->nbins - 1 ?
                    0.5 * (grav_mono->vr_avg[idx] + grav_mono->vr_avg[idx + 1]) :
                    grav_mono->vr_avg[idx];
                double rho = grav_mono->rho_avg[idx];
                double u = grav_mono->uavg_int[idx];
                double erad = grav_mono->erad_avg[idx];
                double vdF = grav_mono->vdotF_avg[idx];
                double integrand;
                double gamma_avg;
                double gsq;

                /* Cell-centered gamma_avg using current gamma_face. */
                gamma_avg = 0.5 * (gamma_face[idx] + gamma_face[idx + 1]);
                integrand = (rho + u / c2 + erad) * gamma_avg + vdF;
                gsq = 1.0 + (vedge / PRJ_CLIGHT) * (vedge / PRJ_CLIGHT) -
                    2.0 * PRJ_GNEWT * (m_tov + shell_vol * integrand) / (r1 * c2);
                if (gsq < 1.0e-12) {
                    gsq = 1.0e-12;
                }
                gamma_face[idx + 1] = sqrt(gsq);

                /* Re-evaluate integrand with updated gamma_face[idx+1]. */
                gamma_avg = 0.5 * (gamma_face[idx] + gamma_face[idx + 1]);
                integrand = (rho + u / c2 + erad) * gamma_avg + vdF;
                m_tov += shell_vol * integrand;
                enclosed_face[idx + 1] = m_tov;
            }
            if (iter > 0) {
                double denom = m_tov != 0.0 ? prj_gravity_abs_double(m_tov) : 1.0;

                if (prj_gravity_abs_double(mlast - m_tov) / denom < 1.0e-10) {
                    break;
                }
            }
            mlast = m_tov;
        }
    }
    for (idx = 0; idx < grav_mono->nbins; ++idx) {
        double r = grav_mono->rf[idx + 1];
        double pgas_edge = idx < grav_mono->nbins - 1 ?
            0.5 * (grav_mono->pgas_avg[idx] + grav_mono->pgas_avg[idx + 1]) :
            grav_mono->pgas_avg[idx];
        double prad_edge = idx < grav_mono->nbins - 1 ?
            0.5 * (grav_mono->prad_avg[idx] + grav_mono->prad_avg[idx + 1]) :
            grav_mono->prad_avg[idx];
        double rho_edge = idx < grav_mono->nbins - 1 ?
            0.5 * (grav_mono->rho_avg[idx] + grav_mono->rho_avg[idx + 1]) :
            grav_mono->rho_avg[idx];
        double uedge = idx < grav_mono->nbins - 1 ?
            0.5 * (grav_mono->uavg_int[idx] + grav_mono->uavg_int[idx + 1]) :
            grav_mono->uavg_int[idx];
        double numerator;
        double gamma_term;
        double enthalpy_term;

        if (rho_edge <= 0.0) {
            grav_mono->accel[idx] = 0.0;
            continue;
        }
        numerator = -PRJ_GNEWT * (enclosed_face[idx + 1] +
            4.0 * M_PI * r * r * r * (pgas_edge + prad_edge) /
            (PRJ_CLIGHT * PRJ_CLIGHT));
        gamma_term = r * gamma_face[idx + 1];
        if (prj_gravity_abs_double(gamma_term) < 1.0e-30) {
            gamma_term = gamma_term < 0.0 ? -1.0e-30 : 1.0e-30;
        }
        enthalpy_term = (rho_edge + (uedge + pgas_edge) /
            (PRJ_CLIGHT * PRJ_CLIGHT)) / rho_edge;
        grav_mono->accel[idx] = numerator / (gamma_term * gamma_term) * enthalpy_term;
    }
    grav_mono->phi[grav_mono->nbins] =
        -PRJ_GNEWT * baryon_mass_face[grav_mono->nbins] / prj_gravity_rmax;
    for (idx = grav_mono->nbins - 1; idx >= 0; --idx) {
        grav_mono->phi[idx] = grav_mono->phi[idx + 1] + 0.5 *
            (grav_mono->rf[idx + 1] - grav_mono->rf[idx]) *
            (idx == grav_mono->nbins - 1 ? grav_mono->accel[idx] : grav_mono->accel[idx + 1] + grav_mono->accel[idx]);
    }
    for (idx = 0; idx <= grav_mono->nbins; ++idx) {
        grav_mono->lapse[idx] = exp(grav_mono->phi[idx] / (PRJ_CLIGHT * PRJ_CLIGHT));
    }
#else
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
    for (idx = 0; idx <= grav_mono->nbins; ++idx) {
        grav_mono->lapse[idx] = 1.0;
    }
#endif

    free(gamma_face);
    free(baryon_mass_face);
    free(enclosed_face);
}

int prj_gravity_apply(void)
{
    return 0;
}
