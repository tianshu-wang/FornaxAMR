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

static prj_grav *prj_gravity_active = 0;
static double prj_gravity_rmax = 0.0;

static void prj_gravity_fill_mesh_fields(prj_mesh *mesh);

static double prj_gravity_min_double(double a, double b)
{
    return a < b ? a : b;
}

static double prj_gravity_radius_center(const prj_grav *grav, int idx)
{
    return 0.5 * (grav->rf[idx] + grav->rf[idx + 1]);
}

static double prj_gravity_abs_double(double a)
{
    return a < 0.0 ? -a : a;
}

static double prj_gravity_assoc_legendre(int l, int m, double x)
{
    double pmm = 1.0;
    double pmmp1;
    double pll = 0.0;
    int i;
    int ll;

    if (l < 0 || m < 0 || m > l) {
        return 0.0;
    }
    if (x > 1.0) {
        x = 1.0;
    } else if (x < -1.0) {
        x = -1.0;
    }
    if (m > 0) {
        double somx2 = sqrt(prj_gravity_abs_double((1.0 - x) * (1.0 + x)));
        double fact = 1.0;

        for (i = 1; i <= m; ++i) {
            pmm *= -fact * somx2;
            fact += 2.0;
        }
    }
    if (l == m) {
        return pmm;
    }
    pmmp1 = x * (double)(2 * m + 1) * pmm;
    if (l == m + 1) {
        return pmmp1;
    }
    for (ll = m + 2; ll <= l; ++ll) {
        pll = ((double)(2 * ll - 1) * x * pmmp1 -
            (double)(ll + m - 1) * pmm) / (double)(ll - m);
        pmm = pmmp1;
        pmmp1 = pll;
    }
    return pll;
}

static double prj_gravity_ylm_norm(int l, int m)
{
    double ratio = 1.0;
    int n;

    for (n = l - m + 1; n <= l + m; ++n) {
        ratio /= (double)n;
    }
    return sqrt(((double)(2 * l + 1) / (4.0 * M_PI)) * ratio);
}

double prj_gravity_real_spherical_harmonic(int l, int m, double x1, double x2, double x3)
{
    int abs_m = m < 0 ? -m : m;
    double r;
    double mu;
    double phi;
    double plm;
    double ylm;

    if (l < 0 || abs_m > l) {
        return 0.0;
    }
    r = sqrt(x1 * x1 + x2 * x2 + x3 * x3);
    if (r <= 0.0) {
        return l == 0 && m == 0 ? sqrt(1.0 / (4.0 * M_PI)) : 0.0;
    }
    mu = x3 / r;
    if (mu > 1.0) {
        mu = 1.0;
    } else if (mu < -1.0) {
        mu = -1.0;
    }
    phi = atan2(x2, x1);
    plm = prj_gravity_assoc_legendre(l, abs_m, mu);
    ylm = prj_gravity_ylm_norm(l, abs_m) * plm;
    if (m > 0) {
        ylm *= sqrt(2.0) * cos((double)abs_m * phi);
    } else if (m < 0) {
        ylm *= sqrt(2.0) * sin((double)abs_m * phi);
    }
    return ylm;
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

static int prj_gravity_multipole_arrays_valid(const prj_grav *grav)
{
    int yidx;

    if (grav == 0) {
        return 0;
    }
    for (yidx = 0; yidx < LMAX*LMAX; ++yidx) {
        if (grav->Clm[yidx] == 0 || grav->Dlm[yidx] == 0) {
            return 0;
        }
    }
    return 1;
}

static void prj_gravity_zero_multipole_coefficients(prj_grav *grav)
{
    int yidx;

    if (grav == 0) {
        return;
    }
    for (yidx = 0; yidx < LMAX*LMAX; ++yidx) {
        if (grav->Clm[yidx] != 0) {
            prj_fill(grav->Clm[yidx], (size_t)grav->nbins, 0.0);
        }
        if (grav->Dlm[yidx] != 0) {
            prj_fill(grav->Dlm[yidx], (size_t)grav->nbins, 0.0);
        }
    }
}

static int prj_gravity_block_is_local_active(const prj_block *block)
{
    prj_mpi *mpi = prj_mpi_current();

    return block != 0 && block->id >= 0 && block->active == 1 &&
        (mpi == 0 || block->rank == mpi->rank);
}

static int prj_gravity_bin_index(const prj_grav *grav, double r)
{
    double log_span;
    int idx;

    if (grav == 0 || grav->nbins <= 0 || r < 0.0 || r >= prj_gravity_rmax) {
        return -1;
    }
    if (r < grav->dr_min) {
        return 0;
    }
    log_span = log(prj_gravity_rmax / grav->dr_min);
    if (log_span <= 0.0 || grav->nbins < 2) {
        return -1;
    }
    idx = 1 + (int)((double)(grav->nbins - 1) * log(r / grav->dr_min) / log_span);
    if (idx < 1) {
        idx = 1;
    }
    if (idx >= grav->nbins) {
        return -1;
    }
    return idx;
}

static void prj_gravity_cache_entry(const prj_grav *grav, double r, int *ridx, double *fr)
{
    int idx;

    if (ridx == 0 || fr == 0) {
        return;
    }
    *ridx = PRJ_GRAVITY_CACHE_INVALID;
    *fr = 0.0;
    if (grav == 0 || grav->nbins <= 0 || grav->rf == 0 || r < 0.0 ||
        prj_gravity_rmax <= 0.0) {
        return;
    }
    if (r >= prj_gravity_rmax) {
        *ridx = grav->nbins - 1;
        *fr = PRJ_GRAVITY_CACHE_SKIP_REDUCE;
        return;
    }

    idx = prj_gravity_bin_index(grav, r);
    if (idx < 0) {
        return;
    }

    *ridx = idx;
    if (grav->nbins == 1) {
        double r0 = prj_gravity_radius_center(grav, 0);

        if (r <= r0 && r0 > 0.0) {
            *fr = -prj_gravity_clamp_double(r / r0, 0.0, 1.0);
        } else {
            *fr = PRJ_GRAVITY_CACHE_LAST_VALUE;
        }
        return;
    }

    if (idx == 0) {
        double r0 = prj_gravity_radius_center(grav, 0);
        double r1 = prj_gravity_radius_center(grav, 1);

        if (r <= r0) {
            *fr = r0 > 0.0 ? -prj_gravity_clamp_double(r / r0, 0.0, 1.0) : 0.0;
        } else {
            *fr = 1.0 + prj_gravity_clamp_double((r - r0) / (r1 - r0), 0.0, 1.0);
        }
        return;
    }

    if (idx == grav->nbins - 1) {
        double r0 = prj_gravity_radius_center(grav, idx - 1);
        double r1 = prj_gravity_radius_center(grav, idx);

        if (r <= r1) {
            *fr = prj_gravity_clamp_double((r - r0) / (r1 - r0), 0.0, 1.0);
        } else {
            *fr = PRJ_GRAVITY_CACHE_LAST_VALUE;
        }
        return;
    }

    {
        double rc = prj_gravity_radius_center(grav, idx);

        if (r <= rc) {
            double rl = prj_gravity_radius_center(grav, idx - 1);

            *fr = prj_gravity_clamp_double((r - rl) / (rc - rl), 0.0, 1.0);
        } else {
            double rr = prj_gravity_radius_center(grav, idx + 1);

            *fr = 1.0 + prj_gravity_clamp_double((r - rc) / (rr - rc), 0.0, 1.0);
        }
    }
}

static void prj_gravity_cache_block_clear(prj_block *block)
{
    int i;
    int n;

    if (block == 0) {
        return;
    }
    for (i = 0; i < PRJ_BLOCK_NCELLS; ++i) {
        if (block->ridx != 0) {
            block->ridx[i] = PRJ_GRAVITY_CACHE_INVALID;
        }
        if (block->fr != 0) {
            block->fr[i] = 0.0;
        }
        if (block->r_com != 0) {
            block->r_com[i] = 0.0;
        }
        for (n = 0; n < LMAX*LMAX; ++n) {
            if (block->Ylm[n] != 0) {
                block->Ylm[n][i] = 0.0;
            }
        }
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

static void prj_gravity_build_rf(prj_grav *grav, const prj_mesh *mesh)
{
    double min_cell;
    double log_span;
    int i;

    if (grav == 0 || grav->nbins <= 0 || grav->rf == 0 ||
        prj_gravity_rmax <= 0.0) {
        return;
    }
    min_cell = prj_gravity_min_cell_size(mesh);
    grav->dr_min = 1.5 * min_cell;
    if (grav->dr_min <= 0.0 || prj_gravity_rmax <= grav->dr_min) {
        grav->dr_min = prj_gravity_rmax / (double)grav->nbins;
    }
    log_span = log(prj_gravity_rmax / grav->dr_min);
    grav->rf[0] = 0.0;
    if (grav->nbins == 1) {
        grav->rf[1] = prj_gravity_rmax;
        return;
    }
    for (i = 1; i <= grav->nbins; ++i) {
        grav->rf[i] = grav->dr_min *
            exp(log_span * (double)(i - 1) / (double)(grav->nbins - 1));
    }
}

void prj_gravity_cache_block(prj_block *block)
{
    const prj_grav *grav = prj_gravity_active;
    int i;
    int j;
    int k;

    if (block == 0 || block->ridx == 0 || block->fr == 0 || block->r_com == 0) {
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
                double dx1 = x1 - (grav != 0 ? grav->x_com[0] : 0.0);
                double dx2 = x2 - (grav != 0 ? grav->x_com[1] : 0.0);
                double dx3 = x3 - (grav != 0 ? grav->x_com[2] : 0.0);
                int cache_idx = prj_block_cache_index(i, j, k);
                int l;

                block->r_com[cache_idx] = sqrt(dx1 * dx1 + dx2 * dx2 + dx3 * dx3);
                prj_gravity_cache_entry(prj_gravity_active, block->r_com[cache_idx],
                    &block->ridx[cache_idx], &block->fr[cache_idx]);
                for (l = 0; l < LMAX; ++l) {
                    int m;

                    for (m = -l; m <= l; ++m) {
                        int yidx = PRJ_YLM_INDEX(l, m);

                        if (block->Ylm[yidx] != 0) {
                            block->Ylm[yidx][cache_idx] =
                                prj_gravity_real_spherical_harmonic(l, m, dx1, dx2, dx3);
                        }
                    }
                }
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
    prj_gravity_build_rf(&sim->grav, &sim->mesh);
    if (sim->grav.accel != 0) {
        for (i = 0; i < sim->grav.nbins; ++i) {
            sim->grav.accel[i] = 0.0;
        }
    }
    if (sim->grav.lapse != 0) {
        for (i = 0; i <= sim->grav.nbins; ++i) {
            sim->grav.lapse[i] = 1.0;
        }
    }
    prj_gravity_zero_multipole_coefficients(&sim->grav);
    prj_gravity_cache_mesh(&sim->mesh);
    prj_gravity_fill_mesh_fields(&sim->mesh);
}

void prj_gravity_free(prj_grav *grav)
{
    int yidx;

    if (grav == 0) {
        return;
    }

    free(grav->rf);
    free(grav->ms);
    free(grav->phi);
    free(grav->accel);
    free(grav->lapse);
    for (yidx = 0; yidx < LMAX*LMAX; ++yidx) {
        free(grav->Clm[yidx]);
        free(grav->Dlm[yidx]);
    }
    free(grav->vol);
    free(grav->rho_avg);
    free(grav->vr_avg);
    free(grav->pgas_avg);
    free(grav->uavg_int);
    free(grav->erad_avg);
    free(grav->prad_avg);
    free(grav->vdotF_avg);
    grav->rf = 0;
    grav->ms = 0;
    grav->phi = 0;
    grav->accel = 0;
    grav->lapse = 0;
    for (yidx = 0; yidx < LMAX*LMAX; ++yidx) {
        grav->Clm[yidx] = 0;
        grav->Dlm[yidx] = 0;
    }
    grav->vol = 0;
    grav->rho_avg = 0;
    grav->vr_avg = 0;
    grav->pgas_avg = 0;
    grav->uavg_int = 0;
    grav->erad_avg = 0;
    grav->prad_avg = 0;
    grav->vdotF_avg = 0;
    grav->nbins = 0;
    grav->dr_min = 0.0;
    grav->x_com[0] = 0.0;
    grav->x_com[1] = 0.0;
    grav->x_com[2] = 0.0;
    grav->x_com_new[0] = 0.0;
    grav->x_com_new[1] = 0.0;
    grav->x_com_new[2] = 0.0;
    if (prj_gravity_active == grav) {
        prj_gravity_active = 0;
        prj_gravity_rmax = 0.0;
    }
}

void prj_gravity_init(prj_sim *sim)
{
    prj_grav *grav;
    double span1;
    double span2;
    double span3;
    int i;
    int yidx;

    if (sim == 0) {
        return;
    }

    grav = &sim->grav;
    prj_gravity_free(grav);

    span1 = fabs(sim->mesh.coord.x1max - sim->mesh.coord.x1min);
    span2 = fabs(sim->mesh.coord.x2max - sim->mesh.coord.x2min);
    span3 = fabs(sim->mesh.coord.x3max - sim->mesh.coord.x3min);
    prj_gravity_rmax = 0.5 * prj_gravity_min_double(span1, prj_gravity_min_double(span2, span3));

    grav->nbins = PRJ_GRAVITY_DEFAULT_NBINS;
    for (i = 0; i < 3; ++i) {
        grav->x_com[i] = 0.0;
        grav->x_com_new[i] = 0.0;
    }

    grav->rf = (double *)calloc((size_t)grav->nbins + 1U, sizeof(*grav->rf));
    grav->ms = (double *)calloc((size_t)grav->nbins, sizeof(*grav->ms));
    grav->phi = (double *)calloc((size_t)grav->nbins + 1U, sizeof(*grav->phi));
    grav->accel = (double *)calloc((size_t)grav->nbins, sizeof(*grav->accel));
    grav->lapse = (double *)calloc((size_t)grav->nbins + 1U, sizeof(*grav->lapse));
    for (yidx = 0; yidx < LMAX*LMAX; ++yidx) {
        grav->Clm[yidx] = (double *)calloc((size_t)grav->nbins, sizeof(*grav->Clm[yidx]));
        grav->Dlm[yidx] = (double *)calloc((size_t)grav->nbins, sizeof(*grav->Dlm[yidx]));
    }
    grav->vol = (double *)calloc((size_t)grav->nbins, sizeof(*grav->vol));
    grav->rho_avg = (double *)calloc((size_t)grav->nbins, sizeof(*grav->rho_avg));
    grav->vr_avg = (double *)calloc((size_t)grav->nbins, sizeof(*grav->vr_avg));
    grav->pgas_avg = (double *)calloc((size_t)grav->nbins, sizeof(*grav->pgas_avg));
    grav->uavg_int = (double *)calloc((size_t)grav->nbins, sizeof(*grav->uavg_int));
    grav->erad_avg = (double *)calloc((size_t)grav->nbins, sizeof(*grav->erad_avg));
    grav->prad_avg = (double *)calloc((size_t)grav->nbins, sizeof(*grav->prad_avg));
    grav->vdotF_avg = (double *)calloc((size_t)grav->nbins, sizeof(*grav->vdotF_avg));
    if (grav->rf == 0 || grav->ms == 0 || grav->phi == 0 || grav->accel == 0 ||
        grav->lapse == 0 || !prj_gravity_multipole_arrays_valid(grav) ||
        grav->vol == 0 || grav->rho_avg == 0 ||
        grav->vr_avg == 0 || grav->pgas_avg == 0 || grav->uavg_int == 0 ||
        grav->erad_avg == 0 || grav->prad_avg == 0 || grav->vdotF_avg == 0) {
        prj_gravity_free(grav);
        return;
    }

    prj_gravity_build_rf(grav, &sim->mesh);
    for (i = 0; i <= grav->nbins; ++i) {
        grav->lapse[i] = 1.0;
    }

    prj_gravity_active = grav;
    prj_gravity_cache_mesh(&sim->mesh);
    prj_gravity_fill_mesh_fields(&sim->mesh);
}

const prj_grav *prj_gravity_active_monopole(void)
{
    return prj_gravity_active;
}

double prj_gravity_interp_accel(const prj_grav *grav, double r)
{
    int idx;

    if (grav == 0 || grav->nbins <= 0 || grav->accel == 0) {
        return 0.0;
    }
    if (r <= prj_gravity_radius_center(grav, 0)) {
        double r0 = prj_gravity_radius_center(grav, 0);

        if (r0 <= 0.0) {
            return grav->accel[0];
        }
        return grav->accel[0] * (r / r0);
    }
    for (idx = 0; idx < grav->nbins - 1; ++idx) {
        double r0 = prj_gravity_radius_center(grav, idx);
        double r1 = prj_gravity_radius_center(grav, idx + 1);

        if (r <= r1) {
            double weight = (r - r0) / (r1 - r0);

            return (1.0 - weight) * grav->accel[idx] + weight * grav->accel[idx + 1];
        }
    }
    return grav->accel[grav->nbins - 1];
}

double prj_gravity_interp_lapse(const prj_grav *grav, double r)
{
    int idx;

    if (grav == 0 || grav->nbins <= 0 || grav->lapse == 0) {
        return 1.0;
    }
    if (r <= prj_gravity_radius_center(grav, 0)) {
        return grav->lapse[0];
    }
    for (idx = 0; idx < grav->nbins - 1; ++idx) {
        double r0 = prj_gravity_radius_center(grav, idx);
        double r1 = prj_gravity_radius_center(grav, idx + 1);

        if (r <= r1) {
            double weight = (r - r0) / (r1 - r0);

            return (1.0 - weight) * grav->lapse[idx] + weight * grav->lapse[idx + 1];
        }
    }
    return grav->lapse[grav->nbins - 1];
}

static double prj_gravity_block_cached_accel_at(const prj_block *block, int i, int j, int k)
{
    const prj_grav *grav = prj_gravity_active;
    int cache_idx;
    int idx;
    double fr;

    if (block == 0 || grav == 0 || grav->nbins <= 0 ||
        grav->accel == 0 || block->ridx == 0 || block->fr == 0) {
        return 0.0;
    }
    cache_idx = prj_block_cache_index(i, j, k);
    idx = block->ridx[cache_idx];
    fr = block->fr[cache_idx];
    if (idx == PRJ_GRAVITY_CACHE_INVALID) {
        return 0.0;
    }
    if (grav->nbins == 1) {
        return fr < 0.0 ? grav->accel[0] * (-fr) : grav->accel[0];
    }
    if (fr < 0.0) {
        return grav->accel[0] * (-fr);
    }
    if (fr <= 1.0) {
        int left = idx > 0 ? idx - 1 : 0;

        return (1.0 - fr) * grav->accel[left] + fr * grav->accel[left + 1];
    }
    if (fr >= PRJ_GRAVITY_CACHE_LAST_VALUE || idx >= grav->nbins - 1) {
        return grav->accel[grav->nbins - 1];
    }
    fr -= 1.0;
    return (1.0 - fr) * grav->accel[idx] + fr * grav->accel[idx + 1];
}

static double prj_gravity_block_cached_lapse_at(const prj_block *block, int i, int j, int k)
{
    const prj_grav *grav = prj_gravity_active;
    int cache_idx;
    int idx;
    double fr;

    if (block == 0 || grav == 0 || grav->nbins <= 0 ||
        grav->lapse == 0 || block->ridx == 0 || block->fr == 0) {
        return 1.0;
    }
    cache_idx = prj_block_cache_index(i, j, k);
    idx = block->ridx[cache_idx];
    fr = block->fr[cache_idx];
    if (idx == PRJ_GRAVITY_CACHE_INVALID) {
        return 1.0;
    }
    if (grav->nbins == 1) {
        return grav->lapse[0];
    }
    if (fr < 0.0) {
        return grav->lapse[0];
    }
    if (fr <= 1.0) {
        int left = idx > 0 ? idx - 1 : 0;

        return (1.0 - fr) * grav->lapse[left] + fr * grav->lapse[left + 1];
    }
    if (fr >= PRJ_GRAVITY_CACHE_LAST_VALUE || idx >= grav->nbins - 1) {
        return grav->lapse[grav->nbins - 1];
    }
    fr -= 1.0;
    return (1.0 - fr) * grav->lapse[idx] + fr * grav->lapse[idx + 1];
}

static double prj_gravity_pow_int(double x, int n)
{
    double value = 1.0;
    int i;

    for (i = 0; i < n; ++i) {
        value *= x;
    }
    return value;
}

static double prj_gravity_block_multipole_phi_at(const prj_block *block, int i, int j, int k)
{
    const prj_grav *grav = prj_gravity_active;
    int cache_idx;
    int idx;
    double r;
    double fr;
    double phi = 0.0;
    int l;

    if (block == 0 || grav == 0 || grav->nbins <= 0 || block->ridx == 0 ||
        block->fr == 0 || block->r_com == 0) {
        return 0.0;
    }
    cache_idx = prj_block_cache_index(i, j, k);
    idx = block->ridx[cache_idx];
    if (idx == PRJ_GRAVITY_CACHE_INVALID) {
        return 0.0;
    }
    r = block->r_com[cache_idx];
    if (r <= 0.0) {
        return 0.0;
    }
    fr = block->fr[cache_idx];
    for (l = 1; l < LMAX; ++l) {
        double r_l = prj_gravity_pow_int(r, l);
        double inv_r_l1 = 1.0 / (r_l * r);
        int cidx = fr < PRJ_GRAVITY_CACHE_SKIP_REDUCE ? idx - 1 : idx;
        int m;

        for (m = -l; m <= l; ++m) {
            int yidx = PRJ_YLM_INDEX(l, m);
            double ylm = block->Ylm[yidx] != 0 ? block->Ylm[yidx][cache_idx] : 0.0;
            double coeff = cidx >= 0 ? grav->Clm[yidx][cidx] * inv_r_l1 : 0.0;

            if (fr < PRJ_GRAVITY_CACHE_SKIP_REDUCE) {
                coeff += grav->Dlm[yidx][idx] * r_l;
            }
            phi += coeff * ylm;
        }
    }
    return -PRJ_GNEWT * phi;
}

static double prj_gravity_phi_axis_gradient(const double *Phi, int axis, int i, int j, int k,
    const double dx[3])
{
    int lo = -PRJ_NGHOST;
    int hi = PRJ_BLOCK_SIZE + PRJ_NGHOST - 1;

    if (Phi == 0 || dx[axis] <= 0.0) {
        return 0.0;
    }
    if (axis == 0) {
        if (i <= lo) {
            return (Phi[IDX(i + 1, j, k)] - Phi[IDX(i, j, k)]) / dx[0];
        }
        if (i >= hi) {
            return (Phi[IDX(i, j, k)] - Phi[IDX(i - 1, j, k)]) / dx[0];
        }
        return 0.5 * (Phi[IDX(i + 1, j, k)] - Phi[IDX(i - 1, j, k)]) / dx[0];
    }
    if (axis == 1) {
        if (j <= lo) {
            return (Phi[IDX(i, j + 1, k)] - Phi[IDX(i, j, k)]) / dx[1];
        }
        if (j >= hi) {
            return (Phi[IDX(i, j, k)] - Phi[IDX(i, j - 1, k)]) / dx[1];
        }
        return 0.5 * (Phi[IDX(i, j + 1, k)] - Phi[IDX(i, j - 1, k)]) / dx[1];
    }
    if (k <= lo) {
        return (Phi[IDX(i, j, k + 1)] - Phi[IDX(i, j, k)]) / dx[2];
    }
    if (k >= hi) {
        return (Phi[IDX(i, j, k)] - Phi[IDX(i, j, k - 1)]) / dx[2];
    }
    return 0.5 * (Phi[IDX(i, j, k + 1)] - Phi[IDX(i, j, k - 1)]) / dx[2];
}

static void prj_gravity_add_block_multipole_fields(prj_block *block)
{
    double Phi[PRJ_BLOCK_NCELLS];
    int i;
    int j;
    int k;

    if (block == 0 || block->grav[0] == 0 || block->grav[1] == 0 || block->grav[2] == 0 ||
        LMAX <= 1 || !prj_gravity_multipole_arrays_valid(prj_gravity_active)) {
        return;
    }
    for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
        for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
            for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                Phi[IDX(i, j, k)] = prj_gravity_block_multipole_phi_at(block, i, j, k);
            }
        }
    }
    for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
        for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
            for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                int cache_idx = prj_block_cache_index(i, j, k);

                block->grav[0][cache_idx] -=
                    prj_gravity_phi_axis_gradient(Phi, 0, i, j, k, block->dx);
                block->grav[1][cache_idx] -=
                    prj_gravity_phi_axis_gradient(Phi, 1, i, j, k, block->dx);
                block->grav[2][cache_idx] -=
                    prj_gravity_phi_axis_gradient(Phi, 2, i, j, k, block->dx);
            }
        }
    }
}

double prj_gravity_block_accel_at(const prj_block *block, int i, int j, int k)
{
    if (block != 0 && block->grav[0] != 0 && block->grav[1] != 0 && block->grav[2] != 0) {
        const prj_grav *grav = prj_gravity_active;
        int cache_idx = prj_block_cache_index(i, j, k);
        double x1 = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
        double x2 = block->xmin[1] + ((double)j + 0.5) * block->dx[1];
        double x3 = block->xmin[2] + ((double)k + 0.5) * block->dx[2];
        double dx1 = x1 - (grav != 0 ? grav->x_com[0] : 0.0);
        double dx2 = x2 - (grav != 0 ? grav->x_com[1] : 0.0);
        double dx3 = x3 - (grav != 0 ? grav->x_com[2] : 0.0);
        double r = block->r_com != 0 ? block->r_com[cache_idx] :
            sqrt(dx1 * dx1 + dx2 * dx2 + dx3 * dx3);

        if (r > 0.0) {
            return (block->grav[0][cache_idx] * dx1 +
                block->grav[1][cache_idx] * dx2 +
                block->grav[2][cache_idx] * dx3) / r;
        }
        return 0.0;
    }
    return prj_gravity_block_cached_accel_at(block, i, j, k);
}

double prj_gravity_block_lapse_at(const prj_block *block, int i, int j, int k)
{
    if (block != 0 && block->lapse != 0) {
        return block->lapse[prj_block_cache_index(i, j, k)];
    }
    return prj_gravity_block_cached_lapse_at(block, i, j, k);
}

static void prj_gravity_fill_block_field_defaults(prj_block *block)
{
    int d;

    if (block == 0) {
        return;
    }
    prj_fill(block->lapse, (size_t)PRJ_BLOCK_NCELLS, 1.0);
    prj_fill(block->r_com, (size_t)PRJ_BLOCK_NCELLS, 0.0);
    for (d = 0; d < 3; ++d) {
        prj_fill(block->grav[d], (size_t)PRJ_BLOCK_NCELLS, 0.0);
    }
}

static void prj_gravity_fill_block_fields(prj_block *block)
{
    int i;
    int j;
    int k;

    if (block == 0 || block->lapse == 0 || block->r_com == 0 ||
        block->grav[0] == 0 || block->grav[1] == 0 || block->grav[2] == 0) {
        return;
    }
    if (block->id < 0 || block->dx[0] <= 0.0 || block->dx[1] <= 0.0 ||
        block->dx[2] <= 0.0 || prj_gravity_active == 0) {
        prj_gravity_fill_block_field_defaults(block);
        return;
    }

    for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
        for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
            for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                double x1 = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
                double x2 = block->xmin[1] + ((double)j + 0.5) * block->dx[1];
                double x3 = block->xmin[2] + ((double)k + 0.5) * block->dx[2];
                double dx1 = x1 - prj_gravity_active->x_com[0];
                double dx2 = x2 - prj_gravity_active->x_com[1];
                double dx3 = x3 - prj_gravity_active->x_com[2];
                int cache_idx = prj_block_cache_index(i, j, k);
                double r = block->r_com[cache_idx];
                double accel = prj_gravity_block_cached_accel_at(block, i, j, k);

                block->lapse[cache_idx] = prj_gravity_block_cached_lapse_at(block, i, j, k);
                if (r > 0.0) {
                    double accel_over_r = accel / r;

                    block->grav[0][cache_idx] = accel_over_r * dx1;
                    block->grav[1][cache_idx] = accel_over_r * dx2;
                    block->grav[2][cache_idx] = accel_over_r * dx3;
                } else {
                    block->grav[0][cache_idx] = 0.0;
                    block->grav[1][cache_idx] = 0.0;
                    block->grav[2][cache_idx] = 0.0;
                }
            }
        }
    }
    prj_gravity_add_block_multipole_fields(block);
}

static void prj_gravity_fill_mesh_fields(prj_mesh *mesh)
{
    int bidx;

    if (mesh == 0) {
        return;
    }
    PRJ_TIMER_CURRENT_START("gravity_fill_block_fields");
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_gravity_fill_block_fields(&mesh->blocks[bidx]);
    }
    PRJ_TIMER_CURRENT_STOP("gravity_fill_block_fields");
}

int prj_gravity_update_center_of_mass(prj_mesh *mesh, double x_com_err_tol)
{
    prj_grav *grav = prj_gravity_active;
    double local[4] = {0.0, 0.0, 0.0, 0.0};
    double global[4] = {0.0, 0.0, 0.0, 0.0};
    double dx0;
    double dx1;
    double dx2;
    double distance;
    double min_cell;
    double threshold;
    int bidx;
    int d;

    if (mesh == 0 || grav == 0) {
        return 0;
    }

    PRJ_TIMER_CURRENT_START("gravity_x_com_reduce_local");
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        const prj_block *block = &mesh->blocks[bidx];
        int i;
        int j;
        int k;

        if (!prj_gravity_block_is_local_active(block) || block->W == 0) {
            continue;
        }
        for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
            for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    double rho = block->W[VIDX(PRJ_PRIM_RHO, i, j, k)];
                    double dm = rho * block->vol;
                    double x1 = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
                    double x2 = block->xmin[1] + ((double)j + 0.5) * block->dx[1];
                    double x3 = block->xmin[2] + ((double)k + 0.5) * block->dx[2];

                    local[0] += dm;
                    local[1] += dm * x1;
                    local[2] += dm * x2;
                    local[3] += dm * x3;
                }
            }
        }
    }
    PRJ_TIMER_CURRENT_STOP("gravity_x_com_reduce_local");

#if defined(PRJ_ENABLE_MPI)
    MPI_Allreduce(local, global, 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
    for (d = 0; d < 4; ++d) {
        global[d] = local[d];
    }
#endif

    if (global[0] > 0.0) {
        grav->x_com_new[0] = global[1] / global[0];
        grav->x_com_new[1] = global[2] / global[0];
        grav->x_com_new[2] = global[3] / global[0];
    } else {
        for (d = 0; d < 3; ++d) {
            grav->x_com_new[d] = grav->x_com[d];
        }
    }

    dx0 = grav->x_com_new[0] - grav->x_com[0];
    dx1 = grav->x_com_new[1] - grav->x_com[1];
    dx2 = grav->x_com_new[2] - grav->x_com[2];
    distance = sqrt(dx0 * dx0 + dx1 * dx1 + dx2 * dx2);
    min_cell = prj_gravity_min_cell_size(mesh);
    threshold = x_com_err_tol * min_cell;
    if (threshold < 0.0) {
        threshold = 0.0;
    }

    if (min_cell <= 0.0 || distance <= threshold) {
        return 0;
    }

    for (d = 0; d < 3; ++d) {
        grav->x_com[d] = grav->x_com_new[d];
    }
    prj_gravity_cache_mesh(mesh);
    prj_gravity_monopole_reduce(mesh, 1);
    prj_gravity_monopole_integrate(mesh);
    return 1;
}

void prj_gravity_monopole_reduce(prj_mesh *mesh, int stage)
{
    prj_grav *grav = prj_gravity_active;
    int bidx;
    int idx;

    if (mesh == 0 || grav == 0 || grav->ms == 0 || grav->vol == 0 ||
        grav->rho_avg == 0 || grav->vr_avg == 0 || grav->pgas_avg == 0 ||
        grav->uavg_int == 0 || grav->erad_avg == 0 || grav->prad_avg == 0 ||
        grav->vdotF_avg == 0 || !prj_gravity_multipole_arrays_valid(grav)) {
        return;
    }

    PRJ_TIMER_CURRENT_START("gravity_reduce_zero");
    for (idx = 0; idx < grav->nbins; ++idx) {
        grav->ms[idx] = 0.0;
        grav->vol[idx] = 0.0;
        grav->rho_avg[idx] = 0.0;
        grav->vr_avg[idx] = 0.0;
        grav->pgas_avg[idx] = 0.0;
        grav->uavg_int[idx] = 0.0;
        grav->erad_avg[idx] = 0.0;
        grav->prad_avg[idx] = 0.0;
        grav->vdotF_avg[idx] = 0.0;
    }
    prj_gravity_zero_multipole_coefficients(grav);
    PRJ_TIMER_CURRENT_STOP("gravity_reduce_zero");

    PRJ_TIMER_CURRENT_START("gravity_reduce_local");
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
                    double dx1 = x1 - grav->x_com[0];
                    double dx2 = x2 - grav->x_com[1];
                    double dx3 = x3 - grav->x_com[2];
                    double r = block->r_com != 0 ? block->r_com[cache_idx] :
                        sqrt(dx1 * dx1 + dx2 * dx2 + dx3 * dx3);
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
                        vr = r > 0.0 ? (v1 * dx1 + v2 * dx2 + v3 * dx3) / r : 0.0;
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
                                    double fr = r > 0.0 ? (f1 * dx1 + f2 * dx2 + f3 * dx3) / r : 0.0;

                                    erad += e_rad / (PRJ_CLIGHT * PRJ_CLIGHT);
                                    prad += (e_rad / 3.0) / (PRJ_CLIGHT * PRJ_CLIGHT * PRJ_CLIGHT);
                                    vdotF += (v1 * f1 + v2 * f2 + v3 * f3) /
                                        (PRJ_CLIGHT * PRJ_CLIGHT * PRJ_CLIGHT);
                                    (void)fr;
                                }
                            }
                        }
#endif
                        grav->vol[idx] += block->vol;
                        grav->rho_avg[idx] += block->vol * rho;
                        grav->vr_avg[idx] += block->vol * vr;
                        grav->pgas_avg[idx] += block->vol * pgas;
                        grav->uavg_int[idx] += block->vol * rho * eint;
                        grav->erad_avg[idx] += block->vol * erad;
                        grav->prad_avg[idx] += block->vol * prad;
                        grav->vdotF_avg[idx] += block->vol * vdotF;
                        grav->ms[idx] += rho * block->vol;
                        if (r > 0.0) {
                            int l;

                            for (l = 0; l < LMAX; ++l) {
                                double r_l = prj_gravity_pow_int(r, l);
                                double inv_r_l1 = 1.0 / (r_l * r);
                                double ylm_scale = (4.0 * M_PI) / (double)(2 * l + 1);
                                int m;

                                for (m = -l; m <= l; ++m) {
                                    int yidx = PRJ_YLM_INDEX(l, m);
                                    double ylm = block->Ylm[yidx] != 0 ?
                                        block->Ylm[yidx][cache_idx] :
                                        prj_gravity_real_spherical_harmonic(l, m, dx1, dx2, dx3);
                                    double coeff = ylm_scale * rho * block->vol * ylm;

                                    grav->Clm[yidx][idx] += coeff * r_l;
                                    grav->Dlm[yidx][idx] += coeff * inv_r_l1;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    PRJ_TIMER_CURRENT_STOP("gravity_reduce_local");

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
            double *restrict global_lm;
            int yidx;

            PRJ_TIMER_CURRENT_START("gravity_reduce_mpi_allreduce");
            global_ms = (double *)calloc((size_t)grav->nbins, sizeof(*global_ms));
            global_vol = (double *)calloc((size_t)grav->nbins, sizeof(*global_vol));
            global_rho_avg = (double *)calloc((size_t)grav->nbins, sizeof(*global_rho_avg));
            global_vr_avg = (double *)calloc((size_t)grav->nbins, sizeof(*global_vr_avg));
            global_pgas_avg = (double *)calloc((size_t)grav->nbins, sizeof(*global_pgas_avg));
            global_uavg_int = (double *)calloc((size_t)grav->nbins, sizeof(*global_uavg_int));
            global_erad_avg = (double *)calloc((size_t)grav->nbins, sizeof(*global_erad_avg));
            global_prad_avg = (double *)calloc((size_t)grav->nbins, sizeof(*global_prad_avg));
            global_vdotF_avg = (double *)calloc((size_t)grav->nbins, sizeof(*global_vdotF_avg));
            global_lm = (double *)calloc((size_t)grav->nbins, sizeof(*global_lm));
            if (global_ms == 0 || global_vol == 0 || global_rho_avg == 0 || global_vr_avg == 0 ||
                global_pgas_avg == 0 || global_uavg_int == 0 || global_erad_avg == 0 ||
                global_prad_avg == 0 || global_vdotF_avg == 0 || global_lm == 0) {
                free(global_lm);
                free(global_vdotF_avg);
                free(global_prad_avg);
                free(global_erad_avg);
                free(global_uavg_int);
                free(global_pgas_avg);
                free(global_vr_avg);
                free(global_rho_avg);
                free(global_vol);
                free(global_ms);
                PRJ_TIMER_CURRENT_STOP("gravity_reduce_mpi_allreduce");
                return;
            }
            MPI_Allreduce(grav->ms, global_ms, grav->nbins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(grav->vol, global_vol, grav->nbins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(grav->rho_avg, global_rho_avg, grav->nbins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(grav->vr_avg, global_vr_avg, grav->nbins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(grav->pgas_avg, global_pgas_avg, grav->nbins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(grav->uavg_int, global_uavg_int, grav->nbins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(grav->erad_avg, global_erad_avg, grav->nbins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(grav->prad_avg, global_prad_avg, grav->nbins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(grav->vdotF_avg, global_vdotF_avg, grav->nbins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            for (idx = 0; idx < grav->nbins; ++idx) {
                grav->ms[idx] = global_ms[idx];
                grav->vol[idx] = global_vol[idx];
                grav->rho_avg[idx] = global_rho_avg[idx];
                grav->vr_avg[idx] = global_vr_avg[idx];
                grav->pgas_avg[idx] = global_pgas_avg[idx];
                grav->uavg_int[idx] = global_uavg_int[idx];
                grav->erad_avg[idx] = global_erad_avg[idx];
                grav->prad_avg[idx] = global_prad_avg[idx];
                grav->vdotF_avg[idx] = global_vdotF_avg[idx];
            }
            for (yidx = 0; yidx < LMAX*LMAX; ++yidx) {
                MPI_Allreduce(grav->Clm[yidx], global_lm, grav->nbins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                for (idx = 0; idx < grav->nbins; ++idx) {
                    grav->Clm[yidx][idx] = global_lm[idx];
                }
                MPI_Allreduce(grav->Dlm[yidx], global_lm, grav->nbins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                for (idx = 0; idx < grav->nbins; ++idx) {
                    grav->Dlm[yidx][idx] = global_lm[idx];
                }
            }
            free(global_lm);
            free(global_vdotF_avg);
            free(global_prad_avg);
            free(global_erad_avg);
            free(global_uavg_int);
            free(global_pgas_avg);
            free(global_vr_avg);
            free(global_rho_avg);
            free(global_vol);
            free(global_ms);
            PRJ_TIMER_CURRENT_STOP("gravity_reduce_mpi_allreduce");
        }
    }
#endif

    PRJ_TIMER_CURRENT_START("gravity_reduce_multipole_cumsum");
    {
        int yidx;

        for (yidx = 0; yidx < LMAX*LMAX; ++yidx) {
            for (idx = 1; idx < grav->nbins; ++idx) {
                grav->Clm[yidx][idx] += grav->Clm[yidx][idx - 1];
            }
            for (idx = grav->nbins - 2; idx >= 0; --idx) {
                grav->Dlm[yidx][idx] += grav->Dlm[yidx][idx + 1];
            }
        }
    }
    PRJ_TIMER_CURRENT_STOP("gravity_reduce_multipole_cumsum");

    PRJ_TIMER_CURRENT_START("gravity_reduce_normalize");
    for (idx = 0; idx < grav->nbins; ++idx) {
        double r0 = grav->rf[idx];
        double r1 = grav->rf[idx + 1];
        double shell_vol = (4.0 / 3.0) * M_PI * (r1 * r1 * r1 - r0 * r0 * r0);

        if (shell_vol > 0.0) {
            grav->rho_avg[idx] /= shell_vol;
            grav->vr_avg[idx] /= shell_vol;
            grav->pgas_avg[idx] /= shell_vol;
            grav->uavg_int[idx] /= shell_vol;
            grav->erad_avg[idx] /= shell_vol;
            grav->prad_avg[idx] /= shell_vol;
            grav->vdotF_avg[idx] /= shell_vol;
        }
    }
    PRJ_TIMER_CURRENT_STOP("gravity_reduce_normalize");
}

void prj_gravity_monopole_integrate(prj_mesh *mesh)
{
    prj_grav *grav = prj_gravity_active;
    double *restrict enclosed_face;
    double *restrict baryon_mass_face;
    double *restrict gamma_face;
    int idx;

    if (grav == 0 || grav->nbins <= 0) {
        return;
    }

    PRJ_TIMER_CURRENT_START("gravity_integrate");
    enclosed_face = (double *)calloc((size_t)grav->nbins + 1U, sizeof(*enclosed_face));
    baryon_mass_face = (double *)calloc((size_t)grav->nbins + 1U, sizeof(*baryon_mass_face));
    gamma_face = (double *)calloc((size_t)grav->nbins + 1U, sizeof(*gamma_face));
    if (enclosed_face == 0 || baryon_mass_face == 0 || gamma_face == 0) {
        free(gamma_face);
        free(baryon_mass_face);
        free(enclosed_face);
        PRJ_TIMER_CURRENT_STOP("gravity_integrate");
        return;
    }

    for (idx = 1; idx < grav->nbins; ++idx) {
        grav->ms[idx] += grav->ms[idx - 1];
    }
    enclosed_face[0] = 0.0;
    for (idx = 0; idx < grav->nbins; ++idx) {
        enclosed_face[idx + 1] = grav->ms[idx];
    }

#if PRJ_GRAVITY_USE_GR
    {
        double c2 = PRJ_CLIGHT * PRJ_CLIGHT;
        int iter;
        double mlast = 0.0;

        baryon_mass_face[0] = 0.0;
        for (idx = 0; idx < grav->nbins; ++idx) {
            baryon_mass_face[idx + 1] = enclosed_face[idx + 1];
        }

        /* Warm-start gamma from baryon mass only. */
        gamma_face[0] = 1.0;
        for (idx = 0; idx < grav->nbins; ++idx) {
            double r1 = grav->rf[idx + 1];
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

            for (idx = 0; idx < grav->nbins; ++idx) {
                double r0 = grav->rf[idx];
                double r1 = grav->rf[idx + 1];
                double shell_vol = (4.0 / 3.0) * M_PI * (r1 * r1 * r1 - r0 * r0 * r0);
                double vedge = idx < grav->nbins - 1 ?
                    0.5 * (grav->vr_avg[idx] + grav->vr_avg[idx + 1]) :
                    grav->vr_avg[idx];
                double rho = grav->rho_avg[idx];
                double u = grav->uavg_int[idx];
                double erad = grav->erad_avg[idx];
                double vdF = grav->vdotF_avg[idx];
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
    for (idx = 0; idx < grav->nbins; ++idx) {
        double r = grav->rf[idx + 1];
        double pgas_edge = idx < grav->nbins - 1 ?
            0.5 * (grav->pgas_avg[idx] + grav->pgas_avg[idx + 1]) :
            grav->pgas_avg[idx];
        double prad_edge = idx < grav->nbins - 1 ?
            0.5 * (grav->prad_avg[idx] + grav->prad_avg[idx + 1]) :
            grav->prad_avg[idx];
        double rho_edge = idx < grav->nbins - 1 ?
            0.5 * (grav->rho_avg[idx] + grav->rho_avg[idx + 1]) :
            grav->rho_avg[idx];
        double uedge = idx < grav->nbins - 1 ?
            0.5 * (grav->uavg_int[idx] + grav->uavg_int[idx + 1]) :
            grav->uavg_int[idx];
        double numerator;
        double gamma_term;
        double enthalpy_term;

        if (rho_edge <= 0.0) {
            grav->accel[idx] = 0.0;
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
        grav->accel[idx] = numerator / (gamma_term * gamma_term) * enthalpy_term;
    }
    grav->phi[grav->nbins] =
        -PRJ_GNEWT * baryon_mass_face[grav->nbins] / prj_gravity_rmax;
    for (idx = grav->nbins - 1; idx >= 0; --idx) {
        grav->phi[idx] = grav->phi[idx + 1] + 0.5 *
            (grav->rf[idx + 1] - grav->rf[idx]) *
            (idx == grav->nbins - 1 ? grav->accel[idx] : grav->accel[idx + 1] + grav->accel[idx]);
    }
    for (idx = 0; idx <= grav->nbins; ++idx) {
        grav->lapse[idx] = exp(grav->phi[idx] / (PRJ_CLIGHT * PRJ_CLIGHT));
    }
#else
    grav->phi[grav->nbins] =
        -PRJ_GNEWT * enclosed_face[grav->nbins] / prj_gravity_rmax;
    for (idx = grav->nbins - 1; idx >= 0; --idx) {
        double r0 = grav->rf[idx];
        double r1 = grav->rf[idx + 1];
        double f0 = PRJ_GNEWT * enclosed_face[idx] / (r0 * r0);
        double f1 = PRJ_GNEWT * enclosed_face[idx + 1] / (r1 * r1);

        grav->phi[idx] = grav->phi[idx + 1] - 0.5 * (r1 - r0) * (f0 + f1);
    }
    for (idx = 0; idx < grav->nbins; ++idx) {
        grav->accel[idx] =
            -(grav->phi[idx + 1] - grav->phi[idx]) /
            (grav->rf[idx + 1] - grav->rf[idx]);
    }
    for (idx = 0; idx <= grav->nbins; ++idx) {
        grav->lapse[idx] = 1.0;
    }
#endif

    free(gamma_face);
    free(baryon_mass_face);
    free(enclosed_face);
    PRJ_TIMER_CURRENT_STOP("gravity_integrate");
    prj_gravity_fill_mesh_fields(mesh);
}

int prj_gravity_apply(void)
{
    return 0;
}
