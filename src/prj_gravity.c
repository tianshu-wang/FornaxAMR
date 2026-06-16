#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

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
#define PRJ_GRAVITY_SQRT2 1.41421356237309504880168872420969808

static void prj_gravity_fill_mesh_fields(prj_mesh *mesh, const prj_grav *grav, const prj_mpi *mpi);

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

void prj_gravity_real_spherical_harmonics_all(double x1, double x2, double x3, double *out)
{
    double r;
    double mu;
    double somx2_sq;
    double somx2;
    double rho;
    double cphi;
    double sphi;
    double Plm[LMAX][LMAX];
    double cos_m[LMAX];
    double sin_m[LMAX];
    int l;
    int m;

    if (out == 0) {
        return;
    }

    for (l = 0; l < LMAX*LMAX; ++l) {
        out[l] = 0.0;
    }

    r = sqrt(x1 * x1 + x2 * x2 + x3 * x3);
    if (r <= 0.0) {
        out[PRJ_YLM_INDEX(0, 0)] = prj_gravity_ylm_norm(0, 0);
        return;
    }

    mu = x3 / r;
    if (mu > 1.0) {
        mu = 1.0;
    } else if (mu < -1.0) {
        mu = -1.0;
    }
    somx2_sq = 1.0 - mu * mu;
    if (somx2_sq < 0.0) {
        somx2_sq = 0.0;
    }
    somx2 = sqrt(somx2_sq);

    for (l = 0; l < LMAX; ++l) {
        for (m = 0; m < LMAX; ++m) {
            Plm[l][m] = 0.0;
        }
    }
    Plm[0][0] = 1.0;
    for (m = 1; m < LMAX; ++m) {
        Plm[m][m] = -(double)(2 * m - 1) * somx2 * Plm[m - 1][m - 1];
    }
    for (m = 0; m < LMAX - 1; ++m) {
        Plm[m + 1][m] = (double)(2 * m + 1) * mu * Plm[m][m];
        for (l = m + 2; l < LMAX; ++l) {
            Plm[l][m] = ((double)(2 * l - 1) * mu * Plm[l - 1][m] -
                (double)(l + m - 1) * Plm[l - 2][m]) / (double)(l - m);
        }
    }

    rho = sqrt(x1 * x1 + x2 * x2);
    if (rho > 0.0) {
        cphi = x1 / rho;
        sphi = x2 / rho;
    } else {
        cphi = 1.0;
        sphi = 0.0;
    }
    cos_m[0] = 1.0;
    sin_m[0] = 0.0;
#if LMAX > 1
    cos_m[1] = cphi;
    sin_m[1] = sphi;
    for (m = 1; m < LMAX - 1; ++m) {
        cos_m[m + 1] = 2.0 * cphi * cos_m[m] - cos_m[m - 1];
        sin_m[m + 1] = 2.0 * cphi * sin_m[m] - sin_m[m - 1];
    }
#endif

    for (l = 0; l < LMAX; ++l) {
        int yidx = PRJ_YLM_INDEX(l, 0);

        out[yidx] = prj_gravity_ylm_norm(l, 0) * Plm[l][0];
        for (m = 1; m <= l; ++m) {
            double ylm = PRJ_GRAVITY_SQRT2 *
                prj_gravity_ylm_norm(l, m) * Plm[l][m];

            out[PRJ_YLM_INDEX(l, m)] = ylm * cos_m[m];
            out[PRJ_YLM_INDEX(l, -m)] = ylm * sin_m[m];
        }
    }
}

double prj_gravity_real_spherical_harmonic(int l, int m, double x1, double x2, double x3)
{
    int abs_m = m < 0 ? -m : m;

    if (l < 0 || abs_m > l) {
        return 0.0;
    }
    if (l < LMAX) {
        double ylm[LMAX*LMAX];

        prj_gravity_real_spherical_harmonics_all(x1, x2, x3, ylm);
        return ylm[PRJ_YLM_INDEX(l, m)];
    } else {
        double r = sqrt(x1 * x1 + x2 * x2 + x3 * x3);
        double mu;
        double phi;
        double plm;
        double ylm;

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
            ylm *= PRJ_GRAVITY_SQRT2 * cos((double)abs_m * phi);
        } else if (m < 0) {
            ylm *= PRJ_GRAVITY_SQRT2 * sin((double)abs_m * phi);
        }
        return ylm;
    }
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

static int prj_gravity_use_multipole(const prj_grav *grav)
{
    return grav != 0 && grav->use_multipole_gravity != 0 && LMAX > 1;
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

static int prj_gravity_block_is_local_active(const prj_mpi *mpi, const prj_block *block)
{
    return block != 0 && block->id >= 0 && block->active == 1 &&
        (mpi == 0 || block->rank == mpi->rank);
}

static int prj_gravity_bin_index(const prj_grav *grav, double r)
{
    double log_span;
    int idx;

    if (grav == 0 || grav->nbins <= 0 || r < 0.0 || r >= grav->rmax) {
        return -1;
    }
    if (r < grav->dr_min) {
        return 0;
    }
    log_span = log(grav->rmax / grav->dr_min);
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
        grav->rmax <= 0.0) {
        return;
    }
    if (r >= grav->rmax) {
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

static void prj_gravity_cache_block_clear(prj_block *block, const prj_grav *grav)
{
    int use_multipole = prj_gravity_use_multipole(grav);
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
        if (use_multipole) {
            for (n = 0; n < LMAX*LMAX; ++n) {
                if (block->Ylm[n] != 0) {
                    block->Ylm[n][i] = 0.0;
                }
            }
        }
    }
}

static void prj_gravity_build_rf(prj_grav *grav, const prj_mesh *mesh)
{
    double min_cell;
    double log_span;
    int i;

    if (grav == 0 || grav->nbins <= 0 || grav->rf == 0 ||
        grav->rmax <= 0.0) {
        return;
    }
    min_cell = mesh != 0 ? mesh->min_allowable_cell_size : 0.0;
    grav->min_cell = min_cell;
    grav->dr_min = 0.5 * min_cell;
    if (grav->dr_min <= 0.0 || grav->rmax <= grav->dr_min) {
        grav->dr_min = grav->rmax / (double)grav->nbins;
    }
    log_span = log(grav->rmax / grav->dr_min);
    grav->rf[0] = 0.0;
    if (grav->nbins == 1) {
        grav->rf[1] = grav->rmax;
        return;
    }
    for (i = 1; i <= grav->nbins; ++i) {
        grav->rf[i] = grav->dr_min *
            exp(log_span * (double)(i - 1) / (double)(grav->nbins - 1));
    }
}

void prj_gravity_cache_block(prj_block *block, const prj_mesh *mesh, const prj_grav *grav)
{
    int use_multipole = prj_gravity_use_multipole(grav);
    double x_com[3] = {0.0, 0.0, 0.0};
    int i;
    int j;
    int k;

    if (block == 0 || block->ridx == 0 || block->fr == 0 || block->r_com == 0) {
        return;
    }
    if (block->id < 0 || block->dx[0] <= 0.0 || block->dx[1] <= 0.0 || block->dx[2] <= 0.0) {
        prj_gravity_cache_block_clear(block, grav);
        return;
    }
    if (mesh != 0) {
        x_com[0] = mesh->x_com[0];
        x_com[1] = mesh->x_com[1];
        x_com[2] = mesh->x_com[2];
    }

    for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
        for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
            for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                double x1 = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
                double x2 = block->xmin[1] + ((double)j + 0.5) * block->dx[1];
                double x3 = block->xmin[2] + ((double)k + 0.5) * block->dx[2];
                double dx1 = x1 - x_com[0];
                double dx2 = x2 - x_com[1];
                double dx3 = x3 - x_com[2];
                int cache_idx = prj_block_cache_index(i, j, k);

                prj_gravity_cache_entry(grav, block->r_com[cache_idx],
                    &block->ridx[cache_idx], &block->fr[cache_idx]);
                if (use_multipole) {
                    double ylm[LMAX*LMAX];
                    int l;

                    prj_gravity_real_spherical_harmonics_all(dx1, dx2, dx3, ylm);
                    for (l = 0; l < LMAX; ++l) {
                        int m;

                        for (m = -l; m <= l; ++m) {
                            int yidx = PRJ_YLM_INDEX(l, m);

                            if (block->Ylm[yidx] != 0) {
                                block->Ylm[yidx][cache_idx] = ylm[yidx];
                            }
                        }
                    }
                }
            }
        }
    }
}

void prj_gravity_cache_mesh(prj_mesh *mesh, const prj_grav *grav)
{
    int bidx;

    if (mesh == 0) {
        return;
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_gravity_cache_block(&mesh->blocks[bidx], mesh, grav);
    }
}

void prj_gravity_rebuild_grid(prj_sim *sim, const prj_mpi *mpi)
{
    int i;

    if (sim == 0) {
        return;
    }
    prj_gravity_build_rf(&sim->grav, &sim->mesh);
    if (sim->grav.lapse != 0) {
        for (i = 0; i <= sim->grav.nbins; ++i) {
            sim->grav.lapse[i] = 1.0;
        }
    }
    if (prj_gravity_use_multipole(&sim->grav)) {
        prj_gravity_zero_multipole_coefficients(&sim->grav);
    }
    prj_gravity_cache_mesh(&sim->mesh, &sim->grav);
    prj_gravity_fill_mesh_fields(&sim->mesh, &sim->grav, mpi);
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
    free(grav->reduce_avg_buf);
    free(grav->reduce_lm_buf);
    grav->rf = 0;
    grav->ms = 0;
    grav->phi = 0;
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
    grav->reduce_avg_buf = 0;
    grav->reduce_lm_buf = 0;
    grav->nbins = 0;
    grav->dr_min = 0.0;
    grav->rmax = 0.0;
    grav->min_cell = 0.0;
}

void prj_gravity_init(prj_sim *sim, const prj_mpi *mpi)
{
    prj_grav *grav;
    double span1;
    double span2;
    double span3;
    int i;
    int yidx;
    int use_multipole;

    if (sim == 0) {
        return;
    }

    grav = &sim->grav;
    prj_gravity_free(grav);
    use_multipole = prj_gravity_use_multipole(grav);

    span1 = fabs(sim->mesh.coord.x1max - sim->mesh.coord.x1min);
    span2 = fabs(sim->mesh.coord.x2max - sim->mesh.coord.x2min);
    span3 = fabs(sim->mesh.coord.x3max - sim->mesh.coord.x3min);
    grav->rmax = 0.5 * prj_gravity_min_double(span1, prj_gravity_min_double(span2, span3));

    grav->nbins = PRJ_GRAVITY_DEFAULT_NBINS;

    grav->rf = (double *)calloc((size_t)grav->nbins + 1U, sizeof(*grav->rf));
    grav->ms = (double *)calloc((size_t)grav->nbins, sizeof(*grav->ms));
    grav->phi = (double *)calloc((size_t)grav->nbins + 1U, sizeof(*grav->phi));
    grav->lapse = (double *)calloc((size_t)grav->nbins + 1U, sizeof(*grav->lapse));
    if (use_multipole) {
        for (yidx = 0; yidx < LMAX*LMAX; ++yidx) {
            grav->Clm[yidx] = (double *)calloc((size_t)grav->nbins, sizeof(*grav->Clm[yidx]));
            grav->Dlm[yidx] = (double *)calloc((size_t)grav->nbins, sizeof(*grav->Dlm[yidx]));
        }
    }
    grav->vol = (double *)calloc((size_t)grav->nbins, sizeof(*grav->vol));
    grav->rho_avg = (double *)calloc((size_t)grav->nbins, sizeof(*grav->rho_avg));
    grav->vr_avg = (double *)calloc((size_t)grav->nbins, sizeof(*grav->vr_avg));
    grav->pgas_avg = (double *)calloc((size_t)grav->nbins, sizeof(*grav->pgas_avg));
    grav->uavg_int = (double *)calloc((size_t)grav->nbins, sizeof(*grav->uavg_int));
    grav->erad_avg = (double *)calloc((size_t)grav->nbins, sizeof(*grav->erad_avg));
    grav->prad_avg = (double *)calloc((size_t)grav->nbins, sizeof(*grav->prad_avg));
    grav->vdotF_avg = (double *)calloc((size_t)grav->nbins, sizeof(*grav->vdotF_avg));
    grav->reduce_avg_buf = (double *)calloc(9U * (size_t)grav->nbins,
        sizeof(*grav->reduce_avg_buf));
    if (use_multipole) {
        grav->reduce_lm_buf = (double *)calloc(2U * (size_t)(LMAX * LMAX) *
            (size_t)grav->nbins, sizeof(*grav->reduce_lm_buf));
    }
    if (grav->rf == 0 || grav->ms == 0 || grav->phi == 0 ||
        grav->lapse == 0 || (use_multipole && !prj_gravity_multipole_arrays_valid(grav)) ||
        grav->vol == 0 || grav->rho_avg == 0 ||
        grav->vr_avg == 0 || grav->pgas_avg == 0 || grav->uavg_int == 0 ||
        grav->erad_avg == 0 || grav->prad_avg == 0 || grav->vdotF_avg == 0 ||
        grav->reduce_avg_buf == 0 || (use_multipole && grav->reduce_lm_buf == 0)) {
        prj_gravity_free(grav);
        return;
    }

    prj_gravity_build_rf(grav, &sim->mesh);
    for (i = 0; i <= grav->nbins; ++i) {
        grav->lapse[i] = 1.0;
    }

    prj_gravity_cache_mesh(&sim->mesh, grav);
    prj_gravity_fill_mesh_fields(&sim->mesh, grav, mpi);
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

static double prj_gravity_block_cached_phi_at(const prj_grav *grav, const prj_block *block,
    int i, int j, int k)
{
    int cache_idx;
    int idx;
    double fr;

    if (block == 0 || grav == 0 || grav->nbins <= 0 ||
        grav->phi == 0 || block->ridx == 0 || block->fr == 0) {
        return 0.0;
    }
    cache_idx = prj_block_cache_index(i, j, k);
    idx = block->ridx[cache_idx];
    fr = block->fr[cache_idx];
    if (idx == PRJ_GRAVITY_CACHE_INVALID) {
        /* Cell could not be binned (degenerate radial grid): flat asymptotic
           potential, i.e. no gravitational force. */
        return grav->phi[grav->nbins];
    }
    if (fr >= PRJ_GRAVITY_CACHE_SKIP_REDUCE) {
        /* r >= rmax: extrapolate the monopole potential as the exact vacuum
           tail Phi(r) = Phi(rmax) * rmax / r. */
        double phi_rmax = grav->phi[grav->nbins];

        if (block->r_com != 0 && grav->rmax > 0.0) {
            double r = block->r_com[cache_idx];

            if (r > grav->rmax) {
                return phi_rmax * grav->rmax / r;
            }
        }
        return phi_rmax;
    }
    if (grav->nbins == 1) {
        return grav->phi[0];
    }
    if (fr < 0.0) {
        return grav->phi[0];
    }
    if (fr <= 1.0) {
        int left = idx > 0 ? idx - 1 : 0;

        return (1.0 - fr) * grav->phi[left] + fr * grav->phi[left + 1];
    }
    if (fr >= PRJ_GRAVITY_CACHE_LAST_VALUE || idx >= grav->nbins - 1) {
        return grav->phi[grav->nbins - 1];
    }
    fr -= 1.0;
    return (1.0 - fr) * grav->phi[idx] + fr * grav->phi[idx + 1];
}

static double prj_gravity_block_cached_lapse_at(const prj_grav *grav, const prj_block *block,
    int i, int j, int k)
{
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

static void prj_gravity_fill_block_multipole_phi(const prj_grav *grav, const prj_block *block,
    double *Phi)
{
    double r_l[PRJ_BLOCK_NCELLS];
    double inv_r_l1[PRJ_BLOCK_NCELLS];
    double inv_r[PRJ_BLOCK_NCELLS];
    int inside[PRJ_BLOCK_NCELLS];
    int outside[PRJ_BLOCK_NCELLS];
    int n_inside = 0;
    int n_outside = 0;
    int cache_idx;
    int l;

    if (Phi == 0) {
        return;
    }
    for (cache_idx = 0; cache_idx < PRJ_BLOCK_NCELLS; ++cache_idx) {
        Phi[cache_idx] = 0.0;
        r_l[cache_idx] = 1.0;
        inv_r_l1[cache_idx] = 0.0;
        inv_r[cache_idx] = 0.0;
    }
    if (block == 0 || grav == 0 || grav->nbins <= 0 || block->ridx == 0 ||
        block->fr == 0 || block->r_com == 0) {
        return;
    }

    for (cache_idx = 0; cache_idx < PRJ_BLOCK_NCELLS; ++cache_idx) {
        int idx = block->ridx[cache_idx];
        double r = block->r_com[cache_idx];
        double fr = block->fr[cache_idx];

        if (idx == PRJ_GRAVITY_CACHE_INVALID || r <= 0.0) {
            continue;
        }
        inv_r[cache_idx] = 1.0 / r;
        inv_r_l1[cache_idx] = inv_r[cache_idx];
        if (fr < PRJ_GRAVITY_CACHE_SKIP_REDUCE) {
            inside[n_inside] = cache_idx;
            n_inside += 1;
        } else {
            outside[n_outside] = cache_idx;
            n_outside += 1;
        }
    }

    for (l = 1; l < LMAX; ++l) {
        int n;
        int yidx;
        int yidx_end;

        for (n = 0; n < n_inside; ++n) {
            cache_idx = inside[n];
            r_l[cache_idx] *= block->r_com[cache_idx];
            inv_r_l1[cache_idx] *= inv_r[cache_idx];
        }
        for (n = 0; n < n_outside; ++n) {
            cache_idx = outside[n];
            r_l[cache_idx] *= block->r_com[cache_idx];
            inv_r_l1[cache_idx] *= inv_r[cache_idx];
        }

        yidx = PRJ_YLM_INDEX(l, -l);
        yidx_end = PRJ_YLM_INDEX(l, l) + 1;
        for (; yidx < yidx_end; ++yidx) {
            const double *Ylm = block->Ylm[yidx];
            const double *Clm = grav->Clm[yidx];
            const double *Dlm = grav->Dlm[yidx];

            if (Ylm == 0 || Clm == 0 || Dlm == 0) {
                continue;
            }
            for (n = 0; n < n_inside; ++n) {
                int idx;
                int cidx;
                double coeff;

                cache_idx = inside[n];
                idx = block->ridx[cache_idx];
                cidx = idx - 1;
                coeff = Dlm[idx] * r_l[cache_idx];
                if (cidx >= 0) {
                    coeff += Clm[cidx] * inv_r_l1[cache_idx];
                }
                Phi[cache_idx] += coeff * Ylm[cache_idx];
            }
            for (n = 0; n < n_outside; ++n) {
                int idx;

                cache_idx = outside[n];
                idx = block->ridx[cache_idx];
                Phi[cache_idx] +=
                    Clm[idx] * inv_r_l1[cache_idx] * Ylm[cache_idx];
            }
        }
    }
    for (cache_idx = 0; cache_idx < PRJ_BLOCK_NCELLS; ++cache_idx) {
        Phi[cache_idx] *= -PRJ_GNEWT;
    }
}

/* Isotropic 3D Sobel gradient of Phi at cell (i,j,k). Along each axis the
   central derivative kernel [-1, 0, +1] is smoothed in the two perpendicular
   directions with the isotropic weights [1, sqrt(2), 1] (sum 2+sqrt(2)),
   which minimises the angular error of the discrete gradient relative to the
   plain [1, 2, 1] Sobel. The caller (fill_block_fields, over the field range)
   guarantees the full 3x3x3 stencil lies inside the filled Phi region, so no
   edge handling is required. */
static void prj_gravity_phi_sobel_gradient(const double *Phi, int i, int j, int k,
    const double dx[3], double grad[3])
{
    static const double w[3] = {1.0, PRJ_GRAVITY_SQRT2, 1.0};
    const double norm = (2.0 + PRJ_GRAVITY_SQRT2) * (2.0 + PRJ_GRAVITY_SQRT2);
    double gx = 0.0;
    double gy = 0.0;
    double gz = 0.0;
    int a;
    int b;

    grad[0] = 0.0;
    grad[1] = 0.0;
    grad[2] = 0.0;
    if (Phi == 0) {
        return;
    }
    for (a = -1; a <= 1; ++a) {
        for (b = -1; b <= 1; ++b) {
            double wab = w[a + 1] * w[b + 1];

            /* d/dx: derivative along x, smoothing over (y, z) = (a, b). */
            gx += wab * (Phi[IDX(i + 1, j + a, k + b)] - Phi[IDX(i - 1, j + a, k + b)]);
            /* d/dy: derivative along y, smoothing over (x, z) = (a, b). */
            gy += wab * (Phi[IDX(i + a, j + 1, k + b)] - Phi[IDX(i + a, j - 1, k + b)]);
            /* d/dz: derivative along z, smoothing over (x, y) = (a, b). */
            gz += wab * (Phi[IDX(i + a, j + b, k + 1)] - Phi[IDX(i + a, j + b, k - 1)]);
        }
    }
    grad[0] = dx[0] > 0.0 ? gx / (norm * 2.0 * dx[0]) : 0.0;
    grad[1] = dx[1] > 0.0 ? gy / (norm * 2.0 * dx[1]) : 0.0;
    grad[2] = dx[2] > 0.0 ? gz / (norm * 2.0 * dx[2]) : 0.0;
}

double prj_gravity_block_lapse_at(const prj_grav *grav, const prj_block *block, int i, int j, int k)
{
    if (block != 0 && block->lapse != 0) {
        return block->lapse[prj_block_cache_index(i, j, k)];
    }
    return prj_gravity_block_cached_lapse_at(grav, block, i, j, k);
}

static void prj_gravity_fill_block_field_defaults(prj_block *block)
{
    int d;

    if (block == 0) {
        return;
    }
    prj_fill(block->lapse, (size_t)PRJ_BLOCK_NCELLS, 1.0);
    for (d = 0; d < 3; ++d) {
        prj_fill(block->grav[d], (size_t)PRJ_BLOCK_NCELLS, 0.0);
    }
}

static void prj_gravity_fill_block_fields(const prj_mesh *mesh, prj_block *block, const prj_grav *grav)
{
    const int field_lo = PRJ_NGHOST > 0 ? -1 : 0;
    const int field_hi = PRJ_BLOCK_SIZE + (PRJ_NGHOST > 0 ? 1 : 0);
    double Phi[PRJ_BLOCK_NCELLS];
    int use_multipole;
    int cache_idx;
    int i;
    int j;
    int k;

    (void)mesh;
    if (block == 0 || block->lapse == 0 || block->r_com == 0 ||
        block->grav[0] == 0 || block->grav[1] == 0 || block->grav[2] == 0) {
        return;
    }
    if (block->id < 0 || block->dx[0] <= 0.0 || block->dx[1] <= 0.0 ||
        block->dx[2] <= 0.0 || grav == 0) {
        prj_gravity_fill_block_field_defaults(block);
        return;
    }

    /* Build the gravitational potential on the grid, then take grav = -grad(Phi)
       with the same finite-difference stencil used for the multipole terms.
       The multipole pass (l >= 1) owns and scales Phi; the monopole (l = 0)
       potential is interpolated radially and added on top. */
    use_multipole = prj_gravity_use_multipole(grav) &&
        prj_gravity_multipole_arrays_valid(grav);
    if (use_multipole) {
        prj_gravity_fill_block_multipole_phi(grav, block, Phi);
    } else {
        for (cache_idx = 0; cache_idx < PRJ_BLOCK_NCELLS; ++cache_idx) {
            Phi[cache_idx] = 0.0;
        }
    }
    for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
        for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
            for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                cache_idx = prj_block_cache_index(i, j, k);
                Phi[cache_idx] += prj_gravity_block_cached_phi_at(grav, block, i, j, k);
            }
        }
    }

    for (i = field_lo; i < field_hi; ++i) {
        for (j = field_lo; j < field_hi; ++j) {
            for (k = field_lo; k < field_hi; ++k) {
                double grad[3];

                cache_idx = prj_block_cache_index(i, j, k);

                block->lapse[cache_idx] = prj_gravity_block_cached_lapse_at(grav, block, i, j, k);
                prj_gravity_phi_sobel_gradient(Phi, i, j, k, block->dx, grad);
                block->grav[0][cache_idx] = -grad[0];
                block->grav[1][cache_idx] = -grad[1];
                block->grav[2][cache_idx] = -grad[2];
            }
        }
    }
}

static void prj_gravity_fill_mesh_fields(prj_mesh *mesh, const prj_grav *grav, const prj_mpi *mpi)
{
    int bidx;

    if (mesh == 0) {
        return;
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        if (!prj_gravity_block_is_local_active(mpi, &mesh->blocks[bidx])) {
            continue;
        }
        prj_gravity_fill_block_fields(mesh, &mesh->blocks[bidx], grav);
    }
}

static double prj_gravity_shell_volume(const prj_grav *grav, int idx)
{
    double r0 = grav->rf[idx];
    double r1 = grav->rf[idx + 1];

    return (4.0 / 3.0) * M_PI * (r1 * r1 * r1 - r0 * r0 * r0);
}

/* Spread the donor bin j's accumulated content uniformly (per unit shell
   volume) over bins i..j, conserving the total. After this every bin in
   [i, j] holds content proportional to its own shell volume. grav->vol is
   spread the same way, so each filled bin's accumulated volume is likewise
   proportional to its shell volume; dividing the spread content by the spread
   grav->vol (in prj_gravity_monopole_reduce) then yields an identical
   volume-average across the range, and grav->vol * density used in the TOV
   integral recovers the donor mass exactly. */
static void prj_gravity_spread_bins(double *arr, const prj_grav *grav,
                                    int i, int j, double inv_total_vol)
{
    double donor = arr[j];
    int k;

    for (k = i; k <= j; ++k) {
        arr[k] = donor * prj_gravity_shell_volume(grav, k) * inv_total_vol;
    }
}

/* The radial face grid (grav->rf) need not align with cell centers, so a
   shell can contain no cells and reduce to an empty bin (zero mass/volume).
   For each maximal run of empty bins [i, j-1] bounded below by a non-empty
   bin (or i == 0) and above by the non-empty bin j, redistribute bin j's
   accumulated reductions over bins i..j weighted by shell volume, so the
   per-volume averages are uniform across the run and the totals (mass and
   multipole moments) are conserved. */
static void prj_gravity_fill_empty_bins(prj_grav *grav, int use_multipole)
{
    int i;

    if (grav == 0 || grav->nbins <= 0 || grav->rf == 0 || grav->vol == 0) {
        return;
    }

    i = 0;
    while (i < grav->nbins) {
        int j;
        int k;
        double total_vol;
        double inv_total_vol;

        if (grav->vol[i] > 0.0) {
            ++i;
            continue;
        }
        /* Empty run starts at i (by construction i == 0 or bin i-1 is
           non-empty). Find the first non-empty bin j above it. */
        j = i;
        while (j < grav->nbins && grav->vol[j] <= 0.0) {
            ++j;
        }
        if (j >= grav->nbins) {
            /* Trailing empty bins with no donor above: nothing to spread. */
            break;
        }

        total_vol = 0.0;
        for (k = i; k <= j; ++k) {
            total_vol += prj_gravity_shell_volume(grav, k);
        }
        if (total_vol <= 0.0) {
            i = j + 1;
            continue;
        }
        inv_total_vol = 1.0 / total_vol;

        prj_gravity_spread_bins(grav->ms,        grav, i, j, inv_total_vol);
        prj_gravity_spread_bins(grav->vol,       grav, i, j, inv_total_vol);
        prj_gravity_spread_bins(grav->rho_avg,   grav, i, j, inv_total_vol);
        prj_gravity_spread_bins(grav->vr_avg,    grav, i, j, inv_total_vol);
        prj_gravity_spread_bins(grav->pgas_avg,  grav, i, j, inv_total_vol);
        prj_gravity_spread_bins(grav->uavg_int,  grav, i, j, inv_total_vol);
        prj_gravity_spread_bins(grav->erad_avg,  grav, i, j, inv_total_vol);
        prj_gravity_spread_bins(grav->prad_avg,  grav, i, j, inv_total_vol);
        prj_gravity_spread_bins(grav->vdotF_avg, grav, i, j, inv_total_vol);

        if (use_multipole) {
            int yidx;

            for (yidx = 0; yidx < LMAX * LMAX; ++yidx) {
                if (grav->Clm[yidx] != 0) {
                    prj_gravity_spread_bins(grav->Clm[yidx], grav, i, j,
                        inv_total_vol);
                }
                if (grav->Dlm[yidx] != 0) {
                    prj_gravity_spread_bins(grav->Dlm[yidx], grav, i, j,
                        inv_total_vol);
                }
            }
        }

        i = j + 1;
    }
}

void prj_gravity_monopole_reduce(prj_mesh *mesh, prj_grav *grav, const prj_mpi *mpi, int stage)
{
    int bidx;
    int idx;
    int use_multipole;
    double x_com[3];

    if (mesh == 0 || grav == 0 || grav->ms == 0 || grav->vol == 0 ||
        grav->rho_avg == 0 || grav->vr_avg == 0 || grav->pgas_avg == 0 ||
        grav->uavg_int == 0 || grav->erad_avg == 0 || grav->prad_avg == 0 ||
        grav->vdotF_avg == 0) {
        return;
    }
    x_com[0] = mesh->x_com[0];
    x_com[1] = mesh->x_com[1];
    x_com[2] = mesh->x_com[2];
    use_multipole = prj_gravity_use_multipole(grav);
    if (use_multipole && !prj_gravity_multipole_arrays_valid(grav)) {
        return;
    }

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
    if (use_multipole) {
        prj_gravity_zero_multipole_coefficients(grav);
    }

    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        const prj_block *block = &mesh->blocks[bidx];
        double *W;
        int i;
        int j;
        int k;

        if (!prj_gravity_block_is_local_active(mpi, block)) {
            continue;
        }
        W = stage == 2 ? block->W1 : block->W;
        for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
            for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    int cache_idx = prj_block_cache_index(i, j, k);
                    double x1 = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
                    double x2 = block->xmin[1] + ((double)j + 0.5) * block->dx[1];
                    double x3 = block->xmin[2] + ((double)k + 0.5) * block->dx[2];
                    double dx1 = x1 - x_com[0];
                    double dx2 = x2 - x_com[1];
                    double dx3 = x3 - x_com[2];
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
                                    /* Radiation E/F are stored in RAD_SCALE*erg
                                       units; convert back to physical erg for the
                                       gravitational source. */
                                    double e_rad = W[VIDX(PRJ_PRIM_RAD_E(field, group), i, j, k)] * RAD_SCALE;
                                    double f1 = W[VIDX(PRJ_PRIM_RAD_F1(field, group), i, j, k)] * RAD_SCALE;
                                    double f2 = W[VIDX(PRJ_PRIM_RAD_F2(field, group), i, j, k)] * RAD_SCALE;
                                    double f3 = W[VIDX(PRJ_PRIM_RAD_F3(field, group), i, j, k)] * RAD_SCALE;
                                    double fr = r > 0.0 ? (f1 * dx1 + f2 * dx2 + f3 * dx3) / r : 0.0;

                                    erad += e_rad / (PRJ_CLIGHT * PRJ_CLIGHT);
                                    prad += (e_rad / 3.0) / (PRJ_CLIGHT * PRJ_CLIGHT * PRJ_CLIGHT);
                                    /* v.F/c^4 has mass-density units for the TOV integrand. */
                                    vdotF += (v1 * f1 + v2 * f2 + v3 * f3) /
                                        (PRJ_CLIGHT * PRJ_CLIGHT * PRJ_CLIGHT * PRJ_CLIGHT);
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
                    }
                }
            }
        }
        if (use_multipole) {
            int l;

            for (l = 0; l < LMAX; ++l) {
                double ylm_scale = (4.0 * M_PI) / (double)(2 * l + 1);
                int m;

                for (m = -l; m <= l; ++m) {
                    int yidx = PRJ_YLM_INDEX(l, m);
                    const double *Ylm = block->Ylm[yidx];
                    double *Clm = grav->Clm[yidx];
                    double *Dlm = grav->Dlm[yidx];

                    if (Ylm == 0) {
                        continue;
                    }
                    for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
                        for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                            for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                                int cache_idx = prj_block_cache_index(i, j, k);
                                double fr = block->fr != 0 ? block->fr[cache_idx] : 0.0;
                                double r = block->r_com != 0 ? block->r_com[cache_idx] : 0.0;
                                double rho;
                                double r_l;
                                double inv_r_l1;
                                double coeff;

                                idx = block->ridx != 0 ? block->ridx[cache_idx] :
                                    PRJ_GRAVITY_CACHE_INVALID;
                                if (idx < 0 || fr >= PRJ_GRAVITY_CACHE_SKIP_REDUCE) {
                                    continue;
                                }
                                if (block->r_com == 0) {
                                    double x1 = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
                                    double x2 = block->xmin[1] + ((double)j + 0.5) * block->dx[1];
                                    double x3 = block->xmin[2] + ((double)k + 0.5) * block->dx[2];
                                    double dx1 = x1 - x_com[0];
                                    double dx2 = x2 - x_com[1];
                                    double dx3 = x3 - x_com[2];

                                    r = sqrt(dx1 * dx1 + dx2 * dx2 + dx3 * dx3);
                                }
                                if (r <= 0.0) {
                                    continue;
                                }
                                rho = W[(size_t)PRJ_PRIM_RHO * (size_t)PRJ_BLOCK_NCELLS +
                                    (size_t)cache_idx];
                                r_l = prj_gravity_pow_int(r, l);
                                inv_r_l1 = 1.0 / (r_l * r);
                                coeff = ylm_scale * rho * block->vol * Ylm[cache_idx];
                                Clm[idx] += coeff * r_l;
                                Dlm[idx] += coeff * inv_r_l1;
                            }
                        }
                    }
                }
            }
        }
    }

#if defined(PRJ_ENABLE_MPI)
    {
        if (mpi != 0 && mpi->totrank > 1 && grav->reduce_avg_buf != 0 &&
            (!use_multipole || grav->reduce_lm_buf != 0)) {
            double *restrict buf = grav->reduce_avg_buf;
            size_t nb = (size_t)grav->nbins;

            memcpy(buf + 0U * nb, grav->ms,        nb * sizeof(double));
            memcpy(buf + 1U * nb, grav->vol,       nb * sizeof(double));
            memcpy(buf + 2U * nb, grav->rho_avg,   nb * sizeof(double));
            memcpy(buf + 3U * nb, grav->vr_avg,    nb * sizeof(double));
            memcpy(buf + 4U * nb, grav->pgas_avg,  nb * sizeof(double));
            memcpy(buf + 5U * nb, grav->uavg_int,  nb * sizeof(double));
            memcpy(buf + 6U * nb, grav->erad_avg,  nb * sizeof(double));
            memcpy(buf + 7U * nb, grav->prad_avg,  nb * sizeof(double));
            memcpy(buf + 8U * nb, grav->vdotF_avg, nb * sizeof(double));
            MPI_Allreduce(MPI_IN_PLACE, buf, 9 * grav->nbins, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
            memcpy(grav->ms,        buf + 0U * nb, nb * sizeof(double));
            memcpy(grav->vol,       buf + 1U * nb, nb * sizeof(double));
            memcpy(grav->rho_avg,   buf + 2U * nb, nb * sizeof(double));
            memcpy(grav->vr_avg,    buf + 3U * nb, nb * sizeof(double));
            memcpy(grav->pgas_avg,  buf + 4U * nb, nb * sizeof(double));
            memcpy(grav->uavg_int,  buf + 5U * nb, nb * sizeof(double));
            memcpy(grav->erad_avg,  buf + 6U * nb, nb * sizeof(double));
            memcpy(grav->prad_avg,  buf + 7U * nb, nb * sizeof(double));
            memcpy(grav->vdotF_avg, buf + 8U * nb, nb * sizeof(double));
            if (use_multipole) {
                double *restrict lmbuf = grav->reduce_lm_buf;
                int yidx;

                for (yidx = 0; yidx < LMAX * LMAX; ++yidx) {
                    memcpy(lmbuf + (size_t)(2 * yidx + 0) * nb, grav->Clm[yidx],
                        nb * sizeof(double));
                    memcpy(lmbuf + (size_t)(2 * yidx + 1) * nb, grav->Dlm[yidx],
                        nb * sizeof(double));
                }
                MPI_Allreduce(MPI_IN_PLACE, lmbuf, 2 * LMAX * LMAX * grav->nbins,
                    MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                for (yidx = 0; yidx < LMAX * LMAX; ++yidx) {
                    memcpy(grav->Clm[yidx], lmbuf + (size_t)(2 * yidx + 0) * nb,
                        nb * sizeof(double));
                    memcpy(grav->Dlm[yidx], lmbuf + (size_t)(2 * yidx + 1) * nb,
                        nb * sizeof(double));
                }
            }
        }
    }
#endif

    /* Fill radial bins that no cell mapped into by spreading the neighboring
       outer bin's reductions across the gap (volume-weighted), before the
       cumulative-sum and shell-average steps consume them. */
    prj_gravity_fill_empty_bins(grav, use_multipole);

    if (use_multipole) {
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
    }

    for (idx = 0; idx < grav->nbins; ++idx) {
        /* Divide by the ACTUAL accumulated cell volume in the bin, not the
           analytical shell volume 4/3 pi (r1^3 - r0^3).  The gravity log-grid
           is finer than the hydro cells, so the cells whose centers land in a
           bin occupy a volume that differs from the analytical shell by a
           radius-locked factor (~+/-30%); dividing by that analytical volume
           biases these volume-averages and feeds a grid-locked spurious force.
           grav->vol[idx] is the true accumulated cell volume, so this yields a
           genuine volume-weighted average. */
        double bin_vol = grav->vol[idx];

        if (bin_vol > 0.0) {
            grav->rho_avg[idx] /= bin_vol;
            grav->vr_avg[idx] /= bin_vol;
            grav->pgas_avg[idx] /= bin_vol;
            grav->uavg_int[idx] /= bin_vol;
            grav->erad_avg[idx] /= bin_vol;
            grav->prad_avg[idx] /= bin_vol;
            grav->vdotF_avg[idx] /= bin_vol;
        }
    }
}

void prj_gravity_monopole_integrate(prj_mesh *mesh, prj_grav *grav, const prj_mpi *mpi)
{
    double *restrict enclosed_face;
    double *restrict baryon_mass_face;
    double *restrict gamma_face;
#if PRJ_GRAVITY_USE_GR
    double *restrict accel_face;
#endif
    int idx;

    if (grav == 0 || grav->nbins <= 0) {
        return;
    }

    enclosed_face = (double *)calloc((size_t)grav->nbins + 1U, sizeof(*enclosed_face));
    baryon_mass_face = (double *)calloc((size_t)grav->nbins + 1U, sizeof(*baryon_mass_face));
    gamma_face = (double *)calloc((size_t)grav->nbins + 1U, sizeof(*gamma_face));
#if PRJ_GRAVITY_USE_GR
    accel_face = (double *)calloc((size_t)grav->nbins + 1U, sizeof(*accel_face));
#endif
    if (enclosed_face == 0 || baryon_mass_face == 0 || gamma_face == 0
#if PRJ_GRAVITY_USE_GR
        || accel_face == 0
#endif
        ) {
#if PRJ_GRAVITY_USE_GR
        free(accel_face);
#endif
        free(gamma_face);
        free(baryon_mass_face);
        free(enclosed_face);
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
                double r1 = grav->rf[idx + 1];
                /* Convert per-bin densities to extensive contributions with the
                   ACTUAL accumulated cell volume (matches the volume-averaged
                   rho_avg/uavg/erad above), so bin_vol*density recovers the
                   deposited mass-energy and total mass is conserved exactly.
                   Geometric factors below (r1, 4 pi r^3 P) keep the physical
                   bin-face radius. */
                double bin_vol = grav->vol[idx];
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
                    2.0 * PRJ_GNEWT * (m_tov + bin_vol * integrand) / (r1 * c2);
                if (gsq < 1.0e-12) {
                    gsq = 1.0e-12;
                }
                gamma_face[idx + 1] = sqrt(gsq);

                /* Re-evaluate integrand with updated gamma_face[idx+1]. */
                gamma_avg = 0.5 * (gamma_face[idx] + gamma_face[idx + 1]);
                integrand = (rho + u / c2 + erad) * gamma_avg + vdF;
                m_tov += bin_vol * integrand;
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
            accel_face[idx + 1] = 0.0;
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
        accel_face[idx + 1] = numerator / (gamma_term * gamma_term) * enthalpy_term;
    }
    grav->phi[grav->nbins] =
        -PRJ_GNEWT * baryon_mass_face[grav->nbins] / grav->rmax;
    for (idx = grav->nbins - 1; idx >= 0; --idx) {
        grav->phi[idx] = grav->phi[idx + 1] + 0.5 *
            (grav->rf[idx + 1] - grav->rf[idx]) *
            (accel_face[idx] + accel_face[idx + 1]);
    }
    for (idx = 0; idx <= grav->nbins; ++idx) {
        grav->lapse[idx] = exp(grav->phi[idx] / (PRJ_CLIGHT * PRJ_CLIGHT));
    }
#else
    grav->phi[grav->nbins] =
        -PRJ_GNEWT * enclosed_face[grav->nbins] / grav->rmax;
    for (idx = grav->nbins - 1; idx >= 0; --idx) {
        double r0 = grav->rf[idx];
        double r1 = grav->rf[idx + 1];
        double f0 = PRJ_GNEWT * enclosed_face[idx] / (r0 * r0);
        double f1 = PRJ_GNEWT * enclosed_face[idx + 1] / (r1 * r1);

        grav->phi[idx] = grav->phi[idx + 1] - 0.5 * (r1 - r0) * (f0 + f1);
    }
    for (idx = 0; idx <= grav->nbins; ++idx) {
        grav->lapse[idx] = 1.0;
    }
#endif

#if PRJ_GRAVITY_USE_GR
    free(accel_face);
#endif
    free(gamma_face);
    free(baryon_mass_face);
    free(enclosed_face);
    prj_gravity_fill_mesh_fields(mesh, grav, mpi);
}

int prj_gravity_apply(void)
{
    return 0;
}
