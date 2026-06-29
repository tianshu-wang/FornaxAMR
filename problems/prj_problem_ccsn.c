#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "prj.h"

#define PRJ_CCSN_KELVIN_PER_MEV 1.160451812e10

typedef struct prj_ccsn_profile {
    int npts;
    double *radius;
    double *rho;
    double *temp;
    double *ye;
    double *vr;
} prj_ccsn_profile;

typedef struct prj_ccsn_init_amr_ctx {
    int npts;
    const double *radius;
    double *LP;
    double scale_factor;
} prj_ccsn_init_amr_ctx;

static int prj_ccsn_init_amr_ctx_build(prj_ccsn_init_amr_ctx *ctx, const prj_ccsn_profile *profile,
    double scale_factor, prj_eos *eos, double r_min_domain, double r_max_domain)
{
    double *pressure;
    int i;
    int i_first;
    int i_last;
    int i_lo;
    int i_hi;

    if (ctx == 0 || profile == 0 || profile->npts < 2 || eos == 0) {
        return 1;
    }
    ctx->npts = profile->npts;
    ctx->radius = profile->radius;
    ctx->scale_factor = scale_factor;
    ctx->LP = (double *)malloc((size_t)profile->npts * sizeof(*ctx->LP));
    if (ctx->LP == 0) {
        return 1;
    }
    for (i = 0; i < profile->npts; ++i) {
        ctx->LP[i] = HUGE_VAL;
    }

    i_first = 0;
    while (i_first < profile->npts && profile->radius[i_first] < r_min_domain) {
        i_first += 1;
    }
    i_last = profile->npts - 1;
    while (i_last >= 0 && profile->radius[i_last] > r_max_domain) {
        i_last -= 1;
    }
    if (i_first > i_last) {
        return 0;
    }
    i_lo = i_first > 0 ? i_first - 1 : 0;
    i_hi = i_last < profile->npts - 1 ? i_last + 1 : profile->npts - 1;
    if (i_hi <= i_lo) {
        return 0;
    }

    pressure = (double *)malloc((size_t)profile->npts * sizeof(*pressure));
    if (pressure == 0) {
        free(ctx->LP);
        ctx->LP = 0;
        return 1;
    }
    for (i = 0; i < profile->npts; ++i) {
        pressure[i] = 0.0;
    }
    for (i = i_lo; i <= i_hi; ++i) {
        double eos_q[PRJ_EOS_NQUANT];

        if (profile->rho[i] <= 0.0) {
            continue;
        }
        prj_eos_rty(eos, profile->rho[i],
            profile->temp[i] / PRJ_CCSN_KELVIN_PER_MEV, profile->ye[i], eos_q, PRJ_EOS_CTX_MAIN);
        pressure[i] = eos_q[PRJ_EOS_PRESSURE];
    }
    for (i = i_lo; i <= i_hi; ++i) {
        double P0;
        double P1;
        double dr;
        double dlnP;

        if (i == i_lo) {
            P0 = pressure[i];
            P1 = pressure[i + 1];
            dr = profile->radius[i + 1] - profile->radius[i];
        } else if (i == i_hi) {
            P0 = pressure[i - 1];
            P1 = pressure[i];
            dr = profile->radius[i] - profile->radius[i - 1];
        } else {
            P0 = pressure[i - 1];
            P1 = pressure[i + 1];
            dr = profile->radius[i + 1] - profile->radius[i - 1];
        }
        if (P0 <= 0.0 || P1 <= 0.0 || dr <= 0.0) {
            ctx->LP[i] = HUGE_VAL;
            continue;
        }
        dlnP = log(P1) - log(P0);
        if (dlnP == 0.0) {
            ctx->LP[i] = HUGE_VAL;
        } else {
            ctx->LP[i] = fabs(dr / dlnP);
        }
    }
    free(pressure);
    return 0;
}

static void prj_ccsn_init_amr_ctx_free(prj_ccsn_init_amr_ctx *ctx)
{
    if (ctx == 0) {
        return;
    }
    free(ctx->LP);
    ctx->LP = 0;
    ctx->radius = 0;
    ctx->npts = 0;
    ctx->scale_factor = 0.0;
}

static double prj_ccsn_init_amr_LP_interp(const prj_ccsn_init_amr_ctx *ctx, double r)
{
    int lo;
    int hi;

    if (ctx == 0 || ctx->npts <= 0) {
        return HUGE_VAL;
    }
    if (r <= ctx->radius[0]) {
        return ctx->LP[0];
    }
    if (r >= ctx->radius[ctx->npts - 1]) {
        return ctx->LP[ctx->npts - 1];
    }
    lo = 0;
    hi = ctx->npts - 1;
    while (hi - lo > 1) {
        int mid = lo + (hi - lo) / 2;

        if (ctx->radius[mid] <= r) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    {
        double w = (r - ctx->radius[lo]) / (ctx->radius[hi] - ctx->radius[lo]);

        return (1.0 - w) * ctx->LP[lo] + w * ctx->LP[hi];
    }
}

static double prj_ccsn_init_amr_LP_min(const prj_ccsn_init_amr_ctx *ctx, double rmin, double rmax)
{
    double Lmin;
    double Lend;
    int i;

    if (ctx == 0 || ctx->npts <= 0) {
        return HUGE_VAL;
    }
    if (rmax < rmin) {
        double t = rmin;

        rmin = rmax;
        rmax = t;
    }
    Lmin = prj_ccsn_init_amr_LP_interp(ctx, rmin);
    Lend = prj_ccsn_init_amr_LP_interp(ctx, rmax);
    if (Lend < Lmin) {
        Lmin = Lend;
    }
    for (i = 0; i < ctx->npts; ++i) {
        if (ctx->radius[i] > rmin && ctx->radius[i] < rmax && ctx->LP[i] < Lmin) {
            Lmin = ctx->LP[i];
        }
    }
    return Lmin;
}

static int prj_ccsn_init_amr_cell_offset(int idx)
{
    if (idx == 0) {
        return 0;
    }
    if (idx == 1) {
        return 1;
    }
    if (idx == 2) {
        return PRJ_BLOCK_SIZE - 1;
    }
    return PRJ_BLOCK_SIZE;
}

static int prj_ccsn_init_amr_refine_block(const prj_block *block, void *userdata)
{
    const prj_ccsn_init_amr_ctx *ctx = (const prj_ccsn_init_amr_ctx *)userdata;
    double rmin = HUGE_VAL;
    double rmax = 0.0;
    double dx;
    double Lmin;
    int ix;
    int iy;
    int iz;

    if (block == 0 || ctx == 0 || ctx->npts <= 0) {
        return 0;
    }
    dx = block->dx[0];
    if (block->dx[1] > dx) {
        dx = block->dx[1];
    }
    if (block->dx[2] > dx) {
        dx = block->dx[2];
    }
    for (ix = 0; ix < 4; ++ix) {
        const int xoff = prj_ccsn_init_amr_cell_offset(ix);
        double x = block->xmin[0] + (double)xoff * block->dx[0];

        for (iy = 0; iy < 4; ++iy) {
            const int yoff = prj_ccsn_init_amr_cell_offset(iy);
            double y = block->xmin[1] + (double)yoff * block->dx[1];

            for (iz = 0; iz < 4; ++iz) {
                const int zoff = prj_ccsn_init_amr_cell_offset(iz);
                double z = block->xmin[2] + (double)zoff * block->dx[2];
                double r = sqrt(x * x + y * y + z * z);

                if (r < rmin) {
                    rmin = r;
                }
                if (r > rmax) {
                    rmax = r;
                }
            }
        }
    }
    Lmin = prj_ccsn_init_amr_LP_min(ctx, rmin, rmax);
    return dx > ctx->scale_factor * Lmin ? 1 : 0;
}

static int prj_problem_local_block(const prj_mpi *mpi, const prj_block *block)
{
    return block != 0 && block->id >= 0 && block->active == 1 &&
        (mpi == 0 || block->rank == mpi->rank);
}

static unsigned long long prj_problem_mesh_signature(const prj_mesh *mesh)
{
    unsigned long long sig = 1469598103934665603ULL;
    int bidx;
    int oct;

    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        const prj_block *block = &mesh->blocks[bidx];

        sig ^= (unsigned long long)(unsigned int)(block->id + 3);
        sig *= 1099511628211ULL;
        sig ^= (unsigned long long)(unsigned int)(block->level + 7);
        sig *= 1099511628211ULL;
        sig ^= (unsigned long long)(unsigned int)(block->active + 11);
        sig *= 1099511628211ULL;
        sig ^= (unsigned long long)(unsigned int)(block->parent + 13);
        sig *= 1099511628211ULL;
        for (oct = 0; oct < 8; ++oct) {
            sig ^= (unsigned long long)(unsigned int)(block->children[oct] + 17 + oct);
            sig *= 1099511628211ULL;
        }
    }
    return sig;
}

static void prj_problem_store_cell(prj_block *block, int i, int j, int k, const double *W, const double *U)
{
    int v;

    for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
        block->W[WIDX(v, i, j, k)] = W[v];
        prj_block_prim_stage(block, 1)[WIDX(v, i, j, k)] = W[v];
    }
    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        block->U[VIDX(v, i, j, k)] = U[v];
    }
}

static double prj_ccsn_kelvin_to_mev(double temperature_kelvin)
{
    return temperature_kelvin / PRJ_CCSN_KELVIN_PER_MEV;
}

static void prj_ccsn_profile_free(prj_ccsn_profile *profile)
{
    if (profile == 0) {
        return;
    }
    free(profile->radius);
    free(profile->rho);
    free(profile->temp);
    free(profile->ye);
    free(profile->vr);
    memset(profile, 0, sizeof(*profile));
}

static int prj_ccsn_profile_load(prj_ccsn_profile *profile, const char *filename)
{
    FILE *fp;
    int capacity;

    if (profile == 0 || filename == 0) {
        return 1;
    }

    memset(profile, 0, sizeof(*profile));
    fp = fopen(filename, "r");
    if (fp == 0) {
        return 1;
    }

    capacity = 1024;
    profile->radius = (double *)malloc((size_t)capacity * sizeof(*profile->radius));
    profile->rho = (double *)malloc((size_t)capacity * sizeof(*profile->rho));
    profile->temp = (double *)malloc((size_t)capacity * sizeof(*profile->temp));
    profile->ye = (double *)malloc((size_t)capacity * sizeof(*profile->ye));
    profile->vr = (double *)malloc((size_t)capacity * sizeof(*profile->vr));
    if (profile->radius == 0 || profile->rho == 0 || profile->temp == 0 ||
        profile->ye == 0 || profile->vr == 0) {
        fclose(fp);
        prj_ccsn_profile_free(profile);
        return 1;
    }

    while (!feof(fp)) {
        double index_col;
        double mass_col;
        double radius_col;
        double rho_col;
        double temp_col;
        double ye_col;
        double vr_col;
        double entropy_col;
        int scanned;

        scanned = fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf",
            &index_col, &mass_col, &radius_col, &rho_col, &temp_col, &ye_col, &vr_col, &entropy_col);
        if (scanned == EOF) {
            break;
        }
        if (scanned != 8) {
            fclose(fp);
            prj_ccsn_profile_free(profile);
            return 1;
        }
        if (profile->npts >= capacity) {
            int next_capacity = capacity * 2;
            double *next_radius;
            double *next_rho;
            double *next_temp;
            double *next_ye;
            double *next_vr;

            next_radius = (double *)malloc((size_t)next_capacity * sizeof(*next_radius));
            next_rho = (double *)malloc((size_t)next_capacity * sizeof(*next_rho));
            next_temp = (double *)malloc((size_t)next_capacity * sizeof(*next_temp));
            next_ye = (double *)malloc((size_t)next_capacity * sizeof(*next_ye));
            next_vr = (double *)malloc((size_t)next_capacity * sizeof(*next_vr));
            if (next_radius == 0 || next_rho == 0 || next_temp == 0 || next_ye == 0 || next_vr == 0) {
                fclose(fp);
                free(next_radius);
                free(next_rho);
                free(next_temp);
                free(next_ye);
                free(next_vr);
                prj_ccsn_profile_free(profile);
                return 1;
            }
            memcpy(next_radius, profile->radius, (size_t)profile->npts * sizeof(*next_radius));
            memcpy(next_rho, profile->rho, (size_t)profile->npts * sizeof(*next_rho));
            memcpy(next_temp, profile->temp, (size_t)profile->npts * sizeof(*next_temp));
            memcpy(next_ye, profile->ye, (size_t)profile->npts * sizeof(*next_ye));
            memcpy(next_vr, profile->vr, (size_t)profile->npts * sizeof(*next_vr));
            free(profile->radius);
            free(profile->rho);
            free(profile->temp);
            free(profile->ye);
            free(profile->vr);
            profile->radius = next_radius;
            profile->rho = next_rho;
            profile->temp = next_temp;
            profile->ye = next_ye;
            profile->vr = next_vr;
            capacity = next_capacity;
        }

        (void)index_col;
        (void)mass_col;
        (void)entropy_col;
        profile->radius[profile->npts] = radius_col;
        profile->rho[profile->npts] = rho_col;
        profile->temp[profile->npts] = temp_col;
        profile->ye[profile->npts] = ye_col;
        profile->vr[profile->npts] = vr_col;
        profile->npts += 1;
    }

    fclose(fp);
    return profile->npts > 1 ? 0 : 1;
}

static void prj_ccsn_profile_sample(const prj_ccsn_profile *profile, double r,
    double *rho, double *temp, double *ye, double *vr)
{
    int lo;
    int hi;

    if (profile->npts <= 0) {
        *rho = 1.0;
        *temp = 1.0e9;
        *ye = 0.5;
        *vr = 0.0;
        return;
    }
    if (r <= profile->radius[0]) {
        *rho = profile->rho[0];
        *temp = profile->temp[0];
        *ye = profile->ye[0];
        *vr = profile->vr[0];
        return;
    }
    if (r >= profile->radius[profile->npts - 1]) {
        *rho = profile->rho[profile->npts - 1];
        *temp = profile->temp[profile->npts - 1];
        *ye = profile->ye[profile->npts - 1];
        *vr = profile->vr[profile->npts - 1];
        return;
    }

    lo = 0;
    hi = profile->npts - 1;
    while (hi - lo > 1) {
        int mid = lo + (hi - lo) / 2;

        if (profile->radius[mid] <= r) {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    {
        double w = (r - profile->radius[lo]) / (profile->radius[hi] - profile->radius[lo]);

        *rho = (1.0 - w) * profile->rho[lo] + w * profile->rho[hi];
        *temp = (1.0 - w) * profile->temp[lo] + w * profile->temp[hi];
        *ye = (1.0 - w) * profile->ye[lo] + w * profile->ye[hi];
        *vr = (1.0 - w) * profile->vr[lo] + w * profile->vr[hi];
    }
}

/* Volume-average the primitive state over a cell using the same 3x3x3
   Gauss-Legendre quadrature (about the mesh center of mass) that
   prj_mesh_update_block_r_com uses for r_com. At every quadrature node we
   compute the radius, sample the progenitor, and convert to primitives
   (density, radially-projected velocity, internal energy, electron fraction),
   then accumulate the volume average. This is more accurate than evaluating
   the primitives once at the cell-averaged radius r_com. */
static void prj_ccsn_quadrature_primitives(prj_sim *sim, const prj_ccsn_profile *profile,
    const prj_block *block, int i, int j, int k, double *W)
{
    static const double gq_node[3] = {
        -0.77459666924148337704, 0.0, 0.77459666924148337704
    };
    static const double gq_wnorm[3] = {5.0 / 18.0, 8.0 / 18.0, 5.0 / 18.0};
    double x_com[3] = {0.0, 0.0, 0.0};
    double xc;
    double yc;
    double zc;
    double hx;
    double hy;
    double hz;
    int a;
    int b;
    int c;
    int v;

    for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
        W[v] = 0.0;
    }

    x_com[0] = sim->mesh.x_com[0];
    x_com[1] = sim->mesh.x_com[1];
    x_com[2] = sim->mesh.x_com[2];

    xc = block->xmin[0] + ((double)i + 0.5) * block->dx[0] - x_com[0];
    yc = block->xmin[1] + ((double)j + 0.5) * block->dx[1] - x_com[1];
    zc = block->xmin[2] + ((double)k + 0.5) * block->dx[2] - x_com[2];
    hx = 0.5 * block->dx[0];
    hy = 0.5 * block->dx[1];
    hz = 0.5 * block->dx[2];

    for (a = 0; a < 3; ++a) {
        double dx1 = xc + hx * gq_node[a];

        for (b = 0; b < 3; ++b) {
            double dx2 = yc + hy * gq_node[b];
            double wab = gq_wnorm[a] * gq_wnorm[b];

            for (c = 0; c < 3; ++c) {
                double dx3 = zc + hz * gq_node[c];
                double w = wab * gq_wnorm[c];
                double r = sqrt(dx1 * dx1 + dx2 * dx2 + dx3 * dx3);
                double rho;
                double temp;
                double ye;
                double vr;
                double eos_q[PRJ_EOS_NQUANT];

                prj_ccsn_profile_sample(profile, r, &rho, &temp, &ye, &vr);
                if (rho == 0.0) {
                    fprintf(stderr,
                        "prj_ccsn_quadrature_primitives: rho=0 before prj_eos_rty for current_block id=%d current_rank=%d level=%d "
                        "cell=(%d,%d,%d) node=(%d,%d,%d) x=(%.17e, %.17e, %.17e) r=%.17e temp=%.17e ye=%.17e\n",
                        block->id, block->rank, block->level, i, j, k, a, b, c,
                        dx1 + x_com[0], dx2 + x_com[1], dx3 + x_com[2], r, temp, ye);
                    exit(EXIT_FAILURE);
                }
                prj_eos_rty(&sim->eos, rho, prj_ccsn_kelvin_to_mev(temp), ye, eos_q, PRJ_EOS_CTX_MAIN);

                W[PRJ_PRIM_RHO] += w * rho;
                if (r > 0.0) {
                    W[PRJ_PRIM_V1] += w * vr * dx1 / r;
                    W[PRJ_PRIM_V2] += w * vr * dx2 / r;
                    W[PRJ_PRIM_V3] += w * vr * dx3 / r;
                }
                W[PRJ_PRIM_EINT] += w * eos_q[PRJ_EOS_EINT];
                W[PRJ_PRIM_YE] += w * ye;
            }
        }
    }
}

static void prj_ccsn_fill_mesh(prj_sim *sim, const prj_mpi *mpi, const prj_ccsn_profile *profile)
{
    int bidx;

    for (bidx = 0; bidx < sim->mesh.nblocks; ++bidx) {
        prj_block *block = &sim->mesh.blocks[bidx];
        int i;
        int j;
        int k;

        if (!prj_problem_local_block(mpi, block)) {
            continue;
        }
        for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
            for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
                for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                    double W[PRJ_NVAR_PRIM] = {0.0};
                    double U[PRJ_NVAR_CONS] = {0.0};

                    /* Volume-averaged primitives over the cell, using the same
                       quadrature as the r_com calculation. */
                    prj_ccsn_quadrature_primitives(sim, profile, block, i, j, k, W);
                    prj_eos_prim2cons(&sim->eos, W, U);
                    prj_problem_store_cell(block, i, j, k, W, U);
                }
            }
        }
    }
}

static void prj_ccsn_initialize_amr(prj_sim *sim, prj_mpi *mpi, const prj_ccsn_profile *profile)
{
    unsigned long long prev_sig;
    unsigned long long next_sig;
    prj_ccsn_init_amr_ctx init_ctx;
    int init_ctx_ok;

    init_ctx.npts = 0;
    init_ctx.radius = 0;
    init_ctx.LP = 0;
    init_ctx.scale_factor = sim->mesh.amr_init_scale_factor;
    {
        double ax = fabs(sim->coord.x1min) > fabs(sim->coord.x1max) ?
            fabs(sim->coord.x1min) : fabs(sim->coord.x1max);
        double ay = fabs(sim->coord.x2min) > fabs(sim->coord.x2max) ?
            fabs(sim->coord.x2min) : fabs(sim->coord.x2max);
        double az = fabs(sim->coord.x3min) > fabs(sim->coord.x3max) ?
            fabs(sim->coord.x3min) : fabs(sim->coord.x3max);
        double cx = sim->coord.x1min > 0.0 ? sim->coord.x1min :
            (sim->coord.x1max < 0.0 ? sim->coord.x1max : 0.0);
        double cy = sim->coord.x2min > 0.0 ? sim->coord.x2min :
            (sim->coord.x2max < 0.0 ? sim->coord.x2max : 0.0);
        double cz = sim->coord.x3min > 0.0 ? sim->coord.x3min :
            (sim->coord.x3max < 0.0 ? sim->coord.x3max : 0.0);
        double r_max_domain = sqrt(ax * ax + ay * ay + az * az);
        double r_min_domain = sqrt(cx * cx + cy * cy + cz * cz);

        init_ctx_ok = (prj_ccsn_init_amr_ctx_build(&init_ctx, profile,
            sim->mesh.amr_init_scale_factor, &sim->eos,
            r_min_domain, r_max_domain) == 0);
    }

    prj_ccsn_fill_mesh(sim, mpi, profile);
    prj_eos_fill_mesh(&sim->mesh, &sim->eos, mpi, 1, PRJ_EOS_CTX_MAIN);
#if PRJ_USE_GRAVITY
    prj_gravity_init(sim, mpi);
    prj_gravity_monopole_reduce(&sim->mesh, &sim->grav, mpi, 1);
    prj_gravity_monopole_integrate(&sim->mesh, &sim->grav, mpi);
#endif
    if (sim->mesh.max_level == 0) {
        prj_ccsn_init_amr_ctx_free(&init_ctx);
        return;
    }


    prj_eos_fill_active_cells(&sim->mesh, &sim->eos, mpi, 1, PRJ_EOS_CTX_MAIN);
    prj_boundary_fill_ghosts(&sim->mesh, mpi, &sim->bc, 1);
    prj_eos_fill_mesh(&sim->mesh, &sim->eos, mpi, 1, PRJ_EOS_CTX_MAIN);
#if PRJ_USE_GRAVITY
    prj_gravity_monopole_reduce(&sim->mesh, &sim->grav, mpi, 1);
    prj_gravity_monopole_integrate(&sim->mesh, &sim->grav, mpi);
#endif

    if (init_ctx_ok) {
        sim->mesh.amr_init_refine_fn = prj_ccsn_init_amr_refine_block;
        sim->mesh.amr_init_refine_userdata = &init_ctx;
    }

    do {
        prev_sig = prj_problem_mesh_signature(&sim->mesh);
        prj_amr_adapt(&sim->mesh, &sim->eos, mpi);
        prj_mpi_rebalance(&sim->mesh, mpi);
    #if PRJ_USE_GRAVITY
        prj_gravity_rebuild_grid(sim, mpi);
    #endif
        prj_ccsn_fill_mesh(sim, mpi, profile);

        prj_eos_fill_active_cells(&sim->mesh, &sim->eos, mpi, 1, PRJ_EOS_CTX_MAIN);
        prj_boundary_fill_ghosts(&sim->mesh, mpi, &sim->bc, 1);
        prj_eos_fill_mesh(&sim->mesh, &sim->eos, mpi, 1, PRJ_EOS_CTX_MAIN);
    #if PRJ_USE_GRAVITY
        prj_gravity_monopole_reduce(&sim->mesh, &sim->grav, mpi, 1);
        prj_gravity_monopole_integrate(&sim->mesh, &sim->grav, mpi);
    #endif

        next_sig = prj_problem_mesh_signature(&sim->mesh);
    } while ((int)prj_mpi_global_sum(mpi, (double)(next_sig != prev_sig ? 1 : 0)) != 0);

    sim->mesh.amr_init_refine_fn = 0;
    sim->mesh.amr_init_refine_userdata = 0;
    prj_ccsn_init_amr_ctx_free(&init_ctx);

    prj_amr_init_neighbors(&sim->mesh);
    prj_mpi_prepare(&sim->mesh, mpi);
}

void prj_problem_ccsn(prj_sim *sim, prj_mpi *mpi)
{
    prj_ccsn_profile profile;

    if (sim->progenitor_file[0] == '\0') {
        return;
    }
    if (prj_mesh_init(&sim->mesh, sim->mesh.root_nx[0], sim->mesh.root_nx[1], sim->mesh.root_nx[2],
        sim->mesh.max_level, &sim->coord) != 0) {
        return;
    }
    prj_mpi_decompose(&sim->mesh, mpi);
    prj_mpi_prepare(&sim->mesh, mpi);

    if (prj_ccsn_profile_load(&profile, sim->progenitor_file) != 0) {
        return;
    }
    prj_ccsn_initialize_amr(sim, mpi, &profile);
    prj_mhd_init(sim, mpi);
    prj_ccsn_profile_free(&profile);
}
