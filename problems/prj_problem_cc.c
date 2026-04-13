#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "prj.h"

#define PRJ_CC_EOS_PATH "../eos_tmp/SFHoEOS__ye__0.035_0.56_50__logT_-4.793_2.176_500__logrho_-8.699_15.5_500_extend.dat"
#define PRJ_CC_KELVIN_PER_MEV 1.160451812e10

static const double prj_cc_Hp_scale = 0.5;

typedef struct prj_cc_profile {
    int npts;
    double *radius;
    double *rho;
    double *temp;
    double *ye;
    double *vr;
    double *pressure;
    double *Hp;
} prj_cc_profile;

static int prj_problem_local_block(const prj_block *block)
{
    prj_mpi *mpi = prj_mpi_current();

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
        block->W[VIDX(v, i, j, k)] = W[v];
        block->W1[VIDX(v, i, j, k)] = W[v];
    }
    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        block->U[VIDX(v, i, j, k)] = U[v];
    }
}

static double prj_cc_kelvin_to_mev(double temperature_kelvin)
{
    return temperature_kelvin / PRJ_CC_KELVIN_PER_MEV;
}

static void prj_cc_profile_free(prj_cc_profile *profile)
{
    if (profile == 0) {
        return;
    }
    free(profile->radius);
    free(profile->rho);
    free(profile->temp);
    free(profile->ye);
    free(profile->vr);
    free(profile->pressure);
    free(profile->Hp);
    memset(profile, 0, sizeof(*profile));
}

static int prj_cc_profile_load(prj_cc_profile *profile, const char *filename)
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
    profile->pressure = (double *)malloc((size_t)capacity * sizeof(*profile->pressure));
    profile->Hp = (double *)malloc((size_t)capacity * sizeof(*profile->Hp));
    if (profile->radius == 0 || profile->rho == 0 || profile->temp == 0 ||
        profile->ye == 0 || profile->vr == 0 || profile->pressure == 0 || profile->Hp == 0) {
        fclose(fp);
        prj_cc_profile_free(profile);
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
            prj_cc_profile_free(profile);
            return 1;
        }
        if (profile->npts >= capacity) {
            int next_capacity = capacity * 2;
            double *next_radius;
            double *next_rho;
            double *next_temp;
            double *next_ye;
            double *next_vr;
            double *next_pressure;
            double *next_Hp;

            next_radius = (double *)malloc((size_t)next_capacity * sizeof(*next_radius));
            next_rho = (double *)malloc((size_t)next_capacity * sizeof(*next_rho));
            next_temp = (double *)malloc((size_t)next_capacity * sizeof(*next_temp));
            next_ye = (double *)malloc((size_t)next_capacity * sizeof(*next_ye));
            next_vr = (double *)malloc((size_t)next_capacity * sizeof(*next_vr));
            next_pressure = (double *)malloc((size_t)next_capacity * sizeof(*next_pressure));
            next_Hp = (double *)malloc((size_t)next_capacity * sizeof(*next_Hp));
            if (next_radius == 0 || next_rho == 0 || next_temp == 0 || next_ye == 0 || next_vr == 0 ||
                next_pressure == 0 || next_Hp == 0) {
                fclose(fp);
                free(next_radius);
                free(next_rho);
                free(next_temp);
                free(next_ye);
                free(next_vr);
                free(next_pressure);
                free(next_Hp);
                prj_cc_profile_free(profile);
                return 1;
            }
            memcpy(next_radius, profile->radius, (size_t)profile->npts * sizeof(*next_radius));
            memcpy(next_rho, profile->rho, (size_t)profile->npts * sizeof(*next_rho));
            memcpy(next_temp, profile->temp, (size_t)profile->npts * sizeof(*next_temp));
            memcpy(next_ye, profile->ye, (size_t)profile->npts * sizeof(*next_ye));
            memcpy(next_vr, profile->vr, (size_t)profile->npts * sizeof(*next_vr));
            memcpy(next_pressure, profile->pressure, (size_t)profile->npts * sizeof(*next_pressure));
            memcpy(next_Hp, profile->Hp, (size_t)profile->npts * sizeof(*next_Hp));
            free(profile->radius);
            free(profile->rho);
            free(profile->temp);
            free(profile->ye);
            free(profile->vr);
            free(profile->pressure);
            free(profile->Hp);
            profile->radius = next_radius;
            profile->rho = next_rho;
            profile->temp = next_temp;
            profile->ye = next_ye;
            profile->vr = next_vr;
            profile->pressure = next_pressure;
            profile->Hp = next_Hp;
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

static int prj_cc_profile_prepare_pressure_scale_height(prj_cc_profile *profile, prj_eos *eos)
{
    int i;

    if (profile == 0 || eos == 0 || profile->npts <= 1) {
        return 1;
    }

    for (i = 0; i < profile->npts; ++i) {
        double eos_q[PRJ_EOS_NQUANT];

        prj_eos_rty(eos,
            profile->rho[i],
            prj_cc_kelvin_to_mev(profile->temp[i]),
            profile->ye[i],
            eos_q);
        profile->pressure[i] = eos_q[PRJ_EOS_PRESSURE];
    }

    for (i = 0; i < profile->npts; ++i) {
        int im = i > 0 ? i - 1 : i;
        int ip = i + 1 < profile->npts ? i + 1 : i;
        double dr;
        double dPdr;
        double P = profile->pressure[i];

        if (ip == im) {
            profile->Hp[i] = 1.0e99;
            continue;
        }
        dr = profile->radius[ip] - profile->radius[im];
        if (fabs(dr) <= 0.0) {
            profile->Hp[i] = 1.0e99;
            continue;
        }
        dPdr = (profile->pressure[ip] - profile->pressure[im]) / dr;
        if (fabs(dPdr) <= 1.0e-300) {
            profile->Hp[i] = 1.0e99;
        } else {
            profile->Hp[i] = fabs(P / dPdr);
        }
    }

    return 0;
}

static void prj_cc_profile_sample(const prj_cc_profile *profile, double r,
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

static void prj_cc_fill_mesh(prj_sim *sim, const prj_cc_profile *profile)
{
    int bidx;

    for (bidx = 0; bidx < sim->mesh.nblocks; ++bidx) {
        prj_block *block = &sim->mesh.blocks[bidx];
        int i;
        int j;
        int k;

        if (!prj_problem_local_block(block)) {
            continue;
        }
        for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
            for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
                for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                    double x1 = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
                    double x2 = block->xmin[1] + ((double)j + 0.5) * block->dx[1];
                    double x3 = block->xmin[2] + ((double)k + 0.5) * block->dx[2];
                    double r = sqrt(x1 * x1 + x2 * x2 + x3 * x3);
                    double rho;
                    double temp;
                    double ye;
                    double vr;
                    double eos_q[PRJ_EOS_NQUANT];
                    double W[PRJ_NVAR_PRIM];
                    double U[PRJ_NVAR_CONS];

                    prj_cc_profile_sample(profile, r, &rho, &temp, &ye, &vr);
                    prj_eos_rty(&sim->eos, rho, prj_cc_kelvin_to_mev(temp), ye, eos_q);
                    W[PRJ_PRIM_RHO] = rho;
                    if (r > 0.0) {
                        W[PRJ_PRIM_V1] = vr * x1 / r;
                        W[PRJ_PRIM_V2] = vr * x2 / r;
                        W[PRJ_PRIM_V3] = vr * x3 / r;
                    } else {
                        W[PRJ_PRIM_V1] = 0.0;
                        W[PRJ_PRIM_V2] = 0.0;
                        W[PRJ_PRIM_V3] = 0.0;
                    }
                    W[PRJ_PRIM_EINT] = eos_q[PRJ_EOS_EINT];
                    W[PRJ_PRIM_YE] = ye;
                    prj_eos_prim2cons(&sim->eos, W, U);
                    prj_problem_store_cell(block, i, j, k, W, U);
                }
            }
        }
    }
}

static void prj_cc_block_radius_range(const prj_block *block, double *rmin, double *rmax)
{
    double local_rmin = 1.0e99;
    double local_rmax = 0.0;
    int ix;
    int iy;
    int iz;

    for (ix = 0; ix < 2; ++ix) {
        for (iy = 0; iy < 2; ++iy) {
            for (iz = 0; iz < 2; ++iz) {
                double x = ix == 0 ? block->xmin[0] : block->xmax[0];
                double y = iy == 0 ? block->xmin[1] : block->xmax[1];
                double z = iz == 0 ? block->xmin[2] : block->xmax[2];
                double r = sqrt(x * x + y * y + z * z);

                if (r < local_rmin) {
                    local_rmin = r;
                }
                if (r > local_rmax) {
                    local_rmax = r;
                }
            }
        }
    }

    if (block->xmin[0] <= 0.0 && block->xmax[0] >= 0.0 &&
        block->xmin[1] <= 0.0 && block->xmax[1] >= 0.0 &&
        block->xmin[2] <= 0.0 && block->xmax[2] >= 0.0) {
        local_rmin = 0.0;
    }

    *rmin = local_rmin;
    *rmax = local_rmax;
}

static double prj_cc_profile_min_Hp_in_range(const prj_cc_profile *profile, double rmin, double rmax)
{
    double min_Hp = 1.0e99;
    int i;

    if (profile == 0 || profile->npts <= 0) {
        return min_Hp;
    }
    if (rmax <= profile->radius[0]) {
        return profile->Hp[0];
    }
    if (rmin >= profile->radius[profile->npts - 1]) {
        return profile->Hp[profile->npts - 1];
    }

    for (i = 0; i < profile->npts; ++i) {
        if (profile->radius[i] >= rmin && profile->radius[i] <= rmax &&
            profile->Hp[i] < min_Hp) {
            min_Hp = profile->Hp[i];
        }
    }

    if (min_Hp >= 1.0e98) {
        int idx = 0;

        while (idx + 1 < profile->npts && profile->radius[idx + 1] < rmin) {
            idx += 1;
        }
        min_Hp = profile->Hp[idx];
        if (idx + 1 < profile->npts && profile->Hp[idx + 1] < min_Hp) {
            min_Hp = profile->Hp[idx + 1];
        }
    }
    return min_Hp;
}

static void prj_cc_initialize_amr(prj_sim *sim, const prj_cc_profile *profile)
{
    int level;

    prj_cc_fill_mesh(sim, profile);
    prj_eos_fill_mesh(&sim->mesh, &sim->eos, 1);
#if PRJ_USE_GRAVITY
    prj_gravity_init(sim);
    prj_gravity_monopole_reduce(&sim->mesh, 1);
    prj_gravity_monopole_integrate(&sim->mesh);
#endif
    if (sim->mesh.max_level <= 0) {
        return;
    }
    for (level = 0; level < sim->mesh.max_level; ++level) {
        int bidx;
        unsigned long long prev_sig = prj_problem_mesh_signature(&sim->mesh);

        for (bidx = 0; bidx < sim->mesh.nblocks; ++bidx) {
            prj_block *block = &sim->mesh.blocks[bidx];

            if (block->id >= 0 && block->active == 1 && block->level == level) {
                double rmin;
                double rmax;
                double min_Hp;
                double cell_size = block->dx[0];

                if (block->dx[1] > cell_size) {
                    cell_size = block->dx[1];
                }
                if (block->dx[2] > cell_size) {
                    cell_size = block->dx[2];
                }
                prj_cc_block_radius_range(block, &rmin, &rmax);
                min_Hp = prj_cc_profile_min_Hp_in_range(profile, rmin, rmax);
                if (cell_size > prj_cc_Hp_scale * min_Hp) {
                    prj_amr_refine_block(&sim->mesh, bidx);
                }
            }
        }
        prj_mpi_rebalance(&sim->mesh);
        prj_cc_fill_mesh(sim, profile);
        prj_eos_fill_mesh(&sim->mesh, &sim->eos, 1);
#if PRJ_USE_GRAVITY
        prj_gravity_monopole_reduce(&sim->mesh, 1);
        prj_gravity_monopole_integrate(&sim->mesh);
#endif
        if (prj_problem_mesh_signature(&sim->mesh) == prev_sig) {
            break;
        }
    }
    prj_amr_init_neighbors(&sim->mesh);
}

void prj_problem_cc(prj_sim *sim)
{
    prj_cc_profile profile;

    if (sim->progenitor_file[0] == '\0') {
        return;
    }
    if (sim->eos.kind == PRJ_EOS_KIND_TABLE && sim->eos.filename[0] == '\0') {
        strncpy(sim->eos.filename, PRJ_CC_EOS_PATH, sizeof(sim->eos.filename) - 1);
        sim->eos.filename[sizeof(sim->eos.filename) - 1] = '\0';
    }
    prj_eos_init(&sim->eos);
    if (prj_mesh_init(&sim->mesh, sim->mesh.root_nx[0], sim->mesh.root_nx[1], sim->mesh.root_nx[2],
        sim->mesh.max_level, &sim->coord) != 0) {
        return;
    }
    prj_mpi_decompose(&sim->mesh);
    prj_mpi_prepare(&sim->mesh, prj_mpi_current());

    if (prj_cc_profile_load(&profile, sim->progenitor_file) != 0) {
        return;
    }
    if (prj_cc_profile_prepare_pressure_scale_height(&profile, &sim->eos) != 0) {
        prj_cc_profile_free(&profile);
        return;
    }
    prj_cc_initialize_amr(sim, &profile);
    prj_cc_profile_free(&profile);
}
