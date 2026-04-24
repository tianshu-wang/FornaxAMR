#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "prj.h"

#define PRJ_CCSN_EOS_PATH "../eos_tmp/SFHoEOS__ye__0.035_0.56_50__logT_-4.793_2.176_500__logrho_-8.699_15.5_500_extend.dat"
#define PRJ_CCSN_KELVIN_PER_MEV 1.160451812e10

typedef struct prj_ccsn_profile {
    int npts;
    double *radius;
    double *rho;
    double *temp;
    double *ye;
    double *vr;
} prj_ccsn_profile;

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

static void prj_ccsn_fill_mesh(prj_sim *sim, const prj_ccsn_profile *profile)
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
                    int v;

                    for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
                        W[v] = 0.0;
                    }
                    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
                        U[v] = 0.0;
                    }

                    prj_ccsn_profile_sample(profile, r, &rho, &temp, &ye, &vr);
                    if (rho == 0.0) {
                        fprintf(stderr,
                            "prj_ccsn_fill_mesh: rho=0 before prj_eos_rty for current_block id=%d current_rank=%d level=%d "
                            "cell=(%d,%d,%d) x=(%.17e, %.17e, %.17e) r=%.17e temp=%.17e ye=%.17e\n",
                            block->id, block->rank, block->level, i, j, k, x1, x2, x3, r, temp, ye);
                        exit(EXIT_FAILURE);
                    }
                    prj_eos_rty(&sim->eos, rho, prj_ccsn_kelvin_to_mev(temp), ye, eos_q);
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

static void prj_ccsn_initialize_amr(prj_sim *sim, const prj_ccsn_profile *profile)
{
    unsigned long long prev_sig;
    unsigned long long next_sig;

    prj_ccsn_fill_mesh(sim, profile);
    prj_eos_fill_mesh(&sim->mesh, &sim->eos, 1);
#if PRJ_USE_GRAVITY
    prj_gravity_init(sim);
    prj_gravity_monopole_reduce(&sim->mesh, 1);
    prj_gravity_monopole_integrate(&sim->mesh);
#endif
    if (sim->mesh.max_level == 0) {
        return;
    }

    do {
        prev_sig = prj_problem_mesh_signature(&sim->mesh);
        prj_eos_fill_active_cells(&sim->mesh, &sim->eos, 1);
        prj_boundary_fill_ghosts(&sim->mesh, &sim->bc, 1);
        prj_eos_fill_mesh(&sim->mesh, &sim->eos, 1);
        prj_amr_adapt(&sim->mesh, &sim->eos);
        prj_mpi_rebalance(&sim->mesh);
        prj_ccsn_fill_mesh(sim, profile);
        prj_eos_fill_active_cells(&sim->mesh, &sim->eos, 1);
        prj_boundary_fill_ghosts(&sim->mesh, &sim->bc, 1);
        prj_eos_fill_mesh(&sim->mesh, &sim->eos, 1);
#if PRJ_USE_GRAVITY
        prj_gravity_monopole_reduce(&sim->mesh, 1);
        prj_gravity_monopole_integrate(&sim->mesh);
#endif
        next_sig = prj_problem_mesh_signature(&sim->mesh);
    } while ((int)prj_mpi_global_sum((double)(next_sig != prev_sig ? 1 : 0)) != 0);

    prj_amr_init_neighbors(&sim->mesh);
    prj_mpi_prepare(&sim->mesh, prj_mpi_current());
}

void prj_problem_ccsn(prj_sim *sim)
{
    prj_ccsn_profile profile;

    if (sim->progenitor_file[0] == '\0') {
        return;
    }
    if (sim->eos.kind == PRJ_EOS_KIND_TABLE && sim->eos.filename[0] == '\0') {
        strncpy(sim->eos.filename, PRJ_CCSN_EOS_PATH, sizeof(sim->eos.filename) - 1);
        sim->eos.filename[sizeof(sim->eos.filename) - 1] = '\0';
    }
    prj_eos_init(&sim->eos);
    if (prj_mesh_init(&sim->mesh, sim->mesh.root_nx[0], sim->mesh.root_nx[1], sim->mesh.root_nx[2],
        sim->mesh.max_level, &sim->coord) != 0) {
        return;
    }
    prj_mpi_decompose(&sim->mesh);
    prj_mpi_prepare(&sim->mesh, prj_mpi_current());

    if (prj_ccsn_profile_load(&profile, sim->progenitor_file) != 0) {
        return;
    }
    prj_ccsn_initialize_amr(sim, &profile);
    prj_mhd_init(sim);
    prj_ccsn_profile_free(&profile);
}
