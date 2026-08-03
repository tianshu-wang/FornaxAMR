#include <math.h>
#include "prj_problem_rad_common.h"

#define ARAD 7.5657e-15
#define KB_MEV 8.617333262145e-11
#define EOS_SCALE 0.95655684e18
#ifndef RAD_SCALE
#define RAD_SCALE 1.0e25
#endif
#ifndef RAD_TEST_VARIANT
#define RAD_TEST_VARIANT 1
#endif
#ifndef RAD_SHOCK_CASE
#define RAD_SHOCK_CASE 1
#endif

typedef struct rad_shock_state {
    double rho, p, ux, erad;
} rad_shock_state;

static double ideal_eint_from_k(double temp_k)
{
    return EOS_SCALE * (KB_MEV * temp_k) / (PRJ_IDEAL_GAMMA - 1.0);
}

static void store_cell(prj_sim *sim, prj_block *block, int i, int j, int k,
    double *W, double *U)
{
    prj_eos_cell_prim2cons(&sim->eos, &sim->mesh, block, 0, i, j, k,
        W, U, PRJ_EOS_CTX_MAIN);
    prj_block_store_prim_cell(block, 0, i, j, k, W);
    prj_block_store_prim_cell(block, 1, i, j, k, W);
    prj_block_store_cons_cell(block, i, j, k, U);
}

static double accretion_velocity(double r)
{
    if (r <= 7.0e6) return 0.0;
    if (r <= 8.0e6) return -0.2 * PRJ_CLIGHT * (r - 7.0e6) / 1.0e6;
    return -0.2 * PRJ_CLIGHT * (8.0e6 * 8.0e6) / (r * r);
}

static void fill_sphere_like(prj_sim *sim, int mode)
{
    int bidx;
    prj_mesh_init(&sim->mesh, sim->mesh.root_nx[0], sim->mesh.root_nx[1],
        sim->mesh.root_nx[2], 0, &sim->coord, 0);
    for (bidx = 0; bidx < sim->mesh.nblocks; ++bidx) {
        prj_block *block = &sim->mesh.blocks[bidx];
        int i, j, k;
        if (block->id < 0 || block->active != 1) continue;
        for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i)
            for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j)
                for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                    double r = fabs(block->xmin[0] + ((double)i + 0.5) * block->dx[0]);
                    double W[PRJ_NVAR_PRIM] = {0.0};
                    double U[PRJ_NVAR_CONS] = {0.0};
                    double rho0 = 9.0e14;
                    double temp_k;
                    double etot;
                    int inside = r <= 1.0e6;
                    int g;

                    if (mode == 4) {
                        etot = RAD_TEST_VARIANT == 1 ? 0.8 : 10.0;
                        temp_k = pow(etot / ARAD, 0.25);
                        if (!inside) {
                            temp_k = 1.0;
                            etot = ARAD;
                        }
                    } else {
                        temp_k = 5.0 / KB_MEV;
                        etot = (7.0 / 8.0) * ARAD * pow(temp_k, 4.0);
                    }
                    W[PRJ_PRIM_RHO] = inside ? rho0 : 1.0e-10 * rho0;
                    W[PRJ_PRIM_EINT] = ideal_eint_from_k(temp_k);
                    W[PRJ_PRIM_YE] = 0.5;
                    if (mode == 5 || (mode == 6 && RAD_TEST_VARIANT == 2))
                        W[PRJ_PRIM_V1] = accretion_velocity(r);
                    for (g = 0; g < PRJ_NEGROUP; ++g) {
                        double E = inside ? etot / PRJ_NEGROUP :
                            etot * 1.0e-20 / PRJ_NEGROUP;
                        W[PRJ_PRIM_RAD_E(0, g)] = E / RAD_SCALE;
                    }
                    store_cell(sim, block, i, j, k, W, U);
                }
    }
}

static void shock_table(int c, rad_shock_state *left, rad_shock_state *right)
{
    static const rad_shock_state l[4] = {
        {1, 3e-5, 0.0015, 1e-8}, {1, 4e-3, 0.25, 2e-5},
        {1, 60, 10, 2}, {1, 6e-3, 0.69, 0.18}};
    static const rad_shock_state r[4] = {
        {2.4, 1.61e-4, 6.25e-3, 2.51e-7},
        {3.11, 0.04512, 0.0804, 3.46e-3},
        {8, 2.34e3, 1.25, 1.14e3},
        {3.65, 3.59e-2, 0.189, 1.30}};
    *left = l[c - 1]; *right = r[c - 1];
}

static void fill_shock(prj_sim *sim)
{
    rad_shock_state left, right;
    double temp_mev_l, temp_k_l, rho_unit;
    int bidx;
    shock_table(RAD_SHOCK_CASE, &left, &right);
    temp_mev_l = (left.p / left.rho) * PRJ_CLIGHT * PRJ_CLIGHT / EOS_SCALE;
    temp_k_l = temp_mev_l / KB_MEV;
    rho_unit = ARAD * pow(temp_k_l, 4.0) /
        (left.erad * PRJ_CLIGHT * PRJ_CLIGHT);
    prj_mesh_init(&sim->mesh, sim->mesh.root_nx[0], sim->mesh.root_nx[1],
        sim->mesh.root_nx[2], 0, &sim->coord, 0);
#if PRJ_DYNAMIC_GR
    prj_z4c_init_mesh_flat(&sim->mesh, 0);
#endif
    for (bidx = 0; bidx < sim->mesh.nblocks; ++bidx) {
        prj_block *block = &sim->mesh.blocks[bidx];
        int i, j, k;
        if (block->id < 0 || block->active != 1) continue;
        for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i)
            for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j)
                for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                    double x = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
                    const rad_shock_state *s = x < 0.0 ? &left : &right;
                    double W[PRJ_NVAR_PRIM] = {0.0};
                    double U[PRJ_NVAR_CONS] = {0.0};
                    double gamma_lor = sqrt(1.0 + s->ux * s->ux);
                    double v = PRJ_CLIGHT * s->ux / gamma_lor;
                    double Ecom = s->erad * rho_unit * PRJ_CLIGHT * PRJ_CLIGHT;
                    double Elab = Ecom * (4.0 * gamma_lor * gamma_lor - 1.0) / 3.0;
                    double Flab = (4.0 / 3.0) * gamma_lor * gamma_lor * v * Ecom;
                    int g;
                    W[PRJ_PRIM_RHO] = s->rho * rho_unit;
                    W[PRJ_PRIM_V1] = v;
                    W[PRJ_PRIM_EINT] = s->p * PRJ_CLIGHT * PRJ_CLIGHT /
                        (s->rho * (PRJ_IDEAL_GAMMA - 1.0));
                    W[PRJ_PRIM_YE] = 0.5;
                    for (g = 0; g < PRJ_NEGROUP; ++g) {
                        W[PRJ_PRIM_RAD_E(0, g)] = Elab / (RAD_SCALE * PRJ_NEGROUP);
                        W[PRJ_PRIM_RAD_F1(0, g)] = Flab / (RAD_SCALE * PRJ_NEGROUP);
                    }
                    store_cell(sim, block, i, j, k, W, U);
                }
    }
}

void prj_problem_rad_sphere(prj_sim *sim, prj_mpi *mpi)
{ (void)mpi; fill_sphere_like(sim, 4); }
void prj_problem_rad_doppler(prj_sim *sim, prj_mpi *mpi)
{ (void)mpi; fill_sphere_like(sim, 5); }
void prj_problem_rad_grav_redshift(prj_sim *sim, prj_mpi *mpi)
{ (void)mpi; fill_sphere_like(sim, 6); }
void prj_problem_rad_shock(prj_sim *sim, prj_mpi *mpi)
{ (void)mpi; fill_shock(sim); }
