/* Anninos & Fragile (2020) radiation benchmarks, sections 4.1--4.3. */
#include <math.h>
#include <string.h>

#include "prj.h"

#define RAD_ARAD 7.5657e-15
#define RAD_KB_MEV 8.617333262145e-11
#define RAD_EOS_SCALE 0.95655684e18
#define RAD_R 2.0e9
#define RAD_KS 2.5e-6
#define RAD_GNEWT 6.67430e-8
#ifndef RAD_SCALE
#define RAD_SCALE 1.0e25
#endif

static double rad_ideal_eint(double temp_k)
{
    return RAD_EOS_SCALE * (RAD_KB_MEV * temp_k) / (2.0 / 3.0);
}

static double rad_t3_eint(double rho, double temp_k)
{
    return RAD_ARAD * pow(temp_k, 4.0) / rho;
}

static double rad_planck_kernel(double x)
{
    if (x < 1.0e-5) return x * x;
    if (x > 80.0) return 0.0;
    return x * x * x / expm1(x);
}

static void rad_picket_planck_weights(double temp_k, double *weight)
{
    const double emin = 1.0e-14;
    const double emax = 1.0e-5;
    const int nquad = 128;
    double sum = 0.0;
    int g;
    for (g = 0; g < PRJ_NEGROUP; ++g) {
        double e0 = emin * pow(emax / emin, (double)g / PRJ_NEGROUP);
        double e1 = emin * pow(emax / emin, (double)(g + 1) / PRJ_NEGROUP);
        double x0 = e0 / (RAD_KB_MEV * temp_k);
        double x1 = e1 / (RAD_KB_MEV * temp_k);
        double h = (x1 - x0) / nquad;
        double integ = rad_planck_kernel(x0) + rad_planck_kernel(x1);
        int n;
        for (n = 1; n < nquad; ++n)
            integ += (n & 1 ? 4.0 : 2.0) * rad_planck_kernel(x0 + n * h);
        weight[g] = integ * h / 3.0;
        sum += weight[g];
    }
    for (g = 0; g < PRJ_NEGROUP; ++g) weight[g] /= sum;
}

static void rad_store(prj_block *block, int i, int j, int k, double *W, double *U)
{
    prj_block_store_prim_cell(block, 0, i, j, k, W);
    prj_block_store_prim_cell(block, 1, i, j, k, W);
    prj_block_store_cons_cell(block, i, j, k, U);
}

void prj_rad_test_fill_problem(prj_sim *sim, int which)
{
    int bidx;
    double picket_weight[PRJ_NEGROUP];
    if (which == 3) rad_picket_planck_weights(1.0, picket_weight);
    prj_mesh_init(&sim->mesh, sim->mesh.root_nx[0], sim->mesh.root_nx[1],
        sim->mesh.root_nx[2], 0, &sim->coord, 0);
    for (bidx = 0; bidx < sim->mesh.nblocks; ++bidx) {
        prj_block *block = &sim->mesh.blocks[bidx];
        int i, j, k;
        if (block->id < 0 || block->active != 1) continue;
        for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
            for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
                for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                    double x = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
                    double W[PRJ_NVAR_PRIM] = {0.0};
                    double U[PRJ_NVAR_CONS] = {0.0};
                    double rho = which == 2 ? 9.0e14 : 1.0;
                    double temp_k = which == 1 ? 1.0e4 : (which == 2 ?
                        5.0e6 / RAD_KB_MEV : 1.0);
                    int g;
                    W[PRJ_PRIM_RHO] = rho;
                    W[PRJ_PRIM_EINT] = which == 3 ? rad_t3_eint(rho, temp_k) :
                        rad_ideal_eint(temp_k);
                    W[PRJ_PRIM_YE] = 0.5;
                    for (g = 0; g < PRJ_NEGROUP; ++g) {
                        double E;
                        double F = 0.0;
                        if (which == 1) {
                            E = RAD_ARAD * pow(temp_k, 4.0) /
                                (RAD_SCALE * (double)PRJ_NEGROUP);
                        } else if (which == 2) {
                            double that = PRJ_CLIGHT * (200.0 * RAD_R / PRJ_CLIGHT);
                            double profile = pow(RAD_KS / that, 1.5) *
                                exp(-3.0 * RAD_KS * x * x / (4.0 * that));
                            E = rho * PRJ_CLIGHT * PRJ_CLIGHT * profile /
                                (RAD_SCALE * (double)PRJ_NEGROUP);
                            F = x / (2.0 * (200.0 * RAD_R / PRJ_CLIGHT)) * E;
                        } else {
                            E = RAD_ARAD * picket_weight[g] / RAD_SCALE;
                        }
                        W[PRJ_PRIM_RAD_E(0, g)] = E;
                        W[PRJ_PRIM_RAD_F1(0, g)] = F;
                    }
                    prj_eos_prim2cons(&sim->eos, W, U);
                    rad_store(block, i, j, k, W, U);
                }
            }
        }
    }
}

void prj_problem_user_boundary(const prj_mesh *mesh, const prj_block *block,
    double *W, int axis, int side, int i, int j, int k,
    const double position[3], double time_seconds)
{
    (void)mesh;
    (void)block;
    (void)position;
    (void)time_seconds;
#if PRJ_RAD_TEST_PROBLEM == 1 && PRJ_USE_RADIATION_M1
    if (axis == 0 && side == 0) {
        int g;
        double E = RAD_ARAD * 1.0e24 / (RAD_SCALE * (double)PRJ_NEGROUP);
        W[WIDX(PRJ_PRIM_EINT, i, j, k)] = rad_ideal_eint(1.0e6);
        for (g = 0; g < PRJ_NEGROUP; ++g) {
            W[WIDX(PRJ_PRIM_RAD_E(0, g), i, j, k)] = E;
            W[WIDX(PRJ_PRIM_RAD_F1(0, g), i, j, k)] = 0.999 * PRJ_CLIGHT * E;
            W[WIDX(PRJ_PRIM_RAD_F2(0, g), i, j, k)] = 0.0;
            W[WIDX(PRJ_PRIM_RAD_F3(0, g), i, j, k)] = 0.0;
        }
    }
#else
    (void)W; (void)axis; (void)side; (void)i; (void)j; (void)k;
#endif
}

void prj_problem_user_source(const prj_mesh *mesh, const prj_block *block,
    double *W_mhd, double *W_rad, double *mhd_rhs, double *rad_rhs)
{
    (void)W_mhd; (void)W_rad; (void)mhd_rhs;
#if PRJ_RAD_TEST_PROBLEM == 3 && PRJ_USE_RADIATION_M1
    if (mesh != 0 && block != 0 && rad_rhs != 0 &&
        mesh->time_seconds * PRJ_CLIGHT * 11.0 <= 10.0) {
        int i, j, k, g;
        double weight[PRJ_NEGROUP];
        double total = PRJ_CLIGHT * RAD_ARAD * 11.0 * 1.0e12 / RAD_SCALE;
        rad_picket_planck_weights(1.0e3, weight);
        for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
            double x = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
            if (x > 0.5 / 11.0) continue;
            for (j = 0; j < PRJ_BLOCK_SIZE; ++j) for (k = 0; k < PRJ_BLOCK_SIZE; ++k)
                for (g = 0; g < PRJ_NEGROUP; ++g)
                    rad_rhs[RADVIDX(PRJ_RAD_CONS_E(0, g), i, j, k)] +=
                        total * weight[g];
        }
    }
#else
    (void)mesh; (void)block; (void)rad_rhs;
#endif
}

void prj_problem_static_metric(const double position[3], double *lapse,
    double *gamma_rr, double *dphi_dr)
{
#if PRJ_USER_STATIC_METRIC
    const double rstar = 1.0e6;
    const double rho = 9.0e14;
    const double mass = (4.0 / 3.0) * M_PI * rho * rstar * rstar * rstar;
    double r = fabs(position[0]);
    double phi;
    double grad;

    if (r <= rstar) {
        phi = -RAD_GNEWT * mass * (3.0 * rstar * rstar - r * r) /
            (2.0 * rstar * rstar * rstar);
        grad = RAD_GNEWT * mass * r / (rstar * rstar * rstar);
    } else {
        phi = -RAD_GNEWT * mass / r;
        grad = RAD_GNEWT * mass / (r * r);
    }
    *lapse = sqrt(1.0 + 2.0 * phi / (PRJ_CLIGHT * PRJ_CLIGHT));
    *gamma_rr = 1.0 - 2.0 * phi / (PRJ_CLIGHT * PRJ_CLIGHT);
    *dphi_dr = grad;
#else
    (void)position;
    *lapse = 1.0;
    *gamma_rr = 1.0;
    *dphi_dr = 0.0;
#endif
}
