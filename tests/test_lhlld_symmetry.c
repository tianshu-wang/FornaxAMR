#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "prj.h"

#if PRJ_MHD
static void die(const char *msg)
{
    fprintf(stderr, "test_lhlld_symmetry: %s\n", msg);
    exit(1);
}

static double test_pressure_from_state(const double *W, double gamma)
{
    return (gamma - 1.0) * W[PRJ_PRIM_RHO] * W[PRJ_PRIM_EINT];
}

static void test_set_state(double *W, double rho, double vx, double vy, double vz,
    double p, double ye, double bx, double by, double bz, double gamma)
{
    int v;

    for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
        W[v] = 0.0;
    }
    W[PRJ_PRIM_RHO] = rho;
    W[PRJ_PRIM_V1] = vx;
    W[PRJ_PRIM_V2] = vy;
    W[PRJ_PRIM_V3] = vz;
    W[PRJ_PRIM_EINT] = p / ((gamma - 1.0) * rho);
    W[PRJ_PRIM_YE] = ye;
    W[PRJ_PRIM_B1] = bx;
    W[PRJ_PRIM_B2] = by;
    W[PRJ_PRIM_B3] = bz;
}

static void test_mirror_state(const double *src, double *dst)
{
    int v;

    for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
        dst[v] = src[v];
    }
    dst[PRJ_PRIM_V1] = -src[PRJ_PRIM_V1];
    dst[PRJ_PRIM_B1] = -src[PRJ_PRIM_B1];
}

static void test_assert_close(const char *label, double a, double b)
{
    double scale = fmax(1.0, fmax(fabs(a), fabs(b)));
    double tol = 1.0e-10 * scale;

    if (!isfinite(a) || !isfinite(b) || fabs(a - b) > tol) {
        fprintf(stderr, "test_lhlld_symmetry: %s mismatch: %.17e vs %.17e (tol %.3e)\n",
            label, a, b, tol);
        exit(1);
    }
}

static void test_assert_flux_parity(const double *F, const double *Fm)
{
    test_assert_close("rho flux parity", Fm[PRJ_CONS_RHO], -F[PRJ_CONS_RHO]);
    test_assert_close("mom1 flux parity", Fm[PRJ_CONS_MOM1], F[PRJ_CONS_MOM1]);
    test_assert_close("mom2 flux parity", Fm[PRJ_CONS_MOM2], -F[PRJ_CONS_MOM2]);
    test_assert_close("mom3 flux parity", Fm[PRJ_CONS_MOM3], -F[PRJ_CONS_MOM3]);
    test_assert_close("etot flux parity", Fm[PRJ_CONS_ETOT], -F[PRJ_CONS_ETOT]);
    test_assert_close("ye flux parity", Fm[PRJ_CONS_YE], -F[PRJ_CONS_YE]);
    test_assert_close("B1 flux parity", Fm[PRJ_CONS_B1], F[PRJ_CONS_B1]);
    test_assert_close("B2 flux parity", Fm[PRJ_CONS_B2], -F[PRJ_CONS_B2]);
    test_assert_close("B3 flux parity", Fm[PRJ_CONS_B3], -F[PRJ_CONS_B3]);
}

static void test_assert_flux_close(const char *label, const double *F,
    const double *Fref)
{
    int v;

    for (v = 0; v < PRJ_NHYDRO; ++v) {
        char name[96];
        snprintf(name, sizeof(name), "%s flux[%d]", label, v);
        test_assert_close(name, F[v], Fref[v]);
    }
}

static void test_assert_face_parity(const double *v_face, const double *v_face_m,
    double bv1, double bv2, double bv1_m, double bv2_m)
{
    test_assert_close("vface1 parity", v_face_m[0], -v_face[0]);
    test_assert_close("vface2 parity", v_face_m[1], v_face[1]);
    test_assert_close("vface3 parity", v_face_m[2], v_face[2]);
    test_assert_close("Bv1 parity", bv1_m, -bv1);
    test_assert_close("Bv2 parity", bv2_m, -bv2);
}

typedef void (*test_mhd_riemann_solver)(const double *WL, const double *WR,
    double pL, double pR, double gL, double gR, const prj_eos *eos,
    double bn, double *flux, double v_face[3], double *Bv1, double *Bv2,
    double deltau, double deltav, double deltaw);

static void test_direct_case_solver(const char *name, const char *solver_name,
    const double *WL, const double *WR, double gamma, double bn,
    double deltau, double deltav, double deltaw, test_mhd_riemann_solver solver)
{
    double WML[PRJ_NVAR_PRIM];
    double WMR[PRJ_NVAR_PRIM];
    double F[PRJ_NVAR_CONS] = {0.0};
    double FM[PRJ_NVAR_CONS] = {0.0};
    double v_face[3] = {0.0, 0.0, 0.0};
    double v_face_m[3] = {0.0, 0.0, 0.0};
    double bv1 = 0.0;
    double bv2 = 0.0;
    double bv1_m = 0.0;
    double bv2_m = 0.0;
    double pL = test_pressure_from_state(WL, gamma);
    double pR = test_pressure_from_state(WR, gamma);

    test_mirror_state(WR, WML);
    test_mirror_state(WL, WMR);

    solver(WL, WR, pL, pR, gamma, gamma, 0, bn, F,
        v_face, &bv1, &bv2, deltau, deltav, deltaw);
    solver(WML, WMR, pR, pL, gamma, gamma, 0, -bn, FM,
        v_face_m, &bv1_m, &bv2_m, deltau, deltav, deltaw);

    test_assert_flux_parity(F, FM);
    test_assert_face_parity(v_face, v_face_m, bv1, bv2, bv1_m, bv2_m);
    printf("test_lhlld_symmetry: direct %s %s ok\n", solver_name, name);
}

static void test_direct_case(const char *name, const double *WL, const double *WR,
    double gamma, double bn, double deltau, double deltav, double deltaw)
{
    test_direct_case_solver(name, "LHLLD", WL, WR, gamma, bn,
        deltau, deltav, deltaw, prj_riemann_lhlld);
    test_direct_case_solver(name, "HLL", WL, WR, gamma, bn,
        deltau, deltav, deltaw, prj_riemann_hll);
}

static void test_direct_solver_symmetry(void)
{
    const double gamma = 5.0 / 3.0;
    double WL[PRJ_NVAR_PRIM];
    double WR[PRJ_NVAR_PRIM];

    test_set_state(WL, 1.1, 0.23, -0.17, 0.11, 0.95, 0.12,
        0.05, 0.07, -0.04, gamma);
    test_set_state(WR, 0.82, -0.12, 0.08, -0.05, 0.70, 0.19,
        0.05, -0.03, 0.06, gamma);
    test_direct_case("generic-oblique", WL, WR, gamma, 0.05, -0.35, -0.02, 0.03);

    test_set_state(WL, 1.0, 2.0e-3, 1.0e-3, -1.5e-3, 1.0e4, 0.10,
        2.0e-4, -1.5e-4, 1.0e-4, gamma);
    test_set_state(WR, 0.98, -1.0e-3, -8.0e-4, 1.2e-3, 0.9999e4, 0.11,
        2.0e-4, 1.0e-4, -2.0e-4, gamma);
    test_direct_case("low-mach", WL, WR, gamma, 2.0e-4, -3.0e-3, -1.0e-3, -5.0e-4);

    test_set_state(WL, 1.0, 5.0, 0.20, -0.10, 1.0, 0.15,
        0.10, 0.02, 0.01, gamma);
    test_set_state(WR, 0.9, 4.2, -0.10, 0.05, 0.8, 0.13,
        0.10, -0.01, 0.03, gamma);
    test_direct_case("supersonic-upwind", WL, WR, gamma, 0.10, -0.8, 0.02, 0.01);

    test_set_state(WL, 1.3, 0.08, -0.05, 0.04, 0.6, 0.16,
        0.0, 0.12, -0.08, gamma);
    test_set_state(WR, 0.7, -0.04, 0.03, -0.02, 0.5, 0.18,
        0.0, -0.06, 0.09, gamma);
    test_direct_case("zero-normal-field", WL, WR, gamma, 0.0, -0.12, -0.01, 0.02);
}

static void test_fill_eosvar(prj_block *block, int i, int j, int k,
    double p, double gamma)
{
    block->eosvar[EIDX(PRJ_EOSVAR_PRESSURE, i, j, k)] = p;
    block->eosvar[EIDX(PRJ_EOSVAR_TEMPERATURE, i, j, k)] = p;
    block->eosvar[EIDX(PRJ_EOSVAR_GAMMA, i, j, k)] = gamma;
}

static void test_fill_cell(prj_block *block, int i, int j, int k,
    const double *W, double p, double gamma)
{
    int v;

    for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
        block->W[WIDX(v, i, j, k)] = W[v];
    }
    test_fill_eosvar(block, i, j, k, p, gamma);
}

static double test_original_bface(int dir, int i)
{
    if (dir == X1DIR) {
        return 0.04 + 0.002 * (double)(i - PRJ_BLOCK_SIZE / 2);
    }
    if (dir == X2DIR) {
        return -0.03 + 0.001 * (double)i;
    }
    return 0.02 - 0.0015 * (double)i;
}

static void test_fill_block_pair(prj_block *block, prj_block *mirror)
{
    const double gamma = 5.0 / 3.0;
    const int mid = PRJ_BLOCK_SIZE / 2;
    int i;
    int j;
    int k;
    int dir;

    for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
        for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
            for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                double W[PRJ_NVAR_PRIM];
                double WM[PRJ_NVAR_PRIM];
                double p;
                double pm;
                int io = PRJ_BLOCK_SIZE - 1 - i;

                if (i < mid) {
                    test_set_state(W, 1.15, 0.14, -0.06, 0.04, 0.90, 0.17,
                        0.0, 0.05, -0.03, gamma);
                } else {
                    test_set_state(W, 0.85, -0.09, 0.04, -0.02, 0.72, 0.20,
                        0.0, -0.04, 0.06, gamma);
                }
                p = test_pressure_from_state(W, gamma);
                test_fill_cell(block, i, j, k, W, p, gamma);

                if (io < mid) {
                    test_set_state(WM, 1.15, 0.14, -0.06, 0.04, 0.90, 0.17,
                        0.0, 0.05, -0.03, gamma);
                } else {
                    test_set_state(WM, 0.85, -0.09, 0.04, -0.02, 0.72, 0.20,
                        0.0, -0.04, 0.06, gamma);
                }
                test_mirror_state(WM, WM);
                pm = test_pressure_from_state(WM, gamma);
                test_fill_cell(mirror, i, j, k, WM, pm, gamma);
            }
        }
    }

    for (dir = 0; dir < 3; ++dir) {
        int imax = PRJ_BLOCK_SIZE + PRJ_NGHOST - 1;
        int jmax = PRJ_BLOCK_SIZE + PRJ_NGHOST - 1;
        int kmax = PRJ_BLOCK_SIZE + PRJ_NGHOST - 1;

        if (dir == X1DIR) {
            imax = PRJ_BLOCK_SIZE + PRJ_NGHOST;
        } else if (dir == X2DIR) {
            jmax = PRJ_BLOCK_SIZE + PRJ_NGHOST;
        } else {
            kmax = PRJ_BLOCK_SIZE + PRJ_NGHOST;
        }
        for (i = -PRJ_NGHOST; i <= imax; ++i) {
            for (j = -PRJ_NGHOST; j <= jmax; ++j) {
                for (k = -PRJ_NGHOST; k <= kmax; ++k) {
                    int io = (dir == X1DIR) ? PRJ_BLOCK_SIZE - i : PRJ_BLOCK_SIZE - 1 - i;
                    double bf = test_original_bface(dir, i);
                    double bfm = test_original_bface(dir, io);

                    if (dir == X1DIR) {
                        block->Bf[dir][FACE_IDX(dir, i, j, k)] = bf;
                        mirror->Bf[dir][FACE_IDX(dir, i, j, k)] = -bfm;
                    } else {
                        block->Bf[dir][FACE_IDX(dir, i, j, k)] = bf;
                        mirror->Bf[dir][FACE_IDX(dir, i, j, k)] = bfm;
                    }
                }
            }
        }
    }
}

static void test_flux_update_symmetry(void)
{
    prj_block block;
    prj_block mirror;
    prj_eos eos;
    double *flux[3];
    double *flux_m[3];
    int iface = PRJ_BLOCK_SIZE / 2;
    int j = 2;
    int k = 3;

    memset(&block, 0, sizeof(block));
    memset(&mirror, 0, sizeof(mirror));
    memset(&eos, 0, sizeof(eos));

    if (prj_block_alloc_data(&block) != 0 || prj_block_alloc_data(&mirror) != 0) {
        die("block allocation failed");
    }
    test_fill_block_pair(&block, &mirror);
    for (int dir = 0; dir < 3; ++dir) {
        flux[dir] = block.flux[dir];
        flux_m[dir] = mirror.flux[dir];
    }

    prj_flux_update(&eos, 0, &block, block.W, block.eosvar, flux, 0);
    prj_flux_update(&eos, 0, &mirror, mirror.W, mirror.eosvar, flux_m, 0);

    {
        double F[PRJ_NHYDRO] = {0.0};
        double FM[PRJ_NHYDRO] = {0.0};
        int v;

        for (v = 0; v < PRJ_NHYDRO; ++v) {
            F[v] = block.flux[X1DIR][VIDX(v, iface, j, k)];
            FM[v] = mirror.flux[X1DIR][VIDX(v, iface, j, k)];
        }
        test_assert_flux_parity(F, FM);
    }
    test_assert_close("integrated vface1 parity",
        mirror.v_riemann[X1DIR][VRIDX(0, iface, j, k)],
        -block.v_riemann[X1DIR][VRIDX(0, iface, j, k)]);
    test_assert_close("integrated vface2 parity",
        mirror.v_riemann[X1DIR][VRIDX(1, iface, j, k)],
        block.v_riemann[X1DIR][VRIDX(1, iface, j, k)]);
    test_assert_close("integrated vface3 parity",
        mirror.v_riemann[X1DIR][VRIDX(2, iface, j, k)],
        block.v_riemann[X1DIR][VRIDX(2, iface, j, k)]);
    test_assert_close("integrated Bv1 parity",
        mirror.Bv1[X1DIR][IDX(iface, j, k)], -block.Bv1[X1DIR][IDX(iface, j, k)]);
    test_assert_close("integrated Bv2 parity",
        mirror.Bv2[X1DIR][IDX(iface, j, k)], -block.Bv2[X1DIR][IDX(iface, j, k)]);

    prj_block_free_data(&mirror);
    prj_block_free_data(&block);
    printf("test_lhlld_symmetry: flux update ok\n");
}

static void test_dense_face_hll_override(void)
{
    prj_block block;
    prj_eos eos;
    double *flux[3];
    const double gamma = 5.0 / 3.0;
    const int iface = PRJ_BLOCK_SIZE / 2;
    const int jface = 3;
    const int kface = 4;
    double WL[PRJ_NVAR_PRIM];
    double WR[PRJ_NVAR_PRIM];
    double F[PRJ_NHYDRO] = {0.0};
    double Fhll[PRJ_NHYDRO] = {0.0};
    double v_face[3] = {0.0, 0.0, 0.0};
    double v_hll[3] = {0.0, 0.0, 0.0};
    double bv1_hll = 0.0;
    double bv2_hll = 0.0;
    double bn;
    double pL;
    double pR;
    int i;
    int j;
    int k;
    int dir;
    int v;

    memset(&block, 0, sizeof(block));
    memset(&eos, 0, sizeof(eos));

    if (prj_block_alloc_data(&block) != 0) {
        die("dense override block allocation failed");
    }

    test_set_state(WL, 1.1 * PRJ_HLL_RIEMANN_DENSITY_THRESHOLD,
        0.14, -0.06, 0.04, 1.0e27, 0.17, 0.0, 0.05, -0.03, gamma);
    test_set_state(WR, 0.9 * PRJ_HLL_RIEMANN_DENSITY_THRESHOLD,
        -0.09, 0.04, -0.02, 0.8e27, 0.20, 0.0, -0.04, 0.06, gamma);
    pL = test_pressure_from_state(WL, gamma);
    pR = test_pressure_from_state(WR, gamma);

    for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
        for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
            for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                const double *W = i < iface ? WL : WR;
                double p = i < iface ? pL : pR;
                test_fill_cell(&block, i, j, k, W, p, gamma);
            }
        }
    }

    for (dir = 0; dir < 3; ++dir) {
        int imax = PRJ_BLOCK_SIZE + PRJ_NGHOST - 1;
        int jmax = PRJ_BLOCK_SIZE + PRJ_NGHOST - 1;
        int kmax = PRJ_BLOCK_SIZE + PRJ_NGHOST - 1;

        if (dir == X1DIR) {
            imax = PRJ_BLOCK_SIZE + PRJ_NGHOST;
        } else if (dir == X2DIR) {
            jmax = PRJ_BLOCK_SIZE + PRJ_NGHOST;
        } else {
            kmax = PRJ_BLOCK_SIZE + PRJ_NGHOST;
        }
        for (i = -PRJ_NGHOST; i <= imax; ++i) {
            for (j = -PRJ_NGHOST; j <= jmax; ++j) {
                for (k = -PRJ_NGHOST; k <= kmax; ++k) {
                    block.Bf[dir][FACE_IDX(dir, i, j, k)] =
                        dir == X1DIR ? 0.04 : (dir == X2DIR ? -0.03 : 0.02);
                }
            }
        }
    }

    for (dir = 0; dir < 3; ++dir) {
        flux[dir] = block.flux[dir];
    }
    prj_flux_update(&eos, 0, &block, block.W, block.eosvar, flux, 0);

    bn = block.Bf[X1DIR][FACE_IDX(X1DIR, iface, jface, kface)];
    WL[PRJ_PRIM_B1] = bn;
    WR[PRJ_PRIM_B1] = bn;
    prj_riemann_hll(WL, WR, pL, pR, gamma, gamma, &eos, bn, Fhll,
        v_hll, &bv1_hll, &bv2_hll, 0.0, 0.0, 0.0);

    for (v = 0; v < PRJ_NHYDRO; ++v) {
        F[v] = block.flux[X1DIR][VIDX(v, iface, jface, kface)];
    }
    v_face[0] = block.v_riemann[X1DIR][VRIDX(0, iface, jface, kface)];
    v_face[1] = block.v_riemann[X1DIR][VRIDX(1, iface, jface, kface)];
    v_face[2] = block.v_riemann[X1DIR][VRIDX(2, iface, jface, kface)];

    test_assert_flux_close("dense HLL override", F, Fhll);
    test_assert_close("dense HLL vface1", v_face[0], v_hll[0]);
    test_assert_close("dense HLL vface2", v_face[1], v_hll[1]);
    test_assert_close("dense HLL vface3", v_face[2], v_hll[2]);
    test_assert_close("dense HLL Bv1", block.Bv1[X1DIR][IDX(iface, jface, kface)],
        bv1_hll);
    test_assert_close("dense HLL Bv2", block.Bv2[X1DIR][IDX(iface, jface, kface)],
        bv2_hll);

    prj_block_free_data(&block);
    printf("test_lhlld_symmetry: dense HLL override ok\n");
}
#endif

int main(int argc, char *argv[])
{
    (void)argc;
    (void)argv;

#if PRJ_MHD
    test_direct_solver_symmetry();
    test_flux_update_symmetry();
    test_dense_face_hll_override();
#else
    fprintf(stderr, "test_lhlld_symmetry: built without MHD (PRJ_MHD=0)\n");
#endif

    return 0;
}
