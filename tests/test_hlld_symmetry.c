#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(PRJ_ENABLE_MPI)
#include <mpi.h>
#endif

#include "prj.h"

#if PRJ_MHD
static void die(const char *msg)
{
    fprintf(stderr, "test_hlld_symmetry: %s\n", msg);
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

    for (v = 0; v < PRJ_NVAR_MHD_PRIM; ++v) {
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

    for (v = 0; v < PRJ_NVAR_MHD_PRIM; ++v) {
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
        fprintf(stderr, "test_hlld_symmetry: %s mismatch: %.17e vs %.17e (tol %.3e)\n",
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

static void test_assert_face_parity(const double *v_face, const double *v_face_m,
    double bv1, double bv2, double bv1_m, double bv2_m)
{
    test_assert_close("vface1 parity", v_face_m[0], -v_face[0]);
    test_assert_close("vface2 parity", v_face_m[1], v_face[1]);
    test_assert_close("vface3 parity", v_face_m[2], v_face[2]);
    test_assert_close("Bv1 parity", bv1_m, -bv1);
    test_assert_close("Bv2 parity", bv2_m, -bv2);
}

static void test_direct_case(const char *name, const double *WL, const double *WR,
    double gamma, double bn, double deltau, double deltav, double deltaw)
{
    double WML[PRJ_NVAR_MHD_PRIM];
    double WMR[PRJ_NVAR_MHD_PRIM];
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

    prj_riemann_hlld(WL, WR, pL, pR, gamma, gamma, 0, bn, F,
        v_face, &bv1, &bv2, deltau, deltav, deltaw);
    prj_riemann_hlld(WML, WMR, pR, pL, gamma, gamma, 0, -bn, FM,
        v_face_m, &bv1_m, &bv2_m, deltau, deltav, deltaw);

    test_assert_flux_parity(F, FM);
    test_assert_face_parity(v_face, v_face_m, bv1, bv2, bv1_m, bv2_m);
    printf("test_hlld_symmetry: direct %s ok\n", name);
}

static void test_direct_solver_symmetry(void)
{
    const double gamma = 5.0 / 3.0;
    double WL[PRJ_NVAR_MHD_PRIM];
    double WR[PRJ_NVAR_MHD_PRIM];

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

#if PRJ_DYNAMIC_GR
static double test_det3(const double g[3][3])
{
    return g[0][0] * (g[1][1] * g[2][2] - g[1][2] * g[2][1]) -
        g[0][1] * (g[1][0] * g[2][2] - g[1][2] * g[2][0]) +
        g[0][2] * (g[1][0] * g[2][1] - g[1][1] * g[2][0]);
}

static void test_mirror_gr_geom(const double gamma_metric[3][3],
    const double beta[3], double gamma_m[3][3], double beta_m[3])
{
    int a;
    int b;

    for (a = 0; a < 3; ++a) {
        double sa = a == 0 ? -1.0 : 1.0;

        beta_m[a] = sa * beta[a];
        for (b = 0; b < 3; ++b) {
            double sb = b == 0 ? -1.0 : 1.0;

            gamma_m[a][b] = sa * sb * gamma_metric[a][b];
        }
    }
}

static void test_gr_direct_case(const char *name,
    const double *WL, const double *WR, double gamma_gas,
    const double gamma_metric[3][3], double alpha, const double beta[3],
    double bn_tilde, double deltau, double deltav, double deltaw)
{
    double gamma_m[3][3];
    double beta_m[3];
    double WML[PRJ_NVAR_MHD_PRIM];
    double WMR[PRJ_NVAR_MHD_PRIM];
    double F[PRJ_NVAR_CONS] = {0.0};
    double FM[PRJ_NVAR_CONS] = {0.0};
    double v_face[3] = {0.0, 0.0, 0.0};
    double v_face_m[3] = {0.0, 0.0, 0.0};
    double bv1 = 0.0;
    double bv2 = 0.0;
    double bv1_m = 0.0;
    double bv2_m = 0.0;
    double pL = test_pressure_from_state(WL, gamma_gas);
    double pR = test_pressure_from_state(WR, gamma_gas);
    double detg = test_det3(gamma_metric);
    double sqrt_gamma;
    prj_eos eos;

    if (detg <= 0.0 || !isfinite(detg)) {
        die("bad GR direct-test metric");
    }
    sqrt_gamma = sqrt(detg);
    memset(&eos, 0, sizeof(eos));
    eos.kind = PRJ_EOS_KIND_IDEAL;
    prj_eos_init(&eos, 0);
    test_mirror_state(WR, WML);
    test_mirror_state(WL, WMR);
    test_mirror_gr_geom(gamma_metric, beta, gamma_m, beta_m);

    prj_riemann_gr_hlld(WL, WR, pL, pR, gamma_gas, gamma_gas,
        &eos, gamma_metric, sqrt_gamma, alpha, beta, bn_tilde, F,
        v_face, &bv1, &bv2, deltau, deltav, deltaw);
    prj_riemann_gr_hlld(WML, WMR, pR, pL, gamma_gas, gamma_gas,
        &eos, gamma_m, sqrt_gamma, alpha, beta_m, -bn_tilde, FM,
        v_face_m, &bv1_m, &bv2_m, deltau, deltav, deltaw);

    test_assert_flux_parity(F, FM);
    test_assert_face_parity(v_face, v_face_m, bv1, bv2, bv1_m, bv2_m);
    printf("test_hlld_symmetry: GR direct HLLD %s ok\n", name);
}

static void test_gr_direct_solver_symmetry(void)
{
    const double gamma_gas = 5.0 / 3.0;
    const double flat_metric[3][3] = {
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0}
    };
    const double curved_metric[3][3] = {
        {1.36, 0.08, -0.04},
        {0.08, 1.18, 0.05},
        {-0.04, 0.05, 0.92}
    };
    const double beta0[3] = {0.0, 0.0, 0.0};
    const double beta_shift[3] = {0.04, -0.015, 0.02};
    double WL[PRJ_NVAR_MHD_PRIM];
    double WR[PRJ_NVAR_MHD_PRIM];

    test_set_state(WL, 1.1, 2.0e7, -1.0e7, 0.6e7, 2.5e20, 0.12,
        3.0e9, 2.0e9, -1.0e9, gamma_gas);
    test_set_state(WR, 0.9, -1.5e7, 0.8e7, -0.4e7, 2.0e20, 0.19,
        3.0e9, -1.5e9, 0.7e9, gamma_gas);
    test_gr_direct_case("flat-metric", WL, WR, gamma_gas,
        flat_metric, 1.0, beta0, 3.0e9,
        -3.5e7, -1.0e7, 0.5e7);
    test_gr_direct_case("curved-shifted", WL, WR, gamma_gas,
        curved_metric, 0.82, beta_shift, 3.0e9,
        -3.5e7, -1.0e7, 0.5e7);
    test_gr_direct_case("small-normal-field-fallback", WL, WR, gamma_gas,
        flat_metric, 1.0, beta0, 0.0,
        -3.5e7, -1.0e7, 0.5e7);
}
#endif

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

    for (v = 0; v < PRJ_NVAR_MHD_PRIM; ++v) {
        block->W_mhd[WIDX(v, i, j, k)] = W[v];
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
                double W[PRJ_NVAR_MHD_PRIM];
                double WM[PRJ_NVAR_MHD_PRIM];
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

    prj_flux_update(&eos, 0, 0, &block, block.W_mhd, block.eosvar, flux, 0);
    prj_flux_update(&eos, 0, 0, &mirror, mirror.W_mhd, mirror.eosvar, flux_m, 0);

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
    test_assert_close("integrated Bv1 parity",
        mirror.Bv1[X1DIR][IDX(iface, j, k)], -block.Bv1[X1DIR][IDX(iface, j, k)]);
    test_assert_close("integrated Bv2 parity",
        mirror.Bv2[X1DIR][IDX(iface, j, k)], -block.Bv2[X1DIR][IDX(iface, j, k)]);

    prj_block_free_data(&mirror);
    prj_block_free_data(&block);
    printf("test_hlld_symmetry: flux update ok\n");
}

#if PRJ_DYNAMIC_GR
static void test_gr_flux_update_symmetry(void)
{
    prj_mesh mesh;
    prj_mesh mirror_mesh;
    prj_coord coord;
    prj_coord mirror_coord;
    prj_block *block;
    prj_block *mirror;
    prj_eos eos;
    double *flux[3];
    double *flux_m[3];
    int iface = PRJ_BLOCK_SIZE / 2;
    int j = 2;
    int k = 3;
    int dir;

    memset(&mesh, 0, sizeof(mesh));
    memset(&mirror_mesh, 0, sizeof(mirror_mesh));
    memset(&coord, 0, sizeof(coord));
    memset(&mirror_coord, 0, sizeof(mirror_coord));
    memset(&eos, 0, sizeof(eos));
    coord.x1max = (double)PRJ_BLOCK_SIZE;
    coord.x2max = (double)PRJ_BLOCK_SIZE;
    coord.x3max = (double)PRJ_BLOCK_SIZE;
    mirror_coord = coord;
    prj_z4c_init_params(&mesh.z4c_params);
    prj_z4c_init_params(&mirror_mesh.z4c_params);
    mesh.use_full_dynamic_gr = 1;
    mirror_mesh.use_full_dynamic_gr = 1;
    if (prj_mesh_init(&mesh, 1, 1, 1, 0, &coord, 0) != 0 ||
        prj_mesh_init(&mirror_mesh, 1, 1, 1, 0, &mirror_coord, 0) != 0) {
        die("GR mesh allocation failed");
    }
    prj_z4c_init_mesh_flat(&mesh, 0);
    prj_z4c_init_mesh_flat(&mirror_mesh, 0);
    block = &mesh.blocks[0];
    mirror = &mirror_mesh.blocks[0];
    test_fill_block_pair(block, mirror);
    for (dir = 0; dir < 3; ++dir) {
        flux[dir] = block->flux[dir];
        flux_m[dir] = mirror->flux[dir];
    }

    eos.kind = PRJ_EOS_KIND_IDEAL;
    prj_eos_init(&eos, 0);
    prj_flux_update(&eos, 0, &mesh, block, block->W_mhd, block->eosvar, flux, 0);
    prj_flux_update(&eos, 0, &mirror_mesh, mirror, mirror->W_mhd,
        mirror->eosvar, flux_m, 0);

    {
        double F[PRJ_NHYDRO] = {0.0};
        double FM[PRJ_NHYDRO] = {0.0};
        int v;

        for (v = 0; v < PRJ_NHYDRO; ++v) {
            F[v] = block->flux[X1DIR][VIDX(v, iface, j, k)];
            FM[v] = mirror->flux[X1DIR][VIDX(v, iface, j, k)];
        }
        test_assert_flux_parity(F, FM);
    }
    test_assert_close("GR integrated Bv1 parity",
        mirror->Bv1[X1DIR][IDX(iface, j, k)], -block->Bv1[X1DIR][IDX(iface, j, k)]);
    test_assert_close("GR integrated Bv2 parity",
        mirror->Bv2[X1DIR][IDX(iface, j, k)], -block->Bv2[X1DIR][IDX(iface, j, k)]);

    prj_mesh_destroy(&mirror_mesh);
    prj_mesh_destroy(&mesh);
    printf("test_hlld_symmetry: GR flux update ok\n");
}
#endif
#endif

int main(int argc, char *argv[])
{
#if defined(PRJ_ENABLE_MPI)
    MPI_Init(&argc, &argv);
#else
    (void)argc;
    (void)argv;
#endif

#if PRJ_MHD
    test_direct_solver_symmetry();
#if PRJ_DYNAMIC_GR
    test_gr_direct_solver_symmetry();
#endif
    test_flux_update_symmetry();
#if PRJ_DYNAMIC_GR
    test_gr_flux_update_symmetry();
#endif
#else
    fprintf(stderr, "test_hlld_symmetry: built without MHD (PRJ_MHD=0)\n");
#endif

#if defined(PRJ_ENABLE_MPI)
    MPI_Finalize();
#endif
    return 0;
}
