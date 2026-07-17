#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(PRJ_ENABLE_MPI)
#include <mpi.h>
#endif

#include "prj.h"

static void die(const char *msg)
{
    fprintf(stderr, "test_gr_prim_cons: %s\n", msg);
    exit(1);
}

static void assert_close(const char *name, double got, double expected, double rel, double abs_tol)
{
    double scale = fmax(1.0, fmax(fabs(got), fabs(expected)));
    double tol = fmax(abs_tol, rel * scale);

    if (!isfinite(got) || !isfinite(expected) || fabs(got - expected) > tol) {
        fprintf(stderr,
            "test_gr_prim_cons: %s got %.17e expected %.17e tol %.3e\n",
            name, got, expected, tol);
        exit(1);
    }
}

static double det3(const double g[3][3])
{
    return g[0][0] * (g[1][1] * g[2][2] - g[1][2] * g[2][1])
        - g[0][1] * (g[1][0] * g[2][2] - g[1][2] * g[2][0])
        + g[0][2] * (g[1][0] * g[2][1] - g[1][1] * g[2][0]);
}

static double geom_sqrt_det(const prj_eos_gr_geom *geom)
{
    return sqrt(det3(geom->gamma));
}

static void set_flat_geom(prj_eos_gr_geom *geom)
{
    memset(geom, 0, sizeof(*geom));
    geom->gamma[0][0] = 1.0;
    geom->gamma[1][1] = 1.0;
    geom->gamma[2][2] = 1.0;
}

static void set_curved_geom(prj_eos_gr_geom *geom)
{
    geom->gamma[0][0] = 1.30;
    geom->gamma[0][1] = 0.08;
    geom->gamma[0][2] = -0.04;
    geom->gamma[1][0] = 0.08;
    geom->gamma[1][1] = 0.95;
    geom->gamma[1][2] = 0.06;
    geom->gamma[2][0] = -0.04;
    geom->gamma[2][1] = 0.06;
    geom->gamma[2][2] = 1.12;
}

static void init_ideal_eos(prj_eos *eos)
{
    memset(eos, 0, sizeof(*eos));
    eos->kind = PRJ_EOS_KIND_IDEAL;
    prj_eos_init(eos, 0);
}

static void zero_prim(double *W)
{
    int v;

    for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
        W[v] = 0.0;
    }
}

static void set_state(double *W, double rho, double v1, double v2, double v3,
    double eint, double ye)
{
    zero_prim(W);
    W[PRJ_PRIM_RHO] = rho;
    W[PRJ_PRIM_V1] = v1;
    W[PRJ_PRIM_V2] = v2;
    W[PRJ_PRIM_V3] = v3;
    W[PRJ_PRIM_EINT] = eint;
    W[PRJ_PRIM_YE] = ye;
#if PRJ_MHD
    W[PRJ_PRIM_B1] = 1.1e12;
    W[PRJ_PRIM_B2] = -0.7e12;
    W[PRJ_PRIM_B3] = 0.4e12;
#endif
}

static void fill_radiation_state(double *W)
{
#if PRJ_USE_RADIATION_FSA
    int field;
    int group;
    int angle;

    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            for (angle = 0; angle < PRJ_NANGLE; ++angle) {
                W[PRJ_PRIM_RAD_I(field, group, angle)] =
                    10.0 + 100.0 * field + 3.0 * group + 0.01 * angle;
            }
        }
    }
#elif PRJ_USE_RADIATION_M1
    int field;
    int group;

    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            W[PRJ_PRIM_RAD_E(field, group)] = 10.0 + 100.0 * field + group;
            W[PRJ_PRIM_RAD_F1(field, group)] = 0.1 + field + 0.01 * group;
            W[PRJ_PRIM_RAD_F2(field, group)] = -0.2 - field - 0.02 * group;
            W[PRJ_PRIM_RAD_F3(field, group)] = 0.3 + field + 0.03 * group;
        }
    }
#else
    (void)W;
#endif
}

static void assert_radiation_densitized(const double *W, const double *U, double sqrt_det)
{
#if PRJ_USE_RADIATION_FSA
    int field;
    int group;
    int angle;

    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            for (angle = 0; angle < PRJ_NANGLE; ++angle) {
                int idx = PRJ_CONS_RAD_I(field, group, angle);

                assert_close("radiation I densitized", U[idx],
                    sqrt_det * W[PRJ_PRIM_RAD_I(field, group, angle)], 1.0e-13, 0.0);
            }
        }
    }
#elif PRJ_USE_RADIATION_M1
    int field;
    int group;

    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            assert_close("radiation E densitized", U[PRJ_CONS_RAD_E(field, group)],
                sqrt_det * W[PRJ_PRIM_RAD_E(field, group)], 1.0e-13, 0.0);
            assert_close("radiation F1 densitized", U[PRJ_CONS_RAD_F1(field, group)],
                sqrt_det * W[PRJ_PRIM_RAD_F1(field, group)], 1.0e-13, 0.0);
            assert_close("radiation F2 densitized", U[PRJ_CONS_RAD_F2(field, group)],
                sqrt_det * W[PRJ_PRIM_RAD_F2(field, group)], 1.0e-13, 0.0);
            assert_close("radiation F3 densitized", U[PRJ_CONS_RAD_F3(field, group)],
                sqrt_det * W[PRJ_PRIM_RAD_F3(field, group)], 1.0e-13, 0.0);
        }
    }
#else
    (void)W;
    (void)U;
    (void)sqrt_det;
#endif
}

static void assert_radiation_roundtrip(const double *W, const double *Wout)
{
#if PRJ_NRAD > 0
    int v;

    for (v = PRJ_NHYDRO; v < PRJ_NVAR_PRIM; ++v) {
        assert_close("radiation roundtrip", Wout[v], W[v], 1.0e-12, 0.0);
    }
#else
    (void)W;
    (void)Wout;
#endif
}

static void assert_mhd_densitized(const double *W, const double *U, double sqrt_det)
{
#if PRJ_MHD
    assert_close("B1 densitized", U[PRJ_CONS_B1],
        sqrt_det * W[PRJ_PRIM_B1], 1.0e-13, 0.0);
    assert_close("B2 densitized", U[PRJ_CONS_B2],
        sqrt_det * W[PRJ_PRIM_B2], 1.0e-13, 0.0);
    assert_close("B3 densitized", U[PRJ_CONS_B3],
        sqrt_det * W[PRJ_PRIM_B3], 1.0e-13, 0.0);
#else
    (void)W;
    (void)U;
    (void)sqrt_det;
#endif
}

static void assert_grmhd_state_helper(prj_eos *eos, const prj_eos_gr_geom *geom,
    const double *W, double sqrt_det)
{
#if PRJ_MHD
    prj_eos_grmhd_state state;
    int status = prj_eos_grmhd_state_from_prim(eos, geom, W, NAN, &state,
        PRJ_EOS_CTX_MAIN);

    if (status != PRJ_EOS_GR_OK) {
        die("GRMHD state helper failed");
    }
    assert_close("state sqrt_gamma", state.sqrt_gamma, sqrt_det, 1.0e-13, 0.0);
    assert_close("state B1", state.Bcon[0], W[PRJ_PRIM_B1], 1.0e-13, 0.0);
    assert_close("state B2", state.Bcon[1], W[PRJ_PRIM_B2], 1.0e-13, 0.0);
    assert_close("state B3", state.Bcon[2], W[PRJ_PRIM_B3], 1.0e-13, 0.0);
    assert_close("state Uloc B1", state.Uloc[PRJ_CONS_B1],
        W[PRJ_PRIM_B1], 1.0e-13, 0.0);
    assert_close("state Uloc B2", state.Uloc[PRJ_CONS_B2],
        W[PRJ_PRIM_B2], 1.0e-13, 0.0);
    assert_close("state Uloc B3", state.Uloc[PRJ_CONS_B3],
        W[PRJ_PRIM_B3], 1.0e-13, 0.0);
    if (!(state.Bsq > 0.0) || !(state.ptot > state.pressure) ||
        !isfinite(state.E) || !isfinite(state.S_cov[0])) {
        die("GRMHD state helper missing magnetic contribution");
    }
#else
    (void)eos;
    (void)geom;
    (void)W;
    (void)sqrt_det;
#endif
}

static void assert_hydro_roundtrip(const double *W, const double *Wout)
{
    assert_close("rho roundtrip", Wout[PRJ_PRIM_RHO], W[PRJ_PRIM_RHO], 1.0e-8, 0.0);
    assert_close("v1 roundtrip", Wout[PRJ_PRIM_V1], W[PRJ_PRIM_V1], 1.0e-6, 1.0e-3);
    assert_close("v2 roundtrip", Wout[PRJ_PRIM_V2], W[PRJ_PRIM_V2], 1.0e-6, 1.0e-3);
    assert_close("v3 roundtrip", Wout[PRJ_PRIM_V3], W[PRJ_PRIM_V3], 1.0e-6, 1.0e-3);
    assert_close("eint roundtrip", Wout[PRJ_PRIM_EINT], W[PRJ_PRIM_EINT], 1.0e-6, 0.0);
    assert_close("Ye roundtrip", Wout[PRJ_PRIM_YE], W[PRJ_PRIM_YE], 1.0e-12, 0.0);
#if PRJ_MHD
    assert_close("B1 roundtrip", Wout[PRJ_PRIM_B1], W[PRJ_PRIM_B1], 1.0e-12, 0.0);
    assert_close("B2 roundtrip", Wout[PRJ_PRIM_B2], W[PRJ_PRIM_B2], 1.0e-12, 0.0);
    assert_close("B3 roundtrip", Wout[PRJ_PRIM_B3], W[PRJ_PRIM_B3], 1.0e-12, 0.0);
#endif
}

#if PRJ_DYNAMIC_GR && PRJ_NRAD > 0
static void init_cell_wrapper_mesh(prj_mesh *mesh, prj_coord *coord)
{
    memset(mesh, 0, sizeof(*mesh));
    memset(coord, 0, sizeof(*coord));
    prj_z4c_init_params(&mesh->z4c_params);
    mesh->use_full_dynamic_gr = 1;
    coord->x1min = 0.0;
    coord->x1max = (double)PRJ_BLOCK_SIZE;
    coord->x2min = 0.0;
    coord->x2max = (double)PRJ_BLOCK_SIZE;
    coord->x3min = 0.0;
    coord->x3max = (double)PRJ_BLOCK_SIZE;
    if (prj_mesh_init(mesh, 1, 1, 1, 0, coord, 0) != 0) {
        die("cell-wrapper mesh init failed");
    }
}

static void set_diagonal_z4c_metric(prj_block *block,
    double chi, double gxx, double gyy, double gzz)
{
    double *z;
    int i;
    int j;
    int k;

    z = prj_block_z4c_stage(block, 0);
    if (z == 0) {
        die("missing z4c storage for cell-wrapper test");
    }
    for (i = -PRJ_NGHOST_Z4C; i < PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C; ++i) {
        for (j = -PRJ_NGHOST_Z4C; j < PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C; ++j) {
            for (k = -PRJ_NGHOST_Z4C; k < PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C; ++k) {
                z[Z4CIDX(PRJ_Z4C_CHI, i, j, k)] = chi;
                z[Z4CIDX(PRJ_Z4C_GXX, i, j, k)] = gxx;
                z[Z4CIDX(PRJ_Z4C_GXY, i, j, k)] = 0.0;
                z[Z4CIDX(PRJ_Z4C_GXZ, i, j, k)] = 0.0;
                z[Z4CIDX(PRJ_Z4C_GYY, i, j, k)] = gyy;
                z[Z4CIDX(PRJ_Z4C_GYZ, i, j, k)] = 0.0;
                z[Z4CIDX(PRJ_Z4C_GZZ, i, j, k)] = gzz;
                z[Z4CIDX(PRJ_Z4C_KHAT, i, j, k)] = 0.0;
                z[Z4CIDX(PRJ_Z4C_AXX, i, j, k)] = 0.0;
                z[Z4CIDX(PRJ_Z4C_AXY, i, j, k)] = 0.0;
                z[Z4CIDX(PRJ_Z4C_AXZ, i, j, k)] = 0.0;
                z[Z4CIDX(PRJ_Z4C_AYY, i, j, k)] = 0.0;
                z[Z4CIDX(PRJ_Z4C_AYZ, i, j, k)] = 0.0;
                z[Z4CIDX(PRJ_Z4C_AZZ, i, j, k)] = 0.0;
                z[Z4CIDX(PRJ_Z4C_GAMX, i, j, k)] = 0.0;
                z[Z4CIDX(PRJ_Z4C_GAMY, i, j, k)] = 0.0;
                z[Z4CIDX(PRJ_Z4C_GAMZ, i, j, k)] = 0.0;
                z[Z4CIDX(PRJ_Z4C_THETA, i, j, k)] = 0.0;
                z[Z4CIDX(PRJ_Z4C_ALPHA, i, j, k)] = 1.0;
                z[Z4CIDX(PRJ_Z4C_BETAX, i, j, k)] = 0.0;
                z[Z4CIDX(PRJ_Z4C_BETAY, i, j, k)] = 0.0;
                z[Z4CIDX(PRJ_Z4C_BETAZ, i, j, k)] = 0.0;
            }
        }
    }
}

static void test_cell_wrappers_preserve_radiation_slots(prj_eos *eos)
{
    prj_mesh mesh;
    prj_coord coord;
    prj_block *block;
    double W[PRJ_NVAR_PRIM];
    double Wout[PRJ_NVAR_PRIM];
    double U[PRJ_NVAR_CONS];
    const int i = 1;
    const int j = 2;
    const int k = 3;
    const double chi = 0.25;
    const double gxx = 2.0;
    const double gyy = 0.5;
    const double gzz = 1.5;
    double sqrt_det;

    init_cell_wrapper_mesh(&mesh, &coord);
    block = &mesh.blocks[0];
    set_diagonal_z4c_metric(block, chi, gxx, gyy, gzz);
    sqrt_det = sqrt(gxx * gyy * gzz / (chi * chi * chi));

    set_state(W, 1.4e8, 0.01 * PRJ_CLIGHT, -0.02 * PRJ_CLIGHT,
        0.015 * PRJ_CLIGHT, 3.0e18, 0.28);
    fill_radiation_state(W);
    prj_eos_cell_prim2cons(eos, &mesh, block, 0, i, j, k, W, U,
        PRJ_EOS_CTX_MAIN);
    assert_radiation_densitized(W, U, sqrt_det);
    prj_eos_cell_cons2prim(eos, &mesh, block, 0, i, j, k, U, Wout,
        PRJ_EOS_CTX_MAIN);
    assert_hydro_roundtrip(W, Wout);
    assert_radiation_roundtrip(W, Wout);
    prj_mesh_destroy(&mesh);
}
#endif

static void test_low_velocity_limit(prj_eos *eos)
{
    prj_eos_gr_geom geom;
    double W[PRJ_NVAR_PRIM];
    double U[PRJ_NVAR_CONS];
    double rho = 1.2e8;
    double v1 = 1.0e5;
    double v2 = -2.0e5;
    double v3 = 3.0e5;
    double eint = 2.0e15;
    double expected_tau;
    int status;

    set_flat_geom(&geom);
    set_state(W, rho, v1, v2, v3, eint, 0.21);
#if PRJ_MHD
    W[PRJ_PRIM_B1] = 2.0e10;
    W[PRJ_PRIM_B2] = -1.5e10;
    W[PRJ_PRIM_B3] = 0.5e10;
#endif
    status = prj_eos_gr_prim2cons(eos, &geom, W, U, PRJ_EOS_CTX_MAIN);
    if (status != PRJ_EOS_GR_OK) {
        die("low-velocity prim2cons failed");
    }
    expected_tau = rho * eint + 0.5 * rho * (v1 * v1 + v2 * v2 + v3 * v3);
#if PRJ_MHD
    expected_tau += 0.5 * (W[PRJ_PRIM_B1] * W[PRJ_PRIM_B1] +
        W[PRJ_PRIM_B2] * W[PRJ_PRIM_B2] + W[PRJ_PRIM_B3] * W[PRJ_PRIM_B3]);
#endif
    assert_close("D low-velocity", U[PRJ_CONS_RHO], rho, 1.0e-9, 0.0);
    assert_close("S1 low-velocity", U[PRJ_CONS_MOM1], rho * v1, 1.0e-5, 0.0);
    assert_close("S2 low-velocity", U[PRJ_CONS_MOM2], rho * v2, 1.0e-5, 0.0);
    assert_close("S3 low-velocity", U[PRJ_CONS_MOM3], rho * v3, 1.0e-5, 0.0);
    assert_close("tau low-velocity", U[PRJ_CONS_ETOT], expected_tau, 1.0e-5, 0.0);
}

static void test_roundtrip(prj_eos *eos, const char *name, const prj_eos_gr_geom *geom,
    double rho, double v1, double v2, double v3, double eint)
{
    double W[PRJ_NVAR_PRIM];
    double Wout[PRJ_NVAR_PRIM];
    double U[PRJ_NVAR_CONS];
    double sqrt_det = geom_sqrt_det(geom);
    int status;
    int v;

    set_state(W, rho, v1, v2, v3, eint, 0.32);
    fill_radiation_state(W);
    for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
        Wout[v] = -9.0e99;
    }
    status = prj_eos_gr_prim2cons(eos, geom, W, U, PRJ_EOS_CTX_MAIN);
    if (status != PRJ_EOS_GR_OK) {
        fprintf(stderr, "test_gr_prim_cons: %s prim2cons status %d\n", name, status);
        exit(1);
    }
    assert_mhd_densitized(W, U, sqrt_det);
    assert_grmhd_state_helper(eos, geom, W, sqrt_det);
    assert_radiation_densitized(W, U, sqrt_det);
    status = prj_eos_gr_cons2prim(eos, geom, U, Wout, PRJ_EOS_CTX_MAIN);
    if (status != PRJ_EOS_GR_OK) {
        fprintf(stderr, "test_gr_prim_cons: %s cons2prim status %d\n", name, status);
        exit(1);
    }
    assert_hydro_roundtrip(W, Wout);
    assert_radiation_roundtrip(W, Wout);
}

static void test_invalid_inputs(prj_eos *eos)
{
    prj_eos_gr_geom geom;
    prj_eos_gr_geom bad_geom;
    double W[PRJ_NVAR_PRIM];
    double U[PRJ_NVAR_CONS];
    double out[PRJ_NVAR_CONS];
    double Wout[PRJ_NVAR_PRIM];
    int status;
    int v;

    set_flat_geom(&geom);
    set_flat_geom(&bad_geom);
    bad_geom.gamma[1][1] = -1.0;
    set_state(W, 1.0, 1.0e7, 0.0, 0.0, 1.0e18, 0.2);
    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        out[v] = 123.0 + v;
    }
    status = prj_eos_gr_prim2cons(eos, &bad_geom, W, out, PRJ_EOS_CTX_MAIN);
    if (status != PRJ_EOS_GR_BAD_METRIC) {
        die("bad metric did not return BAD_METRIC");
    }
    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        assert_close("bad metric output unchanged", out[v], 123.0 + v, 0.0, 0.0);
    }
    W[PRJ_PRIM_V1] = 1.1 * PRJ_CLIGHT;
    status = prj_eos_gr_prim2cons(eos, &geom, W, out, PRJ_EOS_CTX_MAIN);
    if (status != PRJ_EOS_GR_BAD_STATE) {
        die("superluminal primitive did not return BAD_STATE");
    }
    set_state(W, 1.0, 1.0e7, 0.0, 0.0, 1.0e18, 0.2);
    status = prj_eos_gr_prim2cons(eos, &geom, W, U, PRJ_EOS_CTX_MAIN);
    if (status != PRJ_EOS_GR_OK) {
        die("valid state failed before invalid conserved test");
    }
    U[PRJ_CONS_RHO] = -U[PRJ_CONS_RHO];
    for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
        Wout[v] = -321.0 - v;
    }
    status = prj_eos_gr_cons2prim(eos, &geom, U, Wout, PRJ_EOS_CTX_MAIN);
    if (status != PRJ_EOS_GR_BAD_STATE) {
        die("negative conserved density did not return BAD_STATE");
    }
    for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
        assert_close("bad conserved output unchanged", Wout[v], -321.0 - v, 0.0, 0.0);
    }
}

int main(int argc, char **argv)
{
    prj_eos eos;
    prj_eos_gr_geom flat;
    prj_eos_gr_geom curved;

#if defined(PRJ_ENABLE_MPI)
    MPI_Init(&argc, &argv);
#else
    (void)argc;
    (void)argv;
#endif

    init_ideal_eos(&eos);
    set_flat_geom(&flat);
    set_curved_geom(&curved);
    test_low_velocity_limit(&eos);
    test_roundtrip(&eos, "flat", &flat, 1.0e8, 0.08 * PRJ_CLIGHT,
        -0.03 * PRJ_CLIGHT, 0.02 * PRJ_CLIGHT, 4.0e18);
    test_roundtrip(&eos, "curved", &curved, 2.0e8, 0.04 * PRJ_CLIGHT,
        -0.02 * PRJ_CLIGHT, 0.03 * PRJ_CLIGHT, 6.0e18);
#if PRJ_DYNAMIC_GR && PRJ_NRAD > 0
    test_cell_wrappers_preserve_radiation_slots(&eos);
#endif
    test_invalid_inputs(&eos);
    printf("test_gr_prim_cons: ok\n");
#if defined(PRJ_ENABLE_MPI)
    MPI_Finalize();
#endif
    return 0;
}
