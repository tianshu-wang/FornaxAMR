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
    fprintf(stderr, "test_z4c: %s\n", msg);
    exit(1);
}

static void check_enum_names(void)
{
    static const char *const expected[PRJ_NZ4C] = {
        "z4c_chi",
        "z4c_gxx", "z4c_gxy", "z4c_gxz", "z4c_gyy", "z4c_gyz", "z4c_gzz",
        "z4c_Khat",
        "z4c_Axx", "z4c_Axy", "z4c_Axz", "z4c_Ayy", "z4c_Ayz", "z4c_Azz",
        "z4c_Gamx", "z4c_Gamy", "z4c_Gamz",
        "z4c_Theta",
        "z4c_alpha",
        "z4c_betax", "z4c_betay", "z4c_betaz"
    };
    static const char *const expected_tmunu[PRJ_NTMUNU] = {
        "z4c_tmunu_Sxx", "z4c_tmunu_Sxy", "z4c_tmunu_Sxz",
        "z4c_tmunu_Syy", "z4c_tmunu_Syz", "z4c_tmunu_Szz",
        "z4c_tmunu_E", "z4c_tmunu_Sx", "z4c_tmunu_Sy", "z4c_tmunu_Sz"
    };
    int v;

    if (PRJ_NZ4C != 22) {
        die("unexpected PRJ_NZ4C");
    }
    for (v = 0; v < PRJ_NZ4C; ++v) {
        if (strcmp(prj_z4c_var_name(v), expected[v]) != 0) {
            die("unexpected Z4c variable name/order");
        }
    }
    if (PRJ_NTMUNU != 10) {
        die("unexpected PRJ_NTMUNU");
    }
    for (v = 0; v < PRJ_NTMUNU; ++v) {
        if (strcmp(prj_z4c_tmunu_var_name(v), expected_tmunu[v]) != 0) {
            die("unexpected Tmunu variable name/order");
        }
    }
}

#if PRJ_DYNAMIC_GR
static void assert_close(const char *name, double got, double expected, double tol)
{
    if (!isfinite(got) || fabs(got - expected) > tol) {
        fprintf(stderr, "test_z4c: %s got %.17e expected %.17e tol %.3e\n",
            name, got, expected, tol);
        exit(1);
    }
}

static void assert_close_rel(const char *name, double got, double expected, double rel)
{
    double scale = expected != 0.0 ? fabs(expected) : 1.0;
    assert_close(name, got, expected, rel * scale);
}

static double geo_factor(void)
{
    double c2 = PRJ_CLIGHT * PRJ_CLIGHT;
    return PRJ_GNEWT / (c2 * c2);
}

static void init_one_block_mesh(prj_mesh *mesh, prj_coord *coord)
{
    memset(mesh, 0, sizeof(*mesh));
    prj_z4c_init_params(&mesh->z4c_params);
    mesh->use_dynamic_gr = 1;
    coord->x1min = 0.0;
    coord->x1max = (double)PRJ_BLOCK_SIZE * 3.0e10;
    coord->x2min = 0.0;
    coord->x2max = (double)PRJ_BLOCK_SIZE * 3.0e10;
    coord->x3min = 0.0;
    coord->x3max = (double)PRJ_BLOCK_SIZE * 3.0e10;
    if (prj_mesh_init(mesh, 1, 1, 1, 0, coord, 0) != 0) {
        die("mesh init failed");
    }
}

static void check_flat_state(void)
{
    prj_mesh mesh;
    prj_coord coord;
    prj_block *block;
    double *z;
    double *rhs;
    double dt;
    int i, j, k, v;

    init_one_block_mesh(&mesh, &coord);
    block = &mesh.blocks[0];
    if (block->z4c == 0 || block->z4c_rhs == 0 || block->z4c_tmunu == 0) {
        die("Z4c block storage missing");
    }

    dt = prj_z4c_calc_dt_seconds(&mesh, 0, 0.5);
    assert_close("dt_z4c", dt, 0.5, 1.0e-3);

    prj_z4c_init_mesh_flat(&mesh, 0);
    z = prj_block_z4c_stage(block, 0);
    for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
        for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
            for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                assert_close("chi", z[Z4CIDX(PRJ_Z4C_CHI, i, j, k)], 1.0, 0.0);
                assert_close("gxx", z[Z4CIDX(PRJ_Z4C_GXX, i, j, k)], 1.0, 0.0);
                assert_close("gyy", z[Z4CIDX(PRJ_Z4C_GYY, i, j, k)], 1.0, 0.0);
                assert_close("gzz", z[Z4CIDX(PRJ_Z4C_GZZ, i, j, k)], 1.0, 0.0);
                assert_close("alpha", z[Z4CIDX(PRJ_Z4C_ALPHA, i, j, k)], 1.0, 0.0);
                for (v = 0; v < PRJ_NZ4C; ++v) {
                    if (v != PRJ_Z4C_CHI && v != PRJ_Z4C_GXX &&
                        v != PRJ_Z4C_GYY && v != PRJ_Z4C_GZZ &&
                        v != PRJ_Z4C_ALPHA) {
                        assert_close("flat zero", z[Z4CIDX(v, i, j, k)], 0.0, 0.0);
                    }
                }
            }
        }
    }

    prj_z4c_compute_rhs(&mesh, 0, 0, 0, 0.0);
    rhs = prj_block_z4c_rhs_stage(block, 0);
    for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
        for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
            for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                for (v = 0; v < PRJ_NZ4C; ++v) {
                    assert_close("flat rhs", rhs[Z4CIDX(v, i, j, k)], 0.0, 1.0e-50);
                }
            }
        }
    }
    prj_mesh_destroy(&mesh);
}

static void check_tmunu_hydro_projection(void)
{
    prj_mesh mesh;
    prj_coord coord;
    prj_block *block;
    double *W;
    double *tmunu;
    const int i = 1, j = 2, k = 3;
    const double rho = 2.0;
    const double v1 = 3.0;
    const double v2 = 4.0;
    const double v3 = 5.0;
    const double eint = 7.0;
    const double pressure = 11.0;
    double fac = geo_factor();
    double B2 = 0.0;

    init_one_block_mesh(&mesh, &coord);
    block = &mesh.blocks[0];
    prj_z4c_init_mesh_flat(&mesh, 0);
    prj_fill(block->W_mhd, (size_t)PRJ_NVAR_MHD_PRIM *
        (size_t)PRJ_BLOCK_NSTAGES * (size_t)PRJ_BLOCK_NCELLS, 0.0);
    prj_fill(block->eosvar, (size_t)PRJ_NVAR_EOSVAR *
        (size_t)PRJ_BLOCK_NCELLS, 0.0);
    W = prj_block_mhd_stage(block, 0);
    W[WIDX(PRJ_PRIM_RHO, i, j, k)] = rho;
    W[WIDX(PRJ_PRIM_V1, i, j, k)] = v1;
    W[WIDX(PRJ_PRIM_V2, i, j, k)] = v2;
    W[WIDX(PRJ_PRIM_V3, i, j, k)] = v3;
    W[WIDX(PRJ_PRIM_EINT, i, j, k)] = eint;
    block->eosvar[EIDX(PRJ_EOSVAR_PRESSURE, i, j, k)] = pressure;
#if PRJ_MHD
    W[WIDX(PRJ_PRIM_B1, i, j, k)] = 13.0;
    W[WIDX(PRJ_PRIM_B2, i, j, k)] = 17.0;
    W[WIDX(PRJ_PRIM_B3, i, j, k)] = 19.0;
    B2 = 13.0 * 13.0 + 17.0 * 17.0 + 19.0 * 19.0;
#endif

    prj_z4c_build_tmunu_from_matter(&mesh, 0, 0, 0);
    tmunu = prj_block_z4c_tmunu_stage(block, 0);
    assert_close_rel("tmunu E",
        tmunu[TMUNUIDX(PRJ_TMUNU_E, i, j, k)],
        fac * (rho * PRJ_CLIGHT * PRJ_CLIGHT + rho * eint +
            0.5 * rho * (v1 * v1 + v2 * v2 + v3 * v3) + 0.5 * B2),
        1.0e-14);
    assert_close_rel("tmunu Sx", tmunu[TMUNUIDX(PRJ_TMUNU_SX, i, j, k)],
        fac * rho * v1 * PRJ_CLIGHT, 1.0e-14);
    assert_close_rel("tmunu Sxy", tmunu[TMUNUIDX(PRJ_TMUNU_SXY, i, j, k)],
        fac * (rho * v1 * v2
#if PRJ_MHD
        - 13.0 * 17.0
#endif
        ), 1.0e-14);
    assert_close_rel("tmunu Sxx", tmunu[TMUNUIDX(PRJ_TMUNU_SXX, i, j, k)],
        fac * (rho * v1 * v1 + pressure
#if PRJ_MHD
        + 0.5 * B2 - 13.0 * 13.0
#endif
        ), 1.0e-14);
    prj_z4c_clear_tmunu(&mesh, 0, 0);
    assert_close("tmunu clear", tmunu[TMUNUIDX(PRJ_TMUNU_E, i, j, k)], 0.0, 0.0);
    prj_mesh_destroy(&mesh);
}

#if PRJ_USE_RADIATION_M1
static double test_m1_chi_exact(double f)
{
    return (3.0 + 4.0 * f * f) / (5.0 + 2.0 * sqrt(4.0 - 3.0 * f * f));
}

static void check_tmunu_m1_projection(void)
{
    prj_mesh mesh;
    prj_coord coord;
    prj_block *block;
    prj_rad rad;
    double *W_rad;
    double *tmunu;
    const int i = 0, j = 0, k = 0;
    const double Erad = 2.0e20;
    const double f = 0.5;
    const double Frad = f * PRJ_CLIGHT * Erad;
    const double chi = test_m1_chi_exact(f);
    const double a_c = 0.5 * (1.0 - chi);
    const double b_c = 0.5 * (3.0 * chi - 1.0);
    double fac = geo_factor();
    int n;

    init_one_block_mesh(&mesh, &coord);
    block = &mesh.blocks[0];
    memset(&rad, 0, sizeof(rad));
    for (n = 0; n <= NCLOSURE; ++n) {
        rad.chi[n] = test_m1_chi_exact((double)n / (double)NCLOSURE);
    }
    prj_z4c_init_mesh_flat(&mesh, 0);
    prj_fill(block->W_mhd, (size_t)PRJ_NVAR_MHD_PRIM *
        (size_t)PRJ_BLOCK_NSTAGES * (size_t)PRJ_BLOCK_NCELLS, 0.0);
    prj_fill(block->W_rad, (size_t)PRJ_NVAR_RAD_PRIM *
        (size_t)PRJ_BLOCK_NSTAGES * (size_t)PRJ_BLOCK_NCELLS, 0.0);
    W_rad = prj_block_rad_stage(block, 0);
    W_rad[WIDX(PRJ_RAD_PRIM_E(0, 0), i, j, k)] = Erad / RAD_SCALE;
    W_rad[WIDX(PRJ_RAD_PRIM_F1(0, 0), i, j, k)] = Frad / RAD_SCALE;
    prj_z4c_build_tmunu_from_matter(&mesh, 0, &rad, 0);
    tmunu = prj_block_z4c_tmunu_stage(block, 0);
    assert_close_rel("m1 E", tmunu[TMUNUIDX(PRJ_TMUNU_E, i, j, k)],
        fac * Erad, 1.0e-14);
    assert_close_rel("m1 Sx", tmunu[TMUNUIDX(PRJ_TMUNU_SX, i, j, k)],
        fac * Frad / PRJ_CLIGHT, 1.0e-14);
    assert_close_rel("m1 Sxx", tmunu[TMUNUIDX(PRJ_TMUNU_SXX, i, j, k)],
        fac * Erad * (a_c + b_c), 1.0e-14);
    assert_close_rel("m1 Syy", tmunu[TMUNUIDX(PRJ_TMUNU_SYY, i, j, k)],
        fac * Erad * a_c, 1.0e-14);
    prj_mesh_destroy(&mesh);
}
#endif

static void check_rhs_matter_terms(void)
{
    prj_mesh mesh;
    prj_coord coord;
    prj_block *block;
    double *tmunu;
    double *rhs;
    const int i = 2, j = 2, k = 2;
    const double pi = acos(-1.0);
    const double E = 2.0e-6;
    const double Sx = 3.0e-6;
    const double Sxx = 5.0e-6;
    const double Sxy = 13.0e-6;
    const double Syy = 7.0e-6;
    const double Szz = 11.0e-6;
    const double S = Sxx + Syy + Szz;

    init_one_block_mesh(&mesh, &coord);
    block = &mesh.blocks[0];
    prj_z4c_init_mesh_flat(&mesh, 0);
    tmunu = prj_block_z4c_tmunu_stage(block, 0);
    tmunu[TMUNUIDX(PRJ_TMUNU_E, i, j, k)] = E;
    tmunu[TMUNUIDX(PRJ_TMUNU_SX, i, j, k)] = Sx;
    tmunu[TMUNUIDX(PRJ_TMUNU_SXX, i, j, k)] = Sxx;
    tmunu[TMUNUIDX(PRJ_TMUNU_SXY, i, j, k)] = Sxy;
    tmunu[TMUNUIDX(PRJ_TMUNU_SYY, i, j, k)] = Syy;
    tmunu[TMUNUIDX(PRJ_TMUNU_SZZ, i, j, k)] = Szz;
    prj_z4c_compute_rhs(&mesh, 0, 0, 0, 0.0);
    rhs = prj_block_z4c_rhs_stage(block, 0);
    assert_close("rhs Khat", rhs[Z4CIDX(PRJ_Z4C_KHAT, i, j, k)],
        4.0 * pi * (S + E), 1.0e-18);
    assert_close("rhs Theta", rhs[Z4CIDX(PRJ_Z4C_THETA, i, j, k)],
        -8.0 * pi * E, 1.0e-18);
    assert_close("rhs Gamx", rhs[Z4CIDX(PRJ_Z4C_GAMX, i, j, k)],
        -16.0 * pi * Sx, 1.0e-18);
    assert_close("rhs Axx", rhs[Z4CIDX(PRJ_Z4C_AXX, i, j, k)],
        -8.0 * pi * (Sxx - S / 3.0), 1.0e-18);
    assert_close("rhs Axy", rhs[Z4CIDX(PRJ_Z4C_AXY, i, j, k)],
        -8.0 * pi * Sxy, 1.0e-18);
    assert_close("rhs Ayy", rhs[Z4CIDX(PRJ_Z4C_AYY, i, j, k)],
        -8.0 * pi * (Syy - S / 3.0), 1.0e-18);
    assert_close("rhs Azz", rhs[Z4CIDX(PRJ_Z4C_AZZ, i, j, k)],
        -8.0 * pi * (Szz - S / 3.0), 1.0e-18);
    prj_mesh_destroy(&mesh);
}
#endif

int main(int argc, char **argv)
{
#if defined(PRJ_ENABLE_MPI)
    MPI_Init(&argc, &argv);
#else
    (void)argc;
    (void)argv;
#endif
    check_enum_names();
#if PRJ_DYNAMIC_GR
    check_flat_state();
    check_tmunu_hydro_projection();
#if PRJ_USE_RADIATION_M1
    check_tmunu_m1_projection();
#endif
    check_rhs_matter_terms();
#endif
    printf("test_z4c: ok\n");
#if defined(PRJ_ENABLE_MPI)
    MPI_Finalize();
#endif
    return 0;
}
