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
    fprintf(stderr, "test_gr_m1_closure_consistency: %s\n", msg);
    exit(1);
}

static void assert_close(const char *name, double got, double expected,
    double rel)
{
    double scale = fmax(1.0, fmax(fabs(got), fabs(expected)));
    double tol = rel * scale;

    if (!isfinite(got) || !isfinite(expected) || fabs(got - expected) > tol) {
        fprintf(stderr,
            "test_gr_m1_closure_consistency: %s got %.17e expected %.17e tol %.3e\n",
            name, got, expected, tol);
        exit(1);
    }
}

#if PRJ_DYNAMIC_GR && PRJ_USE_RADIATION_M1 && PRJ_NRAD > 0
int prj_flux_gr_m1_closure_test_wrapper(const prj_z4c_hydro_geom *geom,
    double E, const double Fcov[3], double U[4], double F[4]);
int prj_src_gr_m1_pressure_contractions_test_wrapper(
    const prj_z4c_hydro_geom *geom, double E, const double Fcov[3],
    double *pK_out, double pDgamma[3]);
int prj_rad_gr_m1_m3_from_build_R_test_wrapper(
    const prj_rad_gr_m1_closure_ctx *ctx, double E, const double Fcov[3],
    double *J_out, double H_out[4], double L_out[4][4],
    double Q_out[4][4][4]);
int prj_rad_gr_m1_freq_m3_test_wrapper(
    const prj_rad_gr_m1_closure_ctx *ctx,
    const prj_z4c_hydro_geom *geom, double E, const double Fcov[3],
    double *J_out, double H_out[4], double L_out[4][4],
    double Q_out[4][4][4]);
int prj_rad_gr_m1_freq_drift_test_wrapper(
    const prj_rad_gr_m1_closure_ctx *ctx,
    const prj_z4c_hydro_geom *geom, double E, const double Fcov[3],
    double fast[4], double ref[4]);
int prj_rad_gr_m1_build_R_jac_test_wrapper(
    const prj_z4c_hydro_geom *geom, double E, const double Fcov[3],
    double Rcon[4][4]);
int prj_rad_gr_m1_projected_jac_test_wrapper(
    const prj_z4c_hydro_geom *geom, double E, const double Fcov[3],
    const double cov[4], double R_u[4], double dR_u[4][4],
    double dRuu[4]);
int prj_rad_gr_m1_full_jac_projection_test_wrapper(
    const prj_z4c_hydro_geom *geom, double E, const double Fcov[3],
    const double cov[4], double R_u[4], double dR_u[4][4],
    double dRuu[4]);

static double det3(const double g[3][3])
{
    return g[0][0] * (g[1][1] * g[2][2] - g[1][2] * g[2][1])
        - g[0][1] * (g[1][0] * g[2][2] - g[1][2] * g[2][0])
        + g[0][2] * (g[1][0] * g[2][1] - g[1][1] * g[2][0]);
}

static void inv3(const double g[3][3], double inv[3][3], double *det_out)
{
    double det = det3(g);
    double odet;

    if (!isfinite(det) || det <= 0.0) {
        die("invalid test metric determinant");
    }
    odet = 1.0 / det;
    inv[0][0] = (g[1][1] * g[2][2] - g[1][2] * g[2][1]) * odet;
    inv[0][1] = (g[0][2] * g[2][1] - g[0][1] * g[2][2]) * odet;
    inv[0][2] = (g[0][1] * g[1][2] - g[0][2] * g[1][1]) * odet;
    inv[1][0] = inv[0][1];
    inv[1][1] = (g[0][0] * g[2][2] - g[0][2] * g[2][0]) * odet;
    inv[1][2] = (g[0][2] * g[1][0] - g[0][0] * g[1][2]) * odet;
    inv[2][0] = inv[0][2];
    inv[2][1] = inv[1][2];
    inv[2][2] = (g[0][0] * g[1][1] - g[0][1] * g[1][0]) * odet;
    if (det_out != 0) {
        *det_out = det;
    }
}

static void metric4_from_geom(const prj_z4c_hydro_geom *geom,
    double g_cov[4][4], double g_con[4][4])
{
    double beta_cov[3] = {0.0, 0.0, 0.0};
    double beta2 = 0.0;
    double inv_alpha2;
    int a;
    int b;

    memset(g_cov, 0, 16 * sizeof(double));
    memset(g_con, 0, 16 * sizeof(double));
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            beta_cov[a] += geom->gamma[a][b] * geom->beta[b];
        }
        beta2 += beta_cov[a] * geom->beta[a];
    }
    g_cov[0][0] = -geom->alpha * geom->alpha + beta2;
    for (a = 0; a < 3; ++a) {
        g_cov[0][a + 1] = beta_cov[a];
        g_cov[a + 1][0] = beta_cov[a];
        for (b = 0; b < 3; ++b) {
            g_cov[a + 1][b + 1] = geom->gamma[a][b];
        }
    }
    inv_alpha2 = 1.0 / (geom->alpha * geom->alpha);
    g_con[0][0] = -inv_alpha2;
    for (a = 0; a < 3; ++a) {
        g_con[0][a + 1] = geom->beta[a] * inv_alpha2;
        g_con[a + 1][0] = g_con[0][a + 1];
        for (b = 0; b < 3; ++b) {
            g_con[a + 1][b + 1] = geom->gamma_inv[a][b] -
                geom->beta[a] * geom->beta[b] * inv_alpha2;
        }
    }
}

static void normal_metric_from_geom(const prj_z4c_hydro_geom *geom,
    double g_cov[4][4], double g_con[4][4])
{
    int a;
    int b;

    memset(g_cov, 0, 16 * sizeof(double));
    memset(g_con, 0, 16 * sizeof(double));
    g_cov[0][0] = -1.0;
    g_con[0][0] = -1.0;
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            g_cov[a + 1][b + 1] = geom->gamma[a][b];
            g_con[a + 1][b + 1] = geom->gamma_inv[a][b];
        }
    }
}

static void raise_vec(const prj_z4c_hydro_geom *geom, const double vcov[3],
    double vcon[3])
{
    int a;
    int b;

    for (a = 0; a < 3; ++a) {
        vcon[a] = 0.0;
        for (b = 0; b < 3; ++b) {
            vcon[a] += geom->gamma_inv[a][b] * vcov[b];
        }
    }
}

static void init_geom(prj_z4c_hydro_geom *geom)
{
    static const double gamma[3][3] = {
        {1.21, 0.08, -0.04},
        {0.08, 0.92, 0.03},
        {-0.04, 0.03, 1.37}
    };
    static const double K[3][3] = {
        {2.0e-6, -0.7e-6, 0.4e-6},
        {-0.7e-6, -1.3e-6, 0.5e-6},
        {0.4e-6, 0.5e-6, 0.8e-6}
    };
    double det;
    int a;
    int b;
    int d;

    memset(geom, 0, sizeof(*geom));
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            geom->gamma[a][b] = gamma[a][b];
            geom->K_dd[a][b] = K[a][b];
        }
    }
    inv3(geom->gamma, geom->gamma_inv, &det);
    geom->sqrt_gamma = sqrt(det);
    geom->alpha = 0.83;
    geom->beta[0] = 0.035;
    geom->beta[1] = -0.022;
    geom->beta[2] = 0.014;
    geom->dalpha[0] = 1.1e-7;
    geom->dalpha[1] = -0.6e-7;
    geom->dalpha[2] = 0.8e-7;
    for (d = 0; d < 3; ++d) {
        for (a = 0; a < 3; ++a) {
            geom->dbeta[d][a] =
                (0.11 + 0.03 * (double)d - 0.02 * (double)a) * 1.0e-7;
            for (b = a; b < 3; ++b) {
                double value = (0.2 + 0.17 * (double)(d + 1) +
                    0.05 * (double)(a + 1) - 0.04 * (double)(b + 1)) *
                    1.0e-7;

                geom->dgamma[d][a][b] = value;
                geom->dgamma[d][b][a] = value;
            }
        }
    }
}

static void init_closure_ctx(const prj_z4c_hydro_geom *geom,
    prj_rad_gr_m1_closure_ctx *ctx)
{
    static const double vcon[3] = {0.052, -0.031, 0.019};
    int a;
    int b;
    int d;

    memset(ctx, 0, sizeof(*ctx));
    for (a = 0; a < 3; ++a) {
        ctx->vcon[a] = vcon[a];
        for (b = 0; b < 3; ++b) {
            ctx->gamma[a][b] = geom->gamma[a][b];
            ctx->gamma_inv[a][b] = geom->gamma_inv[a][b];
            ctx->K_dd[a][b] = geom->K_dd[a][b];
            for (d = 0; d < 3; ++d) {
                ctx->dgamma[d][a][b] = geom->dgamma[d][a][b];
            }
        }
    }
    for (d = 0; d < 3; ++d) {
        for (a = 0; a < 3; ++a) {
            ctx->dvdx[d][a] =
                (0.09 + 0.04 * (double)d - 0.03 * (double)a) * 1.0e-7;
        }
    }
}

static void set_flux(const prj_z4c_hydro_geom *geom, double E, double f,
    const double dir[3], double Fcov[3])
{
    double dir_cov[3] = {0.0, 0.0, 0.0};
    double norm2 = 0.0;
    double scale;
    int a;
    int b;

    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            dir_cov[a] += geom->gamma[a][b] * dir[b];
        }
        norm2 += dir_cov[a] * dir[a];
    }
    if (!isfinite(norm2) || norm2 <= 0.0) {
        die("invalid flux direction");
    }
    scale = PRJ_CLIGHT * E * f / sqrt(norm2);
    for (a = 0; a < 3; ++a) {
        Fcov[a] = scale * dir_cov[a];
    }
}

static void assert_matrix4_close(const char *name, const double got[4][4],
    const double expected[4][4], double rel)
{
    char label[128];
    int a;
    int b;

    for (a = 0; a < 4; ++a) {
        for (b = 0; b < 4; ++b) {
            snprintf(label, sizeof(label), "%s[%d][%d]", name, a, b);
            assert_close(label, got[a][b], expected[a][b], rel);
        }
    }
}

static void check_flux_fastpath(const prj_z4c_hydro_geom *geom, double E,
    const double Fcov[3], const double Rfull[4][4],
    const double g_cov[4][4])
{
    double U_fast[4];
    double F_fast[4];
    double U_ref[4];
    double F_ref[4];
    int a;
    int b;

    if (!prj_flux_gr_m1_closure_test_wrapper(geom, E, Fcov, U_fast, F_fast)) {
        die("flux closure wrapper failed");
    }
    U_ref[0] = geom->sqrt_gamma * geom->alpha * geom->alpha * Rfull[0][0];
    F_ref[0] = PRJ_CLIGHT * geom->sqrt_gamma * geom->alpha *
        geom->alpha * Rfull[0][1];
    for (a = 0; a < 3; ++a) {
        double mixed0 = 0.0;
        double mixed1 = 0.0;

        for (b = 0; b < 4; ++b) {
            mixed0 += Rfull[0][b] * g_cov[b][a + 1];
            mixed1 += Rfull[1][b] * g_cov[b][a + 1];
        }
        U_ref[a + 1] = PRJ_CLIGHT * geom->alpha *
            geom->sqrt_gamma * mixed0;
        F_ref[a + 1] = PRJ_CLIGHT * PRJ_CLIGHT * geom->alpha *
            geom->sqrt_gamma * mixed1;
    }
    for (a = 0; a < 4; ++a) {
        assert_close("flux closure U", U_fast[a], U_ref[a], 2.0e-13);
        assert_close("flux closure F", F_fast[a], F_ref[a], 2.0e-13);
    }
}

static void check_source_fastpath(const prj_z4c_hydro_geom *geom, double E,
    const double Fcov[3], const double Rnormal[4][4])
{
    double pK_fast;
    double pD_fast[3];
    double pK_ref = 0.0;
    double pD_ref[3] = {0.0, 0.0, 0.0};
    int a;
    int b;
    int d;

    if (!prj_src_gr_m1_pressure_contractions_test_wrapper(geom, E, Fcov,
            &pK_fast, pD_fast)) {
        die("source pressure contraction wrapper failed");
    }
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            double Pab = Rnormal[a + 1][b + 1];

            pK_ref += Pab * geom->K_dd[a][b];
            for (d = 0; d < 3; ++d) {
                pD_ref[d] += Pab * geom->dgamma[d][a][b];
            }
        }
    }
    assert_close("source pK", pK_fast, pK_ref, 2.0e-13);
    for (d = 0; d < 3; ++d) {
        assert_close("source pDgamma", pD_fast[d], pD_ref[d], 2.0e-13);
    }
}

static void check_freq_fastpath(const prj_z4c_hydro_geom *geom,
    const prj_rad_gr_m1_closure_ctx *ctx, double E, const double Fcov[3])
{
    double J_ref;
    double J_fast;
    double H_ref[4];
    double H_fast[4];
    double L_ref[4][4];
    double L_fast[4][4];
    double Q_ref[4][4][4];
    double Q_fast[4][4][4];
    double drift_fast[4];
    double drift_ref[4];
    int a;
    int b;
    int c;

    if (!prj_rad_gr_m1_m3_from_build_R_test_wrapper(ctx, E, Fcov, &J_ref,
            H_ref, L_ref, Q_ref)) {
        die("generic M3 wrapper failed");
    }
    if (!prj_rad_gr_m1_freq_m3_test_wrapper(ctx, geom, E, Fcov, &J_fast,
            H_fast, L_fast, Q_fast)) {
        die("freq M3 wrapper failed");
    }
    assert_close("freq J", J_fast, J_ref, 2.0e-13);
    for (a = 0; a < 4; ++a) {
        assert_close("freq H", H_fast[a], H_ref[a], 2.0e-13);
        for (b = 0; b < 4; ++b) {
            assert_close("freq L", L_fast[a][b], L_ref[a][b], 2.0e-13);
            for (c = 0; c < 4; ++c) {
                assert_close("freq Q", Q_fast[a][b][c], Q_ref[a][b][c],
                    4.0e-13);
            }
        }
    }
    if (!prj_rad_gr_m1_freq_drift_test_wrapper(ctx, geom, E, Fcov,
            drift_fast, drift_ref)) {
        die("freq drift wrapper failed");
    }
    for (a = 0; a < 4; ++a) {
        assert_close("freq drift", drift_fast[a], drift_ref[a], 8.0e-13);
    }
}

static unsigned int freq_rng_state = 0x4d595df4u;

static void check_case(const char *name,
    const prj_z4c_hydro_geom *geom,
    const prj_rad_gr_m1_closure_ctx *ctx, double E, double f,
    const double dir[3]);

static double freq_rand_unit(void)
{
    freq_rng_state = 1664525u * freq_rng_state + 1013904223u;
    return (double)(freq_rng_state >> 8) / 16777216.0;
}

static double freq_rand_signed(double scale)
{
    return scale * (2.0 * freq_rand_unit() - 1.0);
}

static void randomize_freq_case(prj_z4c_hydro_geom *geom,
    prj_rad_gr_m1_closure_ctx *ctx, double dir[3])
{
    double lower[3][3] = {{0.0}};
    double det;
    int a;
    int b;
    int d;
    int m;

    memset(geom, 0, sizeof(*geom));
    for (a = 0; a < 3; ++a) {
        lower[a][a] = 0.8 + 0.4 * freq_rand_unit();
        for (b = 0; b < a; ++b) {
            lower[a][b] = freq_rand_signed(0.12);
        }
    }
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            for (m = 0; m < 3; ++m) {
                geom->gamma[a][b] += lower[a][m] * lower[b][m];
            }
            geom->K_dd[a][b] = freq_rand_signed(2.0e-6);
            for (d = 0; d < 3; ++d) {
                geom->dgamma[d][a][b] = freq_rand_signed(2.0e-7);
            }
        }
    }
    inv3(geom->gamma, geom->gamma_inv, &det);
    geom->sqrt_gamma = sqrt(det);
    geom->alpha = 0.7 + 0.6 * freq_rand_unit();
    for (a = 0; a < 3; ++a) {
        geom->beta[a] = freq_rand_signed(0.04);
        geom->dalpha[a] = freq_rand_signed(2.0e-7);
        dir[a] = freq_rand_signed(1.0);
        for (d = 0; d < 3; ++d) {
            geom->dbeta[d][a] = freq_rand_signed(2.0e-7);
        }
    }
    init_closure_ctx(geom, ctx);
    for (a = 0; a < 3; ++a) {
        ctx->vcon[a] = freq_rand_signed(0.08);
        for (d = 0; d < 3; ++d) {
            ctx->dvdx[d][a] = freq_rand_signed(2.0e-7);
        }
    }
}

static void check_freq_randomized(void)
{
    prj_z4c_hydro_geom geom;
    prj_rad_gr_m1_closure_ctx ctx;
    double dir[3];
    int n;

    for (n = 0; n < 64; ++n) {
        double E = 0.2 + 3.0 * freq_rand_unit();
        double f;

        randomize_freq_case(&geom, &ctx, dir);
        if (n % 4 == 0) {
            f = 0.0;
        } else if (n % 4 == 1) {
            f = 1.0e-10;
        } else if (n % 4 == 2) {
            f = 1.0 - 2.0 * PRJ_RAD_GR_M1_F_MARGIN;
        } else {
            f = 0.98 * freq_rand_unit();
        }
        check_case("randomized", &geom, &ctx, E, f, dir);
    }
}

static void check_implicit_fastpath(const prj_z4c_hydro_geom *geom, double E,
    const double Fcov[3], const double Rfull[4][4])
{
    static const double cov[4] = {-1.07, 0.031, -0.027, 0.019};
    double Rjac[4][4];
    double Ru_fast[4];
    double Ru_ref[4];
    double dRu_fast[4][4];
    double dRu_ref[4][4];
    double dRuu_fast[4];
    double dRuu_ref[4];
    int a;
    int b;

    if (!prj_rad_gr_m1_build_R_jac_test_wrapper(geom, E, Fcov, Rjac)) {
        die("implicit R_jac wrapper failed");
    }
    assert_matrix4_close("implicit R_jac", Rjac, Rfull, 2.0e-13);
    if (!prj_rad_gr_m1_projected_jac_test_wrapper(geom, E, Fcov, cov,
            Ru_fast, dRu_fast, dRuu_fast)) {
        die("implicit projected jac wrapper failed");
    }
    if (!prj_rad_gr_m1_full_jac_projection_test_wrapper(geom, E, Fcov, cov,
            Ru_ref, dRu_ref, dRuu_ref)) {
        die("implicit full jac projection wrapper failed");
    }
    for (a = 0; a < 4; ++a) {
        assert_close("implicit projected R_u", Ru_fast[a], Ru_ref[a],
            2.0e-13);
        assert_close("implicit projected dRuu", dRuu_fast[a], dRuu_ref[a],
            2.0e-13);
        for (b = 0; b < 4; ++b) {
            assert_close("implicit projected dR_u", dRu_fast[a][b],
                dRu_ref[a][b], 2.0e-13);
        }
    }
}

static void check_case(const char *name, const prj_z4c_hydro_geom *geom,
    const prj_rad_gr_m1_closure_ctx *ctx, double E, double f,
    const double dir[3])
{
    double Fcov[3];
    double Fcon[3];
    double g_cov[4][4];
    double g_con[4][4];
    double gn_cov[4][4];
    double gn_con[4][4];
    double Rfull[4][4];
    double Rnormal[4][4];

    (void)name;
    set_flux(geom, E, f, dir, Fcov);
    raise_vec(geom, Fcov, Fcon);
    metric4_from_geom(geom, g_cov, g_con);
    normal_metric_from_geom(geom, gn_cov, gn_con);
    if (!prj_rad_grm1_build_R(g_cov, g_con, geom->alpha, E, Fcon, Rfull)) {
        die("build_R full metric failed");
    }
    if (!prj_rad_grm1_build_R(gn_cov, gn_con, 1.0, E, Fcon, Rnormal)) {
        die("build_R normal metric failed");
    }
    check_flux_fastpath(geom, E, Fcov, Rfull, g_cov);
    check_source_fastpath(geom, E, Fcov, Rnormal);
    check_freq_fastpath(geom, ctx, E, Fcov);
    check_implicit_fastpath(geom, E, Fcov, Rfull);
}

static void check_gr_m1_closure_consistency(void)
{
    prj_z4c_hydro_geom geom;
    prj_rad_gr_m1_closure_ctx ctx;
    const double dir1[3] = {1.0, 0.3, -0.2};
    const double dir2[3] = {-0.4, 0.8, 0.25};
    const double dir3[3] = {0.2, -0.3, 1.0};

    init_geom(&geom);
    init_closure_ctx(&geom, &ctx);
    check_case("moderate", &geom, &ctx, 2.3, 0.37, dir1);
    check_case("nearly isotropic", &geom, &ctx, 0.9, 1.0e-8, dir2);
    check_case("thin", &geom, &ctx, 1.4, 0.92, dir3);
    check_freq_randomized();
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
#if PRJ_DYNAMIC_GR && PRJ_USE_RADIATION_M1 && PRJ_NRAD > 0
    check_gr_m1_closure_consistency();
    printf("test_gr_m1_closure_consistency: ok\n");
#else
    printf("test_gr_m1_closure_consistency: skipped\n");
#endif
#if defined(PRJ_ENABLE_MPI)
    MPI_Finalize();
#endif
    return 0;
}
