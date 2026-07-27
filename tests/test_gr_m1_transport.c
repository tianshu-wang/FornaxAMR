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
    fprintf(stderr, "test_gr_m1_transport: %s\n", msg);
    exit(1);
}

static void assert_close(const char *name, double got, double expected, double rel)
{
    double scale = fmax(1.0, fmax(fabs(got), fabs(expected)));
    double tol = rel * scale;

    if (!isfinite(got) || !isfinite(expected) || fabs(got - expected) > tol) {
        fprintf(stderr,
            "test_gr_m1_transport: %s got %.17e expected %.17e tol %.3e\n",
            name, got, expected, tol);
        exit(1);
    }
}

#if PRJ_DYNAMIC_GR && PRJ_USE_RADIATION_M1 && PRJ_NRAD > 0
static double test_m1_chi_exact(double f)
{
    if (f <= 0.0) {
        return 1.0 / 3.0;
    }
    if (f >= 1.0) {
        return 1.0;
    }
    return (3.0 + 4.0 * f * f) / (5.0 + 2.0 * sqrt(4.0 - 3.0 * f * f));
}

static double test_m1_chi_lookup(const prj_rad *rad, double f)
{
    (void)rad;
    if (f <= 0.0) {
        return test_m1_chi_exact(0.0);
    }
    if (f >= 1.0) {
        return test_m1_chi_exact(1.0);
    }
    return test_m1_chi_exact(f);
}

static void init_test_rad(prj_rad *rad)
{
    int n;

    memset(rad, 0, sizeof(*rad));
    for (n = 0; n <= NCLOSURE; ++n) {
        double f = (double)n / (double)NCLOSURE;

        rad->chi[n] = test_m1_chi_exact(f);
    }
}

/* Flat GR-vs-non-GR comparisons cross the exact GR closure with the legacy
 * tabulated M1 closure, so choose flux factors exactly on chi(f) table nodes. */
static void set_flat_flux_at_chi_node(double E, int node,
    const double dir[3], double Fcov[3])
{
    double norm2 = dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2];
    double norm;
    double f;
    double scaled;
    int a;

    if (node < 0 || node > NCLOSURE || norm2 <= 0.0 || E <= 0.0) {
        die("invalid closure-table flux node");
    }
    norm = sqrt(norm2);
    f = (double)node / (double)NCLOSURE;
    for (a = 0; a < 3; ++a) {
        Fcov[a] = PRJ_CLIGHT * E * f * dir[a] / norm;
    }
    scaled = sqrt(Fcov[0] * Fcov[0] + Fcov[1] * Fcov[1] +
            Fcov[2] * Fcov[2]) / (PRJ_CLIGHT * E) * (double)NCLOSURE;
    if (fabs(scaled - (double)node) > 1.0e-12 * fmax(1.0, scaled)) {
        die("flux factor missed closure-table node");
    }
}

static void make_closure_ctx(const prj_z4c_hydro_geom *geom,
    const double vcon[3], const double dvdx[3][3], int have_shear,
    double opacity, prj_rad_gr_m1_closure_ctx *ctx)
{
    int a;
    int b;
    int d;

    (void)have_shear;
    (void)opacity;
    memset(ctx, 0, sizeof(*ctx));
    for (a = 0; a < 3; ++a) {
        ctx->vcon[a] = vcon != 0 ? vcon[a] : 0.0;
        for (b = 0; b < 3; ++b) {
            ctx->gamma[a][b] = geom->gamma[a][b];
            ctx->gamma_inv[a][b] = geom->gamma_inv[a][b];
            ctx->K_dd[a][b] = geom->K_dd[a][b];
            for (d = 0; d < 3; ++d) {
                ctx->dgamma[d][a][b] = geom->dgamma[d][a][b];
            }
        }
    }
    if (dvdx != 0) {
        for (d = 0; d < 3; ++d) {
            for (a = 0; a < 3; ++a) {
                ctx->dvdx[d][a] = dvdx[d][a];
            }
        }
    }
}

static void check_grm1_build_R_flat(void)
{
    double g_cov[4][4] = {{0.0}};
    double g_inv[4][4] = {{0.0}};
    double ER_expected = 1.7;
    double v[3] = {0.23, -0.11, 0.07};
    double v2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    double w = 1.0 / sqrt(1.0 - v2);
    double uR[4];
    double expected[4][4];
    double R[4][4];
    double Fcon[3];
    double E;
    int ok;
    int a;
    int b;

    g_cov[0][0] = -1.0;
    g_cov[1][1] = 1.0;
    g_cov[2][2] = 1.0;
    g_cov[3][3] = 1.0;
    g_inv[0][0] = -1.0;
    g_inv[1][1] = 1.0;
    g_inv[2][2] = 1.0;
    g_inv[3][3] = 1.0;
    uR[0] = w;
    for (a = 0; a < 3; ++a) {
        uR[a + 1] = w * v[a];
    }

    for (a = 0; a < 4; ++a) {
        for (b = 0; b < 4; ++b) {
            expected[a][b] = ER_expected *
                ((4.0 / 3.0) * uR[a] * uR[b] +
                 (1.0 / 3.0) * g_inv[a][b]);
        }
    }

    E = expected[0][0];
    for (a = 0; a < 3; ++a) {
        Fcon[a] = PRJ_CLIGHT * expected[0][a + 1];
    }

    ok = prj_rad_grm1_build_R(g_cov, g_inv, 1.0, E, Fcon, R);
    if (!ok) {
        die("flat build_R failed");
    }
    for (a = 0; a < 4; ++a) {
        for (b = 0; b < 4; ++b) {
            assert_close("flat build_R", R[a][b], expected[a][b],
                2.0e-13);
        }
    }
}

static void build_test_metric4(double alpha, const double gamma_diag[3],
    const double beta_con[3], double g_cov[4][4], double g_inv[4][4])
{
    double beta_cov[3];
    double beta2 = 0.0;
    int a;
    int b;

    memset(g_cov, 0, 16 * sizeof(double));
    memset(g_inv, 0, 16 * sizeof(double));
    for (a = 0; a < 3; ++a) {
        beta_cov[a] = gamma_diag[a] * beta_con[a];
        beta2 += beta_cov[a] * beta_con[a];
    }
    g_cov[0][0] = -alpha * alpha + beta2;
    g_inv[0][0] = -1.0 / (alpha * alpha);
    for (a = 0; a < 3; ++a) {
        g_cov[0][a + 1] = beta_cov[a];
        g_cov[a + 1][0] = beta_cov[a];
        g_cov[a + 1][a + 1] = gamma_diag[a];
        g_inv[0][a + 1] = beta_con[a] / (alpha * alpha);
        g_inv[a + 1][0] = g_inv[0][a + 1];
        for (b = 0; b < 3; ++b) {
            g_inv[a + 1][b + 1] =
                (a == b ? 1.0 / gamma_diag[a] : 0.0) -
                beta_con[a] * beta_con[b] / (alpha * alpha);
        }
    }
}

static void check_grm1_build_R_shifted_metric(void)
{
    double alpha = 0.83;
    double gamma_diag[3] = {1.2, 0.9, 1.5};
    double beta_con[3] = {0.08, -0.03, 0.04};
    double vcon[3] = {0.18, -0.07, 0.05};
    double g_cov[4][4];
    double g_inv[4][4];
    double ncon[4];
    double ucon[4];
    double expected[4][4];
    double R[4][4];
    double Fcon[3];
    double ER_expected = 0.9;
    double v2 = 0.0;
    double w;
    double E;
    int ok;
    int a;
    int b;

    build_test_metric4(alpha, gamma_diag, beta_con, g_cov, g_inv);
    for (a = 0; a < 3; ++a) {
        v2 += gamma_diag[a] * vcon[a] * vcon[a];
    }
    w = 1.0 / sqrt(1.0 - v2);
    ncon[0] = 1.0 / alpha;
    ucon[0] = w * ncon[0];
    for (a = 0; a < 3; ++a) {
        ncon[a + 1] = -beta_con[a] / alpha;
        ucon[a + 1] = w * (ncon[a + 1] + vcon[a]);
    }

    for (a = 0; a < 4; ++a) {
        for (b = 0; b < 4; ++b) {
            expected[a][b] = ER_expected *
                ((4.0 / 3.0) * ucon[a] * ucon[b] +
                 (1.0 / 3.0) * g_inv[a][b]);
        }
    }
    E = alpha * alpha * expected[0][0];
    for (a = 0; a < 3; ++a) {
        double Fhat = alpha * expected[0][a + 1] - E * ncon[a + 1];

        Fcon[a] = PRJ_CLIGHT * Fhat;
    }

    ok = prj_rad_grm1_build_R(g_cov, g_inv, alpha, E, Fcon, R);
    if (!ok) {
        die("shifted build_R failed");
    }
    for (a = 0; a < 4; ++a) {
        for (b = 0; b < 4; ++b) {
            assert_close("shifted build_R", R[a][b], expected[a][b],
                3.0e-13);
        }
    }
}

static void check_grm1_build_R_free_streaming(void)
{
    double g_cov[4][4] = {{0.0}};
    double g_inv[4][4] = {{0.0}};
    double E = 2.3;
    double Fcon[3] = {2.3 * PRJ_CLIGHT, 0.0, 0.0};
    double R[4][4];
    int ok;
    int a;
    int b;

    g_cov[0][0] = -1.0;
    g_cov[1][1] = 1.0;
    g_cov[2][2] = 1.0;
    g_cov[3][3] = 1.0;
    g_inv[0][0] = -1.0;
    g_inv[1][1] = 1.0;
    g_inv[2][2] = 1.0;
    g_inv[3][3] = 1.0;

    ok = prj_rad_grm1_build_R(g_cov, g_inv, 1.0, E, Fcon, R);
    if (!ok) {
        die("free-streaming build_R failed");
    }
    for (a = 0; a < 4; ++a) {
        for (b = 0; b < 4; ++b) {
            double expected = (a < 2 && b < 2) ? E : 0.0;

            assert_close("free-streaming build_R", R[a][b], expected,
                2.0e-13);
        }
    }
}

static double expected_rest_frame_m3_component(double E, const double H[3],
    double Hmag, const double n[3], double thin_w, double thick_w,
    int a, int b, int c)
{
    int zero_count = (a == 0) + (b == 0) + (c == 0);

    if (zero_count == 3) {
        return E;
    }
    if (zero_count == 2) {
        int s = a != 0 ? a : (b != 0 ? b : c);

        return H[s - 1];
    }
    if (zero_count == 1) {
        int s0;
        int s1;

        if (a == 0) {
            s0 = b;
            s1 = c;
        } else if (b == 0) {
            s0 = a;
            s1 = c;
        } else {
            s0 = a;
            s1 = b;
        }
        {
            double Pthin = Hmag > 0.0 ? E * n[s0 - 1] * n[s1 - 1] : 0.0;
            double Pthick = s0 == s1 ? E / 3.0 : 0.0;

            return thin_w * Pthin + thick_w * Pthick;
        }
    }
    {
        double Nthin = Hmag > 0.0 ? E * n[a - 1] * n[b - 1] * n[c - 1] : 0.0;
        double Nthick = 0.2 *
            (H[a - 1] * (b == c ? 1.0 : 0.0) +
                H[b - 1] * (a == c ? 1.0 : 0.0) +
                H[c - 1] * (a == b ? 1.0 : 0.0));

        return thin_w * Nthin + thick_w * Nthick;
    }
}

static void check_grm1_freq_drift_rest_frame(void)
{
    double g_cov[4][4] = {{0.0}};
    double g_inv[4][4] = {{0.0}};
    double ucon[4] = {1.0, 0.0, 0.0, 0.0};
    double E = 2.0;
    double H[3] = {0.42, -0.24, 0.16};
    double Hmag = sqrt(H[0] * H[0] + H[1] * H[1] + H[2] * H[2]);
    double n[3];
    double chi = test_m1_chi_exact(Hmag / E);
    double thin_w = 0.5 * (3.0 * chi - 1.0);
    double thick_w = 1.5 * (1.0 - chi);
    double R[4][4] = {{0.0}};
    double ducov[4][4];
    double drift[4];
    int ok;
    int a;
    int b;
    int p;
    int q;

    g_cov[0][0] = -1.0;
    g_cov[1][1] = 1.0;
    g_cov[2][2] = 1.0;
    g_cov[3][3] = 1.0;
    g_inv[0][0] = -1.0;
    g_inv[1][1] = 1.0;
    g_inv[2][2] = 1.0;
    g_inv[3][3] = 1.0;
    for (a = 0; a < 3; ++a) {
        n[a] = H[a] / Hmag;
    }
    R[0][0] = E;
    for (a = 0; a < 3; ++a) {
        R[0][a + 1] = H[a];
        R[a + 1][0] = H[a];
        for (b = 0; b < 3; ++b) {
            R[a + 1][b + 1] = expected_rest_frame_m3_component(E, H, Hmag,
                n, thin_w, thick_w, 0, a + 1, b + 1);
        }
    }

    for (p = 0; p < 4; ++p) {
        for (q = 0; q < 4; ++q) {
            memset(ducov, 0, sizeof(ducov));
            ducov[p][q] = 1.0;
            ok = prj_rad_grm1_freq_drift(g_cov, g_inv, ucon, R, ducov,
                drift);
            if (!ok) {
                die("rest-frame freq_drift failed");
            }
            for (a = 0; a < 4; ++a) {
                double expected = expected_rest_frame_m3_component(E, H,
                    Hmag, n, thin_w, thick_w, a, p, q);

                assert_close("rest-frame freq_drift", drift[a], expected,
                    2.0e-13);
            }
        }
    }
}

static void check_grm1_freq_drift_free_streaming(void)
{
    double g_cov[4][4] = {{0.0}};
    double g_inv[4][4] = {{0.0}};
    double ucon[4] = {1.0, 0.0, 0.0, 0.0};
    double R[4][4];
    double ducov[4][4];
    double drift[4];
    double E = 1.9;
    double Fcon[3] = {1.9 * PRJ_CLIGHT, 0.0, 0.0};
    int ok;
    int a;
    int p;
    int q;

    g_cov[0][0] = -1.0;
    g_cov[1][1] = 1.0;
    g_cov[2][2] = 1.0;
    g_cov[3][3] = 1.0;
    g_inv[0][0] = -1.0;
    g_inv[1][1] = 1.0;
    g_inv[2][2] = 1.0;
    g_inv[3][3] = 1.0;

    ok = prj_rad_grm1_build_R(g_cov, g_inv, 1.0, E, Fcon, R);
    if (!ok) {
        die("free-streaming build_R for M3 failed");
    }
    for (p = 0; p < 4; ++p) {
        for (q = 0; q < 4; ++q) {
            memset(ducov, 0, sizeof(ducov));
            ducov[p][q] = 1.0;
            ok = prj_rad_grm1_freq_drift(g_cov, g_inv, ucon, R, ducov,
                drift);
            if (!ok) {
                die("free-streaming freq_drift failed");
            }
            for (a = 0; a < 4; ++a) {
                double expected = (a < 2 && p < 2 && q < 2) ? E : 0.0;

                assert_close("free-streaming freq_drift", drift[a],
                    expected, 2.0e-13);
            }
        }
    }
}

static void check_grm1_freq_drift_u_contraction(void)
{
    double g_cov[4][4] = {{0.0}};
    double g_inv[4][4] = {{0.0}};
    double v[3] = {0.21, -0.09, 0.06};
    double v2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    double w = 1.0 / sqrt(1.0 - v2);
    double ucon[4];
    double ucov[4];
    double E = 1.4;
    double Fcon[3] = {0.32 * 1.4 * PRJ_CLIGHT,
        -0.18 * 1.4 * PRJ_CLIGHT, 0.11 * 1.4 * PRJ_CLIGHT};
    double R[4][4];
    double ducov[4][4];
    double drift[4];
    int ok;
    int a;
    int c;
    int p;

    g_cov[0][0] = -1.0;
    g_cov[1][1] = 1.0;
    g_cov[2][2] = 1.0;
    g_cov[3][3] = 1.0;
    g_inv[0][0] = -1.0;
    g_inv[1][1] = 1.0;
    g_inv[2][2] = 1.0;
    g_inv[3][3] = 1.0;
    ucon[0] = w;
    ucov[0] = -w;
    for (a = 0; a < 3; ++a) {
        ucon[a + 1] = w * v[a];
        ucov[a + 1] = ucon[a + 1];
    }

    ok = prj_rad_grm1_build_R(g_cov, g_inv, 1.0, E, Fcon, R);
    if (!ok) {
        die("contraction build_R for M3 failed");
    }
    for (p = 0; p < 4; ++p) {
        memset(ducov, 0, sizeof(ducov));
        for (c = 0; c < 4; ++c) {
            ducov[p][c] = ucov[c];
        }
        ok = prj_rad_grm1_freq_drift(g_cov, g_inv, ucon, R, ducov, drift);
        if (!ok) {
            die("u-contraction freq_drift failed");
        }
        for (a = 0; a < 4; ++a) {
            assert_close("freq_drift u contraction", drift[a], -R[a][p],
                2.0e-12);
        }
    }
}

static void expected_zero_velocity_gr_pressure(const prj_rad *rad,
    const prj_z4c_hydro_geom *geom, double E, const double Fcov[3],
    double P[3][3])
{
    double g_cov[4][4] = {{0.0}};
    double g_inv[4][4] = {{0.0}};
    double Rcon[4][4];
    double Fcon[3] = {0.0, 0.0, 0.0};
    double F2 = 0.0;
    double Fmag;
    double cE;
    int a;
    int b;

    (void)rad;
    memset(P, 0, 9 * sizeof(double));
    if (!isfinite(E) || E < 0.0) {
        E = 0.0;
    }
    g_cov[0][0] = -1.0;
    g_inv[0][0] = -1.0;
    for (a = 0; a < 3; ++a) {
        double Fcov_a = isfinite(Fcov[a]) ? Fcov[a] : 0.0;

        for (b = 0; b < 3; ++b) {
            g_cov[a + 1][b + 1] = geom->gamma[a][b];
            g_inv[a + 1][b + 1] = geom->gamma_inv[a][b];
            Fcon[a] += geom->gamma_inv[a][b] *
                (isfinite(Fcov[b]) ? Fcov[b] : 0.0);
        }
        F2 += Fcov_a * Fcon[a];
    }
    if (!isfinite(F2) || F2 < 0.0) {
        F2 = 0.0;
    }
    Fmag = sqrt(F2);
    cE = PRJ_CLIGHT * E;
    if (Fmag > cE && Fmag > 0.0) {
        double scale = cE / Fmag;

        for (a = 0; a < 3; ++a) {
            Fcon[a] *= scale;
        }
    }
    if (!prj_rad_grm1_build_R(g_cov, g_inv, 1.0, E, Fcon, Rcon)) {
        return;
    }
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            P[a][b] = Rcon[a + 1][b + 1];
        }
    }
}

#if 0
/* Expected M1 pressure for a boosted state on a flat metric: the explicit
 * blend P = thin_w*Pthin + thick_w*Pthick with the CLOSED-FORM inviscid
 * optically-thick tensor (Shibata et al. 2011; the pure-thick inversion of
 * E, F_i into fluid-frame J, H_i), mirroring the production closure in
 * prj_rad_gr_m1_pressure_thick / prj_rad_gr_m1_pressure_for_fbar. On the
 * flat metric, raising indices is the identity (Fcon = Fcov, V^i = u_i). */
static void expected_flat_boosted_pressure_for_fbar(const prj_rad *rad,
    double E, const double Fcov_in[3], const double vcon_in[3], double fbar,
    double P[3][3])
{
    double Fcov[3];
    double Fhat[3];
    double Pthin[3][3];
    double Pthick[3][3];
    double vcon[3];
    double u_cov[3];
    double Hcov[3];
    double F2 = 0.0;
    double Fmag;
    double cE;
    double beta2 = 0.0;
    double wlor;
    double W2;
    double FdotU = 0.0;
    double J;
    double dH;
    double chi;
    double thin_w;
    double thick_w;
    int a;
    int b;

    memset(Pthin, 0, sizeof(Pthin));
    for (a = 0; a < 3; ++a) {
        Fcov[a] = Fcov_in[a];
        F2 += Fcov[a] * Fcov[a];
    }
    Fmag = sqrt(F2);
    cE = PRJ_CLIGHT * E;
    if (Fmag > cE && Fmag > 0.0) {
        double scale = cE / Fmag;

        for (a = 0; a < 3; ++a) {
            Fcov[a] *= scale;
        }
        F2 *= scale * scale;
    }
    if (E > 0.0 && F2 > 0.0) {
        for (a = 0; a < 3; ++a) {
            for (b = 0; b < 3; ++b) {
                Pthin[a][b] = E * Fcov[a] * Fcov[b] / F2;
            }
        }
    }
    for (a = 0; a < 3; ++a) {
        vcon[a] = vcon_in[a];
        beta2 += vcon[a] * vcon[a];
    }
    if (beta2 >= 1.0) {
        double scale = sqrt((1.0 - 1.0e-12) / beta2);

        for (a = 0; a < 3; ++a) {
            vcon[a] *= scale;
        }
        beta2 = 1.0 - 1.0e-12;
    }
    wlor = 1.0 / sqrt(1.0 - beta2);
    W2 = wlor * wlor;
    for (a = 0; a < 3; ++a) {
        u_cov[a] = wlor * vcon[a];
        Fhat[a] = Fcov[a] / PRJ_CLIGHT;
        FdotU += Fhat[a] * u_cov[a];
    }

    J = 3.0 / (2.0 * W2 + 1.0) * ((2.0 * W2 - 1.0) * E - 2.0 * wlor * FdotU);
    dH = wlor * (2.0 * W2 + 1.0);
    for (a = 0; a < 3; ++a) {
        double ui = u_cov[a];

        Hcov[a] = Fhat[a] / wlor +
            (-(4.0 * W2 * wlor * ui) * E + ((4.0 * W2 + 1.0) * ui) * FdotU) /
                dH;
    }
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            Pthick[a][b] = J * ((a == b ? 1.0 : 0.0) +
                    4.0 * u_cov[a] * u_cov[b]) / 3.0 +
                Hcov[a] * u_cov[b] + Hcov[b] * u_cov[a];
        }
    }

    chi = test_m1_chi_lookup(rad, fbar);
    thin_w = 0.5 * (3.0 * chi - 1.0);
    thick_w = 1.5 * (1.0 - chi);
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            P[a][b] = thin_w * Pthin[a][b] + thick_w * Pthick[a][b];
        }
    }
}
#endif

static double flat_boosted_fbar(double E, const double Fcov[3],
    const double vcon[3], const double P[3][3])
{
    double beta2 = 0.0;
    double wlor;
    double u_cov[3];
    double Fhat[3];
    double R0;
    double Rcon[3];
    double J;
    double numerator;
    int a;
    int b;

    for (a = 0; a < 3; ++a) {
        beta2 += vcon[a] * vcon[a];
    }
    wlor = 1.0 / sqrt(1.0 - beta2);
    R0 = -E * wlor;
    J = E * wlor * wlor;
    for (a = 0; a < 3; ++a) {
        u_cov[a] = wlor * vcon[a];
        Fhat[a] = Fcov[a] / PRJ_CLIGHT;
        R0 += Fhat[a] * u_cov[a];
        J -= 2.0 * wlor * Fhat[a] * u_cov[a];
    }
    for (a = 0; a < 3; ++a) {
        Rcon[a] = -wlor * Fhat[a];
        for (b = 0; b < 3; ++b) {
            Rcon[a] += P[a][b] * u_cov[b];
            J += P[a][b] * u_cov[a] * u_cov[b];
        }
    }
    numerator = (wlor * wlor - 1.0) * R0 * R0;
    for (a = 0; a < 3; ++a) {
        numerator -= 2.0 * wlor * u_cov[a] * R0 * Rcon[a];
        for (b = 0; b < 3; ++b) {
            numerator += ((a == b ? 1.0 : 0.0) + u_cov[a] * u_cov[b]) *
                Rcon[a] * Rcon[b];
        }
    }
    if (numerator < 0.0) {
        numerator = 0.0;
    }
    return sqrt(numerator) / fabs(J);
}

static void init_test_eos(prj_eos *eos)
{
    memset(eos, 0, sizeof(*eos));
    eos->kind = PRJ_EOS_KIND_IDEAL;
    prj_eos_init(eos, 0);
}

static void init_test_mesh(prj_mesh *mesh, prj_coord *coord)
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
        die("mesh init failed");
    }
    prj_fill(mesh->blocks[0].W_mhd, prj_block_data_count(), 0.0);
}

static void set_uniform_z4c(prj_block *block, double alpha,
    const double beta[3], const double gamma_diag[3])
{
    double *z = prj_block_z4c_stage(block, 0);
    int i;
    int j;
    int k;

    if (z == 0) {
        die("missing z4c storage");
    }
    for (i = -PRJ_NGHOST_Z4C; i < PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C; ++i) {
        for (j = -PRJ_NGHOST_Z4C; j < PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C; ++j) {
            for (k = -PRJ_NGHOST_Z4C; k < PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C; ++k) {
                z[Z4CIDX(PRJ_Z4C_CHI, i, j, k)] = 1.0;
                z[Z4CIDX(PRJ_Z4C_GXX, i, j, k)] = gamma_diag[0];
                z[Z4CIDX(PRJ_Z4C_GXY, i, j, k)] = 0.0;
                z[Z4CIDX(PRJ_Z4C_GXZ, i, j, k)] = 0.0;
                z[Z4CIDX(PRJ_Z4C_GYY, i, j, k)] = gamma_diag[1];
                z[Z4CIDX(PRJ_Z4C_GYZ, i, j, k)] = 0.0;
                z[Z4CIDX(PRJ_Z4C_GZZ, i, j, k)] = gamma_diag[2];
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
                z[Z4CIDX(PRJ_Z4C_ALPHA, i, j, k)] = alpha;
                z[Z4CIDX(PRJ_Z4C_BETAX, i, j, k)] = beta[0];
                z[Z4CIDX(PRJ_Z4C_BETAY, i, j, k)] = beta[1];
                z[Z4CIDX(PRJ_Z4C_BETAZ, i, j, k)] = beta[2];
            }
        }
    }
}

static void fill_constant_state(prj_block *block, double E, const double Fcov[3])
{
    const double rho = 1.0;
    const double pressure = 1.0;
    const double gamma_gas = 5.0 / 3.0;
    int i;
    int j;
    int k;
    int field;
    int group;

    for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
        for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
            for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                block->W_mhd[WIDX(PRJ_PRIM_RHO, i, j, k)] = rho;
                block->W_mhd[WIDX(PRJ_PRIM_V1, i, j, k)] = 0.0;
                block->W_mhd[WIDX(PRJ_PRIM_V2, i, j, k)] = 0.0;
                block->W_mhd[WIDX(PRJ_PRIM_V3, i, j, k)] = 0.0;
                block->W_mhd[WIDX(PRJ_PRIM_EINT, i, j, k)] =
                    pressure / ((gamma_gas - 1.0) * rho);
                block->W_mhd[WIDX(PRJ_PRIM_YE, i, j, k)] = 0.2;
#if PRJ_MHD
                block->W_mhd[WIDX(PRJ_PRIM_B1, i, j, k)] = 0.0;
                block->W_mhd[WIDX(PRJ_PRIM_B2, i, j, k)] = 0.0;
                block->W_mhd[WIDX(PRJ_PRIM_B3, i, j, k)] = 0.0;
#endif
                block->eosvar[EIDX(PRJ_EOSVAR_PRESSURE, i, j, k)] = pressure;
                block->eosvar[EIDX(PRJ_EOSVAR_TEMPERATURE, i, j, k)] = pressure;
                block->eosvar[EIDX(PRJ_EOSVAR_GAMMA, i, j, k)] = gamma_gas;
                for (field = 0; field < PRJ_NRAD; ++field) {
                    for (group = 0; group < PRJ_NEGROUP; ++group) {
                        block->W_rad[WIDX(PRJ_RAD_PRIM_E(field, group), i, j, k)] = 0.0;
                        block->W_rad[WIDX(PRJ_RAD_PRIM_F1(field, group), i, j, k)] = 0.0;
                        block->W_rad[WIDX(PRJ_RAD_PRIM_F2(field, group), i, j, k)] = 0.0;
                        block->W_rad[WIDX(PRJ_RAD_PRIM_F3(field, group), i, j, k)] = 0.0;
                    }
                }
                block->W_rad[WIDX(PRJ_RAD_PRIM_E(0, 0), i, j, k)] = E;
                block->W_rad[WIDX(PRJ_RAD_PRIM_F1(0, 0), i, j, k)] = Fcov[0];
                block->W_rad[WIDX(PRJ_RAD_PRIM_F2(0, 0), i, j, k)] = Fcov[1];
                block->W_rad[WIDX(PRJ_RAD_PRIM_F3(0, 0), i, j, k)] = Fcov[2];
            }
        }
    }
#if PRJ_MHD
    for (int dir = 0; dir < 3; ++dir) {
        prj_fill(block->Bf[dir],
            (size_t)PRJ_BLOCK_NSTAGES * (size_t)PRJ_BLOCK_NFACES, 0.0);
    }
#endif
    prj_fill(block->kappa_cell,
        (size_t)PRJ_NRAD * (size_t)PRJ_NEGROUP * (size_t)PRJ_BLOCK_NCELLS, 0.0);
    prj_fill(block->sigma_cell,
        (size_t)PRJ_NRAD * (size_t)PRJ_NEGROUP * (size_t)PRJ_BLOCK_NCELLS, 0.0);
}

static void fill_xface_discontinuity(prj_block *block, int iface,
    double EL, const double FcovL[3], double ER, const double FcovR[3],
    const double vcon[3])
{
    int i;
    int j;
    int k;

    fill_constant_state(block, ER, FcovR);
    for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
        for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
            for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                const double E = i < iface ? EL : ER;
                const double *Fcov = i < iface ? FcovL : FcovR;

                block->W_mhd[WIDX(PRJ_PRIM_V1, i, j, k)] =
                    PRJ_CLIGHT * vcon[0];
                block->W_mhd[WIDX(PRJ_PRIM_V2, i, j, k)] =
                    PRJ_CLIGHT * vcon[1];
                block->W_mhd[WIDX(PRJ_PRIM_V3, i, j, k)] =
                    PRJ_CLIGHT * vcon[2];
                block->W_rad[WIDX(PRJ_RAD_PRIM_E(0, 0), i, j, k)] = E;
                block->W_rad[WIDX(PRJ_RAD_PRIM_F1(0, 0), i, j, k)] = Fcov[0];
                block->W_rad[WIDX(PRJ_RAD_PRIM_F2(0, 0), i, j, k)] = Fcov[1];
                block->W_rad[WIDX(PRJ_RAD_PRIM_F3(0, 0), i, j, k)] = Fcov[2];
            }
        }
    }
}

static void flux_at_xface(prj_eos *eos, prj_rad *rad, prj_mesh *mesh,
    int iface, int j, int k, double out[4])
{
    prj_block *block = &mesh->blocks[0];
    double *flux[3] = {block->flux[0], block->flux[1], block->flux[2]};

    prj_flux_update(eos, rad, mesh, block, block->W_mhd, block->eosvar, flux, 0);
    out[0] = block->flux[X1DIR][VIDX(PRJ_CONS_RAD_E(0, 0), iface, j, k)];
    out[1] = block->flux[X1DIR][VIDX(PRJ_CONS_RAD_F1(0, 0), iface, j, k)];
    out[2] = block->flux[X1DIR][VIDX(PRJ_CONS_RAD_F2(0, 0), iface, j, k)];
    out[3] = block->flux[X1DIR][VIDX(PRJ_CONS_RAD_F3(0, 0), iface, j, k)];
}

static void expected_gr_m1_limit_state(const prj_z4c_hydro_geom *geom,
    double *E, double Fcov[3], double Fcon[3], double *Fmag_out)
{
    double F2 = 0.0;
    double Fmag;
    double cE;
    int a;
    int b;

    if (!isfinite(*E) || *E < PRJ_RAD_GR_M1_E_FLOOR) {
        *E = PRJ_RAD_GR_M1_E_FLOOR;
    }
    for (a = 0; a < 3; ++a) {
        Fcon[a] = 0.0;
        for (b = 0; b < 3; ++b) {
            Fcon[a] += geom->gamma_inv[a][b] * Fcov[b];
        }
        F2 += Fcov[a] * Fcon[a];
    }
    if (!isfinite(F2) || F2 < 0.0) {
        F2 = 0.0;
    }
    Fmag = sqrt(F2);
    cE = (1.0 - PRJ_RAD_GR_M1_F_MARGIN) * PRJ_CLIGHT * *E;
    if (Fmag > cE && Fmag > 0.0) {
        double scale = cE / Fmag;

        for (a = 0; a < 3; ++a) {
            Fcov[a] *= scale;
            Fcon[a] *= scale;
        }
        Fmag = cE;
    }
    *Fmag_out = Fmag;
}

static void expected_gr_m1_speeds(const prj_rad *rad,
    const prj_z4c_hydro_geom *geom, const double vcon[3],
    const double Fcon[3], double Fmag, double zeta,
    double *smin, double *smax)
{
    double beta2 = 0.0;
    double wlor;
    double wlor2;
    double p;
    double r2;
    double r;
    double den;
    double thin_speed = Fmag > 0.0 ? geom->alpha * fabs(Fcon[0]) / Fmag : 0.0;
    double lambda_thin_l = -geom->beta[0] - thin_speed;
    double lambda_thin_r = -geom->beta[0] + thin_speed;
    double lambda_thick_l_a;
    double lambda_thick_r_a;
    double lambda_fluid;
    double lambda_thick_l;
    double lambda_thick_r;
    double chi;
    double thin_w;
    double thick_w;
    double lambda_l;
    double lambda_r;
    int a;
    int b;

    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            beta2 += geom->gamma[a][b] * vcon[a] * vcon[b];
        }
    }
    if (!isfinite(beta2) || beta2 < 0.0) {
        beta2 = 0.0;
    }
    wlor = 1.0 / sqrt(1.0 - beta2);
    wlor2 = wlor * wlor;
    /* p = alpha V^x / w = alpha v^x (Shibata et al. 2011). */
    p = geom->alpha * vcon[0];
    r2 = geom->alpha * geom->alpha * geom->gamma_inv[0][0] *
        (2.0 * wlor2 + 1.0) - 2.0 * wlor2 * p * p;
    if (!isfinite(r2) || r2 < 0.0) {
        r2 = 0.0;
    }
    r = sqrt(r2);
    den = 2.0 * wlor2 + 1.0;
    lambda_fluid = -geom->beta[0] + p;
    lambda_thick_l_a = -geom->beta[0] + (2.0 * p * wlor2 - r) / den;
    lambda_thick_r_a = -geom->beta[0] + (2.0 * p * wlor2 + r) / den;
    lambda_thick_l = lambda_thick_l_a < lambda_fluid ?
        lambda_thick_l_a : lambda_fluid;
    lambda_thick_r = lambda_thick_r_a > lambda_fluid ?
        lambda_thick_r_a : lambda_fluid;
    chi = test_m1_chi_lookup(rad, zeta);
    thin_w = 0.5 * (3.0 * chi - 1.0);
    thick_w = 1.5 * (1.0 - chi);
    lambda_l = thin_w * lambda_thin_l + thick_w * lambda_thick_l;
    lambda_r = thin_w * lambda_thin_r + thick_w * lambda_thick_r;
    *smin = PRJ_CLIGHT * lambda_l;
    *smax = PRJ_CLIGHT * lambda_r;
}

static void expected_gr_m1_side_energy_flux(const prj_rad *rad,
    const prj_z4c_hydro_geom *geom, const double vcon[3],
    double Ein, const double Fcov_in[3],
    double *U, double *Fphys, double *smin, double *smax)
{
    prj_rad_gr_m1_closure_ctx ctx;
    double Fcov[3] = {Fcov_in[0], Fcov_in[1], Fcov_in[2]};
    double Fcon[3];
    double Pcon[3][3];
    double E = Ein;
    double Fmag;
    double fbar = 0.0;

    expected_gr_m1_limit_state(geom, &E, Fcov, Fcon, &Fmag);
    make_closure_ctx(geom, vcon, 0, 0, 0.0, &ctx);
    prj_rad_gr_m1_pressure_fbar(rad, &ctx, E, Fcov, Pcon, &fbar);
    expected_gr_m1_speeds(rad, geom, vcon, Fcon, Fmag, fbar, smin, smax);
    *U = geom->sqrt_gamma * E;
    *Fphys = geom->sqrt_gamma *
        (geom->alpha * Fcon[0] - PRJ_CLIGHT * geom->beta[0] * E);
}

static double expected_gr_m1_hll_energy_flux(const prj_rad *rad,
    const prj_z4c_hydro_geom *geom, const double vcon[3],
    double EL, const double FcovL[3], double ER, const double FcovR[3])
{
    double UL;
    double UR;
    double FL;
    double FR;
    double sminL;
    double smaxL;
    double sminR;
    double smaxR;
    double sL;
    double sR;

    expected_gr_m1_side_energy_flux(rad, geom, vcon, EL, FcovL,
        &UL, &FL, &sminL, &smaxL);
    expected_gr_m1_side_energy_flux(rad, geom, vcon, ER, FcovR,
        &UR, &FR, &sminR, &smaxR);
    sL = sminL < sminR ? sminL : sminR;
    sR = smaxL > smaxR ? smaxL : smaxR;
    if (sL > 0.0) {
        sL = 0.0;
    }
    if (sR < 0.0) {
        sR = 0.0;
    }
    if (sR - sL < 1.0e-30) {
        sL = -PRJ_CLIGHT;
        sR = PRJ_CLIGHT;
    }
    return (sR * FL - sL * FR + sL * sR * (UR - UL)) / (sR - sL);
}

static void set_combined_rad_state(double *W, double E, const double F[3])
{
    W[PRJ_PRIM_RAD_E(0, 0)] = E;
    W[PRJ_PRIM_RAD_F1(0, 0)] = F[0];
    W[PRJ_PRIM_RAD_F2(0, 0)] = F[1];
    W[PRJ_PRIM_RAD_F3(0, 0)] = F[2];
}

static void check_flat_zero_shift_matches_non_gr(void)
{
    prj_mesh mesh;
    prj_coord coord;
    prj_eos eos;
    prj_rad rad;
    double beta[3] = {0.0, 0.0, 0.0};
    double gamma_diag[3] = {1.0, 1.0, 1.0};
    double Fcov[3];
    double Fdir[3] = {12.0, -4.0, 3.0};
    double got[4];
    double expected[PRJ_NVAR_CONS] = {0.0};
    double Wface[PRJ_NVAR_PRIM] = {0.0};
    double chi_face[PRJ_NRAD * PRJ_NEGROUP] = {0.0};

    init_test_mesh(&mesh, &coord);
    init_test_eos(&eos);
    init_test_rad(&rad);
    set_uniform_z4c(&mesh.blocks[0], 1.0, beta, gamma_diag);
    set_flat_flux_at_chi_node(2.0, 13, Fdir, Fcov);
    fill_constant_state(&mesh.blocks[0], 2.0, Fcov);
    flux_at_xface(&eos, &rad, &mesh, PRJ_BLOCK_SIZE / 2, 1, 1, got);
    set_combined_rad_state(Wface, 2.0, Fcov);
    prj_rad_flux(&rad, Wface, Wface, 1.0, chi_face, mesh.blocks[0].dx[X1DIR],
        0.0, expected);
    assert_close("flat E flux", got[0], expected[PRJ_CONS_RAD_E(0, 0)], 2.0e-12);
    assert_close("flat F1 flux", got[1], expected[PRJ_CONS_RAD_F1(0, 0)], 2.0e-12);
    assert_close("flat F2 flux", got[2], expected[PRJ_CONS_RAD_F2(0, 0)], 2.0e-12);
    assert_close("flat F3 flux", got[3], expected[PRJ_CONS_RAD_F3(0, 0)], 2.0e-12);
    prj_mesh_destroy(&mesh);
}

static void check_flat_shift_terms(void)
{
    prj_mesh mesh;
    prj_coord coord;
    prj_eos eos;
    prj_rad rad;
    double beta[3] = {0.025, -0.01, 0.0};
    double gamma_diag[3] = {1.0, 1.0, 1.0};
    double E = 1.7;
    double Fcov[3];
    double Fdir[3] = {2.0, 3.0, -6.0};
    double got[4];
    double base[PRJ_NVAR_CONS] = {0.0};
    double expected[4];
    double Wface[PRJ_NVAR_PRIM] = {0.0};
    double chi_face[PRJ_NRAD * PRJ_NEGROUP] = {0.0};

    init_test_mesh(&mesh, &coord);
    init_test_eos(&eos);
    init_test_rad(&rad);
    set_uniform_z4c(&mesh.blocks[0], 1.0, beta, gamma_diag);
    set_flat_flux_at_chi_node(E, 7, Fdir, Fcov);
    fill_constant_state(&mesh.blocks[0], E, Fcov);
    flux_at_xface(&eos, &rad, &mesh, PRJ_BLOCK_SIZE / 2, 1, 1, got);
    set_combined_rad_state(Wface, E, Fcov);
    prj_rad_flux(&rad, Wface, Wface, 1.0, chi_face, mesh.blocks[0].dx[X1DIR],
        0.0, base);
    expected[0] = base[PRJ_CONS_RAD_E(0, 0)] - PRJ_CLIGHT * beta[0] * E;
    expected[1] = base[PRJ_CONS_RAD_F1(0, 0)] - PRJ_CLIGHT * beta[0] * Fcov[0];
    expected[2] = base[PRJ_CONS_RAD_F2(0, 0)] - PRJ_CLIGHT * beta[0] * Fcov[1];
    expected[3] = base[PRJ_CONS_RAD_F3(0, 0)] - PRJ_CLIGHT * beta[0] * Fcov[2];
    assert_close("shift E flux", got[0], expected[0], 2.0e-12);
    assert_close("shift F1 flux", got[1], expected[1], 2.0e-12);
    assert_close("shift F2 flux", got[2], expected[2], 2.0e-12);
    assert_close("shift F3 flux", got[3], expected[3], 2.0e-12);
    prj_mesh_destroy(&mesh);
}

static void check_curved_diagonal_flux(void)
{
    prj_mesh mesh;
    prj_coord coord;
    prj_eos eos;
    prj_rad rad;
    double beta[3] = {0.018, 0.0, 0.0};
    double gamma_diag[3] = {1.8, 0.7, 1.4};
    double alpha = 0.83;
    double E = 2.3;
    double Fcov[3] = {0.09 * PRJ_CLIGHT, -0.03 * PRJ_CLIGHT, 0.02 * PRJ_CLIGHT};
    prj_z4c_hydro_geom geom;
    double Pcon[3][3];
    double sqrt_gamma = sqrt(gamma_diag[0] * gamma_diag[1] * gamma_diag[2]);
    double got[4];
    double expected[4];
    int d;
    int a;

    init_test_mesh(&mesh, &coord);
    init_test_eos(&eos);
    init_test_rad(&rad);
    set_uniform_z4c(&mesh.blocks[0], alpha, beta, gamma_diag);
    fill_constant_state(&mesh.blocks[0], E, Fcov);
    flux_at_xface(&eos, &rad, &mesh, PRJ_BLOCK_SIZE / 2, 1, 1, got);
    if (!prj_z4c_load_hydro_geom(&mesh, &mesh.blocks[0], 0,
            PRJ_BLOCK_SIZE / 2, 1, 1, &geom)) {
        die("curved flux geometry load failed");
    }
    expected_zero_velocity_gr_pressure(&rad, &geom, E, Fcov, Pcon);
    expected[0] = sqrt_gamma * (alpha * (Fcov[0] / gamma_diag[0]) -
        PRJ_CLIGHT * beta[0] * E);
    for (d = 0; d < 3; ++d) {
        double Pn_i = 0.0;

        for (a = 0; a < 3; ++a) {
            Pn_i += geom.gamma[d][a] * Pcon[0][a];
        }

        expected[1 + d] = sqrt_gamma *
            (alpha * PRJ_CLIGHT * PRJ_CLIGHT * Pn_i -
                PRJ_CLIGHT * beta[0] * Fcov[d]);
    }
    assert_close("curved E flux", got[0], expected[0], 2.0e-12);
    assert_close("curved F1 flux", got[1], expected[1], 2.0e-12);
    assert_close("curved F2 flux", got[2], expected[2], 2.0e-12);
    assert_close("curved F3 flux", got[3], expected[3], 2.0e-12);
    prj_mesh_destroy(&mesh);
}

static void check_gr_m1_characteristic_speed_case(const char *name,
    double alpha, const double beta[3], const double gamma_diag[3],
    const double vcon[3], double EL, const double FcovL[3],
    double ER, const double FcovR[3])
{
    prj_mesh mesh;
    prj_coord coord;
    prj_eos eos;
    prj_rad rad;
    prj_z4c_hydro_geom geom;
    const int iface = PRJ_BLOCK_SIZE / 2;
    double got[4];
    double expected;

    init_test_mesh(&mesh, &coord);
    init_test_eos(&eos);
    init_test_rad(&rad);
    set_uniform_z4c(&mesh.blocks[0], alpha, beta, gamma_diag);
    fill_xface_discontinuity(&mesh.blocks[0], iface, EL, FcovL, ER, FcovR,
        vcon);
    flux_at_xface(&eos, &rad, &mesh, iface, 1, 1, got);
    if (!prj_z4c_load_hydro_geom(&mesh, &mesh.blocks[0], 0, iface, 1, 1,
            &geom)) {
        die("characteristic-speed geometry load failed");
    }
    expected = expected_gr_m1_hll_energy_flux(&rad, &geom, vcon, EL, FcovL,
        ER, FcovR);
    assert_close(name, got[0], expected, 3.0e-12);
    prj_mesh_destroy(&mesh);
}

static void check_gr_m1_characteristic_speeds(void)
{
    {
        double beta[3] = {0.0, 0.0, 0.0};
        double gamma_diag[3] = {1.0, 1.0, 1.0};
        double vcon[3] = {0.22, -0.05, 0.03};
        double Fzero[3] = {0.0, 0.0, 0.0};

        check_gr_m1_characteristic_speed_case("GR M1 thick speed energy flux",
            0.91, beta, gamma_diag, vcon, 2.0, Fzero, 1.3, Fzero);
    }
    {
        double beta[3] = {0.017, -0.006, 0.004};
        double gamma_diag[3] = {1.7, 0.8, 1.3};
        double vcon[3] = {0.12, -0.04, 0.03};
        double FcovL[3] = {0.32 * PRJ_CLIGHT, -0.11 * PRJ_CLIGHT,
            0.08 * PRJ_CLIGHT};
        double FcovR[3] = {0.18 * PRJ_CLIGHT, 0.07 * PRJ_CLIGHT,
            -0.04 * PRJ_CLIGHT};

        check_gr_m1_characteristic_speed_case(
            "GR M1 interpolated speed energy flux", 0.82, beta, gamma_diag,
            vcon, 2.4, FcovL, 1.6, FcovR);
    }
}

static void check_gr_pressure_closure_zero_velocity(void)
{
    prj_mesh mesh;
    prj_coord coord;
    prj_rad rad;
    prj_z4c_hydro_geom geom;
    prj_rad_gr_m1_closure_ctx ctx;
    double beta[3] = {0.0, 0.0, 0.0};
    double gamma_diag[3] = {1.6, 0.8, 1.25};
    double Fcov[3] = {0.11 * PRJ_CLIGHT, -0.02 * PRJ_CLIGHT, 0.04 * PRJ_CLIGHT};
    double vzero[3] = {0.0, 0.0, 0.0};
    double got[3][3];
    double expected[3][3];
    double E = 2.1;
    int a;
    int b;

    init_test_mesh(&mesh, &coord);
    init_test_rad(&rad);
    set_uniform_z4c(&mesh.blocks[0], 1.0, beta, gamma_diag);
    if (!prj_z4c_load_hydro_geom(&mesh, &mesh.blocks[0], 0, 2, 2, 2, &geom)) {
        die("closure geometry load failed");
    }
    make_closure_ctx(&geom, vzero, 0, 0, 0.0, &ctx);
    prj_rad_gr_m1_pressure(&rad, &ctx, E, Fcov, got);
    expected_zero_velocity_gr_pressure(&rad, &geom, E, Fcov, expected);
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            assert_close("GR pressure closure", got[a][b], expected[a][b],
                2.0e-12);
        }
    }
    prj_mesh_destroy(&mesh);
}

static void check_gr_pressure_closure_boosted_fbar(void)
{
    prj_mesh mesh;
    prj_coord coord;
    prj_rad rad;
    prj_z4c_hydro_geom geom;
    prj_rad_gr_m1_closure_ctx ctx;
    double beta[3] = {0.0, 0.0, 0.0};
    double gamma_diag[3] = {1.0, 1.0, 1.0};
    double Fcov[3] = {0.16 * PRJ_CLIGHT, -0.03 * PRJ_CLIGHT,
        0.05 * PRJ_CLIGHT};
    double vcon[3] = {0.24, -0.13, 0.08};
    double got[3][3];
    double expected[3][3];
    double E = 2.4;
    double fbar;
    double feuler;
    int a;
    int b;

    init_test_mesh(&mesh, &coord);
    init_test_rad(&rad);
    set_uniform_z4c(&mesh.blocks[0], 1.0, beta, gamma_diag);
    if (!prj_z4c_load_hydro_geom(&mesh, &mesh.blocks[0], 0, 2, 2, 2, &geom)) {
        die("boosted closure geometry load failed");
    }
    make_closure_ctx(&geom, vcon, 0, 0, 0.0, &ctx);
    prj_rad_gr_m1_pressure(&rad, &ctx, E, Fcov, got);
    fbar = flat_boosted_fbar(E, Fcov, vcon, got);
    feuler = sqrt(Fcov[0] * Fcov[0] + Fcov[1] * Fcov[1] +
        Fcov[2] * Fcov[2]) / (PRJ_CLIGHT * E);
    if (fabs(fbar - feuler) < 1.0e-3) {
        die("boosted closure did not move away from Eulerian flux factor");
    }
    expected_zero_velocity_gr_pressure(&rad, &geom, E, Fcov, expected);
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            assert_close("boosted GR pressure closure", got[a][b],
                expected[a][b], 2.0e-12);
        }
    }
    prj_mesh_destroy(&mesh);
}

static void check_gr_pressure_closure_small_velocity_build_R(void)
{
    prj_mesh mesh;
    prj_coord coord;
    prj_rad rad;
    prj_z4c_hydro_geom geom;
    prj_rad_gr_m1_closure_ctx ctx;
    double beta[3] = {0.0, 0.0, 0.0};
    double gamma_diag[3] = {1.0, 1.0, 1.0};
    double Fcov[3] = {0.13 * PRJ_CLIGHT, -0.025 * PRJ_CLIGHT,
        0.04 * PRJ_CLIGHT};
    double vcon[3] = {6.0e-5, -3.0e-5, 2.0e-5};
    double got[3][3];
    double expected[3][3];
    double E = 2.2;
    int a;
    int b;

    init_test_mesh(&mesh, &coord);
    init_test_rad(&rad);
    set_uniform_z4c(&mesh.blocks[0], 1.0, beta, gamma_diag);
    if (!prj_z4c_load_hydro_geom(&mesh, &mesh.blocks[0], 0, 2, 2, 2, &geom)) {
        die("small-velocity closure geometry load failed");
    }
    make_closure_ctx(&geom, vcon, 0, 0, 0.0, &ctx);
    prj_rad_gr_m1_pressure(&rad, &ctx, E, Fcov, got);
    expected_zero_velocity_gr_pressure(&rad, &geom, E, Fcov, expected);
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            assert_close("small-velocity GR pressure closure", got[a][b],
                expected[a][b], 2.0e-12);
        }
    }
    prj_mesh_destroy(&mesh);
}

/* The geom argument supplies the actually-loaded lapse gradient: even for a
 * uniform lapse the shared Z4c stencil leaves O(eps) coefficient-rounding
 * noise in dalpha, which the (physically correct) gravitational-shift terms
 * of the drift pick up; folding geom->dalpha into the expectation keeps this
 * test tight at 1e-10 instead of at the FD noise floor. */
static void expected_rest_frame_gr_frequency_terms(const prj_rad *rad,
    const prj_z4c_hydro_geom *geom, double E, const double Fcov[3],
    const double grad[3][3], double *energy_drift, double momentum_drift[3])
{
    double H[3];
    double n[3];
    double Hmag = 0.0;
    double fbar;
    double chi;
    double thin_w;
    double thick_w;
    int a;
    int b;
    int d;

    *energy_drift = 0.0;
    momentum_drift[0] = 0.0;
    momentum_drift[1] = 0.0;
    momentum_drift[2] = 0.0;
    for (a = 0; a < 3; ++a) {
        H[a] = Fcov[a] / PRJ_CLIGHT;
        Hmag += H[a] * H[a];
    }
    Hmag = sqrt(Hmag);
    fbar = E > 0.0 ? Hmag / E : 0.0;
    (void)rad;
    chi = test_m1_chi_exact(fbar);
    thin_w = 0.5 * (3.0 * chi - 1.0);
    thick_w = 1.5 * (1.0 - chi);
    for (a = 0; a < 3; ++a) {
        n[a] = Hmag > 0.0 ? H[a] / Hmag : 0.0;
    }
    for (a = 0; a < 3; ++a) {
        *energy_drift += Fcov[a] * geom->dalpha[a] / geom->alpha;
        for (b = 0; b < 3; ++b) {
            double Pthin = E * n[a] * n[b];
            double Pthick = (a == b ? E / 3.0 : 0.0);
            double P = thin_w * Pthin + thick_w * Pthick;

            *energy_drift += P * grad[a][b];
            momentum_drift[a] += PRJ_CLIGHT * P * geom->dalpha[b] /
                geom->alpha;
            for (d = 0; d < 3; ++d) {
                double Nthin = E * n[a] * n[b] * n[d];
                double Nthick = 0.2 *
                    (H[a] * (b == d ? 1.0 : 0.0) +
                        H[b] * (a == d ? 1.0 : 0.0) +
                        H[d] * (a == b ? 1.0 : 0.0));
                double N = thin_w * Nthin + thick_w * Nthick;

                momentum_drift[a] += N * grad[d][b];
            }
        }
    }
    for (a = 0; a < 3; ++a) {
        momentum_drift[a] *= PRJ_CLIGHT;
    }
}

static void check_gr_m1_frequency_third_moment_rest_frame_case(
    double grad_scale, int expect_upper_donor)
{
    prj_mesh mesh;
    prj_coord coord;
    prj_rad rad;
    prj_block *block;
    prj_z4c_hydro_geom geom;
    double beta[3] = {0.0, 0.0, 0.0};
    double gamma_diag[3] = {1.0, 1.0, 1.0};
    double eedge_store[PRJ_NRAD][PRJ_NEGROUP + 1];
    double u[PRJ_NVAR_CONS] = {0.0};
    double E0 = 2.0;
    double flux_factor_vec[3] = {0.09, -0.025, 0.015};
    double Fcov0[3] = {E0 * 0.09 * PRJ_CLIGHT,
        -E0 * 0.025 * PRJ_CLIGHT, E0 * 0.015 * PRJ_CLIGHT};
    double base_grad[3][3] = {
        {0.017, -0.004, 0.006},
        {0.003, -0.011, 0.005},
        {-0.002, 0.007, 0.013}
    };
    double grad[3][3];
    double expected_A[PRJ_NEGROUP];
    double expected_mq[PRJ_NEGROUP][3];
    double expected_energy_face[PRJ_NEGROUP + 1] = {0.0};
    double expected_momentum_face[PRJ_NEGROUP + 1][3] = {{0.0}};
    double dt = 1.0;
    double observer_time_derivative[4] = {0.0, 0.0, 0.0, 0.0};
    int field;
    int group;
    int gf;
    int dir;
    int comp;
    int a;
    const int i = 2;
    const int j = 2;
    const int k = 2;
    const int gtest = 2;
    char name[128];

    init_test_mesh(&mesh, &coord);
    init_test_rad(&rad);
    block = &mesh.blocks[0];
    set_uniform_z4c(block, 1.0, beta, gamma_diag);
    fill_constant_state(block, 0.0, Fcov0);
    if (!prj_z4c_load_hydro_geom(&mesh, block, 0, i, j, k, &geom)) {
        die("rest frame drift geometry load failed");
    }

    for (field = 0; field < PRJ_NRAD; ++field) {
        for (gf = 0; gf <= PRJ_NEGROUP; ++gf) {
            eedge_store[field][gf] = 1.0 + (double)gf;
        }
        rad.eedge[field] = eedge_store[field];
    }
    for (dir = 0; dir < 3; ++dir) {
        for (comp = 0; comp < 3; ++comp) {
            grad[dir][comp] = grad_scale * base_grad[dir][comp];
        }
    }
    for (group = 0; group < PRJ_NEGROUP; ++group) {
        double E = E0 * (1.0 + 0.07 * (double)group);
        double Fcov[3];

        for (a = 0; a < 3; ++a) {
            Fcov[a] = E * flux_factor_vec[a] * PRJ_CLIGHT;
        }
        block->W_rad[WIDX(PRJ_RAD_PRIM_E(0, group), i, j, k)] = E;
        block->W_rad[WIDX(PRJ_RAD_PRIM_F1(0, group), i, j, k)] = Fcov[0];
        block->W_rad[WIDX(PRJ_RAD_PRIM_F2(0, group), i, j, k)] = Fcov[1];
        block->W_rad[WIDX(PRJ_RAD_PRIM_F3(0, group), i, j, k)] = Fcov[2];
        u[PRJ_CONS_RAD_E(0, group)] = 100.0;
        u[PRJ_CONS_RAD_F1(0, group)] = 0.0;
        u[PRJ_CONS_RAD_F2(0, group)] = 0.0;
        u[PRJ_CONS_RAD_F3(0, group)] = 0.0;
        expected_rest_frame_gr_frequency_terms(&rad, &geom, E, Fcov, grad,
            &expected_A[group], expected_mq[group]);
    }

    for (dir = 0; dir < 3; ++dir) {
        prj_fill(block->v_riemann[dir],
            (size_t)PRJ_NDIM * (size_t)PRJ_BLOCK_NCELLS, 0.0);
        for (comp = 0; comp < 3; ++comp) {
            int ir = i;
            int jr = j;
            int kr = k;

            if (dir == X1DIR) {
                ir = i + 1;
            } else if (dir == X2DIR) {
                jr = j + 1;
            } else {
                kr = k + 1;
            }
            block->v_riemann[dir][VRIDX(comp, i, j, k)] = 0.0;
            block->v_riemann[dir][VRIDX(comp, ir, jr, kr)] =
                grad[dir][comp] * block->dx[dir];
        }
    }

    if (expect_upper_donor && expected_A[gtest] <= 0.0) {
        die("positive GR frequency drift test has non-positive drift");
    }
    if (!expect_upper_donor && expected_A[gtest] >= 0.0) {
        die("negative GR frequency drift test has non-negative drift");
    }

    for (gf = 1; gf < PRJ_NEGROUP; ++gf) {
        double face_drift = expected_A[gf - 1] + expected_A[gf];
        int donor = face_drift >= 0.0 ? gf : gf - 1;
        double nu = eedge_store[0][gf];

        if (expect_upper_donor && donor != gf) {
            die("positive GR frequency drift did not select upper donor");
        }
        if (!expect_upper_donor && donor != gf - 1) {
            die("negative GR frequency drift did not select lower donor");
        }
        expected_energy_face[gf] = nu * expected_A[donor];
        for (a = 0; a < 3; ++a) {
            expected_momentum_face[gf][a] = nu * expected_mq[donor][a];
        }
    }

    prj_rad_freq_flux_apply_gr_m1(&rad, &mesh, block, 0, block->W_rad, u,
        i, j, k, dt, observer_time_derivative);
    snprintf(name, sizeof(name), "GR frequency energy drift %s",
        expect_upper_donor ? "positive" : "negative");
    assert_close(name, u[PRJ_CONS_RAD_E(0, gtest)],
        100.0 + dt * (expected_energy_face[gtest + 1] -
            expected_energy_face[gtest]), 1.0e-10);
    for (a = 0; a < 3; ++a) {
        snprintf(name, sizeof(name), "GR frequency third moment F%d %s",
            a + 1, expect_upper_donor ? "positive" : "negative");
        assert_close(name, u[PRJ_CONS_RAD_F1(0, gtest) + a],
            dt * (expected_momentum_face[gtest + 1][a] -
                expected_momentum_face[gtest][a]), 1.0e-10);
    }
    prj_mesh_destroy(&mesh);
}

static void check_gr_m1_frequency_third_moment_rest_frame(void)
{
    check_gr_m1_frequency_third_moment_rest_frame_case(1.0, 1);
    check_gr_m1_frequency_third_moment_rest_frame_case(-1.0, 0);
}

static void set_linear_source_z4c(prj_block *block, int ic, int jc, int kc,
    double alpha0, const double dalpha[3], const double dbeta[3][3],
    const double dgamma_diag[3][3], const double K[3][3]);

/* Gravitational redshift across frequency bins (Cardall, Endeve & Mezzacappa
 * 2013, Eqs. 146/149): a static fluid (v = 0) in a static, flat-3-metric
 * spacetime with a nonuniform lapse must drift spectral energy at the rate
 *   A     = F^i d_i(alpha) / alpha,
 *   Mq_j  = c^2 P_j^i d_i(alpha) / alpha,
 * with the standard M1 closure P of the local (E, F).  Positive A (flux up
 * the lapse gradient) moves energy toward LOWER bins -- redshift. */
static void check_gr_m1_frequency_redshift_case(double slope_scale,
    int expect_upper_donor)
{
    prj_mesh mesh;
    prj_coord coord;
    prj_rad rad;
    prj_block *block;
    prj_z4c_hydro_geom geom;
    double dalpha_idx[3] = {1.2e-11, -0.4e-11, 0.6e-11};
    double dbeta[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    double dgamma_diag[3][3] =
        {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    double Kzero[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    double eedge_store[PRJ_NRAD][PRJ_NEGROUP + 1];
    double u[PRJ_NVAR_CONS] = {0.0};
    double E0 = 2.0;
    double flux_factor_vec[3] = {0.30, -0.08, 0.05};
    double Fcov0[3] = {0.0, 0.0, 0.0};
    double expected_A[PRJ_NEGROUP];
    double expected_mq[PRJ_NEGROUP][3];
    double expected_energy_face[PRJ_NEGROUP + 1] = {0.0};
    double expected_momentum_face[PRJ_NEGROUP + 1][3] = {{0.0}};
    double dt = 1.0;
    double observer_time_derivative[4] = {0.0, 0.0, 0.0, 0.0};
    int field;
    int group;
    int gf;
    int dir;
    int a;
    int b;
    const int i = 2;
    const int j = 2;
    const int k = 2;
    const int gtest = 2;
    char name[128];

    init_test_mesh(&mesh, &coord);
    init_test_rad(&rad);
    block = &mesh.blocks[0];
    for (a = 0; a < 3; ++a) {
        dalpha_idx[a] *= slope_scale;
    }
    set_linear_source_z4c(block, i, j, k, 1.0, dalpha_idx, dbeta,
        dgamma_diag, Kzero);
    fill_constant_state(block, 0.0, Fcov0);
    for (dir = 0; dir < 3; ++dir) {
        prj_fill(block->v_riemann[dir],
            (size_t)PRJ_NDIM * (size_t)PRJ_BLOCK_NCELLS, 0.0);
    }
    if (!prj_z4c_load_hydro_geom(&mesh, block, 0, i, j, k, &geom)) {
        die("redshift test geometry load failed");
    }
    assert_close("redshift test lapse", geom.alpha, 1.0, 1.0e-12);

    for (field = 0; field < PRJ_NRAD; ++field) {
        for (gf = 0; gf <= PRJ_NEGROUP; ++gf) {
            eedge_store[field][gf] = 1.0 + (double)gf;
        }
        rad.eedge[field] = eedge_store[field];
    }
    for (group = 0; group < PRJ_NEGROUP; ++group) {
        double E = E0 * (1.0 + 0.07 * (double)group);
        double Fcov[3];
        double H[3];
        double n[3];
        double Hmag = 0.0;
        double fbar;
        double chi;
        double thin_w;
        double thick_w;

        for (a = 0; a < 3; ++a) {
            Fcov[a] = E * flux_factor_vec[a] * PRJ_CLIGHT;
        }
        block->W_rad[WIDX(PRJ_RAD_PRIM_E(0, group), i, j, k)] = E;
        block->W_rad[WIDX(PRJ_RAD_PRIM_F1(0, group), i, j, k)] = Fcov[0];
        block->W_rad[WIDX(PRJ_RAD_PRIM_F2(0, group), i, j, k)] = Fcov[1];
        block->W_rad[WIDX(PRJ_RAD_PRIM_F3(0, group), i, j, k)] = Fcov[2];
        u[PRJ_CONS_RAD_E(0, group)] = 100.0;
        u[PRJ_CONS_RAD_F1(0, group)] = 0.0;
        u[PRJ_CONS_RAD_F2(0, group)] = 0.0;
        u[PRJ_CONS_RAD_F3(0, group)] = 0.0;

        for (a = 0; a < 3; ++a) {
            H[a] = Fcov[a] / PRJ_CLIGHT;
            Hmag += H[a] * H[a];
        }
        Hmag = sqrt(Hmag);
        fbar = E > 0.0 ? Hmag / E : 0.0;
        chi = test_m1_chi_exact(fbar);
        thin_w = 0.5 * (3.0 * chi - 1.0);
        thick_w = 1.5 * (1.0 - chi);
        for (a = 0; a < 3; ++a) {
            n[a] = Hmag > 0.0 ? H[a] / Hmag : 0.0;
        }
        expected_A[group] = 0.0;
        for (a = 0; a < 3; ++a) {
            expected_A[group] += Fcov[a] * geom.dalpha[a] / geom.alpha;
            expected_mq[group][a] = 0.0;
            for (b = 0; b < 3; ++b) {
                double P = thin_w * E * n[a] * n[b] +
                    thick_w * (a == b ? E / 3.0 : 0.0);

                expected_mq[group][a] += PRJ_CLIGHT * PRJ_CLIGHT * P *
                    geom.dalpha[b] / geom.alpha;
            }
        }
    }

    if (expect_upper_donor && expected_A[gtest] <= 0.0) {
        die("positive redshift test has non-positive drift");
    }
    if (!expect_upper_donor && expected_A[gtest] >= 0.0) {
        die("negative redshift test has non-negative drift");
    }

    for (gf = 1; gf < PRJ_NEGROUP; ++gf) {
        double face_drift = expected_A[gf - 1] + expected_A[gf];
        int donor = face_drift >= 0.0 ? gf : gf - 1;
        double nu = eedge_store[0][gf];

        expected_energy_face[gf] = nu * expected_A[donor];
        for (a = 0; a < 3; ++a) {
            expected_momentum_face[gf][a] = nu * expected_mq[donor][a];
        }
    }

    prj_rad_freq_flux_apply_gr_m1(&rad, &mesh, block, 0, block->W_rad, u,
        i, j, k, dt, observer_time_derivative);
    snprintf(name, sizeof(name), "GR redshift energy drift %s",
        expect_upper_donor ? "positive" : "negative");
    assert_close(name, u[PRJ_CONS_RAD_E(0, gtest)],
        100.0 + dt * (expected_energy_face[gtest + 1] -
            expected_energy_face[gtest]), 1.0e-10);
    for (a = 0; a < 3; ++a) {
        snprintf(name, sizeof(name), "GR redshift momentum drift F%d %s",
            a + 1, expect_upper_donor ? "positive" : "negative");
        assert_close(name, u[PRJ_CONS_RAD_F1(0, gtest) + a],
            dt * (expected_momentum_face[gtest + 1][a] -
                expected_momentum_face[gtest][a]), 1.0e-10);
    }
    {
        double total = 0.0;

        for (group = 0; group < PRJ_NEGROUP; ++group) {
            total += u[PRJ_CONS_RAD_E(0, group)] - 100.0;
        }
        snprintf(name, sizeof(name), "GR redshift bin conservation %s",
            expect_upper_donor ? "positive" : "negative");
        assert_close(name, total, 0.0, 1.0e-10);
    }
    prj_mesh_destroy(&mesh);
}

static void check_gr_m1_frequency_gravitational_redshift(void)
{
    check_gr_m1_frequency_redshift_case(1.0, 1);
    check_gr_m1_frequency_redshift_case(-1.0, 0);
}

static void set_linear_source_z4c(prj_block *block, int ic, int jc, int kc,
    double alpha0, const double dalpha[3], const double dbeta[3][3],
    const double dgamma_diag[3][3], const double K[3][3])
{
    double *z = prj_block_z4c_stage(block, 0);
    int i;
    int j;
    int k;

    if (z == 0) {
        die("missing source z4c storage");
    }
    for (i = -PRJ_NGHOST_Z4C; i < PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C; ++i) {
        for (j = -PRJ_NGHOST_Z4C; j < PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C; ++j) {
            for (k = -PRJ_NGHOST_Z4C; k < PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C; ++k) {
                double dx = (double)(i - ic);
                double dy = (double)(j - jc);
                double dz = (double)(k - kc);
                double x[3] = {dx, dy, dz};
                int a;
                int d;

                z[Z4CIDX(PRJ_Z4C_CHI, i, j, k)] = 1.0;
                z[Z4CIDX(PRJ_Z4C_GXX, i, j, k)] =
                    1.0 + dgamma_diag[0][0] * dx + dgamma_diag[1][0] * dy +
                    dgamma_diag[2][0] * dz;
                z[Z4CIDX(PRJ_Z4C_GXY, i, j, k)] = 0.0;
                z[Z4CIDX(PRJ_Z4C_GXZ, i, j, k)] = 0.0;
                z[Z4CIDX(PRJ_Z4C_GYY, i, j, k)] =
                    1.0 + dgamma_diag[0][1] * dx + dgamma_diag[1][1] * dy +
                    dgamma_diag[2][1] * dz;
                z[Z4CIDX(PRJ_Z4C_GYZ, i, j, k)] = 0.0;
                z[Z4CIDX(PRJ_Z4C_GZZ, i, j, k)] =
                    1.0 + dgamma_diag[0][2] * dx + dgamma_diag[1][2] * dy +
                    dgamma_diag[2][2] * dz;
                z[Z4CIDX(PRJ_Z4C_KHAT, i, j, k)] = 0.0;
                z[Z4CIDX(PRJ_Z4C_AXX, i, j, k)] = K[0][0];
                z[Z4CIDX(PRJ_Z4C_AXY, i, j, k)] = K[0][1];
                z[Z4CIDX(PRJ_Z4C_AXZ, i, j, k)] = K[0][2];
                z[Z4CIDX(PRJ_Z4C_AYY, i, j, k)] = K[1][1];
                z[Z4CIDX(PRJ_Z4C_AYZ, i, j, k)] = K[1][2];
                z[Z4CIDX(PRJ_Z4C_AZZ, i, j, k)] = K[2][2];
                z[Z4CIDX(PRJ_Z4C_GAMX, i, j, k)] = 0.0;
                z[Z4CIDX(PRJ_Z4C_GAMY, i, j, k)] = 0.0;
                z[Z4CIDX(PRJ_Z4C_GAMZ, i, j, k)] = 0.0;
                z[Z4CIDX(PRJ_Z4C_THETA, i, j, k)] = 0.0;
                z[Z4CIDX(PRJ_Z4C_ALPHA, i, j, k)] =
                    alpha0 + dalpha[0] * dx + dalpha[1] * dy + dalpha[2] * dz;
                for (a = 0; a < 3; ++a) {
                    double beta = 0.01 * (double)(a + 1);

                    for (d = 0; d < 3; ++d) {
                        beta += dbeta[d][a] * x[d];
                    }
                    z[Z4CIDX(PRJ_Z4C_BETAX + a, i, j, k)] = beta;
                }
            }
        }
    }
}

static void fill_velocity_faces_with_gradient(prj_block *block)
{
    int dir;
    int i;
    int j;
    int k;
    int c;

    for (dir = 0; dir < 3; ++dir) {
        for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
            for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
                for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                    for (c = 0; c < 3; ++c) {
                        block->v_riemann[dir][VRIDX(c, i, j, k)] =
                            0.01 * (double)(dir + 1) +
                            0.02 * (double)(c + 1) +
                            0.001 * (double)(i + 2 * j - k);
                    }
                }
            }
        }
    }
}

static void cell_dimless_dvdx(const prj_block *block, int i, int j, int k,
    double dvdx[3][3])
{
    int dir;
    int c;

    for (dir = 0; dir < 3; ++dir) {
        for (c = 0; c < 3; ++c) {
            int il = i;
            int jl = j;
            int kl = k;
            int ir = i;
            int jr = j;
            int kr = k;
            double vL;
            double vR;

            if (dir == X1DIR) {
                ir = i + 1;
            } else if (dir == X2DIR) {
                jr = j + 1;
            } else {
                kr = k + 1;
            }
            vL = block->v_riemann[dir][VRIDX(c, il, jl, kl)];
            vR = block->v_riemann[dir][VRIDX(c, ir, jr, kr)];
            dvdx[dir][c] = (vR - vL) / block->dx[dir] / PRJ_CLIGHT;
        }
    }
}

static void check_gr_m1_sources(void)
{
    prj_mesh mesh;
    prj_coord coord;
    prj_eos eos;
    prj_rad rad;
    prj_block *block;
    prj_z4c_hydro_geom geom;
    prj_rad_gr_m1_closure_ctx closure;
    double Fcov[3] = {0.07 * PRJ_CLIGHT, -0.04 * PRJ_CLIGHT, 0.02 * PRJ_CLIGHT};
    double Fcon[3];
    double Pcon[3][3];
    double dvdx[3][3];
    double vcon[3] = {0.0, 0.0, 0.0};
    double expected[4];
    double E = 2.5;
    const int i = 2;
    const int j = 2;
    const int k = 2;
    const double alpha0 = 1.17;
    const double dalpha[3] = {0.031, -0.022, 0.014};
    const double dbeta[3][3] = {
        {0.006, -0.004, 0.003},
        {-0.002, 0.005, -0.001},
        {0.004, 0.002, -0.003}
    };
    const double dgamma_diag[3][3] = {
        {0.021, -0.011, 0.008},
        {-0.014, 0.017, -0.006},
        {0.009, 0.004, -0.012}
    };
    const double K[3][3] = {
        {0.041, 0.007, -0.005},
        {0.007, -0.019, 0.006},
        {-0.005, 0.006, 0.027}
    };
    int a;
    int b;
    int d;

    init_test_mesh(&mesh, &coord);
    init_test_eos(&eos);
    init_test_rad(&rad);
    block = &mesh.blocks[0];
    fill_constant_state(block, 0.0, Fcov);
    block->W_rad[WIDX(PRJ_RAD_PRIM_E(0, 0), i, j, k)] = E;
    block->W_rad[WIDX(PRJ_RAD_PRIM_F1(0, 0), i, j, k)] = Fcov[0];
    block->W_rad[WIDX(PRJ_RAD_PRIM_F2(0, 0), i, j, k)] = Fcov[1];
    block->W_rad[WIDX(PRJ_RAD_PRIM_F3(0, 0), i, j, k)] = Fcov[2];
    set_linear_source_z4c(block, i, j, k, alpha0, dalpha, dbeta,
        dgamma_diag, K);
    fill_velocity_faces_with_gradient(block);

    if (!prj_z4c_load_hydro_geom(&mesh, block, 0, i, j, k, &geom)) {
        die("source geometry load failed");
    }
    for (a = 0; a < 3; ++a) {
        Fcon[a] = 0.0;
        for (b = 0; b < 3; ++b) {
            Fcon[a] += geom.gamma_inv[a][b] * Fcov[b];
        }
    }
    cell_dimless_dvdx(block, i, j, k, dvdx);
    make_closure_ctx(&geom, vcon, dvdx, 1, 0.0, &closure);
    prj_rad_gr_m1_pressure(&rad, &closure, E, Fcov, Pcon);
    expected[0] = 0.0;
    for (a = 0; a < 3; ++a) {
        expected[0] -= Fcon[a] * geom.dalpha[a] / geom.alpha;
        for (b = 0; b < 3; ++b) {
            expected[0] += PRJ_CLIGHT * Pcon[a][b] * geom.K_dd[a][b];
        }
    }
    expected[0] *= geom.alpha * geom.sqrt_gamma;
    for (d = 0; d < 3; ++d) {
        double mom_src = -PRJ_CLIGHT * PRJ_CLIGHT * E * geom.dalpha[d];

        for (a = 0; a < 3; ++a) {
            mom_src += PRJ_CLIGHT * Fcov[a] * geom.dbeta[d][a];
            for (b = 0; b < 3; ++b) {
                mom_src += 0.5 * geom.alpha * PRJ_CLIGHT * PRJ_CLIGHT *
                    Pcon[a][b] *
                    geom.dgamma[d][a][b];
            }
        }
        expected[1 + d] = geom.sqrt_gamma * mom_src;
    }

    prj_src_update(&eos, &rad, 0, &mesh, block, 0, block->W_mhd, block->W_rad,
        block->mhd_rhs, block->rad_rhs);
    assert_close("source E", block->rad_rhs[RADVIDX(PRJ_RAD_CONS_E(0, 0), i, j, k)],
        expected[0], 2.0e-12);
    assert_close("source F1", block->rad_rhs[RADVIDX(PRJ_RAD_CONS_F1(0, 0), i, j, k)],
        expected[1], 2.0e-12);
    assert_close("source F2", block->rad_rhs[RADVIDX(PRJ_RAD_CONS_F2(0, 0), i, j, k)],
        expected[2], 2.0e-12);
    assert_close("source F3", block->rad_rhs[RADVIDX(PRJ_RAD_CONS_F3(0, 0), i, j, k)],
        expected[3], 2.0e-12);
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

#if PRJ_DYNAMIC_GR && PRJ_USE_RADIATION_M1 && PRJ_NRAD > 0
    check_grm1_build_R_flat();
    check_grm1_build_R_shifted_metric();
    check_grm1_build_R_free_streaming();
    check_grm1_freq_drift_rest_frame();
    check_grm1_freq_drift_free_streaming();
    check_grm1_freq_drift_u_contraction();
    check_flat_zero_shift_matches_non_gr();
    check_flat_shift_terms();
    check_curved_diagonal_flux();
    check_gr_m1_characteristic_speeds();
    check_gr_pressure_closure_zero_velocity();
    check_gr_pressure_closure_boosted_fbar();
    check_gr_pressure_closure_small_velocity_build_R();
    check_gr_m1_frequency_third_moment_rest_frame();
    check_gr_m1_frequency_gravitational_redshift();
    check_gr_m1_sources();
    printf("test_gr_m1_transport: ok\n");
#else
    printf("test_gr_m1_transport: skipped\n");
#endif

#if defined(PRJ_ENABLE_MPI)
    MPI_Finalize();
#endif
    return 0;
}
