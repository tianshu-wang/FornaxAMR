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
    double scaled;
    double w;
    int idx;

    if (f <= 0.0) {
        return rad->chi[0];
    }
    if (f >= 1.0) {
        return rad->chi[NCLOSURE];
    }
    scaled = f * (double)NCLOSURE;
    idx = (int)scaled;
    w = scaled - (double)idx;
    return rad->chi[idx] + w * (rad->chi[idx + 1] - rad->chi[idx]);
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

static void make_closure_ctx(const prj_z4c_hydro_geom *geom,
    const double vcon[3], const double dvdx[3][3], int have_shear,
    double opacity, prj_rad_gr_m1_closure_ctx *ctx)
{
    int a;
    int b;
    int d;

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
    ctx->opacity = opacity;
    ctx->have_shear = have_shear;
}

static void expected_zero_velocity_gr_pressure(const prj_rad *rad,
    const prj_z4c_hydro_geom *geom, double E, const double Fcov[3],
    double P[3][3])
{
    double Fcon[3] = {0.0, 0.0, 0.0};
    double F2 = 0.0;
    double f;
    double chi;
    double thin_w;
    double thick_w;
    int a;
    int b;

    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            Fcon[a] += geom->gamma_inv[a][b] * Fcov[b];
        }
        F2 += Fcov[a] * Fcon[a];
    }
    f = F2 > 0.0 && E > 0.0 ? sqrt(F2) / (PRJ_CLIGHT * E) : 0.0;
    if (f > 1.0) {
        f = 1.0;
    }
    chi = test_m1_chi_lookup(rad, f);
    thin_w = 0.5 * (3.0 * chi - 1.0);
    thick_w = 1.5 * (1.0 - chi);
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            double thin = F2 > 0.0 ? E * Fcon[a] * Fcon[b] / F2 : 0.0;
            double thick = (E / 3.0) * geom->gamma_inv[a][b];

            P[a][b] = thin_w * thin + thick_w * thick;
        }
    }
}

static int test_solve6(double A[6][7], double x[6])
{
    int col;
    int row;

    for (col = 0; col < 6; ++col) {
        int pivot = col;
        double pivabs = fabs(A[col][col]);

        for (row = col + 1; row < 6; ++row) {
            double v = fabs(A[row][col]);

            if (v > pivabs) {
                pivabs = v;
                pivot = row;
            }
        }
        if (!isfinite(pivabs) || pivabs < 1.0e-14) {
            return 0;
        }
        if (pivot != col) {
            int j;

            for (j = col; j < 7; ++j) {
                double tmp = A[col][j];

                A[col][j] = A[pivot][j];
                A[pivot][j] = tmp;
            }
        }
        {
            double inv_piv = 1.0 / A[col][col];
            int j;

            for (j = col; j < 7; ++j) {
                A[col][j] *= inv_piv;
            }
        }
        for (row = 0; row < 6; ++row) {
            double factor;
            int j;

            if (row == col) {
                continue;
            }
            factor = A[row][col];
            if (factor == 0.0) {
                continue;
            }
            for (j = col; j < 7; ++j) {
                A[row][j] -= factor * A[col][j];
            }
        }
    }
    for (row = 0; row < 6; ++row) {
        x[row] = A[row][6];
        if (!isfinite(x[row])) {
            return 0;
        }
    }
    return 1;
}

static void expected_flat_boosted_pressure_for_fbar(const prj_rad *rad,
    double E, const double Fcov_in[3], const double vcon_in[3], double fbar,
    double P[3][3])
{
    static const int ia[6] = {0, 1, 2, 0, 0, 1};
    static const int ib[6] = {0, 1, 2, 1, 2, 2};
    double Fcov[3];
    double Fhat[3];
    double Pthin[3][3];
    double A0[3][3];
    double vcon[3];
    double u_cov[3];
    double h_mix[3][3];
    double H0[3];
    double Hcoef[3][6];
    double Jcoef[6];
    double amat[6][7];
    double sol[6];
    double F2 = 0.0;
    double Fmag;
    double cE;
    double beta2 = 0.0;
    double wlor;
    double FdotV = 0.0;
    double Fdotu;
    double J0;
    double Q0;
    double chi;
    double thin_w;
    double thick_w;
    int a;
    int b;
    int m;

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
    for (a = 0; a < 3; ++a) {
        u_cov[a] = wlor * vcon[a];
        Fhat[a] = Fcov[a] / PRJ_CLIGHT;
        FdotV += Fhat[a] * vcon[a];
    }
    Fdotu = wlor * FdotV;
    J0 = E * wlor * wlor - 2.0 * wlor * Fdotu;
    Q0 = E * wlor - Fdotu;
    for (a = 0; a < 3; ++a) {
        double h_n = -wlor * wlor * vcon[a];
        double h_F = Fhat[a] + wlor * wlor * vcon[a] * FdotV;

        H0[a] = Q0 * h_n + wlor * h_F;
        for (b = 0; b < 3; ++b) {
            h_mix[a][b] = (a == b ? 1.0 : 0.0) +
                wlor * wlor * vcon[a] * vcon[b];
            A0[a][b] = ((a == b ? 1.0 : 0.0) +
                4.0 * vcon[a] * vcon[b]) / 3.0;
        }
    }
    for (m = 0; m < 6; ++m) {
        int p = ia[m];
        int q = ib[m];

        Jcoef[m] = u_cov[p] * u_cov[q] * (p == q ? 1.0 : 2.0);
        for (a = 0; a < 3; ++a) {
            Hcoef[a][m] = -h_mix[a][p] * u_cov[q];
            if (p != q) {
                Hcoef[a][m] -= h_mix[a][q] * u_cov[p];
            }
        }
    }
    chi = test_m1_chi_lookup(rad, fbar);
    thin_w = 0.5 * (3.0 * chi - 1.0);
    thick_w = 1.5 * (1.0 - chi);
    memset(amat, 0, sizeof(amat));
    for (m = 0; m < 6; ++m) {
        int i = ia[m];
        int j = ib[m];
        int n;

        amat[m][6] = thin_w * Pthin[i][j] + thick_w *
            (A0[i][j] * J0 + H0[i] * vcon[j] + H0[j] * vcon[i]);
        for (n = 0; n < 6; ++n) {
            double coef = thick_w * (A0[i][j] * Jcoef[n] +
                Hcoef[i][n] * vcon[j] + Hcoef[j][n] * vcon[i]);

            amat[m][n] = (m == n ? 1.0 : 0.0) - coef;
        }
    }
    if (!test_solve6(amat, sol)) {
        die("boosted pressure solve failed");
    }
    memset(P, 0, 9 * sizeof(double));
    for (m = 0; m < 6; ++m) {
        int i = ia[m];
        int j = ib[m];

        P[i][j] = sol[m];
        P[j][i] = sol[m];
    }
}

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
    double Fcov[3] = {0.12 * PRJ_CLIGHT, -0.04 * PRJ_CLIGHT, 0.03 * PRJ_CLIGHT};
    double got[4];
    double expected[PRJ_NVAR_CONS] = {0.0};
    double Wface[PRJ_NVAR_PRIM] = {0.0};
    double chi_face[PRJ_NRAD * PRJ_NEGROUP] = {0.0};

    init_test_mesh(&mesh, &coord);
    init_test_eos(&eos);
    init_test_rad(&rad);
    set_uniform_z4c(&mesh.blocks[0], 1.0, beta, gamma_diag);
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
    double Fcov[3] = {0.08 * PRJ_CLIGHT, 0.05 * PRJ_CLIGHT, -0.02 * PRJ_CLIGHT};
    double got[4];
    double base[PRJ_NVAR_CONS] = {0.0};
    double expected[4];
    double Wface[PRJ_NVAR_PRIM] = {0.0};
    double chi_face[PRJ_NRAD * PRJ_NEGROUP] = {0.0};

    init_test_mesh(&mesh, &coord);
    init_test_eos(&eos);
    init_test_rad(&rad);
    set_uniform_z4c(&mesh.blocks[0], 1.0, beta, gamma_diag);
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
    expected_flat_boosted_pressure_for_fbar(&rad, E, Fcov, vcon, fbar,
        expected);
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            assert_close("boosted GR pressure closure", got[a][b],
                expected[a][b], 2.0e-11);
        }
    }
    prj_mesh_destroy(&mesh);
}

static void check_gr_pressure_closure_small_velocity_uses_eulerian_fbar(void)
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
    double feuler;
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
    feuler = sqrt(Fcov[0] * Fcov[0] + Fcov[1] * Fcov[1] +
        Fcov[2] * Fcov[2]) / (PRJ_CLIGHT * E);
    expected_flat_boosted_pressure_for_fbar(&rad, E, Fcov, vcon, feuler,
        expected);
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            assert_close("small-velocity GR pressure closure", got[a][b],
                expected[a][b], 2.0e-12);
        }
    }
    prj_mesh_destroy(&mesh);
}

static void expected_rest_frame_gr_frequency_terms(const prj_rad *rad,
    double E, const double Fcov[3], const double grad[3][3],
    double *energy_drift, double momentum_drift[3])
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
    chi = test_m1_chi_lookup(rad, fbar);
    thin_w = 0.5 * (3.0 * chi - 1.0);
    thick_w = 1.5 * (1.0 - chi);
    for (a = 0; a < 3; ++a) {
        n[a] = Hmag > 0.0 ? H[a] / Hmag : 0.0;
    }
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            double Pthin = E * n[a] * n[b];
            double Pthick = (a == b ? E / 3.0 : 0.0);
            double P = thin_w * Pthin + thick_w * Pthick;

            *energy_drift += P * grad[a][b];
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
}

static void check_gr_m1_frequency_third_moment_rest_frame_case(
    double grad_scale, int expect_upper_donor)
{
    prj_mesh mesh;
    prj_coord coord;
    prj_rad rad;
    prj_block *block;
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
        expected_rest_frame_gr_frequency_terms(&rad, E, Fcov, grad,
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
        i, j, k, dt);
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
            expected[0] += Pcon[a][b] * geom.K_dd[a][b];
        }
    }
    expected[0] *= geom.alpha * geom.sqrt_gamma;
    for (d = 0; d < 3; ++d) {
        double mom_src = -E * geom.dalpha[d];

        for (a = 0; a < 3; ++a) {
            mom_src += Fcov[a] * geom.dbeta[d][a];
            for (b = 0; b < 3; ++b) {
                mom_src += 0.5 * geom.alpha * Pcon[a][b] *
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
    check_flat_zero_shift_matches_non_gr();
    check_flat_shift_terms();
    check_curved_diagonal_flux();
    check_gr_pressure_closure_zero_velocity();
    check_gr_pressure_closure_boosted_fbar();
    check_gr_pressure_closure_small_velocity_uses_eulerian_fbar();
    check_gr_m1_frequency_third_moment_rest_frame();
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
