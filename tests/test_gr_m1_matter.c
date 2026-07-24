#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(PRJ_ENABLE_MPI)
#include <mpi.h>
#endif

#include "prj.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#if PRJ_DYNAMIC_GR && PRJ_USE_RADIATION_M1 && PRJ_NRAD > 0
int prj_rad_gr_m1_residual_test_wrapper(const prj_rad *rad, prj_eos *eos,
    const prj_z4c_hydro_geom *geom, const double *u_old, const double *P,
    double dt, double *resid, double *u_new_out);
int prj_rad_gr_m1_jacobian_test_wrapper(const prj_rad *rad, prj_eos *eos,
    const prj_z4c_hydro_geom *geom, const double *u_old, const double *P,
    double dt, double *resid, double *jac, double *u_new_out);

#define TEST_GR_M1_RESIDUAL_NP (6 + 4 * PRJ_NRAD * PRJ_NEGROUP)

static void die(const char *msg)
{
    fprintf(stderr, "test_gr_m1_matter: %s\n", msg);
    exit(1);
}

static void assert_close(const char *name, double got, double expected,
    double rel)
{
    double scale = fmax(1.0, fmax(fabs(got), fabs(expected)));
    double tol = rel * scale;

    if (!isfinite(got) || !isfinite(expected) || fabs(got - expected) > tol) {
        fprintf(stderr,
            "test_gr_m1_matter: %s got %.17e expected %.17e tol %.3e\n",
            name, got, expected, tol);
        exit(1);
    }
}

static double test_m1_chi_exact(double f)
{
    if (f <= 0.0) {
        return 1.0 / 3.0;
    }
    if (f >= 1.0) {
        return 1.0;
    }
    return (3.0 + 4.0 * f * f) /
        (5.0 + 2.0 * sqrt(4.0 - 3.0 * f * f));
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

static void init_source_test_rad_full(prj_rad *rad, double kappa_value,
    double sigma_value, double delta_value, double eta_value, double xe_value)
{
    const size_t ncorners = 2u * 2u * 2u * (size_t)PRJ_NEGROUP;
    const double tiny = 1.0e-300;
    const double eta_factor = 4.0 * M_PI / RAD_SCALE;
    int nu;
    int g;

    init_test_rad(rad);
    rad->min_inel_density = 1.0e99;
    rad->maxiter = 50;
    rad->implicit_err_tol = 1.0e-10;
    rad->nromax = 2;
    rad->ntmax = 2;
    rad->nyemax = 2;
    rad->romin = 0.25;
    rad->romax = 4.0;
    rad->tmin = 0.25;
    rad->tmax = 4.0;
    rad->yemin = 0.0;
    rad->yemax = 1.0;
    rad->log_romin = log(rad->romin);
    rad->log_romax = log(rad->romax);
    rad->log_tmin = log(rad->tmin);
    rad->log_tmax = log(rad->tmax);
    rad->inv_logrho_span = 1.0 / (rad->log_romax - rad->log_romin);
    rad->inv_logtemp_span = 1.0 / (rad->log_tmax - rad->log_tmin);
    rad->inv_ye_span = 1.0 / (rad->yemax - rad->yemin);

    for (nu = 0; nu < PRJ_NRAD; ++nu) {
        rad->egroup[nu] = (double *)calloc((size_t)PRJ_NEGROUP, sizeof(double));
        rad->eedge[nu] =
            (double *)calloc((size_t)PRJ_NEGROUP + 1u, sizeof(double));
        rad->egroup_erg[nu] =
            (double *)calloc((size_t)PRJ_NEGROUP, sizeof(double));
        rad->degroup_erg[nu] =
            (double *)calloc((size_t)PRJ_NEGROUP, sizeof(double));
        rad->x_e[nu] = (double *)calloc((size_t)PRJ_NEGROUP, sizeof(double));
        rad->log_egroup[nu] =
            (double *)calloc((size_t)PRJ_NEGROUP, sizeof(double));
        rad->spec_factor[nu] =
            (double *)calloc((size_t)PRJ_NEGROUP, sizeof(double));
        rad->absopac[nu] =
            (prj_table_real *)calloc(ncorners, sizeof(prj_table_real));
        rad->scaopac[nu] =
            (prj_table_real *)calloc(ncorners, sizeof(prj_table_real));
        rad->emis[nu] =
            (prj_table_real *)calloc(ncorners, sizeof(prj_table_real));
        rad->sdelta[nu] =
            (prj_table_real *)calloc(ncorners, sizeof(prj_table_real));
        if (rad->egroup[nu] == 0 || rad->eedge[nu] == 0 ||
            rad->egroup_erg[nu] == 0 || rad->degroup_erg[nu] == 0 ||
            rad->x_e[nu] == 0 || rad->log_egroup[nu] == 0 ||
            rad->spec_factor[nu] == 0 || rad->absopac[nu] == 0 ||
            rad->scaopac[nu] == 0 || rad->emis[nu] == 0 ||
            rad->sdelta[nu] == 0) {
            die("source opacity allocation failed");
        }
        for (g = 0; g < PRJ_NEGROUP; ++g) {
            double e = 1.0 + (double)g;

            rad->egroup[nu][g] = e;
            rad->eedge[nu][g] = e;
            rad->egroup_erg[nu][g] = e * PRJ_MEV_TO_ERG;
            rad->degroup_erg[nu][g] = PRJ_MEV_TO_ERG;
            rad->log_egroup[nu][g] = log(e);
            rad->spec_factor[nu][g] = 1.0;
            rad->x_e[nu][g] = xe_value;
        }
        rad->eedge[nu][PRJ_NEGROUP] = 1.0 + (double)PRJ_NEGROUP;
        for (g = 0; g < (int)ncorners; ++g) {
            rad->absopac[nu][g] =
                (prj_table_real)log(kappa_value > 0.0 ? kappa_value : tiny);
            rad->scaopac[nu][g] =
                (prj_table_real)log(sigma_value > 0.0 ? sigma_value : tiny);
            rad->emis[nu][g] = (prj_table_real)log(eta_value > 0.0 ?
                eta_value / eta_factor : tiny);
            rad->sdelta[nu][g] = (prj_table_real)delta_value;
        }
    }
}

static void init_source_test_rad(prj_rad *rad, double sigma_value)
{
    init_source_test_rad_full(rad, 1.0e-300, sigma_value, 0.0, 0.0, 0.0);
}

static size_t test_opac_cell_idx(int ir, int it, int iy, int group)
{
    return (((((size_t)ir * 2u + (size_t)it) * 2u + (size_t)iy) *
        (size_t)PRJ_NEGROUP) + (size_t)group);
}

static void make_source_test_rad_sloped(prj_rad *rad)
{
    const double eta_factor = 4.0 * M_PI / RAD_SCALE;
    int nu;
    int g;
    int ir;
    int it;
    int iy;

    init_source_test_rad_full(rad, 2.0e-37, 1.4e-37, 0.18, 3.0e-27,
        0.08 * RAD_SCALE);
    for (nu = 0; nu < PRJ_NRAD; ++nu) {
        for (g = 0; g < PRJ_NEGROUP; ++g) {
            double group_shift = 0.01 * (double)g + 0.02 * (double)nu;

            for (ir = 0; ir < 2; ++ir) {
                for (it = 0; it < 2; ++it) {
                    for (iy = 0; iy < 2; ++iy) {
                        size_t idx = test_opac_cell_idx(ir, it, iy, g);
                        double r = (double)ir;
                        double t = (double)it;
                        double y = (double)iy;

                        rad->absopac[nu][idx] = (prj_table_real)
                            (log(2.0e-37) + group_shift +
                             0.07 * r - 0.04 * t + 0.03 * y);
                        rad->scaopac[nu][idx] = (prj_table_real)
                            (log(1.4e-37) - 0.5 * group_shift -
                             0.03 * r + 0.05 * t - 0.02 * y);
                        rad->emis[nu][idx] = (prj_table_real)
                            (log(3.0e-27 / eta_factor) +
                             0.3 * group_shift + 0.02 * r +
                             0.06 * t - 0.04 * y);
                        rad->sdelta[nu][idx] = (prj_table_real)
                            (0.18 + 0.01 * group_shift + 0.02 * r -
                             0.015 * t + 0.01 * y);
                    }
                }
            }
        }
    }
}

static void init_test_eos(prj_eos *eos)
{
    memset(eos, 0, sizeof(*eos));
    eos->kind = PRJ_EOS_KIND_IDEAL;
    prj_eos_init(eos, 0);
}

static void set_diag_geom(prj_z4c_hydro_geom *geom, double alpha,
    const double beta[3], const double gamma_diag[3])
{
    int d;

    memset(geom, 0, sizeof(*geom));
    geom->alpha = alpha;
    geom->sqrt_gamma = sqrt(gamma_diag[0] * gamma_diag[1] * gamma_diag[2]);
    for (d = 0; d < 3; ++d) {
        geom->beta[d] = beta[d];
        geom->gamma[d][d] = gamma_diag[d];
        geom->gamma_inv[d][d] = 1.0 / gamma_diag[d];
    }
}

static void set_flat_prim(double *W, double rho, const double v[3],
    double eint, double ye, double E, const double Fcov[3])
{
    int field;
    int group;

    for (int n = 0; n < PRJ_NVAR_PRIM; ++n) {
        W[n] = 0.0;
    }
    W[PRJ_PRIM_RHO] = rho;
    W[PRJ_PRIM_V1] = v[0];
    W[PRJ_PRIM_V2] = v[1];
    W[PRJ_PRIM_V3] = v[2];
    W[PRJ_PRIM_EINT] = eint;
    W[PRJ_PRIM_YE] = ye;
#if PRJ_MHD
    W[PRJ_PRIM_B1] = 0.0;
    W[PRJ_PRIM_B2] = 0.0;
    W[PRJ_PRIM_B3] = 0.0;
#endif
    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            W[PRJ_PRIM_RAD_E(field, group)] = 0.0;
            W[PRJ_PRIM_RAD_F1(field, group)] = 0.0;
            W[PRJ_PRIM_RAD_F2(field, group)] = 0.0;
            W[PRJ_PRIM_RAD_F3(field, group)] = 0.0;
        }
    }
    W[PRJ_PRIM_RAD_E(0, 0)] = E;
    W[PRJ_PRIM_RAD_F1(0, 0)] = Fcov[0];
    W[PRJ_PRIM_RAD_F2(0, 0)] = Fcov[1];
    W[PRJ_PRIM_RAD_F3(0, 0)] = Fcov[2];
}

static void prim2cons_or_die(prj_eos *eos, const prj_z4c_hydro_geom *geom,
    const double *W, double *U)
{
    prj_eos_gr_geom eos_geom;
    int a;
    int b;

    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            eos_geom.gamma[a][b] = geom->gamma[a][b];
        }
    }
    if (prj_eos_gr_prim2cons(eos, &eos_geom, W, U, PRJ_EOS_CTX_MAIN) !=
        PRJ_EOS_GR_OK) {
        die("GR prim2cons failed");
    }
}

static void set_residual_p_from_flat_state(const double *u, double temperature,
    double P[TEST_GR_M1_RESIDUAL_NP])
{
    double rho = u[PRJ_CONS_RHO];
    int field;
    int group;

    P[0] = rho;
    P[1] = u[PRJ_CONS_MOM1] / rho;
    P[2] = u[PRJ_CONS_MOM2] / rho;
    P[3] = u[PRJ_CONS_MOM3] / rho;
    P[4] = temperature;
    P[5] = u[PRJ_CONS_YE] / rho;
    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            int idx = field * PRJ_NEGROUP + group;
            int pidx = 6 + 4 * idx;
            double E = u[PRJ_CONS_RAD_E(field, group)];
            double Fhat[3];
            double Fhat2;

            Fhat[0] = u[PRJ_CONS_RAD_F1(field, group)] / PRJ_CLIGHT;
            Fhat[1] = u[PRJ_CONS_RAD_F2(field, group)] / PRJ_CLIGHT;
            Fhat[2] = u[PRJ_CONS_RAD_F3(field, group)] / PRJ_CLIGHT;
            Fhat2 = Fhat[0] * Fhat[0] + Fhat[1] * Fhat[1] +
                Fhat[2] * Fhat[2];
            if (E > 0.0 && Fhat2 > 0.0) {
                double disc = 9.0 * E * E - 27.0 * Fhat2 / 4.0;
                double y;
                double qmag;
                double inv_Fhat;
                int d;

                if (disc < 0.0 && disc > -1.0e-12 * E * E) {
                    disc = 0.0;
                }
                if (disc < 0.0) {
                    die("flat residual P has superluminal flux");
                }
                y = 0.5 * (3.0 * E - sqrt(disc));
                P[pidx] = E - 4.0 * y / 3.0;
                if (P[pidx] < 0.0 && P[pidx] > -1.0e-12 * E) {
                    P[pidx] = 0.0;
                }
                qmag = sqrt(y);
                inv_Fhat = 1.0 / sqrt(Fhat2);
                for (d = 0; d < 3; ++d) {
                    P[pidx + 1 + d] = qmag * Fhat[d] * inv_Fhat;
                }
            } else {
                P[pidx] = E;
                P[pidx + 1] = 0.0;
                P[pidx + 2] = 0.0;
                P[pidx + 3] = 0.0;
            }
        }
    }
}

static void check_residual_zero(const char *name, const double *resid,
    double tol)
{
    int v;

    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        if (!isfinite(resid[v]) || fabs(resid[v]) > tol) {
            fprintf(stderr,
                "test_gr_m1_matter: %s residual[%d] = %.17e tol %.3e\n",
                name, v, resid[v], tol);
            exit(1);
        }
    }
}

static void set_gr_m1_jacobian_test_p(
    double P[TEST_GR_M1_RESIDUAL_NP])
{
    int field;
    int group;

    P[0] = 1.17;
    P[1] = 2.0e8;
    P[2] = -1.1e8;
    P[3] = 0.7e8;
    P[4] = 1.21;
    P[5] = 0.37;
    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            int idx = field * PRJ_NEGROUP + group;
            int pidx = 6 + 4 * idx;
            double ER = 0.035 + 0.002 * (double)(idx + 1);
            double qmag = sqrt(ER);

            P[pidx] = ER;
            P[pidx + 1] = 0.018 * qmag;
            P[pidx + 2] = -0.011 * qmag;
            P[pidx + 3] = 0.007 * qmag;
        }
    }
}

static void check_gr_m1_residual_jacobian_fd(void)
{
    prj_eos eos;
    prj_rad rad;
    double W[PRJ_NVAR_PRIM];
    double u_old[PRJ_NVAR_CONS];
    double u_ref[PRJ_NVAR_CONS];
    double P[TEST_GR_M1_RESIDUAL_NP];
    double Pp[TEST_GR_M1_RESIDUAL_NP];
    double Pm[TEST_GR_M1_RESIDUAL_NP];
    double resid[PRJ_NVAR_CONS];
    double resp[PRJ_NVAR_CONS];
    double resm[PRJ_NVAR_CONS];
    double jac[PRJ_NVAR_CONS * TEST_GR_M1_RESIDUAL_NP];
    double u_new[PRJ_NVAR_CONS];
    prj_z4c_hydro_geom geom;
    double beta[3] = {0.015, -0.01, 0.006};
    double gamma_diag[3] = {1.08, 0.94, 1.13};
    double vzero[3] = {0.0, 0.0, 0.0};
    double Fcov[3] = {0.0, 0.0, 0.0};
    double dt = 0.13;
    int col;
    int row;
    int v;

    init_test_eos(&eos);
    make_source_test_rad_sloped(&rad);
    set_diag_geom(&geom, 0.94, beta, gamma_diag);
    set_flat_prim(W, 1.0, vzero, 1.5, 0.2, 0.02, Fcov);
    prim2cons_or_die(&eos, &geom, W, u_old);
    set_gr_m1_jacobian_test_p(P);
    if (!prj_rad_gr_m1_residual_test_wrapper(&rad, &eos, &geom, u_old, P,
            dt, resid, u_ref)) {
        die("jacobian reference residual failed");
    }
    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        u_old[v] = u_ref[v];
    }
    if (!prj_rad_gr_m1_jacobian_test_wrapper(&rad, &eos, &geom, u_old, P,
            dt, resid, jac, u_new)) {
        die("analytic jacobian failed");
    }

    for (col = 0; col < TEST_GR_M1_RESIDUAL_NP; ++col) {
        double h = 1.0e-6 * fmax(1.0, fabs(P[col]));

        for (v = 0; v < TEST_GR_M1_RESIDUAL_NP; ++v) {
            Pp[v] = P[v];
            Pm[v] = P[v];
        }
        if (col >= 6 && ((col - 6) % 4) == 0) {
            h = fmin(h, 0.25 * P[col]);
        }
        Pp[col] += h;
        Pm[col] -= h;
        if (!prj_rad_gr_m1_residual_test_wrapper(&rad, &eos, &geom, u_old,
                Pp, dt, resp, 0) ||
            !prj_rad_gr_m1_residual_test_wrapper(&rad, &eos, &geom, u_old,
                Pm, dt, resm, 0)) {
            die("finite-difference jacobian residual failed");
        }
        for (row = 0; row < PRJ_NVAR_CONS; ++row) {
            double fd = (resp[row] - resm[row]) / (2.0 * h);
            double got = jac[row * TEST_GR_M1_RESIDUAL_NP + col];
            double scale = fmax(1.0, fmax(fabs(fd), fabs(got)));
            double tol = 5.0e-4 * scale;

            if (!isfinite(fd) || !isfinite(got) || fabs(fd - got) > tol) {
                fprintf(stderr,
                    "test_gr_m1_matter: GR residual jac[%d,%d] got %.17e "
                    "expected %.17e tol %.3e\n",
                    row, col, got, fd, tol);
                exit(1);
            }
        }
    }

    prj_rad3_opac_free(&rad);
}

static void check_gr_m1_residual_rest_energy_matches_non_gr(void)
{
    prj_eos eos;
    prj_rad rad;
    double W[PRJ_NVAR_PRIM];
    double u_old[PRJ_NVAR_CONS];
    double u_non_gr[PRJ_NVAR_CONS];
    double resid[PRJ_NVAR_CONS];
    double u_new[PRJ_NVAR_CONS];
    double P[TEST_GR_M1_RESIDUAL_NP];
    double kappa[PRJ_NRAD * PRJ_NEGROUP];
    prj_z4c_hydro_geom geom;
    double beta[3] = {0.0, 0.0, 0.0};
    double gamma_diag[3] = {1.0, 1.0, 1.0};
    double vzero[3] = {0.0, 0.0, 0.0};
    double Fcov[3] = {0.0, 0.0, 0.0};
    double T = 0.0;
    double dt = 0.17;
    int ok;
    int v;

    init_test_eos(&eos);
    init_source_test_rad_full(&rad, 0.14 / (PRJ_CLIGHT * RAD_SCALE),
        0.0, 0.0, 0.03 / RAD_SCALE, 0.25 * RAD_SCALE);
    set_diag_geom(&geom, 1.0, beta, gamma_diag);
    set_flat_prim(W, 1.0, vzero, 1.5, 0.2, 0.0, Fcov);
    prim2cons_or_die(&eos, &geom, W, u_old);
    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        u_non_gr[v] = u_old[v];
    }
    prj_rad_energy_update(&rad, &eos, u_non_gr, dt, 1.0, &T, kappa);
    set_residual_p_from_flat_state(u_non_gr, T, P);
    ok = prj_rad_gr_m1_residual_test_wrapper(&rad, &eos, &geom, u_old, P,
        dt, resid, u_new);
    if (!ok) {
        die("rest residual failed");
    }

    check_residual_zero("rest energy non-GR match", resid, 2.0e-10);
    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        assert_close("rest residual u_new", u_new[v], u_non_gr[v], 2.0e-12);
    }
    prj_rad3_opac_free(&rad);
}

static void check_gr_m1_residual_invalid_states(void)
{
    prj_eos eos;
    prj_rad rad;
    double W[PRJ_NVAR_PRIM];
    double u_old[PRJ_NVAR_CONS];
    double resid[PRJ_NVAR_CONS];
    double P[TEST_GR_M1_RESIDUAL_NP];
    prj_z4c_hydro_geom geom;
    double beta[3] = {0.0, 0.0, 0.0};
    double gamma_diag[3] = {1.0, 1.0, 1.0};
    double vzero[3] = {0.0, 0.0, 0.0};
    double Fcov[3] = {0.0, 0.0, 0.0};
    prj_z4c_hydro_geom bad_geom;

    init_test_eos(&eos);
    init_source_test_rad(&rad, 0.0);
    set_diag_geom(&geom, 1.0, beta, gamma_diag);
    set_flat_prim(W, 1.0, vzero, 1.5, 0.2, 0.8, Fcov);
    prim2cons_or_die(&eos, &geom, W, u_old);
    set_residual_p_from_flat_state(u_old, 1.0, P);
    if (!prj_rad_gr_m1_residual_test_wrapper(&rad, &eos, &geom, u_old, P,
            0.1, resid, 0)) {
        die("valid residual rejected");
    }

    P[0] = 0.0;
    if (prj_rad_gr_m1_residual_test_wrapper(&rad, &eos, &geom, u_old, P,
            0.1, resid, 0)) {
        die("invalid rho accepted");
    }
    set_residual_p_from_flat_state(u_old, 1.0, P);
    P[4] = 0.0;
    if (prj_rad_gr_m1_residual_test_wrapper(&rad, &eos, &geom, u_old, P,
            0.1, resid, 0)) {
        die("invalid temperature accepted");
    }
    set_residual_p_from_flat_state(u_old, 1.0, P);
    P[1] = 1.01 * PRJ_CLIGHT;
    if (prj_rad_gr_m1_residual_test_wrapper(&rad, &eos, &geom, u_old, P,
            0.1, resid, 0)) {
        die("superluminal velocity accepted");
    }
    set_residual_p_from_flat_state(u_old, 1.0, P);
    P[6] = -1.0;
    if (prj_rad_gr_m1_residual_test_wrapper(&rad, &eos, &geom, u_old, P,
            0.1, resid, 0)) {
        die("negative ER accepted");
    }
    set_residual_p_from_flat_state(u_old, 1.0, P);
    bad_geom = geom;
    bad_geom.alpha = 0.1;
    bad_geom.beta[0] = 1.0;
    if (prj_rad_gr_m1_residual_test_wrapper(&rad, &eos, &bad_geom, u_old, P,
            0.1, resid, 0)) {
        die("no future q0 accepted");
    }
    prj_rad3_opac_free(&rad);
}

static void check_gr_m1_matter_source_rest_momentum(void)
{
    prj_eos eos;
    prj_rad rad;
    double W[PRJ_NVAR_PRIM];
    double u_old[PRJ_NVAR_CONS];
    double u_non_gr[PRJ_NVAR_CONS];
    double resid[PRJ_NVAR_CONS];
    double u_new[PRJ_NVAR_CONS];
    double P[TEST_GR_M1_RESIDUAL_NP];
    double kappa[PRJ_NRAD * PRJ_NEGROUP];
    prj_z4c_hydro_geom geom;
    double beta[3] = {0.0, 0.0, 0.0};
    double gamma_diag[3] = {1.0, 1.0, 1.0};
    double vzero[3] = {0.0, 0.0, 0.0};
    double flux_factor = 1.0e-13;
    double E = 2.0 * flux_factor;
    double Fcov[3] = {flux_factor * PRJ_CLIGHT * 2.0, 0.0, 0.0};
    double sigma_table = 2.0 / PRJ_CLIGHT;
    double delta_value = 0.6;
    double sigma_eff = sigma_table * (1.0 - delta_value / 3.0);
    double eos_q[PRJ_EOS_NQUANT];
    double T;
    double dt = 0.25;
    int ok;
    int v;

    init_test_eos(&eos);
    init_source_test_rad_full(&rad, 0.0, sigma_table, delta_value, 0.0, 0.0);
    set_diag_geom(&geom, 1.0, beta, gamma_diag);
    set_flat_prim(W, 1.0, vzero, 1.5, 0.2, E, Fcov);
    prim2cons_or_die(&eos, &geom, W, u_old);
    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        u_non_gr[v] = u_old[v];
    }
    for (v = 0; v < PRJ_NRAD * PRJ_NEGROUP; ++v) {
        kappa[v] = 0.0;
    }
    prj_rad_momentum_update(&rad, &eos, u_non_gr, dt, 1.0, 1.0, kappa);
    {
        double rho = u_non_gr[PRJ_CONS_RHO];
        double kinetic = 0.5 *
            (u_non_gr[PRJ_CONS_MOM1] * u_non_gr[PRJ_CONS_MOM1] +
             u_non_gr[PRJ_CONS_MOM2] * u_non_gr[PRJ_CONS_MOM2] +
             u_non_gr[PRJ_CONS_MOM3] * u_non_gr[PRJ_CONS_MOM3]) / rho;
        double eint = (u_non_gr[PRJ_CONS_ETOT] - kinetic) / rho;
        double Ye = u_non_gr[PRJ_CONS_YE] / rho;

        prj_eos_rey(&eos, rho, eint, Ye, eos_q, PRJ_EOS_CTX_MAIN);
        T = eos_q[PRJ_EOS_TEMPERATURE];
    }
    set_residual_p_from_flat_state(u_non_gr, T, P);
    ok = prj_rad_gr_m1_residual_test_wrapper(&rad, &eos, &geom, u_old, P,
        dt, resid, u_new);
    if (!ok) {
        die("momentum residual failed");
    }

    assert_close("momentum residual effective sigma",
        u_non_gr[PRJ_CONS_RAD_F1(0, 0)],
        Fcov[0] / (1.0 + dt * PRJ_CLIGHT * sigma_eff), 5.0e-4);
    assert_close("momentum residual gas MOM1", resid[PRJ_CONS_MOM1],
        0.0, 2.0e-6);
    assert_close("momentum residual gas MOM2", resid[PRJ_CONS_MOM2],
        0.0, 2.0e-6);
    assert_close("momentum residual gas MOM3", resid[PRJ_CONS_MOM3],
        0.0, 2.0e-6);
    assert_close("momentum residual radiation F1",
        resid[PRJ_CONS_RAD_F1(0, 0)], 0.0, 2.0e-6);
    assert_close("momentum residual radiation F2",
        resid[PRJ_CONS_RAD_F2(0, 0)], 0.0, 2.0e-6);
    assert_close("momentum residual radiation F3",
        resid[PRJ_CONS_RAD_F3(0, 0)], 0.0, 2.0e-6);
    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        assert_close("momentum residual u_new", u_new[v], u_non_gr[v],
            2.0e-6);
    }
    prj_rad3_opac_free(&rad);
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
    check_gr_m1_residual_rest_energy_matches_non_gr();
    check_gr_m1_residual_jacobian_fd();
    check_gr_m1_residual_invalid_states();
    check_gr_m1_matter_source_rest_momentum();
    printf("test_gr_m1_matter: ok\n");
#else
    printf("test_gr_m1_matter: skipped\n");
#endif

#if defined(PRJ_ENABLE_MPI)
    MPI_Finalize();
#endif
    return 0;
}
