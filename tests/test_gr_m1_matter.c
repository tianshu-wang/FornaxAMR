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

#define TEST_RAD_TABLE_PARAM_FILE "../opacbin.extendT.param"
#define TEST_RAD_TABLE_FILE       "../opacity.SFHo.juo.horo.brem1.extendedT.bin"
#define TEST_EOS_FILE \
    "../eos_tmp/SFHoEOS__ye__0.035_0.56_50__logT_-4.793_2.176_500__logrho_-8.699_15.5_500_extend.dat"

#if PRJ_DYNAMIC_GR && PRJ_USE_RADIATION_M1 && PRJ_NRAD > 0
int prj_rad_gr_m1_residual_test_wrapper(const prj_rad *rad, prj_eos *eos,
    const prj_z4c_hydro_geom *geom, const double *u_old, const double *P,
    double dt, double *resid, double *u_new_out);
int prj_rad_gr_m1_jacobian_test_wrapper(const prj_rad *rad, prj_eos *eos,
    const prj_z4c_hydro_geom *geom, const double *u_old, const double *P,
    double dt, double *resid, double *jac, double *u_new_out);
int prj_rad_gr_m1_implicit_solve_test_wrapper(const prj_rad *rad,
    prj_eos *eos, const prj_z4c_hydro_geom *geom, const double *u_old,
    double dt, double *P, double *resid_out, double *u_new_out);

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

static void assert_abs_close(const char *name, double got, double expected,
    double tol)
{
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

static void init_real_test_rad(prj_rad *rad)
{
    const double emin_list[3] = {1.0, 1.0, 1.0};
    const double emax_list[3] = {300.0, 100.0, 100.0};
    int nu;

    if (PRJ_NRAD > 3) {
        die("real-table Jacobian test expects PRJ_NRAD <= 3");
    }

    memset(rad, 0, sizeof(*rad));
    rad->maxiter = 50;
    rad->implicit_err_tol = 1.0e-10;
    for (nu = 0; nu < PRJ_NRAD; ++nu) {
        rad->emin[nu] = emin_list[nu];
        rad->emax[nu] = emax_list[nu];
    }
    strncpy(rad->table_param_file, TEST_RAD_TABLE_PARAM_FILE,
        sizeof(rad->table_param_file) - 1);
    strncpy(rad->table_file, TEST_RAD_TABLE_FILE,
        sizeof(rad->table_file) - 1);
    prj_rad_init(rad);
}

static void init_real_test_eos(prj_eos *eos)
{
    memset(eos, 0, sizeof(*eos));
    eos->kind = PRJ_EOS_KIND_TABLE;
    strncpy(eos->filename, TEST_EOS_FILE, sizeof(eos->filename) - 1);
    prj_eos_init(eos, 0);
    if (eos->table_loaded != 1) {
        die("EOS table failed to load");
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

static void init_matter_test_mesh(prj_mesh *mesh, prj_coord *coord)
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
        die("matter mesh init failed");
    }
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
            for (k = -PRJ_NGHOST_Z4C;
                 k < PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C; ++k) {
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
                double sqrt_er;

                y = 0.5 * (3.0 * E - sqrt(disc));
                P[pidx] = E - 4.0 * y / 3.0;
                if (P[pidx] < 0.0 && P[pidx] > -1.0e-12 * E) {
                    P[pidx] = 0.0;
                }
                if (P[pidx] <= 0.0) {
                    die("flat residual P has non-positive ER");
                }
                /* qR^i = sqrt(y) Fhat/|Fhat|; store ur^i = qR^i / sqrt(ER). */
                qmag = sqrt(y);
                sqrt_er = sqrt(P[pidx]);
                inv_Fhat = 1.0 / sqrt(Fhat2);
                for (d = 0; d < 3; ++d) {
                    P[pidx + 1 + d] = qmag * Fhat[d] * inv_Fhat / sqrt_er;
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

static double test_gr_m1_solver_norm(const double *u_old, const double *resid,
    double threshold)
{
    const double vmin = 1.0e5;
    double max_norm = 0.0;
    double rho_scale = fabs(u_old[PRJ_CONS_RHO]);
    double etot_scale = fabs(u_old[PRJ_CONS_ETOT]);
    double mom2 = 0.0;
    double mom_norm;
    int field;
    int group;
    int d;

    if (!isfinite(threshold) || threshold <= 0.0 ||
        !isfinite(u_old[PRJ_CONS_YE]) || u_old[PRJ_CONS_YE] <= 0.0 ||
        !isfinite(rho_scale) || rho_scale <= 0.0 ||
        !isfinite(etot_scale) || etot_scale <= 0.0) {
        return HUGE_VAL;
    }

    max_norm = fmax(max_norm, fabs(resid[PRJ_CONS_RHO]) / rho_scale);
    max_norm = fmax(max_norm, fabs(resid[PRJ_CONS_ETOT]) / etot_scale);
    max_norm = fmax(max_norm,
        fabs(resid[PRJ_CONS_YE]) / u_old[PRJ_CONS_YE]);

    for (d = 0; d < 3; ++d) {
        mom2 += u_old[PRJ_CONS_MOM1 + d] *
            u_old[PRJ_CONS_MOM1 + d];
    }
    if (mom2 == 0.0) {
        mom2 = rho_scale * rho_scale * vmin * vmin;
    }
    mom_norm = sqrt(mom2);
    {
        double rmom2 = 0.0;

        for (d = 0; d < 3; ++d) {
            rmom2 += resid[PRJ_CONS_MOM1 + d] *
                resid[PRJ_CONS_MOM1 + d];
        }
        max_norm = fmax(max_norm, sqrt(rmom2) / mom_norm);
    }

    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            int eidx = PRJ_CONS_RAD_E(field, group);
            int fidx = PRJ_CONS_RAD_F1(field, group);
            double scale_e = fmax(threshold * etot_scale / RAD_SCALE,
                fabs(u_old[eidx]));
            double flux2 = 0.0;
            double rflux2 = 0.0;
            double scale_f;

            max_norm = fmax(max_norm, fabs(resid[eidx]) / scale_e);
            for (d = 0; d < 3; ++d) {
                flux2 += u_old[fidx + d] * u_old[fidx + d];
                rflux2 += resid[fidx + d] * resid[fidx + d];
            }
            scale_f = fmax(threshold * mom_norm * PRJ_CLIGHT *
                    PRJ_CLIGHT / RAD_SCALE, sqrt(flux2));
            max_norm = fmax(max_norm, sqrt(rflux2) / scale_f);
        }
    }

    return max_norm;
}

static void check_gr_m1_solver_converged(const char *name,
    const double *u_old, const double *resid, double threshold)
{
    double norm = test_gr_m1_solver_norm(u_old, resid, threshold);

    if (!isfinite(norm) || norm >= threshold) {
        fprintf(stderr,
            "test_gr_m1_matter: %s solver norm %.17e threshold %.3e\n",
            name, norm, threshold);
        exit(1);
    }
}

static void make_gr_m1_exact_old_from_p(const prj_rad *rad, prj_eos *eos,
    const prj_z4c_hydro_geom *geom, const double *P, double dt,
    double *u_old)
{
    double zero[PRJ_NVAR_CONS];
    double resid[PRJ_NVAR_CONS];
    double u_new[PRJ_NVAR_CONS];
    int v;

    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        zero[v] = 0.0;
    }
    if (!prj_rad_gr_m1_residual_test_wrapper(rad, eos, geom, zero, P, dt,
            resid, u_new)) {
        die("exact old residual construction failed");
    }
    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        u_old[v] = resid[v];
    }
    if (u_old[PRJ_CONS_YE] <= 0.0) {
        die("exact old construction produced nonpositive YE");
    }
}

static void set_gr_m1_jacobian_test_p(
    double P[TEST_GR_M1_RESIDUAL_NP])
{
    int field;
    int group;

    P[0] = 1.0e8;
    P[1] = 1.31e9;
    P[2] = -8.4e8;
    P[3] = 5.7e8;
    P[4] = 1.0;
    P[5] = 0.5;
    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            int idx = field * PRJ_NEGROUP + group;
            int pidx = 6 + 4 * idx;
            double phase = (double)(idx + 1);
            double ER = 1.0 + 0.22 * sin(1.37 * phase) +
                0.11 * cos(0.73 * phase);

            /* Newton variables now carry the radiation four-velocity ur^i
             * directly (ur = qR / sqrt(ER)), so no sqrt(ER) scaling here. */
            P[pidx] = ER;
            P[pidx + 1] = 0.030 * sin(2.11 * phase);
            P[pidx + 2] = 0.024 * cos(1.67 * phase);
            P[pidx + 3] = 0.021 * sin(0.91 * phase + 0.4);
        }
    }
}

static void copy_gr_m1_p(double *dst, const double *src)
{
    int n;

    for (n = 0; n < TEST_GR_M1_RESIDUAL_NP; ++n) {
        dst[n] = src[n];
    }
}

static void perturb_gr_m1_solver_guess(double *P)
{
    int field;
    int group;

    P[0] *= 1.0002;
    P[1] += 2.5e5;
    P[2] -= 1.7e5;
    P[3] += 1.1e5;
    P[4] *= 0.9998;
    P[5] += 2.0e-5;
    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            int idx = field * PRJ_NEGROUP + group;
            int pidx = 6 + 4 * idx;
            double phase = (double)(idx + 1);

            /* ur^i are O(1) Newton variables; perturb them directly. */
            P[pidx] *= 1.0 + 1.5e-4 * sin(0.43 * phase);
            P[pidx + 1] += 1.0e-4 * sin(0.71 * phase);
            P[pidx + 2] -= 8.0e-5 * cos(0.59 * phase);
            P[pidx + 3] += 7.0e-5 * sin(0.37 * phase);
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
    double eos_q[PRJ_EOS_NQUANT];
    prj_z4c_hydro_geom geom;
    double beta[3] = {0.015, -0.01, 0.006};
    double gamma_diag[3] = {1.08, 0.94, 1.13};
    double vinit[3] = {1.31e9, -8.4e8, 5.7e8};
    double Fcov[3] = {0.07 * PRJ_CLIGHT, -0.03 * PRJ_CLIGHT,
        0.04 * PRJ_CLIGHT};
    double dt = 0.13;
    int col;
    int row;
    int v;
    int field;
    int group;

    init_real_test_eos(&eos);
    init_real_test_rad(&rad);
    set_diag_geom(&geom, 0.94, beta, gamma_diag);
    set_gr_m1_jacobian_test_p(P);
    prj_eos_rty(&eos, P[0], P[4], P[5], eos_q, PRJ_EOS_CTX_MAIN);
    set_flat_prim(W, P[0], vinit, eos_q[PRJ_EOS_EINT], P[5], 1.0, Fcov);
    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            int idx = field * PRJ_NEGROUP + group;
            double phase = (double)(idx + 1);
            double E = 0.9 + 0.18 * sin(0.61 * phase) +
                0.08 * cos(1.13 * phase);

            W[PRJ_PRIM_RAD_E(field, group)] = E;
            W[PRJ_PRIM_RAD_F1(field, group)] =
                0.045 * sin(1.19 * phase) * PRJ_CLIGHT * E;
            W[PRJ_PRIM_RAD_F2(field, group)] =
                0.035 * cos(0.97 * phase) * PRJ_CLIGHT * E;
            W[PRJ_PRIM_RAD_F3(field, group)] =
                0.030 * sin(0.53 * phase + 0.2) * PRJ_CLIGHT * E;
        }
    }
    prim2cons_or_die(&eos, &geom, W, u_old);
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
        if (col >= 6) {
            h = 1.0e-5 * fmax(1.0, fabs(P[col]));
            if (((col - 6) % 4) == 0) {
                h = fmin(h, 0.25 * P[col]);
            }
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

static void check_gr_m1_implicit_solve_real_tables(void)
{
    prj_eos eos;
    prj_rad rad;
    prj_z4c_hydro_geom geom;
    double beta[3] = {0.015, -0.01, 0.006};
    double gamma_diag[3] = {1.08, 0.94, 1.13};
    double P_root[TEST_GR_M1_RESIDUAL_NP];
    double P_guess[TEST_GR_M1_RESIDUAL_NP];
    double u_old[PRJ_NVAR_CONS];
    double resid[PRJ_NVAR_CONS];
    double u_new[PRJ_NVAR_CONS];
    double dt = 1.0e-6;

    init_real_test_eos(&eos);
    init_real_test_rad(&rad);
    rad.implicit_err_tol = 1.0e-8;
    rad.maxiter = 30;
    set_diag_geom(&geom, 0.94, beta, gamma_diag);
    set_gr_m1_jacobian_test_p(P_root);
    make_gr_m1_exact_old_from_p(&rad, &eos, &geom, P_root, dt, u_old);
    copy_gr_m1_p(P_guess, P_root);
    perturb_gr_m1_solver_guess(P_guess);

    if (!prj_rad_gr_m1_implicit_solve_test_wrapper(&rad, &eos, &geom,
            u_old, dt, P_guess, resid, u_new)) {
        die("real-table implicit solve failed");
    }
    check_gr_m1_solver_converged("real-table implicit solve", u_old, resid,
        rad.implicit_err_tol);
    prj_rad3_opac_free(&rad);
}

static void check_gr_m1_implicit_solve_zero_momentum_floor(void)
{
    prj_eos eos;
    prj_rad rad;
    prj_z4c_hydro_geom geom;
    double beta[3] = {0.0, 0.0, 0.0};
    double gamma_diag[3] = {1.0, 1.0, 1.0};
    double P_root[TEST_GR_M1_RESIDUAL_NP];
    double P_guess[TEST_GR_M1_RESIDUAL_NP];
    double u_old[PRJ_NVAR_CONS];
    double resid[PRJ_NVAR_CONS];
    double u_new[PRJ_NVAR_CONS];
    double dt = 0.05;
    int field;
    int group;
    int d;

    init_test_eos(&eos);
    init_source_test_rad_full(&rad, 0.0, 0.0, 0.0, 0.0, 0.0);
    rad.implicit_err_tol = 0.0;
    rad.maxiter = 30;
    set_diag_geom(&geom, 1.0, beta, gamma_diag);
    for (d = 0; d < TEST_GR_M1_RESIDUAL_NP; ++d) {
        P_root[d] = 0.0;
    }
    P_root[0] = 1.0;
    P_root[4] = 1.0;
    P_root[5] = 0.2;
    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            int idx = field * PRJ_NEGROUP + group;
            int pidx = 6 + 4 * idx;

            P_root[pidx] = 0.5 + 0.01 * (double)idx;
        }
    }
    make_gr_m1_exact_old_from_p(&rad, &eos, &geom, P_root, dt, u_old);
    if (u_old[PRJ_CONS_MOM1] != 0.0 || u_old[PRJ_CONS_MOM2] != 0.0 ||
        u_old[PRJ_CONS_MOM3] != 0.0) {
        die("zero-momentum floor setup has nonzero momentum");
    }
    copy_gr_m1_p(P_guess, P_root);
    P_guess[1] = 2.0e2;
    P_guess[6 + 1] = 1.0e-5 * sqrt(P_guess[6]);

    if (!prj_rad_gr_m1_implicit_solve_test_wrapper(&rad, &eos, &geom,
            u_old, dt, P_guess, resid, u_new)) {
        die("zero-momentum implicit solve failed");
    }
    check_gr_m1_solver_converged("zero-momentum implicit solve", u_old,
        resid, 1.0e-6);
    prj_rad3_opac_free(&rad);
}

static void check_gr_m1_matter_update_overwrites_solution(void)
{
    prj_mesh mesh;
    prj_coord coord;
    prj_eos eos;
    prj_rad rad;
    prj_z4c_hydro_geom geom;
    double beta[3] = {0.015, -0.01, 0.006};
    double gamma_diag[3] = {1.08, 0.94, 1.13};
    double P_root[TEST_GR_M1_RESIDUAL_NP];
    double u_old[PRJ_NVAR_CONS];
    double u_update[PRJ_NVAR_CONS];
    double resid[PRJ_NVAR_CONS];
    double u_expected[PRJ_NVAR_CONS];
    double prim_update[PRJ_NVAR_PRIM];
    double eos_q[PRJ_EOS_NQUANT];
    double final_temperature = -1.0;
    double dt = 1.0e-6;
    int ic = PRJ_BLOCK_SIZE / 2;
    int v;
    int field;
    int group;

    init_matter_test_mesh(&mesh, &coord);
    set_uniform_z4c(&mesh.blocks[0], 0.94, beta, gamma_diag);
    if (!prj_z4c_load_hydro_geom(&mesh, &mesh.blocks[0], 0, ic, ic, ic,
            &geom)) {
        die("matter update geometry load failed");
    }
    init_real_test_eos(&eos);
    init_real_test_rad(&rad);
    rad.implicit_err_tol = 1.0e-8;
    rad.maxiter = 30;
    set_gr_m1_jacobian_test_p(P_root);
    make_gr_m1_exact_old_from_p(&rad, &eos, &geom, P_root, dt, u_old);
    if (!prj_rad_gr_m1_residual_test_wrapper(&rad, &eos, &geom, u_old,
            P_root, dt, resid, u_expected)) {
        die("matter update expected state failed");
    }
    check_gr_m1_solver_converged("matter update expected state", u_old,
        resid, rad.implicit_err_tol);
    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        u_update[v] = u_old[v];
    }
    for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
        prim_update[v] = -1.0;
    }
    prj_eos_rty(&eos, P_root[0], P_root[4], P_root[5], eos_q,
        PRJ_EOS_CTX_MAIN);

    prj_rad_gr_m1_matter_update(&rad, &eos, &mesh, &mesh.blocks[0], 0,
        u_update, prim_update, ic, ic, ic, dt, &final_temperature);
    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        assert_close("matter update overwritten conserved", u_update[v],
            u_expected[v], 1.0e-6);
    }
    assert_close("matter update primitive rho", prim_update[PRJ_PRIM_RHO],
        P_root[0], 1.0e-10);
    assert_close("matter update primitive v1", prim_update[PRJ_PRIM_V1],
        P_root[1], 1.0e-10);
    assert_close("matter update primitive v2", prim_update[PRJ_PRIM_V2],
        P_root[2], 1.0e-10);
    assert_close("matter update primitive v3", prim_update[PRJ_PRIM_V3],
        P_root[3], 1.0e-10);
    assert_close("matter update primitive eint", prim_update[PRJ_PRIM_EINT],
        eos_q[PRJ_EOS_EINT], 1.0e-10);
    assert_close("matter update primitive Ye", prim_update[PRJ_PRIM_YE],
        P_root[5], 1.0e-10);
    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            int eidx = PRJ_CONS_RAD_E(field, group);
            int fidx = PRJ_CONS_RAD_F1(field, group);

            assert_close("matter update primitive radiation E",
                prim_update[PRJ_PRIM_RAD_E(field, group)],
                u_expected[eidx] / geom.sqrt_gamma, 1.0e-10);
            assert_close("matter update primitive radiation F1",
                prim_update[PRJ_PRIM_RAD_F1(field, group)],
                u_expected[fidx] / geom.sqrt_gamma, 1.0e-10);
            assert_close("matter update primitive radiation F2",
                prim_update[PRJ_PRIM_RAD_F2(field, group)],
                u_expected[fidx + 1] / geom.sqrt_gamma, 1.0e-10);
            assert_close("matter update primitive radiation F3",
                prim_update[PRJ_PRIM_RAD_F3(field, group)],
                u_expected[fidx + 2] / geom.sqrt_gamma, 1.0e-10);
        }
    }
    assert_close("matter update final temperature", final_temperature,
        P_root[4], 1.0e-6);
    prj_rad3_opac_free(&rad);
}

static void check_gr_m1_matter_update_clamps_flux(void)
{
    prj_mesh mesh;
    prj_coord coord;
    prj_eos eos;
    prj_rad rad;
    prj_z4c_hydro_geom geom;
    double beta[3] = {0.0, 0.0, 0.0};
    double gamma_diag[3] = {1.0, 1.0, 1.0};
    double vzero[3] = {0.0, 0.0, 0.0};
    double Fcov[3];
    double W[PRJ_NVAR_PRIM];
    double u[PRJ_NVAR_CONS];
    double u_initial[PRJ_NVAR_CONS];
    double prim_update[PRJ_NVAR_PRIM];
    double eos_q[PRJ_EOS_NQUANT];
    double final_temperature = -1.0;
    double E = 1.2;
    double E_floor = 1.0e-50;
    double F_clamped = (1.0 - 1.0e-6) * PRJ_CLIGHT * E;
    int ic = PRJ_BLOCK_SIZE / 2;
    int v;
    int field;
    int group;

    init_matter_test_mesh(&mesh, &coord);
    set_uniform_z4c(&mesh.blocks[0], 1.0, beta, gamma_diag);
    if (!prj_z4c_load_hydro_geom(&mesh, &mesh.blocks[0], 0, ic, ic, ic,
            &geom)) {
        die("matter clamp geometry load failed");
    }
    init_test_eos(&eos);
    init_source_test_rad_full(&rad, 0.0, 0.0, 0.0, 0.0, 0.0);
    rad.implicit_err_tol = 1.0e-8;
    rad.maxiter = 30;
    Fcov[0] = 2.0 * PRJ_CLIGHT * E;
    Fcov[1] = 0.0;
    Fcov[2] = 0.0;
    set_flat_prim(W, 1.0, vzero, 1.5, 0.2, E, Fcov);
    prim2cons_or_die(&eos, &geom, W, u);
    prj_eos_rey(&eos, W[PRJ_PRIM_RHO], W[PRJ_PRIM_EINT],
        W[PRJ_PRIM_YE], eos_q, PRJ_EOS_CTX_MAIN);
    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        u_initial[v] = u[v];
    }
    for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
        prim_update[v] = -1.0;
    }

    prj_rad_gr_m1_matter_update(&rad, &eos, &mesh, &mesh.blocks[0], 0,
        u, prim_update, ic, ic, ic, 0.0, &final_temperature);
    assert_close("matter clamp radiation E", u[PRJ_CONS_RAD_E(0, 0)],
        u_initial[PRJ_CONS_RAD_E(0, 0)], 1.0e-12);
    assert_close("matter clamp radiation F1", u[PRJ_CONS_RAD_F1(0, 0)],
        F_clamped, 1.0e-11);
    assert_close("matter clamp radiation F2", u[PRJ_CONS_RAD_F2(0, 0)],
        0.0, 1.0e-12);
    assert_close("matter clamp radiation F3", u[PRJ_CONS_RAD_F3(0, 0)],
        0.0, 1.0e-12);
    for (v = 0; v < PRJ_NVAR_MHD_CONS; ++v) {
        assert_close("matter clamp hydro unchanged", u[v], u_initial[v],
            1.0e-12);
    }
    assert_close("matter clamp primitive rho", prim_update[PRJ_PRIM_RHO],
        W[PRJ_PRIM_RHO], 1.0e-12);
    assert_close("matter clamp primitive eint", prim_update[PRJ_PRIM_EINT],
        W[PRJ_PRIM_EINT], 1.0e-12);
    assert_close("matter clamp primitive Ye", prim_update[PRJ_PRIM_YE],
        W[PRJ_PRIM_YE], 1.0e-12);
    assert_close("matter clamp primitive radiation E",
        prim_update[PRJ_PRIM_RAD_E(0, 0)], E, 1.0e-12);
    assert_close("matter clamp primitive radiation F1",
        prim_update[PRJ_PRIM_RAD_F1(0, 0)], F_clamped, 1.0e-11);
    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            if (field == 0 && group == 0) {
                continue;
            }
            assert_abs_close("matter clamp floored radiation E",
                u[PRJ_CONS_RAD_E(field, group)], E_floor, 1.0e-60);
            assert_abs_close("matter clamp floored primitive radiation E",
                prim_update[PRJ_PRIM_RAD_E(field, group)], E_floor,
                1.0e-60);
            assert_abs_close("matter clamp floored primitive radiation F1",
                prim_update[PRJ_PRIM_RAD_F1(field, group)], 0.0, 1.0e-60);
        }
    }
    assert_close("matter clamp final temperature", final_temperature,
        eos_q[PRJ_EOS_TEMPERATURE], 1.0e-12);
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

static void check_gr_m1_implicit_solve_invalid_states(void)
{
    prj_eos eos;
    prj_rad rad;
    prj_z4c_hydro_geom geom;
    prj_z4c_hydro_geom bad_geom;
    double beta[3] = {0.0, 0.0, 0.0};
    double gamma_diag[3] = {1.0, 1.0, 1.0};
    double P_root[TEST_GR_M1_RESIDUAL_NP];
    double P_bad[TEST_GR_M1_RESIDUAL_NP];
    double u_old[PRJ_NVAR_CONS];
    double u_bad[PRJ_NVAR_CONS];
    double resid[PRJ_NVAR_CONS];
    double dt = 0.05;
    int n;
    int field;
    int group;

    init_test_eos(&eos);
    init_source_test_rad_full(&rad, 0.0, 0.0, 0.0, 0.0, 0.0);
    rad.implicit_err_tol = 1.0e-8;
    rad.maxiter = 10;
    set_diag_geom(&geom, 1.0, beta, gamma_diag);
    for (n = 0; n < TEST_GR_M1_RESIDUAL_NP; ++n) {
        P_root[n] = 0.0;
    }
    P_root[0] = 1.0;
    P_root[4] = 1.0;
    P_root[5] = 0.2;
    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            int idx = field * PRJ_NEGROUP + group;

            P_root[6 + 4 * idx] = 0.5 + 0.01 * (double)idx;
        }
    }
    make_gr_m1_exact_old_from_p(&rad, &eos, &geom, P_root, dt, u_old);
    copy_gr_m1_p(P_bad, P_root);
    if (!prj_rad_gr_m1_implicit_solve_test_wrapper(&rad, &eos, &geom,
            u_old, dt, P_bad, resid, 0)) {
        die("valid implicit solve rejected");
    }

    copy_gr_m1_p(P_bad, P_root);
    P_bad[0] = 0.0;
    if (prj_rad_gr_m1_implicit_solve_test_wrapper(&rad, &eos, &geom,
            u_old, dt, P_bad, resid, 0)) {
        die("implicit solve accepted nonpositive rho");
    }
    copy_gr_m1_p(P_bad, P_root);
    P_bad[4] = 0.0;
    if (prj_rad_gr_m1_implicit_solve_test_wrapper(&rad, &eos, &geom,
            u_old, dt, P_bad, resid, 0)) {
        die("implicit solve accepted nonpositive temperature");
    }
    copy_gr_m1_p(P_bad, P_root);
    P_bad[1] = 1.01 * PRJ_CLIGHT;
    if (prj_rad_gr_m1_implicit_solve_test_wrapper(&rad, &eos, &geom,
            u_old, dt, P_bad, resid, 0)) {
        die("implicit solve accepted superluminal velocity");
    }
    copy_gr_m1_p(P_bad, P_root);
    P_bad[6] = -1.0;
    if (prj_rad_gr_m1_implicit_solve_test_wrapper(&rad, &eos, &geom,
            u_old, dt, P_bad, resid, 0)) {
        die("implicit solve accepted negative ER");
    }
    copy_gr_m1_p(P_bad, P_root);
    bad_geom = geom;
    bad_geom.alpha = 0.1;
    bad_geom.beta[0] = 1.0;
    if (prj_rad_gr_m1_implicit_solve_test_wrapper(&rad, &eos, &bad_geom,
            u_old, dt, P_bad, resid, 0)) {
        die("implicit solve accepted invalid q0");
    }
    copy_gr_m1_p(P_bad, P_root);
    for (n = 0; n < PRJ_NVAR_CONS; ++n) {
        u_bad[n] = u_old[n];
    }
    u_bad[PRJ_CONS_YE] = 0.0;
    if (prj_rad_gr_m1_implicit_solve_test_wrapper(&rad, &eos, &geom,
            u_bad, dt, P_bad, resid, 0)) {
        die("implicit solve accepted zero u_old YE");
    }
    for (n = 0; n < PRJ_NVAR_CONS; ++n) {
        u_bad[n] = u_old[n];
    }
    u_bad[PRJ_CONS_YE] = -fabs(u_old[PRJ_CONS_YE]);
    if (prj_rad_gr_m1_implicit_solve_test_wrapper(&rad, &eos, &geom,
            u_bad, dt, P_bad, resid, 0)) {
        die("implicit solve accepted negative u_old YE");
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
    check_gr_m1_implicit_solve_real_tables();
    check_gr_m1_implicit_solve_zero_momentum_floor();
    check_gr_m1_matter_update_overwrites_solution();
    check_gr_m1_matter_update_clamps_flux();
    check_gr_m1_implicit_solve_invalid_states();
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
