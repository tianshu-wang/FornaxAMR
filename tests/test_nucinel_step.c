#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(PRJ_ENABLE_MPI)
#include <mpi.h>
#endif

#include "prj.h"

/*
 * Drives prj_rad_nucinel_step (Kompaneets nucleon inelastic scattering) at a
 * fixed (rho, T, Ye) and per-group u[PRJ_CONS_RAD_E(...)] read from
 * ../test1.txt.
 *
 * test1.txt layout:
 *   row 1:   rho [g/cc]  T [MeV]  Ye
 *   rows 2+: g_flat  u[CONS_RAD_E]  (extra columns ignored)
 * with g_flat = nu * PRJ_NEGROUP + g.
 *
 * Per-species xs[g] = egroup[nu][g] / T and Js[g] = E_g * spec_factor[nu][g]
 * (clamped to [0,1]) are the variables that the solver actually evolves; we
 * print them along with E_g before and after the Kompaneets step.
 */

#define TEST1_TXT_PATH            "../test1.txt"
#define TEST_RAD_TABLE_PARAM_FILE "../opacbin.extendT.param"
#define TEST_RAD_TABLE_FILE       "../opacity.SFHo.juo.horo.brem1.extendedT.bin"
#define TEST_ELEINEL_TABLE_DIR    "../300.mixed/"
#define TEST_EOS_FILE \
    "../eos_tmp/SFHoEOS__ye__0.035_0.56_50__logT_-4.793_2.176_500__logrho_-8.699_15.5_500_extend.dat"

/* dt (seconds) used by the Kompaneets step. */
#ifndef TEST_NUCINEL_DT
#define TEST_NUCINEL_DT 1.2228358071587572e-06
#endif

static void die(const char *msg)
{
    fprintf(stderr, "test_nucinel_step: %s\n", msg);
    exit(1);
}

int main(int argc, char *argv[])
{
#if PRJ_NRAD == 0
    (void)argc;
    (void)argv;
    fprintf(stderr, "test_nucinel_step: built without radiation (PRJ_NRAD=0)\n");
    return 0;
#else
    prj_rad rad;
    prj_eos eos;
    const double emin_list[3] = {1.0, 1.0, 1.0};
    const double emax_list[3] = {300.0, 100.0, 100.0};
    const int ntot = PRJ_NRAD * PRJ_NEGROUP;
    double rho;
    double T;
    double Ye;
    double E_in[PRJ_NRAD * PRJ_NEGROUP];
    double *u = 0;
    double eos_q[PRJ_EOS_NQUANT];
    double xs[PRJ_NRAD][PRJ_NEGROUP];
    double Js_clamped[PRJ_NRAD][PRJ_NEGROUP];
    double E_before[PRJ_NRAD * PRJ_NEGROUP];
    double E_after[PRJ_NRAD * PRJ_NEGROUP];
    double dt = TEST_NUCINEL_DT;
    FILE *fp;
    char line[1024];
    int status;
    int nu;
    int g;
    int i;
    int d;

#if defined(PRJ_ENABLE_MPI)
    MPI_Init(&argc, &argv);
#else
    (void)argc;
    (void)argv;
#endif

    if (PRJ_NRAD != 3) die("expected PRJ_NRAD=3");
    if (PRJ_NEGROUP != 12) die("expected PRJ_NEGROUP=12");

    /* --- Read inputs --- */
    fp = fopen(TEST1_TXT_PATH, "r");
    if (fp == 0) die("cannot open ../test1.txt");
    if (fgets(line, sizeof(line), fp) == 0) die("cannot read header row");
    if (sscanf(line, "%lf %lf %lf", &rho, &T, &Ye) != 3) die("bad header row");

    for (i = 0; i < ntot; ++i) {
        int idx_in;
        double e_val;

        if (fgets(line, sizeof(line), fp) == 0) die("not enough data rows in test1.txt");
        if (sscanf(line, "%d %lf", &idx_in, &e_val) < 2) die("bad data row");
        if (idx_in != i) {
            fprintf(stderr, "test_nucinel_step: expected index %d, got %d\n", i, idx_in);
            exit(1);
        }
        E_in[i] = e_val;
    }
    fclose(fp);

    /* --- Initialize radiation tables (need egroup + spec_factor + kom params). --- */
    memset(&rad, 0, sizeof(rad));
    rad.maxiter = 20;
    rad.implicit_err_tol = 1.0e-6;
    rad.kom_epsilon = 0.1;
    rad.kom_delta = 0.01;
    rad.kom_dtmin = 1.0e-20;
    rad.kom_rhocut = 1.0e13;
    for (nu = 0; nu < PRJ_NRAD; ++nu) {
        rad.emin[nu] = emin_list[nu];
        rad.emax[nu] = emax_list[nu];
        rad.kom_Ecut[nu] = 300.0;
    }
    strncpy(rad.table_param_file, TEST_RAD_TABLE_PARAM_FILE,
        sizeof(rad.table_param_file) - 1);
    strncpy(rad.table_file, TEST_RAD_TABLE_FILE,
        sizeof(rad.table_file) - 1);
    strncpy(rad.eleinel_table_dir, TEST_ELEINEL_TABLE_DIR,
        sizeof(rad.eleinel_table_dir) - 1);
    prj_rad_init(&rad);

    /* --- Initialize EOS (so prj_rad_nucinel_step can compute T internally). --- */
    memset(&eos, 0, sizeof(eos));
    eos.kind = PRJ_EOS_KIND_TABLE;
    strncpy(eos.filename, TEST_EOS_FILE, sizeof(eos.filename) - 1);
    prj_eos_init(&eos);
    if (eos.table_loaded != 1) die("EOS table failed to load");

    /* Get eint at the requested (rho, T, Ye) and stash it into u[ETOT] so that
     * prj_rad_nucinel_step's internal prj_eos_rey roundtrip recovers T. */
    prj_eos_rty(&eos, rho, T, Ye, eos_q);

    u = (double *)calloc((size_t)PRJ_NVAR_CONS, sizeof(double));
    if (u == 0) die("u allocation failed");
    u[PRJ_CONS_RHO] = rho;
    u[PRJ_CONS_MOM1] = 0.0;
    u[PRJ_CONS_MOM2] = 0.0;
    u[PRJ_CONS_MOM3] = 0.0;
    u[PRJ_CONS_ETOT] = rho * eos_q[PRJ_EOS_EINT];
    u[PRJ_CONS_YE] = rho * Ye;
#if PRJ_MHD
    u[PRJ_CONS_B1] = 0.0;
    u[PRJ_CONS_B2] = 0.0;
    u[PRJ_CONS_B3] = 0.0;
#endif
    for (nu = 0; nu < PRJ_NRAD; ++nu) {
        for (g = 0; g < PRJ_NEGROUP; ++g) {
            int idx = nu * PRJ_NEGROUP + g;
            u[PRJ_CONS_RAD_E(nu, g)] = E_in[idx];
            for (d = 0; d < PRJ_NDIM; ++d) {
                switch (d) {
                case 0: u[PRJ_CONS_RAD_F1(nu, g)] = 0.0; break;
                case 1: u[PRJ_CONS_RAD_F2(nu, g)] = 0.0; break;
                default: u[PRJ_CONS_RAD_F3(nu, g)] = 0.0; break;
                }
            }
            E_before[idx] = E_in[idx];
        }
    }

    /* --- Pre-step diagnostics: xs and Js the solver will actually use. --- */
    for (nu = 0; nu < PRJ_NRAD; ++nu) {
        for (g = 0; g < PRJ_NEGROUP; ++g) {
            double sf = rad.spec_factor[nu][g];
            double js = E_in[nu * PRJ_NEGROUP + g] * sf;

            xs[nu][g] = rad.egroup[nu][g] / T;
            if (js > 1.0) js = 1.0;
            if (js < 0.0) js = 0.0;
            Js_clamped[nu][g] = js;
        }
    }

    /* --- Run the Kompaneets step. --- */
    status = prj_rad_nucinel_step(&rad, &eos, u, dt);
    if (status != 1) {
        fprintf(stderr, "test_nucinel_step: prj_rad_nucinel_step reported status=%d\n",
            status);
    }

    for (nu = 0; nu < PRJ_NRAD; ++nu) {
        for (g = 0; g < PRJ_NEGROUP; ++g) {
            int idx = nu * PRJ_NEGROUP + g;
            E_after[idx] = u[PRJ_CONS_RAD_E(nu, g)];
        }
    }

    /* --- Emit results. --- */
    {
        const char *out_path = "nucinel_step.txt";
        FILE *ofp = fopen(out_path, "w");

        if (ofp == 0) die("cannot open nucinel_step.txt for writing");

        fprintf(ofp, "# rho = %.17e  T_MeV = %.17e  Ye = %.17e  dt = %.6e  status = %d\n",
            rho, T, Ye, dt, status);
        fprintf(ofp, "# %4s %4s %22s %22s %22s %22s\n",
            "nu", "g", "xs", "Js_clamped", "E_before", "E_after");
        printf("# rho = %.17e  T_MeV = %.17e  Ye = %.17e  dt = %.6e  status = %d\n",
            rho, T, Ye, dt, status);
        printf("# %4s %4s %22s %22s %22s %22s\n",
            "nu", "g", "xs", "Js_clamped", "E_before", "E_after");

        for (nu = 0; nu < PRJ_NRAD; ++nu) {
            for (g = 0; g < PRJ_NEGROUP; ++g) {
                int idx = nu * PRJ_NEGROUP + g;

                fprintf(ofp, "  %4d %4d %22.14e %22.14e %22.14e %22.14e\n",
                    nu, g, xs[nu][g], Js_clamped[nu][g],
                    E_before[idx], E_after[idx]);
                printf("  %4d %4d %22.14e %22.14e %22.14e %22.14e\n",
                    nu, g, xs[nu][g], Js_clamped[nu][g],
                    E_before[idx], E_after[idx]);
            }
        }
        fclose(ofp);
        printf("wrote %s\n", out_path);
    }

    free(u);
    prj_rad_eleinel_free(&rad);
    prj_rad3_opac_free(&rad);

#if defined(PRJ_ENABLE_MPI)
    MPI_Finalize();
#endif
    return 0;
#endif
}
