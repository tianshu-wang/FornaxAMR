#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(PRJ_ENABLE_MPI)
#include <mpi.h>
#endif

#include "prj.h"

/*
 * Reproduces the table-lookup portion of prj_rad_eleinel_step at a fixed
 * (rho, T, Ye) and per-group (E_g, F1_g) read from ../test.txt.
 *
 * test.txt layout:
 *   row 1:   rho [g/cc]  T [MeV]  Ye
 *   rows 2+: g_flat  u[CONS_RAD_E]  u[CONS_RAD_F1]  ...
 * with g_flat = nu * PRJ_NEGROUP + g running 0..PRJ_NRAD*PRJ_NEGROUP-1.
 * F2 and F3 are taken to be zero (only F1 is in the file).
 *
 * For each (nu, g) we print the same source_phys and CLIGHT * sink_arr
 * that prj_rad_eleinel_step uses to form
 *     E_new = (E_old + dt*source_phys) / (1 + dt * CLIGHT * sink_arr).
 */

#define TEST_TXT_PATH             "../test.txt"
#define TEST_RAD_TABLE_PARAM_FILE "../opacbin.extendT.param"
#define TEST_RAD_TABLE_FILE       "../opacity.SFHo.juo.horo.brem1.extendedT.bin"
#define TEST_ELEINEL_TABLE_DIR    "../300.mixed/"
#define TEST_EOS_FILE \
    "../eos_tmp/SFHoEOS__ye__0.035_0.56_50__logT_-4.793_2.176_500__logrho_-8.699_15.5_500_extend.dat"

static void die(const char *msg)
{
    fprintf(stderr, "test_eleinel_lookup: %s\n", msg);
    exit(1);
}

int main(int argc, char *argv[])
{
#if PRJ_NRAD == 0
    (void)argc;
    (void)argv;
    fprintf(stderr, "test_eleinel_lookup: built without radiation (PRJ_NRAD=0)\n");
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
    double F1_in[PRJ_NRAD * PRJ_NEGROUP];
    double je[PRJ_NRAD * PRJ_NEGROUP];
    double he[PRJ_NRAD * PRJ_NEGROUP * PRJ_NDIM];
    double source_arr[PRJ_NRAD * PRJ_NEGROUP];
    double sink_arr[PRJ_NRAD * PRJ_NEGROUP];
    double scatt_arr[PRJ_NRAD * PRJ_NEGROUP];
    const double eta_factor = 1.0 / (4.0 * 3.14159265358979323846);
    double etael;
    FILE *fp;
    char line[1024];
    int nu;
    int g;
    int d;
    int i;

#if defined(PRJ_ENABLE_MPI)
    MPI_Init(&argc, &argv);
#else
    (void)argc;
    (void)argv;
#endif

    if (PRJ_NRAD != 3) die("expected PRJ_NRAD=3");
    if (PRJ_NEGROUP != 12) die("expected PRJ_NEGROUP=12");
    if (PRJ_NDIM != 3) die("expected PRJ_NDIM=3");

    fp = fopen(TEST_TXT_PATH, "r");
    if (fp == 0) die("cannot open ../test.txt");

    if (fgets(line, sizeof(line), fp) == 0) die("cannot read header row");
    if (sscanf(line, "%lf %lf %lf", &rho, &T, &Ye) != 3) die("bad header row");

    for (i = 0; i < ntot; ++i) {
        int idx_in;
        double e_val;
        double f1_val;
        double col4;
        double col5;

        if (fgets(line, sizeof(line), fp) == 0) die("not enough data rows in test.txt");
        if (sscanf(line, "%d %lf %lf %lf %lf", &idx_in, &e_val, &f1_val, &col4, &col5) < 3)
            die("bad data row");
        (void)col4;
        (void)col5;
        if (idx_in != i) {
            fprintf(stderr, "test_eleinel_lookup: expected index %d, got %d\n", i, idx_in);
            exit(1);
        }
        E_in[i] = e_val;
        F1_in[i] = f1_val;
    }
    fclose(fp);

    /* --- Set up the radiation tables (egroups + electron inelastic phi). --- */
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
    if (!rad.eleinel_table_loaded) die("electron inelastic tables failed to load");

    /* --- Set up the EOS so we can read etael at (rho, T, Ye). --- */
    memset(&eos, 0, sizeof(eos));
    eos.kind = PRJ_EOS_KIND_TABLE;
    strncpy(eos.filename, TEST_EOS_FILE, sizeof(eos.filename) - 1);
    prj_eos_init(&eos, 0);
    if (eos.table_loaded != 1) die("EOS table failed to load");

    etael = prj_eos_rty_geteta(&eos, rho, T, Ye, PRJ_EOS_CTX_MAIN);
    if (etael < -20.0) etael = -20.0;

    /* --- Build je / he exactly the way prj_rad_eleinel_step does, then look up. --- */
    for (nu = 0; nu < PRJ_NRAD; ++nu) {
        for (g = 0; g < PRJ_NEGROUP; ++g) {
            int idx = nu * PRJ_NEGROUP + g;
            double E_g = E_in[idx];
            double F1_g = F1_in[idx];
            double he_mag = 0.0;

            je[idx] = PRJ_CLIGHT * E_g * eta_factor * PRJ_MEV_TO_ERG / rad.degroup_erg[nu][g];
            if (je[idx] < 0.0) je[idx] = 0.0;

            for (d = 0; d < PRJ_NDIM; ++d) {
                double F_gd = (d == 0) ? F1_g : 0.0;
                he[idx * PRJ_NDIM + d] = F_gd * eta_factor * PRJ_MEV_TO_ERG / rad.degroup_erg[nu][g];
                he_mag += (he[idx * PRJ_NDIM + d] / (je[idx] + 1.0e-15))
                    * (he[idx * PRJ_NDIM + d] / (je[idx] + 1.0e-15));
            }
            he_mag = sqrt(he_mag);
            if (he_mag > 1.0) {
                for (d = 0; d < PRJ_NDIM; ++d) {
                    he[idx * PRJ_NDIM + d] /= he_mag;
                }
            }
        }
    }

    prj_rad_eleinel_lookup(&rad, rho, T, Ye, etael, je, he,
        source_arr, sink_arr, scatt_arr);

    /* --- Print source_phys and CLIGHT * sink_arr, the two quantities that
     *     enter E_new = (E_old + dt*source_phys) / (1 + dt * CLIGHT*sink). */
    {
        const char *out_path = "eleinel_lookup.txt";
        FILE *ofp = fopen(out_path, "w");

        if (ofp == 0) die("cannot open eleinel_lookup.txt for writing");

        fprintf(ofp, "# rho = %.17e  T_MeV = %.17e  Ye = %.17e  etael = %.17e\n",
            rho, T, Ye, etael);
        fprintf(ofp, "# %4s %4s %22s %22s %22s\n",
            "nu", "g", "e_center[MeV]", "source_phys", "CLIGHT*sink");
        printf("# rho = %.17e  T_MeV = %.17e  Ye = %.17e  etael = %.17e\n",
            rho, T, Ye, etael);
        printf("# %4s %4s %22s %22s %22s\n",
            "nu", "g", "e_center[MeV]", "source_phys", "CLIGHT*sink");

        for (nu = 0; nu < PRJ_NRAD; ++nu) {
            for (g = 0; g < PRJ_NEGROUP; ++g) {
                int idx = nu * PRJ_NEGROUP + g;
                double source_phys = source_arr[idx] * rad.degroup_erg[nu][g]
                    / (PRJ_MEV_TO_ERG * eta_factor);
                double c_sink = PRJ_CLIGHT * sink_arr[idx];

                fprintf(ofp, "  %4d %4d %22.14e %22.14e %22.14e\n",
                    nu, g, rad.egroup[nu][g], source_phys, c_sink);
                printf("  %4d %4d %22.14e %22.14e %22.14e\n",
                    nu, g, rad.egroup[nu][g], source_phys, c_sink);
            }
        }
        fclose(ofp);
        printf("wrote %s\n", out_path);
    }

    prj_rad_eleinel_free(&rad);
    prj_rad3_opac_free(&rad);

#if defined(PRJ_ENABLE_MPI)
    MPI_Finalize();
#endif
    return 0;
#endif
}
