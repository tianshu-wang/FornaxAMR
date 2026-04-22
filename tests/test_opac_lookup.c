#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(PRJ_ENABLE_MPI)
#include <mpi.h>
#endif

#include "prj.h"

/*
 * Dumps the central energy (MeV), absorption, scattering, emissivity and
 * delta for every radiation field / energy group at a fixed (rho, T, Ye).
 * Expects NRAD=3, NEGROUP=12 and the opacity table files under ../ .
 */

#define TEST_RAD_TABLE_PARAM_FILE "../opacbin.extendT.param"
#define TEST_RAD_TABLE_FILE       "../opacity.SFHo.juo.horo.brem1.extendedT.bin"

static void die(const char *msg)
{
    fprintf(stderr, "test_opac_lookup: %s\n", msg);
    exit(1);
}

int main(int argc, char *argv[])
{
#if PRJ_NRAD == 0
    (void)argc;
    (void)argv;
    fprintf(stderr, "test_opac_lookup: built without radiation (PRJ_NRAD=0)\n");
    return 0;
#else
    prj_rad rad;
    const double emin_list[3] = {1.0, 1.0, 1.0};
    const double emax_list[3] = {300.0, 100.0, 100.0};
    double rho = 16400469286.503883;
    double temp = 0.60886570190175793;
    double ye = 0.42999964603729146;
    double *kappa = 0;
    double *sigma = 0;
    double *delta = 0;
    double *eta = 0;
    int nu;
    int g;
    size_t nvals;

#if defined(PRJ_ENABLE_MPI)
    MPI_Init(&argc, &argv);
#else
    (void)argc;
    (void)argv;
#endif

    if (PRJ_NRAD != 3) {
        die("expected PRJ_NRAD=3");
    }
    if (PRJ_NEGROUP != 12) {
        die("expected PRJ_NEGROUP=12");
    }

    memset(&rad, 0, sizeof(rad));
    rad.maxiter = 20;
    rad.implicit_err_tol = 1.0e-6;
    for (nu = 0; nu < PRJ_NRAD; ++nu) {
        rad.emin[nu] = emin_list[nu];
        rad.emax[nu] = emax_list[nu];
    }
    strncpy(rad.table_param_file, TEST_RAD_TABLE_PARAM_FILE,
        sizeof(rad.table_param_file) - 1);
    strncpy(rad.table_file, TEST_RAD_TABLE_FILE,
        sizeof(rad.table_file) - 1);

    prj_rad_init(&rad);

    nvals = (size_t)PRJ_NRAD * (size_t)PRJ_NEGROUP;
    kappa = (double *)malloc(nvals * sizeof(double));
    sigma = (double *)malloc(nvals * sizeof(double));
    delta = (double *)malloc(nvals * sizeof(double));
    eta = (double *)malloc(nvals * sizeof(double));
    if (kappa == 0 || sigma == 0 || delta == 0 || eta == 0) {
        die("allocation failed");
    }

    prj_rad3_opac_lookup(&rad, rho, temp, ye, kappa, sigma, delta, eta);

    {
        const char *out_path = "opac_emis.txt";
        FILE *ofp = fopen(out_path, "w");

        if (ofp == 0) {
            die("cannot open opac_emis.txt for writing");
        }
        fprintf(ofp, "# rho = %.17e  T = %.17e  Ye = %.17e\n", rho, temp, ye);
        fprintf(ofp, "# NRAD = %d  NEGROUP = %d\n", PRJ_NRAD, PRJ_NEGROUP);
        fprintf(ofp, "# emin/emax (MeV):");
        for (nu = 0; nu < PRJ_NRAD; ++nu) {
            fprintf(ofp, " [%g, %g]", rad.emin[nu], rad.emax[nu]);
        }
        fprintf(ofp, "\n");
        fprintf(ofp, "# %4s %4s %22s %22s %22s %22s %22s\n",
            "field", "g", "e_center[MeV]", "kappa_abs", "kappa_sca",
            "eta_emis", "delta");
        for (nu = 0; nu < PRJ_NRAD; ++nu) {
            for (g = 0; g < PRJ_NEGROUP; ++g) {
                size_t idx = (size_t)nu * (size_t)PRJ_NEGROUP + (size_t)g;

                fprintf(ofp, "  %4d %4d %22.14e %22.14e %22.14e %22.14e %22.14e\n",
                    nu, g, rad.egroup[nu][g], kappa[idx], sigma[idx],
                    eta[idx], delta[idx]);
            }
        }
        fclose(ofp);
        printf("wrote %s (%d x %d entries)\n", out_path, PRJ_NRAD, PRJ_NEGROUP);
    }

    free(kappa);
    free(sigma);
    free(delta);
    free(eta);
    prj_rad3_opac_free(&rad);

#if defined(PRJ_ENABLE_MPI)
    MPI_Finalize();
#endif
    return 0;
#endif
}
