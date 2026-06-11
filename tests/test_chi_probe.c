#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#if defined(PRJ_ENABLE_MPI)
#include <mpi.h>
#endif

#include "prj.h"

/* Probe chi = kappa + sigma*(1 - delta/3) (the opacity entering the radiation
 * momentum update) at a given (rho, T, Ye), and report the per-step stiffness
 * dt*c*chi and the implicit flux-attenuation factor.  Usage:
 *   test_chi_probe rho T Ye [dt]
 */

#define TEST_RAD_TABLE_PARAM_FILE "../opacbin.extendT.param"
#define TEST_RAD_TABLE_FILE       "../opacity.SFHo.juo.horo.brem1.extendedT.bin"
#define CLIGHT 2.99792458e10

int main(int argc, char *argv[])
{
#if PRJ_NRAD == 0
    (void)argc; (void)argv;
    fprintf(stderr, "built without radiation\n");
    return 0;
#else
    prj_rad rad;
    const double emin_list[3] = {1.0, 1.0, 1.0};
    const double emax_list[3] = {300.0, 100.0, 100.0};
    double rho = (argc > 1) ? atof(argv[1]) : 3.1115297280e+09;
    double temp = (argc > 2) ? atof(argv[2]) : 5.6545656919e-01;
    double ye  = (argc > 3) ? atof(argv[3]) : 4.4519278407e-01;
    double dt  = (argc > 4) ? atof(argv[4]) : 8.523255e-06;
    double *kappa, *sigma, *delta, *eta;
    int nu, g;
    size_t nvals;

#if defined(PRJ_ENABLE_MPI)
    MPI_Init(&argc, &argv);
#endif
    memset(&rad, 0, sizeof(rad));
    rad.maxiter = 20;
    rad.implicit_err_tol = 1.0e-6;
    for (nu = 0; nu < PRJ_NRAD; ++nu) { rad.emin[nu] = emin_list[nu]; rad.emax[nu] = emax_list[nu]; }
    strncpy(rad.table_param_file, TEST_RAD_TABLE_PARAM_FILE, sizeof(rad.table_param_file) - 1);
    strncpy(rad.table_file, TEST_RAD_TABLE_FILE, sizeof(rad.table_file) - 1);
    prj_rad_init(&rad);

    nvals = (size_t)PRJ_NRAD * (size_t)PRJ_NEGROUP;
    kappa = malloc(nvals*sizeof(double)); sigma = malloc(nvals*sizeof(double));
    delta = malloc(nvals*sizeof(double)); eta = malloc(nvals*sizeof(double));
    prj_rad3_opac_lookup(&rad, rho, temp, ye, kappa, sigma, delta, eta);

    printf("# rho=%.6e  T=%.6e MeV  Ye=%.6e  dt=%.4e s  dt*c=%.4e cm\n", rho, temp, ye, dt, dt*CLIGHT);
    /* machine-readable: nu g kappa eta chi  (eta = emissivity from table) */
    for (nu = 0; nu < PRJ_NRAD; ++nu) {
        for (g = 0; g < PRJ_NEGROUP; ++g) {
            size_t idx = (size_t)nu*PRJ_NEGROUP + g;
            double chi = kappa[idx] + sigma[idx]*(1.0 - delta[idx]/3.0);
            printf("DAT %d %d %.8e %.8e %.8e\n", nu, g, kappa[idx], eta[idx], chi);
        }
    }
    double chi_max = 0.0;
    for (nu = 0; nu < PRJ_NRAD; ++nu) for (g = 0; g < PRJ_NEGROUP; ++g) {
        size_t idx=(size_t)nu*PRJ_NEGROUP+g; double chi=kappa[idx]+sigma[idx]*(1.0-delta[idx]/3.0);
        if (chi>chi_max) chi_max=chi;
    }
    printf("# chi_max=%.4e 1/cm  -> dt*c*chi_max=%.4e ; implicit factor 1/(1+dtcchi)-1 = %.4e\n",
           chi_max, dt*CLIGHT*chi_max, 1.0/(1.0+dt*CLIGHT*chi_max)-1.0);
    printf("# (force is at dt-cap when dt*c*chi >> 1, i.e. chi >> %.3e 1/cm)\n", 1.0/(dt*CLIGHT));
    free(kappa); free(sigma); free(delta); free(eta);
    prj_rad3_opac_free(&rad);
#if defined(PRJ_ENABLE_MPI)
    MPI_Finalize();
#endif
    return 0;
#endif
}
