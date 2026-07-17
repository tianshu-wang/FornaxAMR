#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(PRJ_ENABLE_MPI)
#include <mpi.h>
#endif

#include "prj.h"

#define EOS_FILE "../eos_tmp/SFHoEOS__ye__0.035_0.56_50__logT_-4.793_2.176_500__logrho_-8.699_15.5_500_extend.dat"
#define OPAC_PARAM_FILE "../opacbin.extendT.param"
#define OPAC_FILE "../opacity.SFHo.juo.horo.brem1.extendedT.bin"

#if PRJ_NRAD > 0
static void die(const char *message)
{
    fprintf(stderr, "reproduce_rad_newton: %s\n", message);
    exit(2);
}

static void read_failure_input(const char *filename, double *u, double *dt, double *lapse)
{
    FILE *fp = fopen(filename, "r");
    char line[512];
    int found_dt = 0;
    int found[PRJ_NVAR_CONS] = {0};
    int nfound = 0;

    if (fp == 0) die("cannot open .err file");
    while (fgets(line, sizeof(line), fp) != 0) {
        double dt_in;
        double lapse_in;
        int idx;
        double value;

        if (sscanf(line, " raw input: dt=%lf lapse=%lf", &dt_in, &lapse_in) == 2) {
            *dt = dt_in;
            *lapse = lapse_in;
            found_dt = 1;
        }
        if (sscanf(line, " [%d] = %lf,", &idx, &value) == 2 &&
            idx >= 0 && idx < PRJ_NVAR_CONS) {
            u[idx] = value;
            if (!found[idx]) {
                found[idx] = 1;
                nfound += 1;
            }
        }
    }
    fclose(fp);
    if (!found_dt) die("did not find dt/lapse in .err file");
    if (nfound != PRJ_NVAR_CONS) die("did not find a complete u array in .err file");
}

int main(int argc, char **argv)
{
    prj_eos eos;
    prj_rad rad;
    double u[PRJ_NVAR_CONS] = {0};
    double eos_q[PRJ_EOS_NQUANT];
    double dt = 0.0;
    double lapse = 0.0;
    double final_temperature = 0.0;
    double u_before[PRJ_NVAR_CONS];
    double rho;
    double kinetic;
    double magnetic = 0.0;
    double eint;
    int maxiter = 20;
    int clamp_negative = 0;
    int nu;
    int g;

#if defined(PRJ_ENABLE_MPI)
    MPI_Init(&argc, &argv);
#endif
    if (argc < 2) {
        printf("reproduce_rad_newton: skipped (pass crash.err to run reproducer)\n");
#if defined(PRJ_ENABLE_MPI)
        MPI_Finalize();
#endif
        return 0;
    }
    if (argc > 4) {
        die("usage: reproduce_rad_newton crash.err [maxiter] [clamp_negative]");
    }
    if (argc >= 3) maxiter = atoi(argv[2]);
    if (argc >= 4) clamp_negative = atoi(argv[3]) != 0;
    if (PRJ_NRAD != 3 || PRJ_NEGROUP != 12 || !PRJ_MHD) {
        die("requires the magnetized_ccsn build (MHD=1, NRAD=3, NEGROUP=12)");
    }

    read_failure_input(argv[1], u, &dt, &lapse);
    if (clamp_negative) {
        for (nu = 0; nu < PRJ_NRAD; ++nu) {
            for (g = 0; g < PRJ_NEGROUP; ++g) {
                int idx = PRJ_CONS_RAD_E(nu, g);

                if (u[idx] < 0.0) u[idx] = 0.0;
            }
        }
    }
    memcpy(u_before, u, sizeof(u_before));

    memset(&rad, 0, sizeof(rad));
    rad.maxiter = maxiter;
    rad.implicit_err_tol = 1.0e-6;
    rad.emin[0] = 1.0; rad.emax[0] = 300.0;
    rad.emin[1] = 1.0; rad.emax[1] = 100.0;
    rad.emin[2] = 1.0; rad.emax[2] = 100.0;
    strncpy(rad.table_param_file, OPAC_PARAM_FILE, sizeof(rad.table_param_file) - 1);
    strncpy(rad.table_file, OPAC_FILE, sizeof(rad.table_file) - 1);
    prj_rad3_opac_init(&rad);

    memset(&eos, 0, sizeof(eos));
    eos.kind = PRJ_EOS_KIND_TABLE;
    strncpy(eos.filename, EOS_FILE, sizeof(eos.filename) - 1);
    prj_eos_init(&eos, 0);

    rho = u[PRJ_CONS_RHO];
    kinetic = 0.5 * (u[PRJ_CONS_MOM1] * u[PRJ_CONS_MOM1] +
        u[PRJ_CONS_MOM2] * u[PRJ_CONS_MOM2] +
        u[PRJ_CONS_MOM3] * u[PRJ_CONS_MOM3]) / rho;
#if PRJ_MHD
    magnetic = 0.5 * (u[PRJ_CONS_B1] * u[PRJ_CONS_B1] +
        u[PRJ_CONS_B2] * u[PRJ_CONS_B2] +
        u[PRJ_CONS_B3] * u[PRJ_CONS_B3]);
#endif
    eint = (u[PRJ_CONS_ETOT] - kinetic - magnetic) / rho;
    prj_eos_rey(&eos, rho, eint, u[PRJ_CONS_YE] / rho, eos_q, PRJ_EOS_CTX_MAIN);

    printf("input: dt=%.17e lapse=%.17e rho=%.17e T0=%.17e Ye0=%.17e "
           "maxiter=%d clamp_negative=%d\n",
        dt, lapse, rho, eos_q[PRJ_EOS_TEMPERATURE],
        u[PRJ_CONS_YE] / rho, maxiter, clamp_negative);
    for (nu = 0; nu < PRJ_NRAD; ++nu) {
        printf("E[%d]:", nu);
        for (g = 0; g < PRJ_NEGROUP; ++g) {
            printf(" %.9e", u[PRJ_CONS_RAD_E(nu, g)]);
        }
        printf("\n");
    }
    fflush(stdout);

    prj_rad_energy_update(&rad, &eos, u, dt, lapse, &final_temperature, 0);
    printf("converged: T=%.17e Ye=%.17e\n",
        final_temperature, u[PRJ_CONS_YE] / u[PRJ_CONS_RHO]);
    {
        double rad_before = 0.0;
        double rad_after = 0.0;
        double rad_negative_before = 0.0;
        double gas_before;
        double gas_after;

        for (nu = 0; nu < PRJ_NRAD; ++nu) {
            for (g = 0; g < PRJ_NEGROUP; ++g) {
                double eb = u_before[PRJ_CONS_RAD_E(nu, g)];
                double ea = u[PRJ_CONS_RAD_E(nu, g)];

                rad_before += eb * RAD_SCALE;
                rad_after += ea * RAD_SCALE;
                if (eb < 0.0) rad_negative_before += eb * RAD_SCALE;
            }
        }
        gas_before = u_before[PRJ_CONS_ETOT] - kinetic - magnetic;
        gas_after = u[PRJ_CONS_ETOT] - kinetic - magnetic;
        printf("energy: gas_before=%.17e gas_after=%.17e "
               "rad_before=%.17e rad_after=%.17e negative_rad_before=%.17e "
               "total_delta=%.17e\n",
            gas_before, gas_after, rad_before, rad_after, rad_negative_before,
            (gas_after + rad_after) - (gas_before + rad_before));
        for (nu = 0; nu < PRJ_NRAD; ++nu) {
            printf("E_final[%d]:", nu);
            for (g = 0; g < PRJ_NEGROUP; ++g) {
                printf(" %.9e", u[PRJ_CONS_RAD_E(nu, g)]);
            }
            printf("\n");
        }
    }

    prj_rad3_opac_free(&rad);
    free(eos.table);
#if defined(PRJ_ENABLE_MPI)
    MPI_Finalize();
#endif
    return 0;
}
#else
int main(int argc, char **argv)
{
#if defined(PRJ_ENABLE_MPI)
    MPI_Init(&argc, &argv);
#else
    (void)argc;
    (void)argv;
#endif
    printf("reproduce_rad_newton: skipped (radiation disabled)\n");
#if defined(PRJ_ENABLE_MPI)
    MPI_Finalize();
#endif
    return 0;
}
#endif
