#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(PRJ_ENABLE_MPI)
#include <mpi.h>
#endif

#include "prj.h"

#define TEST_ELEINEL_TABLE_DIR "../300.mixed/"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static void die(const char *msg)
{
    fprintf(stderr, "test_eleinel_detailed_balance: %s\n", msg);
    exit(1);
}

#if PRJ_NRAD > 0
static void build_test_egroups(prj_rad *rad)
{
    const double emin_list[3] = {1.0, 1.0, 1.0};
    const double emax_list[3] = {300.0, 100.0, 100.0};
    int nu;

    for (nu = 0; nu < PRJ_NRAD; ++nu) {
        double emin = emin_list[nu < 3 ? nu : 2];
        double emax = emax_list[nu < 3 ? nu : 2];
        double log_min = log(emin);
        double log_max = log(emax);
        double dlog = (log_max - log_min) / (double)PRJ_NEGROUP;
        int g;

        rad->egroup[nu] = (double *)malloc((size_t)PRJ_NEGROUP * sizeof(double));
        rad->eedge[nu] = (double *)malloc((size_t)(PRJ_NEGROUP + 1) * sizeof(double));
        rad->egroup_erg[nu] = (double *)malloc((size_t)PRJ_NEGROUP * sizeof(double));
        rad->degroup_erg[nu] = (double *)malloc((size_t)PRJ_NEGROUP * sizeof(double));
        rad->x_e[nu] = (double *)malloc((size_t)PRJ_NEGROUP * sizeof(double));
        if (rad->egroup[nu] == 0 || rad->eedge[nu] == 0 ||
            rad->egroup_erg[nu] == 0 || rad->degroup_erg[nu] == 0 ||
            rad->x_e[nu] == 0) {
            die("energy-group allocation failed");
        }

        for (g = 0; g <= PRJ_NEGROUP; ++g) {
            rad->eedge[nu][g] = exp(log_min + (double)g * dlog);
        }
        rad->eedge[nu][0] = emin;
        rad->eedge[nu][PRJ_NEGROUP] = emax;

        for (g = 0; g < PRJ_NEGROUP; ++g) {
            double ec = exp(log_min + ((double)g + 0.5) * dlog);
            double erg = ec * PRJ_MEV_TO_ERG;
            double derg = (rad->eedge[nu][g + 1] - rad->eedge[nu][g]) *
                PRJ_MEV_TO_ERG;

            rad->egroup[nu][g] = ec;
            rad->egroup_erg[nu][g] = erg;
            rad->degroup_erg[nu][g] = derg;
            if (nu == 0) {
                rad->x_e[nu][g] = RAD_SCALE / (PRJ_AVOGADRO * erg);
            } else if (nu == 1) {
                rad->x_e[nu][g] = -RAD_SCALE / (PRJ_AVOGADRO * erg);
            } else {
                rad->x_e[nu][g] = 0.0;
            }
        }
    }
}

static void free_test_egroups(prj_rad *rad)
{
    int nu;

    for (nu = 0; nu < PRJ_NRAD; ++nu) {
        free(rad->egroup[nu]);
        free(rad->eedge[nu]);
        free(rad->egroup_erg[nu]);
        free(rad->degroup_erg[nu]);
        free(rad->x_e[nu]);
        rad->egroup[nu] = 0;
        rad->eedge[nu] = 0;
        rad->egroup_erg[nu] = 0;
        rad->degroup_erg[nu] = 0;
        rad->x_e[nu] = 0;
    }
}

static double occupation_for_group(int nu, int g)
{
    double phase = (double)((3 * nu + 5 * g + 7) % 17) / 17.0;
    return 0.05 + 0.70 * phase;
}

static void check_state(const prj_rad *rad, double T, double etael,
    double *max_rel_out)
{
    const double rho = 1.0e12;
    const double Ye = 0.3;
    const double eta_factor = 1.0 / (4.0 * M_PI);
    double je[PRJ_NRAD * PRJ_NEGROUP];
    double he[PRJ_NRAD * PRJ_NEGROUP * PRJ_NDIM];
    double source[PRJ_NRAD * PRJ_NEGROUP];
    double sink[PRJ_NRAD * PRJ_NEGROUP];
    double max_rel = 0.0;
    double ye_resid = 0.0;
    double ye_abs = 0.0;
    int nu;
    int g;
    int d;

    memset(he, 0, sizeof(he));

    for (nu = 0; nu < PRJ_NRAD; ++nu) {
        double species_cut = (nu == 2) ? 4.0 : 1.0;

        for (g = 0; g < PRJ_NEGROUP; ++g) {
            int idx = nu * PRJ_NEGROUP + g;
            double occ = occupation_for_group(nu, g);
            double flux_occ = 0.03 * occ;

            je[idx] = occ * species_cut /
                rad->eleinel_factf_over_freqe3[idx];
            he[idx * PRJ_NDIM + 0] = flux_occ * species_cut /
                rad->eleinel_factf_over_freqe3[idx];
            for (d = 1; d < PRJ_NDIM; ++d) {
                he[idx * PRJ_NDIM + d] = 0.0;
            }
        }
    }

    prj_rad_eleinel_lookup(rad, rho, T, Ye, etael, je, he, source, sink, 0);

    for (nu = 0; nu < PRJ_NRAD; ++nu) {
        double species_resid = 0.0;
        double species_abs = 0.0;

        for (g = 0; g < PRJ_NEGROUP; ++g) {
            int idx = nu * PRJ_NEGROUP + g;
            double E_old = je[idx] * rad->degroup_erg[nu][g] /
                (PRJ_CLIGHT * RAD_SCALE * eta_factor * PRJ_MEV_TO_ERG);
            double source_phys = source[idx] * rad->degroup_erg[nu][g] /
                (PRJ_MEV_TO_ERG * eta_factor);
            double dEdt = source_phys - PRJ_CLIGHT * sink[idx] * E_old;
            double number_weight = RAD_SCALE /
                (PRJ_AVOGADRO * rad->egroup_erg[nu][g]);
            double species_term = number_weight * dEdt;
            double ye_term = rad->x_e[nu][g] * dEdt;

            species_resid += species_term;
            species_abs += fabs(species_term);
            ye_resid += ye_term;
            ye_abs += fabs(ye_term);
        }

        if (species_abs > 0.0) {
            double rel = fabs(species_resid) / species_abs;

            if (rel > max_rel) {
                max_rel = rel;
            }
        }
    }

    if (ye_abs > 0.0) {
        double rel = fabs(ye_resid) / ye_abs;

        if (rel > max_rel) {
            max_rel = rel;
        }
    }

    if (max_rel > *max_rel_out) {
        *max_rel_out = max_rel;
    }

    printf("T=%.8e etael=%.8e relative_number_residual=%.8e\n",
        T, etael, max_rel);
}
#endif

int main(int argc, char *argv[])
{
#if PRJ_NRAD == 0
    (void)argc;
    (void)argv;
    fprintf(stderr, "test_eleinel_detailed_balance: built without radiation (PRJ_NRAD=0)\n");
    return 0;
#else
    prj_rad rad;
    double max_rel = 0.0;
    const double tol = 5.0e-11;

#if defined(PRJ_ENABLE_MPI)
    MPI_Init(&argc, &argv);
#else
    (void)argc;
    (void)argv;
#endif

    if (PRJ_NRAD != 3) die("expected PRJ_NRAD=3");
    if (PRJ_NEGROUP != 12) die("expected PRJ_NEGROUP=12");
    if (PRJ_NDIM != 3) die("expected PRJ_NDIM=3");

    memset(&rad, 0, sizeof(rad));
    strncpy(rad.eleinel_table_dir, TEST_ELEINEL_TABLE_DIR,
        sizeof(rad.eleinel_table_dir) - 1);
    build_test_egroups(&rad);
    prj_rad_eleinel_init(&rad);

    check_state(&rad, 0.35, 12.0, &max_rel);
    check_state(&rad, 0.90, 22.0, &max_rel);
    check_state(&rad, 1.60, 35.0, &max_rel);

    prj_rad_eleinel_free(&rad);
    free_test_egroups(&rad);

    if (max_rel > tol) {
        fprintf(stderr,
            "test_eleinel_detailed_balance: max relative residual %.17e exceeds %.17e\n",
            max_rel, tol);
#if defined(PRJ_ENABLE_MPI)
        MPI_Finalize();
#endif
        return 1;
    }

#if defined(PRJ_ENABLE_MPI)
    MPI_Finalize();
#endif
    return 0;
#endif
}
