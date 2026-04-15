#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "prj.h"
#include "prj_rad3_opac.h"

#if PRJ_NRAD > 0

#if defined(PRJ_ENABLE_MPI)
#include <mpi.h>
#endif

#define OPAC_IDX(ng, i, j, k, NG, NR, NT, NYE) \
    ((((ng) * (NR) + (i)) * (NT) + (j)) * (NYE) + (k))

#define TEMP_IDX(i, j, k, ng, NR, NT, NYE) \
    ((((ng) * (NYE) + (k)) * (NT) + (j)) * (NR) + (i))

static int prj_rad3_mpi_rank(void)
{
#if defined(PRJ_ENABLE_MPI)
    int flag = 0;
    int rank = 0;

    MPI_Initialized(&flag);
    if (flag) {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    }
    return rank;
#else
    return 0;
#endif
}

static void prj_rad3_build_egroups(prj_rad *rad)
{
    int nu;
    int g;

    for (nu = 0; nu < PRJ_NRAD; ++nu) {
        double emin = rad->emin[nu];
        double emax = rad->emax[nu];
        double log_min;
        double log_max;
        double dlog;

        if (emin <= 0.0 || emax <= emin) {
            fprintf(stderr, "prj_rad3_opac: invalid emin/emax for field %d\n", nu);
            exit(1);
        }
        rad->egroup[nu] = (double *)malloc((size_t)PRJ_NEGROUP * sizeof(double));
        rad->eedge[nu] = (double *)malloc((size_t)(PRJ_NEGROUP + 1) * sizeof(double));
        rad->egroup_erg[nu] = (double *)malloc((size_t)PRJ_NEGROUP * sizeof(double));
        rad->degroup_erg[nu] = (double *)malloc((size_t)PRJ_NEGROUP * sizeof(double));
        rad->x_e[nu] = (double *)malloc((size_t)PRJ_NEGROUP * sizeof(double));
        if (rad->egroup[nu] == 0 || rad->eedge[nu] == 0 || rad->egroup_erg[nu] == 0 ||
            rad->degroup_erg[nu] == 0 || rad->x_e[nu] == 0) {
            fprintf(stderr, "prj_rad3_opac: allocation failed for egroup field %d\n", nu);
            exit(1);
        }
        log_min = log(emin);
        log_max = log(emax);
        dlog = (log_max - log_min) / (double)PRJ_NEGROUP;
        for (g = 0; g <= PRJ_NEGROUP; ++g) {
            rad->eedge[nu][g] = exp(log_min + (double)g * dlog);
        }
        rad->eedge[nu][0] = emin;
        rad->eedge[nu][PRJ_NEGROUP] = emax;
        for (g = 0; g < PRJ_NEGROUP; ++g) {
            double ec = exp(log_min + ((double)g + 0.5) * dlog);
            double erg = ec * PRJ_MEV_TO_ERG;
            double derg = (rad->eedge[nu][g + 1] - rad->eedge[nu][g]) * PRJ_MEV_TO_ERG;

            rad->egroup[nu][g] = ec;
            rad->egroup_erg[nu][g] = erg;
            rad->degroup_erg[nu][g] = derg;
            if (nu == 0) {
                rad->x_e[nu][g] = 1.0 / (PRJ_AVOGADRO * erg);
            } else if (nu == 1) {
                rad->x_e[nu][g] = -1.0 / (PRJ_AVOGADRO * erg);
            } else {
                rad->x_e[nu][g] = 0.0;
            }
        }
    }
}

static void prj_rad3_read_param(prj_rad *rad)
{
    int rank = prj_rad3_mpi_rank();
    int i;

    if (rank == 0) {
        FILE *fp = fopen(rad->table_param_file, "r");
        int ngmax_check;

        if (fp == 0) {
            fprintf(stderr, "prj_rad3_opac: cannot open %s\n", rad->table_param_file);
            exit(1);
        }
        if (fscanf(fp, "%d", &rad->ngmax) != 1) {
            fprintf(stderr, "prj_rad3_opac: failed to read ngmax\n");
            exit(1);
        }
        if (fscanf(fp, "%lf %lf %d", &rad->romin, &rad->romax, &rad->nromax) != 3) {
            fprintf(stderr, "prj_rad3_opac: failed to read rho grid\n");
            exit(1);
        }
        if (fscanf(fp, "%lf %lf %d", &rad->tmin, &rad->tmax, &rad->ntmax) != 3) {
            fprintf(stderr, "prj_rad3_opac: failed to read T grid\n");
            exit(1);
        }
        if (fscanf(fp, "%lf %lf %d", &rad->yemin, &rad->yemax, &rad->nyemax) != 3) {
            fprintf(stderr, "prj_rad3_opac: failed to read Ye grid\n");
            exit(1);
        }
        if (fscanf(fp, "%d", &ngmax_check) != 1 || ngmax_check != rad->ngmax) {
            fprintf(stderr, "prj_rad3_opac: ngmax inconsistency in param file\n");
            exit(1);
        }
        rad->enuk = (double *)malloc((size_t)rad->ngmax * sizeof(double));
        if (rad->enuk == 0) {
            fprintf(stderr, "prj_rad3_opac: enuk allocation failed\n");
            exit(1);
        }
        for (i = 0; i < rad->ngmax; ++i) {
            if (fscanf(fp, "%lf", &rad->enuk[i]) != 1) {
                fprintf(stderr, "prj_rad3_opac: failed to read enuk[%d]\n", i);
                exit(1);
            }
        }
        fclose(fp);
    }

#if defined(PRJ_ENABLE_MPI)
    MPI_Bcast(&rad->ngmax, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&rad->nromax, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&rad->ntmax, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&rad->nyemax, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&rad->romin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&rad->romax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&rad->tmin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&rad->tmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&rad->yemin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&rad->yemax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (rank != 0) {
        rad->enuk = (double *)malloc((size_t)rad->ngmax * sizeof(double));
        if (rad->enuk == 0) {
            fprintf(stderr, "prj_rad3_opac: enuk allocation failed\n");
            exit(1);
        }
    }
    MPI_Bcast(rad->enuk, rad->ngmax, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
}

static void prj_rad3_load_block(prj_rad *rad, FILE *fp, int record_base,
    double **dest, float *tmp)
{
    int rank = prj_rad3_mpi_rank();
    int nromax = rad->nromax;
    int ntmax = rad->ntmax;
    int nyemax = rad->nyemax;
    int ngmax = rad->ngmax;
    size_t field_bytes = (size_t)nromax * (size_t)ntmax * (size_t)nyemax * (size_t)ngmax;
    size_t opac_bytes = (size_t)PRJ_NEGROUP * (size_t)nromax * (size_t)ntmax * (size_t)nyemax;
    int irecl = 4 * nromax * ntmax * nyemax * ngmax;
    double frmin;
    double frmax;
    double dfr;
    int nu;

    frmin = log(rad->enuk[0]);
    frmax = log(rad->enuk[ngmax - 1]);
    dfr = (frmax - frmin) / (double)(ngmax - 1);

    for (nu = 0; nu < PRJ_NRAD; ++nu) {
        dest[nu] = (double *)malloc(opac_bytes * sizeof(double));
        if (dest[nu] == 0) {
            fprintf(stderr, "prj_rad3_opac: opacity allocation failed for field %d\n", nu);
            exit(1);
        }
    }

    if (rank == 0) {
        int nu_r;
        int i;
        int j;
        int k;
        int ng;

        for (nu_r = 0; nu_r < 3; ++nu_r) {
            long offset = (long)(record_base + nu_r) * (long)irecl;
            size_t read_count;

            if (fseek(fp, offset, SEEK_SET) != 0) {
                fprintf(stderr, "prj_rad3_opac: fseek failed\n");
                exit(1);
            }
            read_count = fread(&tmp[(size_t)nu_r * field_bytes], sizeof(float),
                field_bytes, fp);
            if (read_count != field_bytes) {
                fprintf(stderr, "prj_rad3_opac: fread failed\n");
                exit(1);
            }
        }
        for (nu_r = 0; nu_r < PRJ_NRAD; ++nu_r) {
            int src_nu = nu_r < 3 ? nu_r : 2;
            double *restrict out = dest[nu_r];

            for (k = 0; k < nyemax; ++k) {
                for (j = 0; j < ntmax; ++j) {
                    for (i = 0; i < nromax; ++i) {
                        for (ng = 0; ng < PRJ_NEGROUP; ++ng) {
                            double el = log(rad->egroup[nu_r][ng]);
                            int jf = (int)((el - frmin) / dfr);
                            double dfi;
                            float v0;
                            float v1;

                            if (jf < 0) {
                                jf = 0;
                                dfi = 0.0;
                            } else if (jf > ngmax - 2) {
                                jf = ngmax - 2;
                                dfi = 1.0;
                            } else {
                                double f1i = frmin + (double)jf * dfr;

                                dfi = (el - f1i) / dfr;
                            }
                            v0 = tmp[TEMP_IDX(i, j, k, jf, nromax, ntmax, nyemax) +
                                (size_t)src_nu * field_bytes];
                            v1 = tmp[TEMP_IDX(i, j, k, jf + 1, nromax, ntmax, nyemax) +
                                (size_t)src_nu * field_bytes];
                            out[OPAC_IDX(ng, i, j, k, PRJ_NEGROUP, nromax, ntmax, nyemax)] =
                                (1.0 - dfi) * (double)v0 + dfi * (double)v1;
                        }
                    }
                }
            }
        }
    }

#if defined(PRJ_ENABLE_MPI)
    for (nu = 0; nu < PRJ_NRAD; ++nu) {
        MPI_Bcast(dest[nu], (int)opac_bytes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
#endif
}

void prj_rad3_opac_init(prj_rad *rad)
{
    int rank = prj_rad3_mpi_rank();
    FILE *fp = 0;
    float *tmp = 0;
    size_t tmp_count;

    if (rad == 0) {
        fprintf(stderr, "prj_rad3_opac: null rad\n");
        exit(1);
    }

    prj_rad3_build_egroups(rad);
    prj_rad3_read_param(rad);

    tmp_count = (size_t)rad->nromax * (size_t)rad->ntmax * (size_t)rad->nyemax *
        (size_t)rad->ngmax * 3U;
    if (rank == 0) {
        fp = fopen(rad->table_file, "rb");
        if (fp == 0) {
            fprintf(stderr, "prj_rad3_opac: cannot open %s\n", rad->table_file);
            exit(1);
        }
        tmp = (float *)malloc(tmp_count * sizeof(float));
        if (tmp == 0) {
            fprintf(stderr, "prj_rad3_opac: temp allocation failed\n");
            exit(1);
        }
    }

    prj_rad3_load_block(rad, fp, 0, rad->absopac, tmp);
    prj_rad3_load_block(rad, fp, 3, rad->scaopac, tmp);
    prj_rad3_load_block(rad, fp, 6, rad->emis, tmp);
    prj_rad3_load_block(rad, fp, 9, rad->sdelta, tmp);

    if (rank == 0) {
        fclose(fp);
        free(tmp);
    }
}

void prj_rad3_opac_free(prj_rad *rad)
{
    int nu;

    if (rad == 0) {
        return;
    }
    for (nu = 0; nu < PRJ_NRAD; ++nu) {
        free(rad->egroup[nu]);
        free(rad->eedge[nu]);
        free(rad->egroup_erg[nu]);
        free(rad->degroup_erg[nu]);
        free(rad->x_e[nu]);
        free(rad->absopac[nu]);
        free(rad->scaopac[nu]);
        free(rad->emis[nu]);
        free(rad->sdelta[nu]);
        rad->egroup[nu] = 0;
        rad->eedge[nu] = 0;
        rad->egroup_erg[nu] = 0;
        rad->degroup_erg[nu] = 0;
        rad->x_e[nu] = 0;
        rad->absopac[nu] = 0;
        rad->scaopac[nu] = 0;
        rad->emis[nu] = 0;
        rad->sdelta[nu] = 0;
    }
    free(rad->enuk);
    rad->enuk = 0;
}

void prj_rad3_opac_lookup(const prj_rad *rad, double rho, double temp, double ye,
    double *kappa, double *sigma, double *delta, double *eta)
{
    int nromax = rad->nromax;
    int ntmax = rad->ntmax;
    int nyemax = rad->nyemax;
    double rl = log(rho);
    double tl = log(temp);
    double ye_c = ye;
    double tl_c = tl;
    double deltar;
    double deltat;
    double deltaye;
    int jr;
    int jt;
    int jye;
    double r1i;
    double r2i;
    double dri;
    double t1i;
    double t2i;
    double dti;
    double ye1i;
    double ye2i;
    double dyei;
    double coff[8];
    int nu;
    int ng;

    if (ye_c > rad->yemax) ye_c = rad->yemax;
    if (ye_c < rad->yemin) ye_c = rad->yemin;
    if (tl_c > log(rad->tmax)) tl_c = log(rad->tmax);
    if (tl_c < log(rad->tmin)) tl_c = log(rad->tmin);

    deltar = (rl - log(rad->romin)) / (log(rad->romax) - log(rad->romin)) *
        (double)(nromax - 1);
    deltat = (tl_c - log(rad->tmin)) / (log(rad->tmax) - log(rad->tmin)) *
        (double)(ntmax - 1);
    deltaye = (ye_c - rad->yemin) / (rad->yemax - rad->yemin) * (double)(nyemax - 1);

    jr = (int)deltar;
    jt = (int)deltat;
    jye = (int)deltaye;
    if (jr < 0) jr = 0;
    if (jr > nromax - 2) jr = nromax - 2;
    if (jt < 0) jt = 0;
    if (jt > ntmax - 2) jt = ntmax - 2;
    if (jye < 0) jye = 0;
    if (jye > nyemax - 2) jye = nyemax - 2;

    r1i = log(rad->romin) + (log(rad->romax) - log(rad->romin)) * (double)jr / (double)(nromax - 1);
    r2i = log(rad->romin) + (log(rad->romax) - log(rad->romin)) * (double)(jr + 1) / (double)(nromax - 1);
    dri = (rl - r1i) / (r2i - r1i);
    if (dri < 0.0) dri = 0.0;

    t1i = log(rad->tmin) + (log(rad->tmax) - log(rad->tmin)) * (double)jt / (double)(ntmax - 1);
    t2i = log(rad->tmin) + (log(rad->tmax) - log(rad->tmin)) * (double)(jt + 1) / (double)(ntmax - 1);
    dti = (tl_c - t1i) / (t2i - t1i);
    if (dti < 0.0) dti = 0.0;

    ye1i = rad->yemin + (rad->yemax - rad->yemin) * (double)jye / (double)(nyemax - 1);
    ye2i = rad->yemin + (rad->yemax - rad->yemin) * (double)(jye + 1) / (double)(nyemax - 1);
    dyei = (ye_c - ye1i) / (ye2i - ye1i);

    coff[0] = (1.0 - dri) * (1.0 - dti) * (1.0 - dyei);
    coff[1] = dri * (1.0 - dti) * (1.0 - dyei);
    coff[2] = (1.0 - dri) * dti * (1.0 - dyei);
    coff[3] = dri * dti * (1.0 - dyei);
    coff[4] = (1.0 - dri) * (1.0 - dti) * dyei;
    coff[5] = dri * (1.0 - dti) * dyei;
    coff[6] = (1.0 - dri) * dti * dyei;
    coff[7] = dri * dti * dyei;

    for (nu = 0; nu < PRJ_NRAD; ++nu) {
        const double *absopac = rad->absopac[nu];
        const double *scaopac = rad->scaopac[nu];
        const double *emis = rad->emis[nu];
        const double *sdelta = rad->sdelta[nu];

        for (ng = 0; ng < PRJ_NEGROUP; ++ng) {
            size_t base = (size_t)nu * (size_t)PRJ_NEGROUP + (size_t)ng;
            double k_val;
            double s_val;
            double d_val;
            double j_val;

#define SAMPLE(ARR) \
            (coff[0] * (ARR)[OPAC_IDX(ng, jr,     jt,     jye,     PRJ_NEGROUP, nromax, ntmax, nyemax)] + \
             coff[1] * (ARR)[OPAC_IDX(ng, jr + 1, jt,     jye,     PRJ_NEGROUP, nromax, ntmax, nyemax)] + \
             coff[2] * (ARR)[OPAC_IDX(ng, jr,     jt + 1, jye,     PRJ_NEGROUP, nromax, ntmax, nyemax)] + \
             coff[3] * (ARR)[OPAC_IDX(ng, jr + 1, jt + 1, jye,     PRJ_NEGROUP, nromax, ntmax, nyemax)] + \
             coff[4] * (ARR)[OPAC_IDX(ng, jr,     jt,     jye + 1, PRJ_NEGROUP, nromax, ntmax, nyemax)] + \
             coff[5] * (ARR)[OPAC_IDX(ng, jr + 1, jt,     jye + 1, PRJ_NEGROUP, nromax, ntmax, nyemax)] + \
             coff[6] * (ARR)[OPAC_IDX(ng, jr,     jt + 1, jye + 1, PRJ_NEGROUP, nromax, ntmax, nyemax)] + \
             coff[7] * (ARR)[OPAC_IDX(ng, jr + 1, jt + 1, jye + 1, PRJ_NEGROUP, nromax, ntmax, nyemax)])

            k_val = SAMPLE(absopac);
            s_val = SAMPLE(scaopac);
            d_val = SAMPLE(sdelta);
            j_val = SAMPLE(emis);

#undef SAMPLE

            kappa[base] = exp(k_val) * rho;
            sigma[base] = exp(s_val) * rho;
            delta[base] = d_val;
            eta[base] = exp(j_val) * rho;
        }
    }
}

#endif /* PRJ_NRAD > 0 */
