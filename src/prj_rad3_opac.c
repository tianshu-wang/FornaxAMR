#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "prj.h"
#include "prj_rad3_opac.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#if PRJ_NRAD > 0

#if defined(PRJ_ENABLE_MPI)
#include <mpi.h>
#if PRJ_MIXED_PRECISION_TABLE
#define PRJ_TABLE_MPI_TYPE MPI_FLOAT
#else
#define PRJ_TABLE_MPI_TYPE MPI_DOUBLE
#endif
#endif

/* Runtime opacity tables are corner-major with energy groups contiguous:
 * (rho index, temperature index, Ye index) -> group. */
#define OPAC_CELL_IDX(i, j, k, NG, NR, NT, NYE) \
    ((((((size_t)(i) * (size_t)(NT)) + (size_t)(j)) * (size_t)(NYE) + \
        (size_t)(k)) * (size_t)(NG)))

#define OPAC_IDX(ng, i, j, k, NG, NR, NT, NYE) \
    (OPAC_CELL_IDX(i, j, k, NG, NR, NT, NYE) + (size_t)(ng))

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
        rad->egroup[nu] = (double *)prj_malloc((size_t)PRJ_NEGROUP * sizeof(double));
        rad->eedge[nu] = (double *)prj_malloc((size_t)(PRJ_NEGROUP + 1) * sizeof(double));
        rad->egroup_erg[nu] = (double *)prj_malloc((size_t)PRJ_NEGROUP * sizeof(double));
        rad->degroup_erg[nu] = (double *)prj_malloc((size_t)PRJ_NEGROUP * sizeof(double));
        rad->x_e[nu] = (double *)prj_malloc((size_t)PRJ_NEGROUP * sizeof(double));
        rad->log_egroup[nu] = (double *)prj_malloc((size_t)PRJ_NEGROUP * sizeof(double));
        rad->spec_factor[nu] = (double *)prj_malloc((size_t)PRJ_NEGROUP * sizeof(double));
        if (rad->egroup[nu] == 0 || rad->eedge[nu] == 0 || rad->egroup_erg[nu] == 0 ||
            rad->degroup_erg[nu] == 0 || rad->x_e[nu] == 0 || rad->log_egroup[nu] == 0 ||
            rad->spec_factor[nu] == 0) {
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
            rad->log_egroup[nu][g] = log(ec);
            rad->egroup_erg[nu][g] = erg;
            rad->degroup_erg[nu][g] = derg;
            if (nu == 2) {
                rad->spec_factor[nu][g] = 0.25 * pow(PRJ_CLIGHT * PRJ_HPLANCK, 3)
                    / (4.0 * M_PI * derg * erg * erg * erg);
            } else {
                rad->spec_factor[nu][g] = pow(PRJ_CLIGHT * PRJ_HPLANCK, 3)
                    / (4.0 * M_PI * derg * erg * erg * erg);
            }
            /* spec_factor converts radiation energy density to occupation
               number; with E stored in RAD_SCALE*erg units the conversion
               carries an extra RAD_SCALE so the occupation stays physical. */
            rad->spec_factor[nu][g] *= RAD_SCALE;
            /* x_e converts a radiation-energy change to a lepton-number change;
               RAD_SCALE is baked in here so all Ye couplings need no further
               rescaling. */
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
        rad->enuk = (double *)prj_malloc((size_t)rad->ngmax * sizeof(double));
        rad->log_enuk = (double *)prj_malloc((size_t)rad->ngmax * sizeof(double));
        if (rad->enuk == 0 || rad->log_enuk == 0) {
            fprintf(stderr, "prj_rad3_opac: enuk allocation failed\n");
            exit(1);
        }
        for (i = 0; i < rad->ngmax; ++i) {
            if (fscanf(fp, "%lf", &rad->enuk[i]) != 1) {
                fprintf(stderr, "prj_rad3_opac: failed to read enuk[%d]\n", i);
                exit(1);
            }
            rad->log_enuk[i] = log(rad->enuk[i]);
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
        rad->enuk = (double *)prj_malloc((size_t)rad->ngmax * sizeof(double));
        rad->log_enuk = (double *)prj_malloc((size_t)rad->ngmax * sizeof(double));
        if (rad->enuk == 0 || rad->log_enuk == 0) {
            fprintf(stderr, "prj_rad3_opac: enuk allocation failed\n");
            exit(1);
        }
    }
    MPI_Bcast(rad->enuk, rad->ngmax, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

    for (i = 0; i < rad->ngmax; ++i) {
        rad->log_enuk[i] = log(rad->enuk[i]);
    }

    rad->log_romin = log(rad->romin);
    rad->log_romax = log(rad->romax);
    rad->log_tmin = log(rad->tmin);
    rad->log_tmax = log(rad->tmax);
    rad->inv_logrho_span = (rad->log_romax > rad->log_romin) ?
        1.0 / (rad->log_romax - rad->log_romin) : 0.0;
    rad->inv_logtemp_span = (rad->log_tmax > rad->log_tmin) ?
        1.0 / (rad->log_tmax - rad->log_tmin) : 0.0;
    rad->inv_ye_span = (rad->yemax > rad->yemin) ?
        1.0 / (rad->yemax - rad->yemin) : 0.0;
}

static void prj_rad3_load_block(prj_rad *rad, FILE *fp, int record_base,
    prj_table_real **dest, float *tmp)
{
    int rank = prj_rad3_mpi_rank();
    int nromax = rad->nromax;
    int ntmax = rad->ntmax;
    int nyemax = rad->nyemax;
    int ngmax = rad->ngmax;
    size_t field_bytes = (size_t)nromax * (size_t)ntmax * (size_t)nyemax * (size_t)ngmax;
    size_t opac_bytes = (size_t)PRJ_NEGROUP * (size_t)nromax * (size_t)ntmax * (size_t)nyemax;
    int irecl = 4 * nromax * ntmax * nyemax * ngmax + 10;
    double frmin;
    double frmax;
    double dfr;
    int nu;

    frmin = rad->log_enuk[0];
    frmax = rad->log_enuk[ngmax - 1];
    dfr = (frmax - frmin) / (double)(ngmax - 1);

    for (nu = 0; nu < PRJ_NRAD; ++nu) {
        dest[nu] = (prj_table_real *)prj_malloc(opac_bytes * sizeof(**dest));
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
            prj_table_real *restrict out = dest[nu_r];

            for (k = 0; k < nyemax; ++k) {
                for (j = 0; j < ntmax; ++j) {
                    for (i = 0; i < nromax; ++i) {
                        for (ng = 0; ng < PRJ_NEGROUP; ++ng) {
                            double el = rad->log_egroup[nu_r][ng];
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
                                (prj_table_real)((1.0 - dfi) * (double)v0 + dfi * (double)v1);
                        }
                    }
                }
            }
        }
    }

#if defined(PRJ_ENABLE_MPI)
    for (nu = 0; nu < PRJ_NRAD; ++nu) {
        MPI_Bcast(dest[nu], (int)opac_bytes, PRJ_TABLE_MPI_TYPE, 0, MPI_COMM_WORLD);
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
        tmp = (float *)prj_malloc(tmp_count * sizeof(float));
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
        free(rad->log_egroup[nu]);
        free(rad->spec_factor[nu]);
        free(rad->absopac[nu]);
        free(rad->scaopac[nu]);
        free(rad->emis[nu]);
        free(rad->sdelta[nu]);
        rad->egroup[nu] = 0;
        rad->eedge[nu] = 0;
        rad->egroup_erg[nu] = 0;
        rad->degroup_erg[nu] = 0;
        rad->x_e[nu] = 0;
        rad->log_egroup[nu] = 0;
        rad->spec_factor[nu] = 0;
        rad->absopac[nu] = 0;
        rad->scaopac[nu] = 0;
        rad->emis[nu] = 0;
        rad->sdelta[nu] = 0;
    }
    free(rad->enuk);
    free(rad->log_enuk);
    rad->enuk = 0;
    rad->log_enuk = 0;
}

/* Specialized opacity lookup for the radiation implicit solver: always returns
 * both kappa and eta (no scattering/delta, no per-output NULL checks).  The
 * interpolation is identical to prj_rad3_opac_lookup(); this variant just drops
 * the branches the hot path doesn't need. */
void prj_rad3_opac_lookup_ke(const prj_rad *restrict rad, double rho, double temp, double ye,
    double *restrict kappa, double *restrict eta,
    double *restrict dlnkappa_dlnT, double *restrict dlnkappa_dYe,
    double *restrict dlneta_dlnT, double *restrict dlneta_dYe)
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
    double dri;
    double dti;
    double dyei;
    double o_r;
    double o_t;
    double o_y;
    double value_weight[8];
    double dtemp_weight[4];
    double dye_weight[4];
    double inv_dlnT;
    double inv_dYe;
    size_t corner[8];
    double factor;
    int nu;
    int ng;

    if (ye_c > rad->yemax) ye_c = rad->yemax;
    if (ye_c < rad->yemin) ye_c = rad->yemin;
    if (tl_c > rad->log_tmax) tl_c = rad->log_tmax;
    if (tl_c < rad->log_tmin) tl_c = rad->log_tmin;

    deltar = (rl - rad->log_romin) * rad->inv_logrho_span * (double)(nromax - 1);
    deltat = (tl_c - rad->log_tmin) * rad->inv_logtemp_span * (double)(ntmax - 1);
    deltaye = (ye_c - rad->yemin) * rad->inv_ye_span * (double)(nyemax - 1);

    jr = (int)deltar;
    jt = (int)deltat;
    jye = (int)deltaye;
    if (jr < 0) jr = 0;
    if (jr > nromax - 2) jr = nromax - 2;
    if (jt < 0) jt = 0;
    if (jt > ntmax - 2) jt = ntmax - 2;
    if (jye < 0) jye = 0;
    if (jye > nyemax - 2) jye = nyemax - 2;

    dri = deltar - (double)jr;
    if (dri < 0.0) dri = 0.0;
    dti = deltat - (double)jt;
    if (dti < 0.0) dti = 0.0;
    dyei = deltaye - (double)jye;

    o_r = 1.0 - dri;
    o_t = 1.0 - dti;
    o_y = 1.0 - dyei;
    value_weight[0] = o_r * o_t * o_y;
    value_weight[1] = dri * o_t * o_y;
    value_weight[2] = o_r * dti * o_y;
    value_weight[3] = dri * dti * o_y;
    value_weight[4] = o_r * o_t * dyei;
    value_weight[5] = dri * o_t * dyei;
    value_weight[6] = o_r * dti * dyei;
    value_weight[7] = dri * dti * dyei;

    dtemp_weight[0] = o_r * o_y;
    dtemp_weight[1] = dri * o_y;
    dtemp_weight[2] = o_r * dyei;
    dtemp_weight[3] = dri * dyei;
    dye_weight[0] = o_r * o_t;
    dye_weight[1] = dri * o_t;
    dye_weight[2] = o_r * dti;
    dye_weight[3] = dri * dti;
    inv_dlnT = rad->inv_logtemp_span * (double)(ntmax - 1);
    inv_dYe = rad->inv_ye_span * (double)(nyemax - 1);

    corner[0] = OPAC_CELL_IDX(jr,     jt,     jye,     PRJ_NEGROUP, nromax, ntmax, nyemax);
    corner[1] = OPAC_CELL_IDX(jr + 1, jt,     jye,     PRJ_NEGROUP, nromax, ntmax, nyemax);
    corner[2] = OPAC_CELL_IDX(jr,     jt + 1, jye,     PRJ_NEGROUP, nromax, ntmax, nyemax);
    corner[3] = OPAC_CELL_IDX(jr + 1, jt + 1, jye,     PRJ_NEGROUP, nromax, ntmax, nyemax);
    corner[4] = OPAC_CELL_IDX(jr,     jt,     jye + 1, PRJ_NEGROUP, nromax, ntmax, nyemax);
    corner[5] = OPAC_CELL_IDX(jr + 1, jt,     jye + 1, PRJ_NEGROUP, nromax, ntmax, nyemax);
    corner[6] = OPAC_CELL_IDX(jr,     jt + 1, jye + 1, PRJ_NEGROUP, nromax, ntmax, nyemax);
    corner[7] = OPAC_CELL_IDX(jr + 1, jt + 1, jye + 1, PRJ_NEGROUP, nromax, ntmax, nyemax);

    /* Emissivity is divided by RAD_SCALE so it adds radiation energy in the
       internal RAD_SCALE*erg units. */
    factor = 4.0 * M_PI / PRJ_MEV_TO_ERG / RAD_SCALE;

    for (nu = 0; nu < PRJ_NRAD; ++nu) {
        const prj_table_real *restrict absopac = rad->absopac[nu];
        const prj_table_real *restrict emis = rad->emis[nu];
        const double *restrict degroup_erg = rad->degroup_erg[nu];

        for (ng = 0; ng < PRJ_NEGROUP; ++ng) {
            size_t base = (size_t)nu * (size_t)PRJ_NEGROUP + (size_t)ng;
            double k0 = absopac[corner[0] + (size_t)ng];
            double k1 = absopac[corner[1] + (size_t)ng];
            double k2 = absopac[corner[2] + (size_t)ng];
            double k3 = absopac[corner[3] + (size_t)ng];
            double k4 = absopac[corner[4] + (size_t)ng];
            double k5 = absopac[corner[5] + (size_t)ng];
            double k6 = absopac[corner[6] + (size_t)ng];
            double k7 = absopac[corner[7] + (size_t)ng];
            double j0 = emis[corner[0] + (size_t)ng];
            double j1 = emis[corner[1] + (size_t)ng];
            double j2 = emis[corner[2] + (size_t)ng];
            double j3 = emis[corner[3] + (size_t)ng];
            double j4 = emis[corner[4] + (size_t)ng];
            double j5 = emis[corner[5] + (size_t)ng];
            double j6 = emis[corner[6] + (size_t)ng];
            double j7 = emis[corner[7] + (size_t)ng];
            double kv =
                value_weight[0] * k0 + value_weight[1] * k1 +
                value_weight[2] * k2 + value_weight[3] * k3 +
                value_weight[4] * k4 + value_weight[5] * k5 +
                value_weight[6] * k6 + value_weight[7] * k7;
            double jv =
                value_weight[0] * j0 + value_weight[1] * j1 +
                value_weight[2] * j2 + value_weight[3] * j3 +
                value_weight[4] * j4 + value_weight[5] * j5 +
                value_weight[6] * j6 + value_weight[7] * j7;
            double dk_dt =
                dtemp_weight[0] * (k2 - k0) + dtemp_weight[1] * (k3 - k1) +
                dtemp_weight[2] * (k6 - k4) + dtemp_weight[3] * (k7 - k5);
            double dk_dye =
                dye_weight[0] * (k4 - k0) + dye_weight[1] * (k5 - k1) +
                dye_weight[2] * (k6 - k2) + dye_weight[3] * (k7 - k3);
            double dj_dt =
                dtemp_weight[0] * (j2 - j0) + dtemp_weight[1] * (j3 - j1) +
                dtemp_weight[2] * (j6 - j4) + dtemp_weight[3] * (j7 - j5);
            double dj_dye =
                dye_weight[0] * (j4 - j0) + dye_weight[1] * (j5 - j1) +
                dye_weight[2] * (j6 - j2) + dye_weight[3] * (j7 - j3);

            kappa[base] = exp(kv + rl);
            eta[base] = exp(jv + rl) * factor * degroup_erg[ng];
            dlnkappa_dlnT[base] = dk_dt * inv_dlnT;
            dlnkappa_dYe[base] = dk_dye * inv_dYe;
            dlneta_dlnT[base] = dj_dt * inv_dlnT;
            dlneta_dYe[base] = dj_dye * inv_dYe;
        }
    }
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
    size_t corner[8];
    int nu;
    int ng;

    if (ye_c > rad->yemax) ye_c = rad->yemax;
    if (ye_c < rad->yemin) ye_c = rad->yemin;
    if (tl_c > rad->log_tmax) tl_c = rad->log_tmax;
    if (tl_c < rad->log_tmin) tl_c = rad->log_tmin;

    deltar = (rl - rad->log_romin) * rad->inv_logrho_span * (double)(nromax - 1);
    deltat = (tl_c - rad->log_tmin) * rad->inv_logtemp_span * (double)(ntmax - 1);
    deltaye = (ye_c - rad->yemin) * rad->inv_ye_span * (double)(nyemax - 1);

    jr = (int)deltar;
    jt = (int)deltat;
    jye = (int)deltaye;
    if (jr < 0) jr = 0;
    if (jr > nromax - 2) jr = nromax - 2;
    if (jt < 0) jt = 0;
    if (jt > ntmax - 2) jt = ntmax - 2;
    if (jye < 0) jye = 0;
    if (jye > nyemax - 2) jye = nyemax - 2;

    r1i = rad->log_romin + (rad->log_romax - rad->log_romin) * (double)jr / (double)(nromax - 1);
    r2i = rad->log_romin + (rad->log_romax - rad->log_romin) * (double)(jr + 1) / (double)(nromax - 1);
    dri = (rl - r1i) / (r2i - r1i);
    if (dri < 0.0) dri = 0.0;

    t1i = rad->log_tmin + (rad->log_tmax - rad->log_tmin) * (double)jt / (double)(ntmax - 1);
    t2i = rad->log_tmin + (rad->log_tmax - rad->log_tmin) * (double)(jt + 1) / (double)(ntmax - 1);
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

    corner[0] = OPAC_CELL_IDX(jr,     jt,     jye,     PRJ_NEGROUP, nromax, ntmax, nyemax);
    corner[1] = OPAC_CELL_IDX(jr + 1, jt,     jye,     PRJ_NEGROUP, nromax, ntmax, nyemax);
    corner[2] = OPAC_CELL_IDX(jr,     jt + 1, jye,     PRJ_NEGROUP, nromax, ntmax, nyemax);
    corner[3] = OPAC_CELL_IDX(jr + 1, jt + 1, jye,     PRJ_NEGROUP, nromax, ntmax, nyemax);
    corner[4] = OPAC_CELL_IDX(jr,     jt,     jye + 1, PRJ_NEGROUP, nromax, ntmax, nyemax);
    corner[5] = OPAC_CELL_IDX(jr + 1, jt,     jye + 1, PRJ_NEGROUP, nromax, ntmax, nyemax);
    corner[6] = OPAC_CELL_IDX(jr,     jt + 1, jye + 1, PRJ_NEGROUP, nromax, ntmax, nyemax);
    corner[7] = OPAC_CELL_IDX(jr + 1, jt + 1, jye + 1, PRJ_NEGROUP, nromax, ntmax, nyemax);

    double factor;
    if (eta != 0) {
        /* Emissivity is divided by RAD_SCALE so it adds radiation energy in
           the internal RAD_SCALE*erg units. */
        factor = 4.0 * M_PI / PRJ_MEV_TO_ERG / RAD_SCALE;
    }

    for (nu = 0; nu < PRJ_NRAD; ++nu) {
        const prj_table_real *absopac = rad->absopac[nu];
        const prj_table_real *scaopac = rad->scaopac[nu];
        const prj_table_real *emis = rad->emis[nu];
        const prj_table_real *sdelta = rad->sdelta[nu];

        for (ng = 0; ng < PRJ_NEGROUP; ++ng) {
            size_t base = (size_t)nu * (size_t)PRJ_NEGROUP + (size_t)ng;
            double k_val;
            double s_val;
            double d_val;
            double j_val;

#define SAMPLE(ARR) \
            (coff[0] * (ARR)[corner[0] + (size_t)ng] + \
             coff[1] * (ARR)[corner[1] + (size_t)ng] + \
             coff[2] * (ARR)[corner[2] + (size_t)ng] + \
             coff[3] * (ARR)[corner[3] + (size_t)ng] + \
             coff[4] * (ARR)[corner[4] + (size_t)ng] + \
             coff[5] * (ARR)[corner[5] + (size_t)ng] + \
             coff[6] * (ARR)[corner[6] + (size_t)ng] + \
             coff[7] * (ARR)[corner[7] + (size_t)ng])

            if (kappa != 0) {
                k_val = SAMPLE(absopac);
                kappa[base] = exp(k_val+rl);
            }
            if (sigma != 0) {
                s_val = SAMPLE(scaopac);
                sigma[base] = exp(s_val+rl);
            }
            if (delta != 0) {
                d_val = SAMPLE(sdelta);
                delta[base] = d_val;
            }
            if (eta != 0) {
                j_val = SAMPLE(emis);
                eta[base] = exp(j_val+rl) * factor * rad->degroup_erg[nu][ng];
            }

#undef SAMPLE
        }
    }
}

#endif /* PRJ_NRAD > 0 */
