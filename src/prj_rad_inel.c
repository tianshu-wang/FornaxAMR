#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "prj.h"
#include "prj_rad_inel.h"

#if PRJ_NRAD > 0

#if defined(PRJ_ENABLE_MPI)
#include <mpi.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ================================================================== */
/* Electron inelastic scattering (phi-table based)                    */
/* ================================================================== */

#define INEL_PHI_NETA 30
#ifndef INEL_PHI_NT
#define INEL_PHI_NT 30
#endif

#define INEL_ELEM_ELE(table, m, ke, le, jeta, jq, ng) \
    (table)[(((((le) * INEL_PHI_NETA + (jeta)) * INEL_PHI_NT + (jq)) * 2 + (m)) * (ng) + (ke))]

#define INEL_ELEM_ELE_FILE(table, m, ke, le, jeta, jq, ng) \
    (table)[((((m) * (ng) + (ke)) * (ng) + (le)) * INEL_PHI_NETA + (jeta)) * INEL_PHI_NT + (jq)]

static int prj_rad_inel_mpi_rank(void)
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

static void prj_rad_eleinel_phi_interp_make(double log10xtemp, double etalep,
    int *jeta_out, int *jq_out,
    double *coeff0_out, double *coeff1_out, double *coeff2_out,
    double *coeff3_out, double *coeff4_out, double *coeff5_out)
{
    double t1 = -1.0;
    double t2 = 1.5;
    double t12 = -1.0;
    double t22 = 1.5;
    double eta1 = 0.0;
    double eta2 = 50.0;
    double etal = etalep;
    double tl = log10xtemp;
    double alpha = t1 + (etal - eta1) / (eta2 - eta1) * (t12 - t1);
    double beta = t2 - t1 + ((t22 - t12) - (t2 - t1)) * (etal - eta1) / (eta2 - eta1);
    double ql = (tl - alpha) / beta;
    double delta = (etal - eta1) / (eta2 - eta1) * (double)INEL_PHI_NETA;
    int jeta = (int)delta;
    int jq = (int)((double)INEL_PHI_NT * ql);
    double sp;
    double sq;

    if (jeta < 1)
        jeta = 1;
    if (jq < 1)
        jq = 1;
    if (jeta > INEL_PHI_NETA - 2)
        jeta = INEL_PHI_NETA - 2;
    if (jq > INEL_PHI_NT - 2)
        jq = INEL_PHI_NT - 2;

    sp = delta - jeta;
    sq = (double)INEL_PHI_NT * ql - jq;

    *jeta_out = jeta;
    *jq_out = jq;
    *coeff0_out = 0.5 * sq * (sq - 1.0);
    *coeff1_out = 0.5 * sp * (sp - 1.0);
    *coeff2_out = 1.0 + sp * sq - sp * sp - sq * sq;
    *coeff3_out = 0.5 * sp * (sp - 2.0 * sq + 1.0);
    *coeff4_out = 0.5 * sq * (sq - 2.0 * sp + 1.0);
    *coeff5_out = sp * sq;
}

static inline void prj_rad_eleinel_phifind_interp(const prj_rad *rad, int nutype,
    int nf, int nfp, int jeta, int jq,
    double coeff0, double coeff1, double coeff2,
    double coeff3, double coeff4, double coeff5,
    double *phi0e, double *phi1e)
{
    const double *table = rad->eleinel_phi_ee[nutype];
    int ng = PRJ_NEGROUP;

    *phi0e = coeff0 * INEL_ELEM_ELE(table, 0, nf, nfp, jeta, jq - 1, ng)
        + coeff1 * INEL_ELEM_ELE(table, 0, nf, nfp, jeta - 1, jq, ng)
        + coeff2 * INEL_ELEM_ELE(table, 0, nf, nfp, jeta, jq, ng)
        + coeff3 * INEL_ELEM_ELE(table, 0, nf, nfp, jeta + 1, jq, ng)
        + coeff4 * INEL_ELEM_ELE(table, 0, nf, nfp, jeta, jq + 1, ng)
        + coeff5 * INEL_ELEM_ELE(table, 0, nf, nfp, jeta + 1, jq + 1, ng);
    *phi1e = coeff0 * INEL_ELEM_ELE(table, 1, nf, nfp, jeta, jq - 1, ng)
        + coeff1 * INEL_ELEM_ELE(table, 1, nf, nfp, jeta - 1, jq, ng)
        + coeff2 * INEL_ELEM_ELE(table, 1, nf, nfp, jeta, jq, ng)
        + coeff3 * INEL_ELEM_ELE(table, 1, nf, nfp, jeta + 1, jq, ng)
        + coeff4 * INEL_ELEM_ELE(table, 1, nf, nfp, jeta, jq + 1, ng)
        + coeff5 * INEL_ELEM_ELE(table, 1, nf, nfp, jeta + 1, jq + 1, ng);
}

static void prj_rad_eleinel_read_table(const prj_rad *rad, int nu,
    double *dest, size_t count)
{
    int rank = prj_rad_inel_mpi_rank();
    char species_char[] = {'e', 'a', 'm'};
    char filename[PRJ_PATH_MAX];
    char nstr[4];
    long irecl;

    snprintf(nstr, sizeof(nstr), "%03d", PRJ_NEGROUP);
    snprintf(filename, sizeof(filename), "%sphi_%ce%s_030x030_a.dat",
        rad->eleinel_table_dir, species_char[nu], nstr);

    irecl = (long)(2 * 8 * PRJ_NEGROUP * PRJ_NEGROUP * INEL_PHI_NETA * INEL_PHI_NT);

    if (rank == 0) {
        FILE *fp = fopen(filename, "rb");

        if (fp == NULL) {
            fprintf(stderr, "prj_rad_inel: cannot open %s\n", filename);
            exit(1);
        }
        fseek(fp, irecl, SEEK_SET);
        if (fread(dest, sizeof(double), count, fp) != count) {
            fprintf(stderr, "prj_rad_inel: fread failed for %s\n", filename);
            fclose(fp);
            exit(1);
        }
        fclose(fp);
    }

#if defined(PRJ_ENABLE_MPI)
    MPI_Bcast(dest, (int)count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
}

static double *prj_rad_eleinel_relayout_table(double *old_table, size_t count, int nu)
{
    double *new_table = (double *)malloc(count * sizeof(*new_table));
    int ng = PRJ_NEGROUP;
    int m;
    int ke;
    int le;
    int jeta;
    int jq;

    if (new_table == 0) {
        fprintf(stderr, "prj_rad_inel: relayout allocation failed for species %d\n", nu);
        free(old_table);
        exit(1);
    }

    for (m = 0; m < 2; m++) {
        for (ke = 0; ke < ng; ke++) {
            for (le = 0; le < ng; le++) {
                for (jeta = 0; jeta < INEL_PHI_NETA; jeta++) {
                    for (jq = 0; jq < INEL_PHI_NT; jq++) {
                        INEL_ELEM_ELE(new_table, m, ke, le, jeta, jq, ng) =
                            INEL_ELEM_ELE_FILE(old_table, m, ke, le, jeta, jq, ng);
                    }
                }
            }
        }
    }

    free(old_table);
    return new_table;
}

void prj_rad_eleinel_init(prj_rad *rad)
{
    int rank = prj_rad_inel_mpi_rank();
    int nu;
    size_t ee_count;

    if (rad == 0) {
        fprintf(stderr, "prj_rad_inel: null rad\n");
        exit(1);
    }
    if (rad->eleinel_table_dir[0] == '\0') {
        rad->eleinel_table_loaded = 0;
        return;
    }

    ee_count = (size_t)2 * (size_t)PRJ_NEGROUP * (size_t)PRJ_NEGROUP
        * (size_t)INEL_PHI_NETA * (size_t)INEL_PHI_NT;

    for (nu = 0; nu < PRJ_NRAD; nu++) {
        double *old_table = (double *)malloc(ee_count * sizeof(*old_table));

        if (old_table == 0) {
            fprintf(stderr, "prj_rad_inel: allocation failed for species %d\n", nu);
            exit(1);
        }
        prj_rad_eleinel_read_table(rad, nu, old_table, ee_count);
        rad->eleinel_phi_ee[nu] = prj_rad_eleinel_relayout_table(old_table, ee_count, nu);
    }

    if (rank == 0) {
        fprintf(stderr, "prj_rad_inel: loaded electron scattering tables from %s\n",
            rad->eleinel_table_dir);
    }
    rad->eleinel_table_loaded = 1;

    {
        double hbar = 6.582122e-22;
        double bigG = 3.937e-17;
        double clt = 2.99792458e10;
        double fourpi = 12.56637061e0;
        double pi = 3.1415926535898e0;
        double hc2pi = 2.0 * pi * hbar * clt;

        rad->eleinel_factf = (hc2pi * hc2pi * hc2pi) / clt / 1.60217733e-6;
        /* constin scales the emission (source) term; dividing by RAD_SCALE
           expresses it in the internal RAD_SCALE*erg units.  constout drives
           the sink/scatter rates, which are scale-invariant and unchanged. */
        rad->eleinel_constin = fourpi * (bigG * bigG) / (hc2pi * hc2pi * hc2pi * hc2pi * hc2pi * hc2pi) * 1.60217733e-6 / RAD_SCALE;
        rad->eleinel_constout = fourpi * (bigG * bigG) / (hc2pi * hc2pi * hc2pi) / clt;
    }

    for (nu = 0; nu < PRJ_NRAD; nu++) {
        int g;
        for (g = 0; g < PRJ_NEGROUP; g++) {
            int idx = nu * PRJ_NEGROUP + g;
            double f3 = rad->egroup[nu][g] * rad->egroup[nu][g] * rad->egroup[nu][g];
            int gp = PRJ_MIN(g + 1, PRJ_NEGROUP - 1);
            int gm = PRJ_MAX(g - 1, 0);
            double dnue = (rad->egroup[nu][gp] - rad->egroup[nu][gm]) / 2.0;

            rad->eleinel_freqe3[idx] = f3;
            rad->eleinel_factf_over_freqe3[idx] = rad->eleinel_factf / f3;
            rad->eleinel_freqe2_dnue[idx] = rad->egroup[nu][g] * rad->egroup[nu][g] * dnue;
        }
    }

    for (nu = 0; nu < PRJ_NRAD; nu++) {
        int nf;
        for (nf = 0; nf < PRJ_NEGROUP; nf++) {
            int nfp;
            for (nfp = 0; nfp < PRJ_NEGROUP; nfp++) {
                double omegae = rad->egroup[nu][nf] - rad->egroup[nu][nfp];
                int jq;
                for (jq = 0; jq < INEL_PHI_NT; jq++) {
                    double log10t_jq = -1.0 + 2.5 * (double)jq / (double)INEL_PHI_NT;
                    double T_jq = pow(10.0, log10t_jq);
                    rad->expe[nu][(nf * PRJ_NEGROUP + nfp) * INEL_PHI_NT + jq] =
                        exp(PRJ_MAX(PRJ_MIN(-omegae / T_jq, 207.0), -207.0));
                }
            }
        }
    }
}

void prj_rad_eleinel_free(prj_rad *rad)
{
    int nu;

    if (rad == 0 || !rad->eleinel_table_loaded) {
        return;
    }
    for (nu = 0; nu < PRJ_NRAD; nu++) {
        free(rad->eleinel_phi_ee[nu]);
        rad->eleinel_phi_ee[nu] = 0;
    }
    rad->eleinel_table_loaded = 0;
}

void prj_rad_eleinel_lookup(const prj_rad *rad,
    double rho, double T, double Ye,
    double etael,
    const double *je, const double *he,
    double *source, double *sink, double *scatt)
{
    int nu;
    int g;
    int jeta;
    int jq;
    double log10t;
    int expe_jq;
    double expe_coeff0;
    double expe_coeff1;
    double rho_cut;
    double coeff0;
    double coeff1;
    double coeff2;
    double coeff3;
    double coeff4;
    double coeff5;
    const size_t total = (size_t)PRJ_NRAD * (size_t)PRJ_NEGROUP;

    (void)Ye;

    if (rad == 0 || !rad->eleinel_table_loaded || rho <= rad->min_inel_density) {
        memset(source, 0, total * sizeof(*source));
        memset(sink, 0, total * sizeof(*sink));
        memset(scatt, 0, total * sizeof(*scatt));
        return;
    }

    log10t = log10(T);
    prj_rad_eleinel_phi_interp_make(log10t, etael,
        &jeta, &jq, &coeff0, &coeff1, &coeff2, &coeff3, &coeff4, &coeff5);
    {
        double ql = (log10t + 1.0) / 2.5;
        double sq;

        expe_jq = (int)((double)INEL_PHI_NT * ql);
        if (expe_jq < 0) expe_jq = 0;
        if (expe_jq > INEL_PHI_NT - 2) expe_jq = INEL_PHI_NT - 2;
        sq = (double)INEL_PHI_NT * ql - (double)expe_jq;
        expe_coeff0 = 1.0 - sq;
        expe_coeff1 = sq;
    }

    rho_cut = 1.0;
    if (rho > 1e13) {
        rho_cut = rho / 1e13;
    }

    for (nu = 0; nu < PRJ_NRAD; nu++) {
        const double *freqe = rad->egroup[nu];
        const double *freqe2_dnue = &rad->eleinel_freqe2_dnue[nu * PRJ_NEGROUP];
        const double *factf_over_freqe3 = &rad->eleinel_factf_over_freqe3[nu * PRJ_NEGROUP];
        const double *je_nu = &je[nu * PRJ_NEGROUP];
        const double *he_nu = &he[nu * PRJ_NEGROUP * PRJ_NDIM];
        double xj[PRJ_NEGROUP];
        double xh[PRJ_NEGROUP][PRJ_NDIM];
        double sumin[PRJ_NEGROUP];
        double sumout[PRJ_NEGROUP];
        double ssum[PRJ_NEGROUP];
        int active[PRJ_NEGROUP];
        double species_cut = (nu == 2) ? 4.0 : 1.0;
        double czero = (nu == 2) ? 0.0 : 1.0;
        double constin = rad->eleinel_constin;
        double constout = rad->eleinel_constout;
        int nfreq = PRJ_NEGROUP;
        int nfp;

        for (g = 0; g < PRJ_NEGROUP; g++) {
            double fac = factf_over_freqe3[g] / species_cut;
            int d;

            xj[g] = PRJ_MAX(0.0, PRJ_MIN(je_nu[g] * fac, 1.0));
            for (d = 0; d < PRJ_NDIM; d++) {
                xh[g][d] = he_nu[g * PRJ_NDIM + d] * fac;
            }
        }

        for (g = 0; g < PRJ_NEGROUP; g++) {
            int idx = nu * PRJ_NEGROUP + g;
            active[g] = (je_nu[g] > 0.0);
            sumin[g] = 0.0;
            sumout[g] = 0.0;
            ssum[g] = 0.0;
            source[idx] = 0.0;
            sink[idx] = 0.0;
            scatt[idx] = 0.0;
        }

        for (nfp = 0; nfp < nfreq; nfp++) {
            const double *xh_nfp = xh[nfp];
            double xjpe = xj[nfp];
            double one_minus_xjpe = 1.0 - xjpe;
            double term = freqe2_dnue[nfp];

            for (g = 0; g < nfreq; g++) {
                const double *xh_g;
                double xje;
                double fdotf;
                double expe;
                double phi0;
                double phi1;

                if (!active[g]) {
                    continue;
                }

                xh_g = xh[g];
                xje = xj[g];
                fdotf = xh_g[0] * xh_nfp[0] + xh_g[1] * xh_nfp[1] + xh_g[2] * xh_nfp[2];
                expe = expe_coeff0 * rad->expe[nu][(g * PRJ_NEGROUP + nfp) * INEL_PHI_NT + expe_jq]
                     + expe_coeff1 * rad->expe[nu][(g * PRJ_NEGROUP + nfp) * INEL_PHI_NT + expe_jq + 1];
                prj_rad_eleinel_phifind_interp(rad, nu, g, nfp,
                    jeta, jq, coeff0, coeff1, coeff2, coeff3, coeff4, coeff5,
                    &phi0, &phi1);

                sumin[g] += term * expe *
                    (0.5 * phi0 * xjpe * (1.0 - xje) - czero * 1.5 * phi1 * fdotf);
                sumout[g] += term *
                    (0.5 * phi0 * one_minus_xjpe - czero * 1.5 * phi1 * fdotf / xje);
                ssum[g] += term * (0.5 * phi0 * (one_minus_xjpe + expe * xjpe));
            }
        }

        for (g = 0; g < PRJ_NEGROUP; g++) {
            int idx = nu * PRJ_NEGROUP + g;
            double enu;
            double enu3;

            if (!active[g]) {
                continue;
            }

            enu = freqe[g];
            enu3 = enu * enu * enu;
            source[idx] = constin * enu3 * species_cut * sumin[g] / rho_cut;
            sink[idx] = constout * sumout[g] / rho_cut;
            scatt[idx] = constout * ssum[g] / rho_cut;

            source[idx] = PRJ_MAX(source[idx], 0.0);
            sink[idx] = PRJ_MAX(sink[idx], 0.0);
        }
    }
}

void prj_rad_eleinel_step(prj_rad *rad, prj_eos *eos, double *u, double dt, double T_cell)
{
    double KE;
    double Emag = 0.0;
    double Uint_old;
    double rho;
    double Ye;
    double je[PRJ_NRAD * PRJ_NEGROUP];
    double he[PRJ_NRAD * PRJ_NEGROUP * PRJ_NDIM];
    double source_arr[PRJ_NRAD * PRJ_NEGROUP];
    double sink_arr[PRJ_NRAD * PRJ_NEGROUP];
    double scatt_arr[PRJ_NRAD * PRJ_NEGROUP];
    double eta_factor = 1.0 / (4.0 * M_PI);
    double etael;
    double du;
    double dy;
    double Uint_new;
    int nu;
    int g;
    int d;

    if (!rad->eleinel_table_loaded) return;

    rho = u[PRJ_CONS_RHO];
    KE = 0.5 * (u[PRJ_CONS_MOM1] * u[PRJ_CONS_MOM1] +
        u[PRJ_CONS_MOM2] * u[PRJ_CONS_MOM2] +
        u[PRJ_CONS_MOM3] * u[PRJ_CONS_MOM3]) / rho;
#if PRJ_MHD
    Emag = 0.5 * (u[PRJ_CONS_B1] * u[PRJ_CONS_B1] +
        u[PRJ_CONS_B2] * u[PRJ_CONS_B2] +
        u[PRJ_CONS_B3] * u[PRJ_CONS_B3]);
#endif
    Uint_old = u[PRJ_CONS_ETOT] - KE - Emag;
    Ye = u[PRJ_CONS_YE] / rho;

    etael = prj_eos_rty_geteta(eos, rho, T_cell, Ye, PRJ_EOS_CTX_MAIN);
    if (etael < -20.0) etael = -20.0;

    for (nu = 0; nu < PRJ_NRAD; nu++) {
        for (g = 0; g < PRJ_NEGROUP; g++) {
            int idx = nu * PRJ_NEGROUP + g;
            double E_g = u[PRJ_CONS_RAD_E(nu, g)];
            double he_mag = 0.0;

            je[idx] = PRJ_CLIGHT * E_g * RAD_SCALE * eta_factor * PRJ_MEV_TO_ERG / rad->degroup_erg[nu][g];
            if (je[idx] < 0.0) je[idx] = 0.0;

            for (d = 0; d < PRJ_NDIM; d++) {
                int fidx = idx * PRJ_NDIM + d;
                double F_gd;
                switch (d) {
                case 0: F_gd = u[PRJ_CONS_RAD_F1(nu, g)]; break;
                case 1: F_gd = u[PRJ_CONS_RAD_F2(nu, g)]; break;
                default: F_gd = u[PRJ_CONS_RAD_F3(nu, g)]; break;
                }
                he[fidx] = F_gd * RAD_SCALE * eta_factor * PRJ_MEV_TO_ERG / rad->degroup_erg[nu][g];
                he_mag += (he[fidx] / (je[idx] + 1.0e-15))
                    * (he[fidx] / (je[idx] + 1.0e-15));
            }
            he_mag = sqrt(he_mag);
            if (he_mag > 1.0) {
                for (d = 0; d < PRJ_NDIM; d++) {
                    he[idx * PRJ_NDIM + d] /= he_mag;
                }
            }
        }
    }

    prj_rad_eleinel_lookup(rad, rho, T_cell, Ye, etael, je, he,
        source_arr, sink_arr, scatt_arr);

    du = 0.0;
    dy = 0.0;
    for (nu = 0; nu < PRJ_NRAD; nu++) {
        for (g = 0; g < PRJ_NEGROUP; g++) {
            int idx = nu * PRJ_NEGROUP + g;
            double source_phys = source_arr[idx] * rad->degroup_erg[nu][g]
                / (PRJ_MEV_TO_ERG * eta_factor);
            double E_old = u[PRJ_CONS_RAD_E(nu, g)];
            double E_new = (E_old + dt * source_phys)
                / (1.0 + PRJ_CLIGHT * dt * sink_arr[idx]);

            u[PRJ_CONS_RAD_E(nu, g)] = E_new;
            du += E_new - E_old;
            dy += rad->x_e[nu][g] * (E_new - E_old);
        }
    }

    /* du is a change in RAD_SCALE*erg units; multiply back to erg for the gas.
       dy already carries RAD_SCALE through x_e. */
    Uint_new = Uint_old - du * RAD_SCALE;
    u[PRJ_CONS_ETOT] = Uint_new + KE + Emag;
    u[PRJ_CONS_YE] += dy;
}

/* ================================================================== */
/* Nucleon inelastic scattering (Kompaneets solver)                   */
/* ================================================================== */

static double prj_nucinel_compute_coeff(const prj_rad *rad, double kT, double rho_N, double Ye)
{
    double beta = 1.0 / kT;
    double beta3 = beta * beta * beta;
    double nN = rad->kom_nucinel_const * rho_N / beta3;
    return nN * (rad->kom_nucinel_prot * Ye + rad->kom_nucinel_neut * (1.0 - Ye));
}

static int prj_nucinel_tdma(double *A, double *B, double *C, double *D, int n)
{
    double temp;
    int i;

    C[0] = C[0] / B[0];
    D[0] = D[0] / B[0];
    for (i = 1; i < n; i++) {
        temp = B[i] - A[i] * C[i - 1];
        if (temp == 0.0) {
            printf("TDMA unstable!\n");
            return 0;
        }
        temp = 1.0 / temp;
        C[i] = C[i] * temp;
        D[i] = (D[i] - A[i] * D[i - 1]) * temp;
    }
    D[n - 1] = D[n - 1];
    for (i = n - 2; i >= 0; i--) {
        D[i] -= C[i] * D[i + 1];
        D[i] = D[i];
    }
    return 1;
}

static int prj_nucinel_compute_step(const prj_rad *rad,
    double kT, double rho, double Ye,
    double *xs, double *Js, int ncut, double dt)
{
    int status = 1;
    double totx3J = 0.0;
    double totx3J_temp;
    double x3s[PRJ_NEGROUP];
    double x3J[PRJ_NEGROUP];
    double TDMA_A[PRJ_NEGROUP];
    double TDMA_B[PRJ_NEGROUP];
    double TDMA_C[PRJ_NEGROUP];
    double TDMA_D[PRJ_NEGROUP];
    double TDMA_D1[PRJ_NEGROUP];
    double exph[PRJ_NEGROUP + 1];
    double coeff;
    double ratio;
    double dlogxs;
    double remain_dt;
    double current_dt;
    int method = 2;
    int i;
    double rho_N = rho;

    if (rho_N > rad->kom_rhocut) {
        rho_N = rad->kom_rhocut;
    }

    for (i = 0; i < ncut; i++) {
        x3s[i] = xs[i] * xs[i] * xs[i];
        x3J[i] = x3s[i] * PRJ_MAX(0.0, PRJ_MIN(1.0, Js[i]));
        totx3J += x3J[i];
    }
    dlogxs = log(xs[1]) - log(xs[0]);
    coeff = prj_nucinel_compute_coeff(rad, kT, rho_N, Ye) / dlogxs;

    ratio = xs[1] / xs[0];
    exph[0] = 1.0 / (exp(xs[0] * (1.0 - 1.0 / ratio)) - 1.0);
    for (i = 0; i < ncut; i++) {
        exph[i + 1] = 1.0 / (exp(xs[i] * (ratio - 1.0)) - 1.0);
        TDMA_D[i] = x3J[i];
    }

    remain_dt = dt;

    while (remain_dt > 0.0) {
        double dydt_val;

        current_dt = remain_dt;

        for (i = 0; i < ncut; i++) {
            TDMA_A[i] = 0.0;
            TDMA_B[i] = 0.0;
            TDMA_C[i] = 0.0;
            if (i == 0) {
                TDMA_B[i] = -coeff * (0.0 + 0.0 - x3s[i + 1] * exph[i + 1] - 0.0);
                TDMA_C[i] = -coeff * (x3s[i] - TDMA_D[i] + x3s[i] * exph[i + 1]);
            } else if (i == ncut - 1) {
                TDMA_A[i] = -coeff * x3s[i] * exph[i];
                TDMA_B[i] = -coeff * (-x3s[i - 1] + TDMA_D[i - 1] - 0.0 - x3s[i - 1] * exph[i]);
            } else {
                TDMA_A[i] = -coeff * x3s[i] * exph[i];
                TDMA_B[i] = -coeff * (-x3s[i - 1] + TDMA_D[i - 1] - x3s[i + 1] * exph[i + 1]
                    - x3s[i - 1] * exph[i]);
                TDMA_C[i] = -coeff * (x3s[i] - TDMA_D[i] + x3s[i] * exph[i + 1]);
            }
        }

        for (i = 0; i < ncut; i++) {
            if (i == 0) {
                dydt_val = -(TDMA_B[i] * TDMA_D[i] + TDMA_C[i] * TDMA_D[i + 1]);
            } else if (i == ncut - 1) {
                dydt_val = -(TDMA_B[i] * TDMA_D[i] + TDMA_A[i] * TDMA_D[i - 1]);
            } else {
                dydt_val = -(TDMA_B[i] * TDMA_D[i] + TDMA_A[i] * TDMA_D[i - 1]
                    + TDMA_C[i] * TDMA_D[i + 1]);
            }
            if (dydt_val < 0.0 && dydt_val * current_dt < -rad->kom_epsilon * (1e-10 + TDMA_D[i])) {
                current_dt = -rad->kom_epsilon * (1e-10 + TDMA_D[i]) / dydt_val;
            }
            if (dydt_val > 0.0 && dydt_val * current_dt > rad->kom_epsilon * (x3s[i] - TDMA_D[i] + 1e-10)) {
                current_dt = rad->kom_epsilon * (x3s[i] - TDMA_D[i] + 1e-10) / dydt_val;
            }
        }
        if (current_dt < PRJ_MIN(rad->kom_dtmin, remain_dt)) {
            printf("Kom time step too small!\n");
            current_dt = PRJ_MIN(rad->kom_dtmin, remain_dt);
        }

        if (method == 0) {
            for (i = 0; i < ncut; i++) {
                TDMA_A[i] *= current_dt;
                TDMA_B[i] *= current_dt;
                TDMA_C[i] *= current_dt;
                TDMA_D1[i] = TDMA_D[i];
            }
            for (i = 0; i < ncut; i++) {
                if (i == 0) {
                    TDMA_D[i] -= (TDMA_B[i] * TDMA_D1[i] + TDMA_C[i] * TDMA_D1[i + 1]);
                } else if (i == ncut - 1) {
                    TDMA_D[i] -= (TDMA_A[i] * TDMA_D1[i - 1] + TDMA_B[i] * TDMA_D1[i]);
                } else {
                    TDMA_D[i] -= (TDMA_A[i] * TDMA_D1[i - 1] + TDMA_B[i] * TDMA_D1[i]
                        + TDMA_C[i] * TDMA_D1[i + 1]);
                }
            }
        } else if (method == 1) {
            for (i = 0; i < ncut; i++) {
                TDMA_A[i] *= current_dt;
                TDMA_B[i] *= current_dt;
                TDMA_C[i] *= current_dt;
            }
            for (i = 0; i < ncut; i++) {
                TDMA_B[i] += 1.0;
            }
            status = prj_nucinel_tdma(TDMA_A, TDMA_B, TDMA_C, TDMA_D, ncut);
        } else {
            for (i = 0; i < ncut; i++) {
                TDMA_A[i] *= 0.5 * current_dt;
                TDMA_B[i] *= 0.5 * current_dt;
                TDMA_C[i] *= 0.5 * current_dt;
                TDMA_D1[i] = TDMA_D[i];
            }
            for (i = 0; i < ncut; i++) {
                if (i == 0) {
                    TDMA_D[i] -= (TDMA_B[i] * TDMA_D1[i] + TDMA_C[i] * TDMA_D1[i + 1]);
                } else if (i == ncut - 1) {
                    TDMA_D[i] -= (TDMA_A[i] * TDMA_D1[i - 1] + TDMA_B[i] * TDMA_D1[i]);
                } else {
                    TDMA_D[i] -= (TDMA_A[i] * TDMA_D1[i - 1] + TDMA_B[i] * TDMA_D1[i]
                        + TDMA_C[i] * TDMA_D1[i + 1]);
                }
            }
            for (i = 0; i < ncut; i++) {
                TDMA_B[i] += 1.0;
            }
            status = prj_nucinel_tdma(TDMA_A, TDMA_B, TDMA_C, TDMA_D, ncut);
        }

        totx3J_temp = 0.0;
        for (i = 0; i < ncut; i++) {
            TDMA_D[i] = PRJ_MAX(0.0, PRJ_MIN(x3s[i], TDMA_D[i]));
            totx3J_temp += TDMA_D[i];
        }
        if (fabs(totx3J_temp - totx3J) > rad->kom_delta * totx3J) {
            printf("Neutrino number changes too much! Current dt = %e\n", current_dt);
            status = 0;
        }

        if (status == 0) {
            break;
        }
        remain_dt -= current_dt;
    }

    if (status == 1) {
        for (i = 0; i < ncut; i++) {
            Js[i] = TDMA_D[i] / x3s[i];
        }
    }

    return status;
}

int prj_rad_nucinel_step(prj_rad *rad, prj_eos *eos, double *u, double dt, double T_cell)
{
    double KE;
    double Emag = 0.0;
    double Uint_old;
    double rho;
    double Ye;
    double xs1[PRJ_NEGROUP];
    double xs2[PRJ_NEGROUP];
    double xs3[PRJ_NEGROUP];
    double Js1[PRJ_NEGROUP];
    double Js2[PRJ_NEGROUP];
    double Js3[PRJ_NEGROUP];
    double u_res[PRJ_NRAD * PRJ_NEGROUP];
    double E_old[PRJ_NRAD * PRJ_NEGROUP];
    int ncut1 = PRJ_NEGROUP;
    int ncut2 = PRJ_NEGROUP;
    int ncut3 = PRJ_NEGROUP;
    int status;
    int nu;
    int g;
    double du;
    double dy;
    double Uint_new;

    rho = u[PRJ_CONS_RHO];
    if (rho < rad->min_inel_density) {
        return 1;
    }
    KE = 0.5 * (u[PRJ_CONS_MOM1] * u[PRJ_CONS_MOM1] +
        u[PRJ_CONS_MOM2] * u[PRJ_CONS_MOM2] +
        u[PRJ_CONS_MOM3] * u[PRJ_CONS_MOM3]) / rho;
#if PRJ_MHD
    Emag = 0.5 * (u[PRJ_CONS_B1] * u[PRJ_CONS_B1] +
        u[PRJ_CONS_B2] * u[PRJ_CONS_B2] +
        u[PRJ_CONS_B3] * u[PRJ_CONS_B3]);
#endif
    Uint_old = u[PRJ_CONS_ETOT] - KE - Emag;
    Ye = u[PRJ_CONS_YE] / rho;

    for (nu = 0; nu < PRJ_NRAD; ++nu) {
        for (g = 0; g < PRJ_NEGROUP; ++g) {
            E_old[nu * PRJ_NEGROUP + g] = u[PRJ_CONS_RAD_E(nu, g)];
        }
    }

    for (g = 0; g < PRJ_NEGROUP; ++g) {
        double sf;

        sf = rad->spec_factor[0][g];
        if (rad->egroup[0][g] > rad->kom_Ecut[0] && ncut1 == PRJ_NEGROUP) {
            ncut1 = g;
        }
        xs1[g] = rad->egroup[0][g] / T_cell;
        Js1[g] = u[PRJ_CONS_RAD_E(0, g)] * sf;
        if (Js1[g] > 1.0) Js1[g] = 1.0;
        if (Js1[g] < 0.0) Js1[g] = 0.0;
        u_res[0 * PRJ_NEGROUP + g] = u[PRJ_CONS_RAD_E(0, g)] - Js1[g] / sf;

        sf = rad->spec_factor[1][g];
        if (rad->egroup[1][g] > rad->kom_Ecut[1] && ncut2 == PRJ_NEGROUP) {
            ncut2 = g;
        }
        xs2[g] = rad->egroup[1][g] / T_cell;
        Js2[g] = u[PRJ_CONS_RAD_E(1, g)] * sf;
        if (Js2[g] > 1.0) Js2[g] = 1.0;
        if (Js2[g] < 0.0) Js2[g] = 0.0;
        u_res[1 * PRJ_NEGROUP + g] = u[PRJ_CONS_RAD_E(1, g)] - Js2[g] / sf;

        sf = rad->spec_factor[2][g];
        if (rad->egroup[2][g] > rad->kom_Ecut[2] && ncut3 == PRJ_NEGROUP) {
            ncut3 = g;
        }
        xs3[g] = rad->egroup[2][g] / T_cell;
        Js3[g] = u[PRJ_CONS_RAD_E(2, g)] * sf;
        if (Js3[g] > 1.0) Js3[g] = 1.0;
        if (Js3[g] < 0.0) Js3[g] = 0.0;
        u_res[2 * PRJ_NEGROUP + g] = u[PRJ_CONS_RAD_E(2, g)] - Js3[g] / sf;
    }

    status = prj_nucinel_compute_step(rad, T_cell, rho, Ye, xs1, Js1, ncut1, dt);
    if (status == 1) {
        status = prj_nucinel_compute_step(rad, T_cell, rho, Ye, xs2, Js2, ncut2, dt);
    }
    if (status == 1) {
        status = prj_nucinel_compute_step(rad, T_cell, rho, Ye, xs3, Js3, ncut3, dt);
    }

    if (status == 1) {
        du = 0.0;
        dy = 0.0;
        for (g = 0; g < PRJ_NEGROUP; ++g) {
            double E_new_g;

            E_new_g = u_res[0 * PRJ_NEGROUP + g] + Js1[g] / rad->spec_factor[0][g];
            u[PRJ_CONS_RAD_E(0, g)] = E_new_g;
            du += E_new_g - E_old[0 * PRJ_NEGROUP + g];
            dy += rad->x_e[0][g] * (E_new_g - E_old[0 * PRJ_NEGROUP + g]);

            E_new_g = u_res[1 * PRJ_NEGROUP + g] + Js2[g] / rad->spec_factor[1][g];
            u[PRJ_CONS_RAD_E(1, g)] = E_new_g;
            du += E_new_g - E_old[1 * PRJ_NEGROUP + g];
            dy += rad->x_e[1][g] * (E_new_g - E_old[1 * PRJ_NEGROUP + g]);

            E_new_g = u_res[2 * PRJ_NEGROUP + g] + Js3[g] / rad->spec_factor[2][g];
            u[PRJ_CONS_RAD_E(2, g)] = E_new_g;
            du += E_new_g - E_old[2 * PRJ_NEGROUP + g];
            dy += rad->x_e[2][g] * (E_new_g - E_old[2 * PRJ_NEGROUP + g]);
        }

        /* du is in RAD_SCALE*erg units; multiply back to erg for the gas.
           dy already carries RAD_SCALE through x_e. */
        Uint_new = Uint_old - du * RAD_SCALE;
        u[PRJ_CONS_ETOT] = Uint_new + KE + Emag;
        u[PRJ_CONS_YE] += dy;
    }

    return status;
}

#endif
