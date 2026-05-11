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
#define INEL_PHI_NT 30

#define INEL_ELEM_ELE(table, m, ke, le, jeta, jq, ng) \
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

static void prj_rad_eleinel_phifind_interp(const prj_rad *rad, int nutype,
    int nf, int nfp, double log10xtemp, double etalep,
    double *phi0e, double *phi1e)
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
    double coeff0;
    double coeff1;
    double coeff2;
    double coeff3;
    double coeff4;
    double coeff5;
    const double *table = rad->eleinel_phi_ee[nutype];
    int ng = PRJ_NEGROUP;

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

    coeff0 = 0.5 * sq * (sq - 1.0);
    coeff1 = 0.5 * sp * (sp - 1.0);
    coeff2 = 1.0 + sp * sq - sp * sp - sq * sq;
    coeff3 = 0.5 * sp * (sp - 2.0 * sq + 1.0);
    coeff4 = 0.5 * sq * (sq - 2.0 * sp + 1.0);
    coeff5 = sp * sq;

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

static void prj_rad_eleinel_resourcesink(const prj_rad *rad,
    int nutype, int nfreq, int nf,
    const double *je, const double *he, const double *freqe,
    double T, double etael,
    double *srce, double *sinke, double *scatte)
{
    double factf = rad->eleinel_factf;
    double constin = rad->eleinel_constin;
    double constout = rad->eleinel_constout;
    double cut;
    double czero;
    double xje;
    double xhe[PRJ_NDIM];
    double sumin;
    double sumout;
    double ssum;
    int nfp;
    int id;
    double enu;
    double log10t;
    int nu_offset = nutype * PRJ_NEGROUP;
    double factf_over_freqe3_nf = rad->eleinel_factf_over_freqe3[nu_offset + nf];

    cut = 1.0;
    czero = 1.0;
    if (nutype == 2) {
        cut = 4.0;
        czero = 0.0;
    }

    xje = PRJ_MAX(0.0, PRJ_MIN(je[nf] * factf_over_freqe3_nf / cut, 1.0));
    for (id = 0; id < PRJ_NDIM; id++) {
        xhe[id] = he[nf * PRJ_NDIM + id] * factf_over_freqe3_nf / cut;
    }

    sumin = 0.0;
    sumout = 0.0;
    ssum = 0.0;
    log10t = log10(T);
    enu = freqe[nf];

    for (nfp = 0; nfp < nfreq; nfp++) {
        int nfpp = PRJ_MIN(nfp + 1, nfreq - 1);
        int nfpm = PRJ_MAX(nfp - 1, 0);
        double dnue = (freqe[nfpp] - freqe[nfpm]) / 2.0;
        double factf_over_freqe3_nfp = rad->eleinel_factf_over_freqe3[nu_offset + nfp];
        double xjpe = PRJ_MAX(0.0, PRJ_MIN(je[nfp] * factf_over_freqe3_nfp / cut, 1.0));
        double fdotf = 0.0;
        double omegae;
        double expe;
        double phi0;
        double phi1;
        double term;
        double phi0ee;
        double phi1ee;
        int did;

        for (did = 0; did < PRJ_NDIM; did++) {
            double xhpe = he[nfp * PRJ_NDIM + did] * factf_over_freqe3_nfp / cut;
            fdotf += xhe[did] * xhpe;
        }
        omegae = enu - freqe[nfp];
        expe = exp(PRJ_MAX(PRJ_MIN(-omegae / T, 207.0), -207.0));
        prj_rad_eleinel_phifind_interp(rad, nutype, nf, nfp, log10t, etael, &phi0ee, &phi1ee);
        phi0 = phi0ee;
        phi1 = phi1ee;
        term = freqe[nfp] * freqe[nfp] * dnue;
        sumin += term * expe * (0.5 * phi0 * xjpe * (1.0 - xje) - czero * 3.0 * 0.5 * phi1 * fdotf);
        sumout += term * (0.5 * phi0 * (1.0 - xjpe) - czero * 3.0 * 0.5 * phi1 * fdotf / xje);
        ssum += term * (0.5 * phi0 * (1.0 - xjpe + expe * xjpe));
    }

    sumin *= cut;
    *srce = constin * (enu * enu * enu) * sumin;
    *sinke = constout * sumout;
    *scatte = constout * ssum;
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
        rad->eleinel_phi_ee[nu] = (double *)malloc(ee_count * sizeof(double));
        if (rad->eleinel_phi_ee[nu] == 0) {
            fprintf(stderr, "prj_rad_inel: allocation failed for species %d\n", nu);
            exit(1);
        }
    }

    for (nu = 0; nu < PRJ_NRAD; nu++) {
        prj_rad_eleinel_read_table(rad, nu, rad->eleinel_phi_ee[nu], ee_count);
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
        rad->eleinel_constin = fourpi * (bigG * bigG) / (hc2pi * hc2pi * hc2pi * hc2pi * hc2pi * hc2pi) * 1.60217733e-6;
        rad->eleinel_constout = fourpi * (bigG * bigG) / (hc2pi * hc2pi * hc2pi) / clt;
    }

    for (nu = 0; nu < PRJ_NRAD; nu++) {
        int g;
        for (g = 0; g < PRJ_NEGROUP; g++) {
            int idx = nu * PRJ_NEGROUP + g;
            double f3 = rad->egroup[nu][g] * rad->egroup[nu][g] * rad->egroup[nu][g];
            rad->eleinel_freqe3[idx] = f3;
            rad->eleinel_factf_over_freqe3[idx] = rad->eleinel_factf / f3;
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

    (void)Ye;

    if (rad == 0 || !rad->eleinel_table_loaded) {
        int total = PRJ_NRAD * PRJ_NEGROUP;
        for (g = 0; g < total; g++) {
            source[g] = 0.0;
            sink[g] = 0.0;
            scatt[g] = 0.0;
        }
        return;
    }

    for (nu = 0; nu < PRJ_NRAD; nu++) {
        const double *freqe = rad->egroup[nu];
        int nfreq = PRJ_NEGROUP;

        for (g = 0; g < PRJ_NEGROUP; g++) {
            int idx = nu * PRJ_NEGROUP + g;
            double xxj = je[idx];
            double cut;

            if (rho > rad->min_inel_density && xxj > 0.0) {
                prj_rad_eleinel_resourcesink(rad, nu, nfreq, g,
                    &je[nu * PRJ_NEGROUP], &he[nu * PRJ_NEGROUP * PRJ_NDIM],
                    freqe, T, etael,
                    &source[idx], &sink[idx], &scatt[idx]);
            } else {
                source[idx] = 0.0;
                sink[idx] = 0.0;
                scatt[idx] = 0.0;
            }

            cut = 1.0;
            if (rho > 1e13)
                cut = rho / 1e13;
            source[idx] /= cut;
            sink[idx] /= cut;
            scatt[idx] /= cut;

            source[idx] = PRJ_MAX(source[idx], 0.0);
            sink[idx] = PRJ_MAX(sink[idx], 0.0);
        }
    }
}

void prj_rad_eleinel_step(prj_rad *rad, prj_eos *eos, double *u, double dt)
{
    double KE;
    double Emag = 0.0;
    double Uint_old;
    double rho;
    double Ye;
    double eint;
    double T;
    double eos_q[PRJ_EOS_NQUANT];
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
    eint = Uint_old / rho;

    prj_eos_rey(eos, rho, eint, Ye, eos_q);
    T = eos_q[PRJ_EOS_TEMPERATURE];

    etael = prj_eos_rty_geteta(eos, rho, T, Ye);
    if (etael < -20.0) etael = -20.0;

    for (nu = 0; nu < PRJ_NRAD; nu++) {
        for (g = 0; g < PRJ_NEGROUP; g++) {
            int idx = nu * PRJ_NEGROUP + g;
            double E_g = u[PRJ_CONS_RAD_E(nu, g)];
            double he_mag = 0.0;

            je[idx] = PRJ_CLIGHT * E_g * eta_factor * PRJ_MEV_TO_ERG / rad->degroup_erg[nu][g];
            if (je[idx] < 0.0) je[idx] = 0.0;

            for (d = 0; d < PRJ_NDIM; d++) {
                int fidx = idx * PRJ_NDIM + d;
                double F_gd;
                switch (d) {
                case 0: F_gd = u[PRJ_CONS_RAD_F1(nu, g)]; break;
                case 1: F_gd = u[PRJ_CONS_RAD_F2(nu, g)]; break;
                default: F_gd = u[PRJ_CONS_RAD_F3(nu, g)]; break;
                }
                he[fidx] = F_gd * eta_factor * PRJ_MEV_TO_ERG / rad->degroup_erg[nu][g];
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

    prj_rad_eleinel_lookup(rad, rho, T, Ye, etael, je, he,
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

    Uint_new = Uint_old - du;
    u[PRJ_CONS_ETOT] = Uint_new + KE + Emag;
    u[PRJ_CONS_YE] += dy;

    Ye = u[PRJ_CONS_YE] / rho;
    eint = Uint_new / rho;
    prj_eos_rey(eos, rho, eint, Ye, eos_q);
}

/* ================================================================== */
/* Nucleon inelastic scattering (Kompaneets solver)                   */
/* ================================================================== */

static double prj_nucinel_compute_coeff(double kT, double rho_N, double Ye)
{
    double kom_hbar = 1.0545718e-34;
    double kom_c = 299792458.0;
    double kom_e = 1.60217662e-19;
    double kom_m_n = 939.0;
    double kom_Vp = -0.5 * (1.0 - 4.0 * 0.23122);
    double kom_Ap = -1.2723 / 2.0;
    double kom_Vn = -0.5;
    double kom_An = 1.2723 / 2.0;
    double kom_e_unit = 1e6 * 1.60217662e-19;
    double kom_t_unit = 1.0545718e-34 / (1e6 * 1.60217662e-19);
    double kom_l_unit = 1.0545718e-34 / (1e6 * 1.60217662e-19) * 299792458.0;
    double kom_G2 = 1.327817e-22;
    double kom_m = 939.0;
    double kom_prot = (-0.5 * (1.0 - 4.0 * 0.23122)) * (-0.5 * (1.0 - 4.0 * 0.23122))
        + 5.0 * (-1.2723 / 2.0) * (-1.2723 / 2.0);
    double kom_neut = (-0.5) * (-0.5) + 5.0 * (1.2723 / 2.0) * (1.2723 / 2.0);

    double beta = 1.0 / kT;
    double nN = rho_N * 1e3 * kom_c * kom_c / (kom_m * kom_e_unit)
        * kom_l_unit * kom_l_unit * kom_l_unit;
    double coeff_p;
    double coeff_n;
    double coeff;

    (void)kom_hbar;
    (void)kom_e;
    (void)kom_m_n;
    (void)kom_Vp;
    (void)kom_Ap;
    (void)kom_Vn;
    (void)kom_An;
    (void)kom_e_unit;
    (void)kom_t_unit;
    (void)kom_l_unit;

    coeff = 2.0 * kom_G2 * nN / (3.0 * M_PI * beta * beta * beta * kom_m) / kom_t_unit;
    coeff_p = kom_prot * Ye;
    coeff_n = kom_neut * (1.0 - Ye);
    coeff *= (coeff_p + coeff_n);
    return coeff;
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
    coeff = prj_nucinel_compute_coeff(kT, rho_N, Ye) / dlogxs;

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

int prj_rad_nucinel_step(prj_rad *rad, prj_eos *eos, double *u, double dt)
{
    double KE;
    double Emag = 0.0;
    double Uint_old;
    double rho;
    double Ye;
    double eint;
    double T;
    double eos_q[PRJ_EOS_NQUANT];
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
    eint = Uint_old / rho;

    prj_eos_rey(eos, rho, eint, Ye, eos_q);
    T = eos_q[PRJ_EOS_TEMPERATURE];

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
        xs1[g] = rad->egroup[0][g] / T;
        Js1[g] = u[PRJ_CONS_RAD_E(0, g)] * sf;
        if (Js1[g] > 1.0) Js1[g] = 1.0;
        if (Js1[g] < 0.0) Js1[g] = 0.0;
        u_res[0 * PRJ_NEGROUP + g] = u[PRJ_CONS_RAD_E(0, g)] - Js1[g] / sf;

        sf = rad->spec_factor[1][g];
        if (rad->egroup[1][g] > rad->kom_Ecut[1] && ncut2 == PRJ_NEGROUP) {
            ncut2 = g;
        }
        xs2[g] = rad->egroup[1][g] / T;
        Js2[g] = u[PRJ_CONS_RAD_E(1, g)] * sf;
        if (Js2[g] > 1.0) Js2[g] = 1.0;
        if (Js2[g] < 0.0) Js2[g] = 0.0;
        u_res[1 * PRJ_NEGROUP + g] = u[PRJ_CONS_RAD_E(1, g)] - Js2[g] / sf;

        sf = rad->spec_factor[2][g];
        if (rad->egroup[2][g] > rad->kom_Ecut[2] && ncut3 == PRJ_NEGROUP) {
            ncut3 = g;
        }
        xs3[g] = rad->egroup[2][g] / T;
        Js3[g] = u[PRJ_CONS_RAD_E(2, g)] * sf;
        if (Js3[g] > 1.0) Js3[g] = 1.0;
        if (Js3[g] < 0.0) Js3[g] = 0.0;
        u_res[2 * PRJ_NEGROUP + g] = u[PRJ_CONS_RAD_E(2, g)] - Js3[g] / sf;
    }

    status = prj_nucinel_compute_step(rad, T, rho, Ye, xs1, Js1, ncut1, dt);
    if (status == 1) {
        status = prj_nucinel_compute_step(rad, T, rho, Ye, xs2, Js2, ncut2, dt);
    }
    if (status == 1) {
        status = prj_nucinel_compute_step(rad, T, rho, Ye, xs3, Js3, ncut3, dt);
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

        Uint_new = Uint_old - du;
        u[PRJ_CONS_ETOT] = Uint_new + KE + Emag;
        u[PRJ_CONS_YE] += dy;

        Ye = u[PRJ_CONS_YE] / rho;
        eint = Uint_new / rho;
        prj_eos_rey(eos, rho, eint, Ye, eos_q);
    }

    return status;
}

#endif
