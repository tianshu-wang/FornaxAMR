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

/* Runtime layout: [jeta][jq][m][ke][le].  For a fixed thermodynamic state the
   interpolation fixes (jeta,jq), so an entire (m,ke,le) plane is contiguous.
   This lets prj_rad_eleinel_lookup build the interpolated kernel as a streaming
   weighted sum over the 6 neighboring planes instead of strided gathers. */
#define INEL_ELE_PLANE(ng) (2 * (ng) * (ng))
#define INEL_ELEM_ELE(table, m, ke, le, jeta, jq, ng) \
    (table)[((((jeta) * INEL_PHI_NT + (jq)) * 2 + (m)) * (ng) + (ke)) * (ng) + (le)]

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

/* Returns the lower-left grid cell index (jeta,jq) and the in-cell fractional
   coordinates (sp,sq) for bilinear interpolation. sp,sq are clamped to [0,1] so
   the bilinear weights stay a convex combination -- the interpolant is then
   bounded by its 4 corner values, hence (essentially) non-negative for phi0 and
   sign-consistent for phi1, which lets the detailed-balance loop run without
   sign/positivity gates. Edge extrapolation is intentionally dropped: those are
   the high-eta / low-T degenerate corners where neutrino radiation is
   negligible, so nearest-edge values are sufficient. */
static void prj_rad_eleinel_phi_interp_make(double log10xtemp, double etalep,
    int *jeta_out, int *jq_out, double *sp_out, double *sq_out)
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
    if (sp < 0.0) sp = 0.0; else if (sp > 1.0) sp = 1.0;
    if (sq < 0.0) sq = 0.0; else if (sq > 1.0) sq = 1.0;

    *jeta_out = jeta;
    *jq_out = jq;
    *sp_out = sp;
    *sq_out = sq;
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

static prj_table_real *prj_rad_eleinel_relayout_table(double *old_table, size_t count, int nu)
{
    prj_table_real *new_table = (prj_table_real *)prj_malloc(count * sizeof(*new_table));
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
                            (prj_table_real)INEL_ELEM_ELE_FILE(old_table, m, ke, le, jeta, jq, ng);
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
        double *old_table = (double *)prj_malloc(ee_count * sizeof(*old_table));

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
            double dnue = (rad->degroup_erg[nu] != 0) ?
                rad->degroup_erg[nu][g] / PRJ_MEV_TO_ERG :
                (rad->egroup[nu][gp] - rad->egroup[nu][gm]) / 2.0;

            rad->eleinel_freqe3[idx] = f3;
            rad->eleinel_factf_over_freqe3[idx] = rad->eleinel_factf / f3;
            rad->eleinel_freqe2_dnue[idx] = rad->egroup[nu][g] * rad->egroup[nu][g] * dnue;
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

/* Inner detailed-balance accumulation over all (nfp,g) pairs for one species.
   do_phi1 and want_scatt are passed as constants from the call sites below so
   the compiler prunes the corresponding branches out of this 144-iteration hot
   loop entirely (czero and want_scatt are loop-invariant). phi0/phi1 are the
   contiguous [g*nfreq+nfp] interpolated kernels; xh is flat [g*PRJ_NDIM+d]. */
static inline void prj_rad_eleinel_accumulate(
    int nfreq, int do_phi1, int want_scatt,
    const double *phi0, const double *phi1,
    const double *xj, const double *one_minus_xj, const double *inv_xj,
    const double *xh, const int *active, const double *freqe2_dnue,
    double *sumin, double *sumout, double *ssum)
{
    int nfp;
    int g;

    for (nfp = 0; nfp < nfreq; nfp++) {
        double xh_nfp0 = xh[nfp * PRJ_NDIM + 0];
        double xh_nfp1 = xh[nfp * PRJ_NDIM + 1];
        double xh_nfp2 = xh[nfp * PRJ_NDIM + 2];
        double xjpe = xj[nfp];
        double one_minus_xjpe = 1.0 - xjpe;
        double term = freqe2_dnue[nfp];

        for (g = 0; g < nfreq; g++) {
            double phi0v = phi0[g * nfreq + nfp];
            double phi0_rev = phi0[nfp * nfreq + g];
            double fdotf = xh[g * PRJ_NDIM + 0] * xh_nfp0
                + xh[g * PRJ_NDIM + 1] * xh_nfp1
                + xh[g * PRJ_NDIM + 2] * xh_nfp2;
            double mask = active[g] ? 1.0 : 0.0;
            /* Detailed balance applied directly: q_l(g,g') Phi_l(g,g') =
               Phi_l(g',g), so q cancels and the source uses Phi_rev directly --
               no division, no ratio helpers. No sign/positivity gates: bilinear
               interpolation keeps the kernel bounded by its cell corners (so
               phi0>=0, phi1 sign-consistent), and any underflowed forward entry
               simply contributes ~0 to the sink while the finite reverse entry
               supplies the physical in-scattering source. */
            double pair_source = 0.5 * phi0_rev * xjpe * one_minus_xj[g];
            double pair_sink = 0.5 * phi0v * one_minus_xjpe;

            if (want_scatt) {
                ssum[g] += mask * (term *
                    (0.5 * (phi0v * one_minus_xjpe + phi0_rev * xjpe)));
            }
            if (do_phi1) {
                double phi1v = phi1[g * nfreq + nfp];
                double phi1_rev = phi1[nfp * nfreq + g];
                pair_source -= 1.5 * phi1_rev * fdotf;
                pair_sink -= 1.5 * phi1v * fdotf * inv_xj[g];
            }

            sumin[g] += mask * (term * pair_source);
            sumout[g] += mask * (term * pair_sink);
        }
    }
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
    double rho_cut;
    double sp;
    double sq;
    double w00;
    double w10;
    double w01;
    double w11;
    const size_t total = (size_t)PRJ_NRAD * (size_t)PRJ_NEGROUP;
    int want_scatt = (scatt != 0);

    (void)Ye;

    if (rad == 0 || !rad->eleinel_table_loaded || rho <= rad->min_inel_density) {
        memset(source, 0, total * sizeof(*source));
        memset(sink, 0, total * sizeof(*sink));
        if (want_scatt) {
            memset(scatt, 0, total * sizeof(*scatt));
        }
        return;
    }

    log10t = log10(T);
    prj_rad_eleinel_phi_interp_make(log10t, etael, &jeta, &jq, &sp, &sq);
    /* bilinear weights over the (jeta,jq) cell corners */
    w00 = (1.0 - sp) * (1.0 - sq);
    w10 = sp * (1.0 - sq);
    w01 = (1.0 - sp) * sq;
    w11 = sp * sq;

    rho_cut = 1.0;
    if (rho > 1e13) {
        rho_cut = rho / 1e13;
    }

    for (nu = 0; nu < PRJ_NRAD; nu++) {
        const double *freqe = rad->egroup[nu];
        const prj_table_real *table = rad->eleinel_phi_ee[nu];
        const double *freqe2_dnue = &rad->eleinel_freqe2_dnue[nu * PRJ_NEGROUP];
        const double *factf_over_freqe3 = &rad->eleinel_factf_over_freqe3[nu * PRJ_NEGROUP];
        const double *je_nu = &je[nu * PRJ_NEGROUP];
        const double *he_nu = &he[nu * PRJ_NEGROUP * PRJ_NDIM];
        double xj[PRJ_NEGROUP];
        double one_minus_xj[PRJ_NEGROUP];
        double inv_xj[PRJ_NEGROUP];
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
        double phi0_interp[PRJ_NEGROUP][PRJ_NEGROUP];
        double phi1_interp[PRJ_NEGROUP][PRJ_NEGROUP];

        /* Build the interpolated kernel as a streaming bilinear blend of the 4
           (jeta,jq) cell-corner planes (see INEL_ELEM_ELE layout note).
           phi0_interp and phi1_interp are contiguous [ke][le], matching the m=0
           and m=1 halves of each plane. */
        {
            const int plane = INEL_ELE_PLANE(nfreq);
            const int half = nfreq * nfreq;
            const prj_table_real *pc = table
                + (size_t)(jeta * INEL_PHI_NT + jq) * plane;
            const prj_table_real *c00 = pc;                              /* jeta,   jq   */
            const prj_table_real *c10 = pc + INEL_PHI_NT * plane;        /* jeta+1, jq   */
            const prj_table_real *c01 = pc + plane;                      /* jeta,   jq+1 */
            const prj_table_real *c11 = pc + (INEL_PHI_NT + 1) * plane;  /* jeta+1, jq+1 */
            double *o0 = &phi0_interp[0][0];
            int i;

            for (i = 0; i < half; i++) {
                o0[i] = w00 * c00[i] + w10 * c10[i] + w01 * c01[i] + w11 * c11[i];
            }
            /* phi1 (m=1) is only consumed when czero != 0 (nu != 2); skip its
               build entirely for the heavy-lepton species. */
            if (czero != 0.0) {
                double *o1 = &phi1_interp[0][0];
                for (i = 0; i < half; i++) {
                    int j = half + i;
                    o1[i] = w00 * c00[j] + w10 * c10[j] + w01 * c01[j] + w11 * c11[j];
                }
            }
        }

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
            one_minus_xj[g] = 1.0 - xj[g];
            inv_xj[g] = active[g] ? 1.0 / xj[g] : 0.0;
            sumin[g] = 0.0;
            sumout[g] = 0.0;
            if (want_scatt) {
                ssum[g] = 0.0;
            }
            source[idx] = 0.0;
            sink[idx] = 0.0;
            if (want_scatt) {
                scatt[idx] = 0.0;
            }
        }

        /* Dispatch on the loop-invariant flags with literal constants so the
           hot loop is specialized (czero/want_scatt branches pruned). */
        if (czero != 0.0) {
            if (want_scatt) {
                prj_rad_eleinel_accumulate(nfreq, 1, 1, &phi0_interp[0][0],
                    &phi1_interp[0][0], xj, one_minus_xj, inv_xj, &xh[0][0],
                    active, freqe2_dnue, sumin, sumout, ssum);
            } else {
                prj_rad_eleinel_accumulate(nfreq, 1, 0, &phi0_interp[0][0],
                    &phi1_interp[0][0], xj, one_minus_xj, inv_xj, &xh[0][0],
                    active, freqe2_dnue, sumin, sumout, ssum);
            }
        } else {
            if (want_scatt) {
                prj_rad_eleinel_accumulate(nfreq, 0, 1, &phi0_interp[0][0],
                    0, xj, one_minus_xj, inv_xj, &xh[0][0],
                    active, freqe2_dnue, sumin, sumout, ssum);
            } else {
                prj_rad_eleinel_accumulate(nfreq, 0, 0, &phi0_interp[0][0],
                    0, xj, one_minus_xj, inv_xj, &xh[0][0],
                    active, freqe2_dnue, sumin, sumout, ssum);
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
            if (want_scatt) {
                scatt[idx] = constout * ssum[g] / rho_cut;
            }

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
    double E_pre[PRJ_NRAD * PRJ_NEGROUP];
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
    if (rho <= rad->min_inel_density) return;

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
        source_arr, sink_arr, 0);

    /* Implicit eleinel update; stash the pre-scatter energy per group. */
    for (nu = 0; nu < PRJ_NRAD; nu++) {
        for (g = 0; g < PRJ_NEGROUP; g++) {
            int idx = nu * PRJ_NEGROUP + g;
            double source_phys = source_arr[idx] * rad->degroup_erg[nu][g]
                / (PRJ_MEV_TO_ERG * eta_factor);
            double E_old = u[PRJ_CONS_RAD_E(nu, g)];

            E_pre[idx] = E_old;
            u[PRJ_CONS_RAD_E(nu, g)] = (E_old + dt * source_phys)
                / (1.0 + PRJ_CLIGHT * dt * sink_arr[idx]);
        }
    }

    /* Project the implicit update onto E_g >= 0 before species-number
       normalization.  The normalization factor remains signed so the original
       signed species number is conserved even when n_old < 0. */
    for (nu = 0; nu < PRJ_NRAD; nu++) {
        for (g = 0; g < PRJ_NEGROUP; g++) {
            int eidx = PRJ_CONS_RAD_E(nu, g);

            if (u[eidx] < 0.0) {
                u[eidx] = 0.0;
            }
        }
    }

    /* eleinel is pure scattering: it must conserve each species' total neutrino
       number. The implicit update drifts it by ~1e-4; restore exact conservation
       by rescaling each species back to its pre-scatter number N = sum_g E/erg_g.
       Since x_e[nu][g] ~ +-1/erg_g, dN=0 per species => dYe=0 exactly. */
    for (nu = 0; nu < PRJ_NRAD; nu++) {
        double n_old = 0.0;
        double n_new = 0.0;
        double scale;

        for (g = 0; g < PRJ_NEGROUP; g++) {
            int idx = nu * PRJ_NEGROUP + g;
            double inv_erg = 1.0 / rad->egroup_erg[nu][g];
            n_old += E_pre[idx] * inv_erg;
            n_new += u[PRJ_CONS_RAD_E(nu, g)] * inv_erg;
        }
        scale = (n_new != 0.0) ? n_old / n_new : 1.0;
        for (g = 0; g < PRJ_NEGROUP; g++) {
            u[PRJ_CONS_RAD_E(nu, g)] *= scale;
        }
    }

    /* Energy/Ye exchange computed from the rescaled (final) radiation field so
       total energy stays conserved; dy is now ~roundoff. */
    du = 0.0;
    dy = 0.0;
    for (nu = 0; nu < PRJ_NRAD; nu++) {
        for (g = 0; g < PRJ_NEGROUP; g++) {
            int idx = nu * PRJ_NEGROUP + g;
            double dE = u[PRJ_CONS_RAD_E(nu, g)] - E_pre[idx];

            du += dE;
            dy += rad->x_e[nu][g] * dE;
        }
    }

    /* du is a change in RAD_SCALE*erg units; multiply back to erg for the gas. */
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

#if PRJ_USE_RADIATION_FSA
static void prj_rad_inel_fsa_build_m1_tmp(const prj_rad *rad,
    const double *u, double *u_tmp)
{
    int v;
    int nu;
    int g;
    int angle;
    int d;

    /* The M1 inelastic step (prj_rad_eleinel_step/prj_rad_nucinel_step) reads
     * only the hydro slots and the per-group E/F moment slots, and apply_m1_tmp
     * reads only u_tmp's E and Ye slots.  The reconstruction below writes every
     * E/F slot, so copying only the hydro block is bit-identical and skips the
     * dead angular-intensity slots of u_tmp. */
    for (v = 0; v < PRJ_NHYDRO; ++v) {
        u_tmp[v] = u[v];
    }

    for (nu = 0; nu < PRJ_NRAD; ++nu) {
        for (g = 0; g < PRJ_NEGROUP; ++g) {
            double E = 0.0;
            double F[PRJ_NDIM] = {0.0, 0.0, 0.0};

            for (angle = 0; angle < PRJ_NANGLE; ++angle) {
                int iv = PRJ_CONS_RAD_I(nu, g, angle);
                double J = u[iv];

                E += J;
                for (d = 0; d < PRJ_NDIM; ++d) {
                    F[d] += PRJ_CLIGHT * J * rad->n0[angle][d];
                }
            }

            u_tmp[PRJ_CONS_RAD_E(nu, g)] = E;
            u_tmp[PRJ_CONS_RAD_F1(nu, g)] = F[0];
            u_tmp[PRJ_CONS_RAD_F2(nu, g)] = F[1];
            u_tmp[PRJ_CONS_RAD_F3(nu, g)] = F[2];
        }
    }
}

static void prj_rad_inel_fsa_apply_m1_tmp(const prj_rad *rad,
    double *u, const double *u_tmp)
{
    const double four_pi = 4.0 * M_PI;
    double rho = u[PRJ_CONS_RHO];
    double e_unchanged;
    double dmom[PRJ_NDIM] = {0.0, 0.0, 0.0};
    int nu;
    int g;
    int angle;
    int d;

    /* The eleinel/nucinel steps already applied the matter energy and Ye
       exchange to u_tmp (internal energy in u_tmp[ETOT], Ye in u_tmp[YE]).  Take
       the internal+magnetic energy from there (it uses the unchanged momentum,
       since the steps never touch it) and re-add only the kinetic energy after
       the radiation-momentum back-reaction below. */
    e_unchanged = u_tmp[PRJ_CONS_ETOT] - 0.5 *
        (u[PRJ_CONS_MOM1] * u[PRJ_CONS_MOM1] +
         u[PRJ_CONS_MOM2] * u[PRJ_CONS_MOM2] +
         u[PRJ_CONS_MOM3] * u[PRJ_CONS_MOM3]) / rho;

    for (nu = 0; nu < PRJ_NRAD; ++nu) {
        for (g = 0; g < PRJ_NEGROUP; ++g) {
            double E_old = 0.0;
            double E_new = u_tmp[PRJ_CONS_RAD_E(nu, g)];
            double scale = 0.0;

            /* Reset any unphysical negative angular intensities to zero before
               rescaling, so the group total and the rescale factor are built
               from non-negative intensities only. */
            for (angle = 0; angle < PRJ_NANGLE; ++angle) {
                int iv = PRJ_CONS_RAD_I(nu, g, angle);

                if (u[iv] < 0.0) {
                    u[iv] = 0.0;
                }
                E_old += u[iv];
            }
            /* Rescale the angular intensities so the group's total (and hence its
               neutrino number, since all angles share the group energy) matches
               the post-eleinel/nucinel moment E_new, keeping the angular shape. */
            if (E_old != 0.0) {
                scale = E_new / E_old;
            }

            for (angle = 0; angle < PRJ_NANGLE; ++angle) {
                int iv = PRJ_CONS_RAD_I(nu, g, angle);
                double J_old = u[iv];
                double J_new;

                if (E_old != 0.0) {
                    J_new = J_old * scale;
                } else if (E_new != 0.0) {
                    J_new = E_new * rad->solid_angle[angle] / four_pi;
                } else {
                    J_new = 0.0;
                }

                u[iv] = J_new;
                for (d = 0; d < PRJ_NDIM; ++d) {
                    dmom[d] += (J_old - J_new) * rad->n0[angle][d] / PRJ_CLIGHT;
                }
            }
        }
    }

    /* Ye already updated on u_tmp by the inelastic steps.  Rescaling J changes
       the radiation flux, so the gas absorbs the momentum change; rebuild the
       total energy as unchanged (internal+magnetic) energy plus the new kinetic
       energy.  The corresponding energy exchange is already in e_unchanged. */
    u[PRJ_CONS_YE] = u_tmp[PRJ_CONS_YE];
    u[PRJ_CONS_MOM1] += dmom[0] * RAD_SCALE;
    u[PRJ_CONS_MOM2] += dmom[1] * RAD_SCALE;
    u[PRJ_CONS_MOM3] += dmom[2] * RAD_SCALE;
    u[PRJ_CONS_ETOT] = e_unchanged + 0.5 *
        (u[PRJ_CONS_MOM1] * u[PRJ_CONS_MOM1] +
         u[PRJ_CONS_MOM2] * u[PRJ_CONS_MOM2] +
         u[PRJ_CONS_MOM3] * u[PRJ_CONS_MOM3]) / rho;
}

/* Combined FSA inelastic update.  The J<->(E,F) conversion is done once:
   J -> (E,F) -> eleinel -> E1 -> nucinel -> E2 -> J, instead of once per
   process.  eleinel_step and nucinel_step run in-place on the M1 moments in
   u_tmp (nucinel sees eleinel's updated state), then the single apply projects
   the combined moment change back onto the angular intensities. */
int prj_rad_inel_fsa(prj_rad *rad, prj_eos *eos, double *u, double dt, double T_cell)
{
    double u_tmp[PRJ_NVAR_CONS];
    int status;

    prj_rad_inel_fsa_build_m1_tmp(rad, u, u_tmp);
    prj_rad_eleinel_step(rad, eos, u_tmp, dt, T_cell);
    status = prj_rad_nucinel_step(rad, eos, u_tmp, dt, T_cell);
    prj_rad_inel_fsa_apply_m1_tmp(rad, u, u_tmp);

    return status;
}
#endif

#endif
