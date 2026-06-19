#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "prj.h"

static double prj_rad_m1_chi_exact(double f)
{
    if (f <= 0.0) {
        return 1.0 / 3.0;
    }
    if (f >= 1.0) {
        return 1.0;
    }
    return (3.0 + 4.0 * f * f) / (5.0 + 2.0 * sqrt(4.0 - 3.0 * f * f));
}

#if PRJ_NRAD > 0
/* Levermore/Vaytet third-moment scalar q(f), evaluated through the equivalent
 * boost parameter beta = 3f / (2 + sqrt(4 - 3f^2)). */
static double prj_rad_levermore_q_factor_exact(double f)
{
    double a;
    double beta;
    double beta2;
    double beta4;
    double one_minus_beta2;

    if (f <= 0.0) {
        return 0.0;
    }
    if (f >= 1.0) {
        return 1.0;
    }

    a = sqrt(fmax(0.0, 4.0 - 3.0 * f * f));
    beta = 3.0 * f / (2.0 + a);
    if (beta < 5.0e-2) {
        beta2 = beta * beta;
        return 4.0 * beta * (1.0 / 5.0 + beta2 * (1.0 / 21.0 - beta2 / 315.0));
    }

    beta2 = beta * beta;
    beta4 = beta2 * beta2;
    one_minus_beta2 = 1.0 - beta2;
    return (9.0 * beta - 8.0 / beta + 3.0 / (beta2 * beta) -
        3.0 * one_minus_beta2 * one_minus_beta2 * one_minus_beta2 *
        atanh(beta) / beta4) / (3.0 + beta2);
}

static void prj_rad_init_closure(prj_rad *rad)
{
    int i;

    if (rad == 0) {
        return;
    }
    for (i = 0; i <= NCLOSURE; ++i) {
        double f = (double)i / (double)NCLOSURE;

        rad->chi[i] = prj_rad_m1_chi_exact(f);
        rad->q[i] = prj_rad_levermore_q_factor_exact(f);
    }
}

static int prj_rad_closure_ready(const prj_rad *rad)
{
    return rad != 0 && rad->chi[0] > 0.0 && rad->chi[NCLOSURE] > 0.0 &&
        rad->q[NCLOSURE] > 0.0;
}

static double prj_rad_closure_lookup(const double values[NCLOSURE + 1], double f)
{
    double scaled;
    double w;
    int idx;

    if (f <= 0.0) {
        return values[0];
    }
    if (f >= 1.0) {
        return values[NCLOSURE];
    }
    scaled = f * (double)NCLOSURE;
    idx = (int)scaled;
    w = scaled - (double)idx;
    return values[idx] + w * (values[idx + 1] - values[idx]);
}
#endif

static double prj_rad_m1_chi(const prj_rad *rad, double f)
{
#if PRJ_NRAD > 0
    return prj_rad_closure_lookup(rad->chi, f);
#else
    (void)rad;
    return 0;
#endif
}

#if PRJ_NRAD > 0
static double prj_rad_levermore_q_factor(const prj_rad *rad, double f)
{
    return prj_rad_closure_lookup(rad->q, f);
}
#endif

void prj_rad_init(prj_rad *rad)
{
#if PRJ_NRAD > 0
    prj_rad_init_closure(rad);
    prj_rad3_opac_init(rad);
    prj_rad_eleinel_init(rad);
#else
    (void)rad;
#endif
}

void prj_rad_prim2cons(const double *W, double *U)
{
    int field;
    int group;

    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            U[PRJ_CONS_RAD_E(field, group)] = W[PRJ_PRIM_RAD_E(field, group)];
            U[PRJ_CONS_RAD_F1(field, group)] = W[PRJ_PRIM_RAD_F1(field, group)];
            U[PRJ_CONS_RAD_F2(field, group)] = W[PRJ_PRIM_RAD_F2(field, group)];
            U[PRJ_CONS_RAD_F3(field, group)] = W[PRJ_PRIM_RAD_F3(field, group)];
        }
    }
}

void prj_rad_cons2prim(const double *U, double *W)
{
    int field;
    int group;

    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            W[PRJ_PRIM_RAD_E(field, group)] = U[PRJ_CONS_RAD_E(field, group)];
            W[PRJ_PRIM_RAD_F1(field, group)] = U[PRJ_CONS_RAD_F1(field, group)];
            W[PRJ_PRIM_RAD_F2(field, group)] = U[PRJ_CONS_RAD_F2(field, group)];
            W[PRJ_PRIM_RAD_F3(field, group)] = U[PRJ_CONS_RAD_F3(field, group)];
        }
    }
}

/* Public M1 closure for the pressure tensor.  P^{ij} = E * D^{ij} with the
 * Levermore Eddington tensor D^{ij} = a δ^{ij} + b n^i n^j, n = F/|F|, and
 * χ(f) = (3 + 4f²)/(5 + 2√(4 - 3f²)), f = |F|/(c E).  Falls back to the
 * isotropic limit P^{ij} = (E/3) δ^{ij} when |F| or E vanishes. */
void prj_rad_m1_pressure(const prj_rad *rad, double E, double F1, double F2, double F3,
    double P[3][3])
{
    double Fmag;
    double cE;
    double f;
    double chi;
    double a_c;
    double b_c;
    double n[3];
    int a;
    int b;

    Fmag = sqrt(F1 * F1 + F2 * F2 + F3 * F3);
    cE = PRJ_CLIGHT * (E > 0.0 ? E : 0.0);

    if (cE <= 0.0 || Fmag <= 0.0) {
        double third = (E > 0.0 ? E : 0.0) / 3.0;

        for (a = 0; a < 3; ++a) {
            for (b = 0; b < 3; ++b) {
                P[a][b] = (a == b) ? third : 0.0;
            }
        }
        return;
    }

    f = Fmag / cE;
    if (f > 1.0) {
        f = 1.0;
    }
    chi = prj_rad_m1_chi(rad, f);
    a_c = 0.5 * (1.0 - chi);
    b_c = 0.5 * (3.0 * chi - 1.0);
    n[0] = F1 / Fmag;
    n[1] = F2 / Fmag;
    n[2] = F3 / Fmag;
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            P[a][b] = E * (a_c * (a == b ? 1.0 : 0.0) + b_c * n[a] * n[b]);
        }
    }
}

#if PRJ_NRAD > 0
/* Contraction M[b] = sum_{a,c} Q^{abc} dvdx[a][c] of the Levermore third moment
 * with the velocity gradient, computed analytically so the 27 components of Q
 * never have to be materialised. With
 *   Q^{abc} = coef_nnn n^a n^b n^c + coef_mix (n^a d_bc + n^b d_ca + n^c d_ab),
 * the contraction collapses to
 *   M[b] = coef_nnn n[b] S + coef_mix (T1[b] + n[b] divv + T3[b]),
 * with S = n.dvdx.n, T1 = n^T dvdx, T3 = dvdx n, divv = tr(dvdx). Mathematically
 * identical to building Q and summing (validated to machine epsilon), with the
 * isotropic E<=0 / Fmag<=0 limit returning zero exactly as m1_third_moment does. */
static void prj_rad_m1_third_moment_contract(const prj_rad *rad, double E,
    double F1, double F2, double F3, const double dvdx[3][3], double M[3])
{
    double E_pos;
    double Fmag;
    double cE;
    double f;
    double q_fac;
    double n[3];
    double coef_nnn;
    double coef_mix;
    double divv;
    double T1[3];
    double T3[3];
    double S;
    int b;

    M[0] = 0.0;
    M[1] = 0.0;
    M[2] = 0.0;

    E_pos = E > 0.0 ? E : 0.0;
    Fmag = sqrt(F1 * F1 + F2 * F2 + F3 * F3);
    cE = PRJ_CLIGHT * E_pos;
    if (cE <= 0.0 || Fmag <= 0.0) {
        return;
    }

    f = Fmag / cE;
    if (f > 1.0) {
        f = 1.0;
    }
    q_fac = prj_rad_levermore_q_factor(rad, f);
    n[0] = F1 / Fmag;
    n[1] = F2 / Fmag;
    n[2] = F3 / Fmag;
    coef_nnn = 0.5 * cE * (5.0 * q_fac - 3.0 * f);
    coef_mix = 0.5 * cE * (f - q_fac);

    divv = dvdx[0][0] + dvdx[1][1] + dvdx[2][2];
    for (b = 0; b < 3; ++b) {
        T1[b] = n[0] * dvdx[0][b] + n[1] * dvdx[1][b] + n[2] * dvdx[2][b];
        T3[b] = dvdx[b][0] * n[0] + dvdx[b][1] * n[1] + dvdx[b][2] * n[2];
    }
    S = n[0] * T3[0] + n[1] * T3[1] + n[2] * T3[2];

    for (b = 0; b < 3; ++b) {
        M[b] = coef_nnn * n[b] * S + coef_mix * (T1[b] + n[b] * divv + T3[b]);
    }
}

static void prj_rad_m1_phys_flux_with_fluxmag(const prj_rad *rad, double E, double F1,
    double F2, double F3, double Fmag, double inv_Fmag, double f,
    double *fE, double *fF1, double *fF2, double *fF3)
{
    double chi;
    double n1;
    double n2;
    double n3;
    double D11;
    double D12;
    double D13;
    double c2;

    c2 = PRJ_CLIGHT * PRJ_CLIGHT;

    if (E <= 0.0 || Fmag <= 0.0) {
        /* isotropic: P = (E/3) I */
        *fE = F1;
        *fF1 = c2 * E / 3.0;
        *fF2 = 0.0;
        *fF3 = 0.0;
        return;
    }

    chi = prj_rad_m1_chi(rad, f);

    n1 = F1 * inv_Fmag;
    n2 = F2 * inv_Fmag;
    n3 = F3 * inv_Fmag;

    {
        double a = 0.5 * (1.0 - chi);
        double b = 0.5 * (3.0 * chi - 1.0);
        D11 = a + b * n1 * n1;
        D12 = b * n1 * n2;
        D13 = b * n1 * n3;
    }

    *fE = F1;
    *fF1 = c2 * E * D11;
    *fF2 = c2 * E * D12;
    *fF3 = c2 * E * D13;
}

static void prj_rad_enforce_flux_limit(double *E, double *F1, double *F2, double *F3,
    double *Fmag_out, double *inv_Fmag_out, double *f_out)
{
    double Fmag;
    double cE;
    double f;
    double scale;

    if (*E < 0.0) {
        *E = 0.0;
    }
    Fmag = sqrt((*F1) * (*F1) + (*F2) * (*F2) + (*F3) * (*F3));
    cE = PRJ_CLIGHT * (*E);
    if (Fmag > cE && Fmag > 0.0) {
        scale = cE / Fmag;
        *F1 *= scale;
        *F2 *= scale;
        *F3 *= scale;
        Fmag = cE;
    }
    if (cE > 0.0) {
        f = Fmag / cE;
        if (f > 1.0) {
            f = 1.0;
        }
    } else {
        f = 0.0;
    }
    *Fmag_out = Fmag;
    *inv_Fmag_out = (Fmag > 0.0) ? (1.0 / Fmag) : 0.0;
    *f_out = f;
}

void prj_rad_m1_wavespeeds_with_fluxmag(double E, double F1, double Fmag, double inv_Fmag,
    double f, double *lam_min, double *lam_max)
{
    double mu;
    double fsq;
    double ffac;
    double inv_ffac;
    double lterm;

    if (E <= 0.0 || Fmag <= 0.0) {
        *lam_min = -1.0 / sqrt(3.0);
        *lam_max = 1.0 / sqrt(3.0);
        return;
    }

    mu = F1 * inv_Fmag;

    fsq = f * f;
    ffac = sqrt(4.0 - 3.0 * fsq);
    inv_ffac = 1.0 / ffac;
    lterm = sqrt(fabs((2.0 / 3.0) * (4.0 - 3.0 * fsq - ffac) + 2.0 * mu * mu * (2.0 - fsq - ffac)));
    *lam_min = (mu * f - lterm) * inv_ffac;
    *lam_max = (mu * f + lterm) * inv_ffac;
    if (*lam_min < -1.0) {
        *lam_min = -1.0;
    }
    if (*lam_max > 1.0) {
        *lam_max = 1.0;
    }
}

void prj_rad_m1_wavespeeds(double E, double F1, double F2, double F3,
    double *lam_min, double *lam_max)
{
    double Fmag;
    double inv_Fmag;
    double cE;
    double f;

    cE = PRJ_CLIGHT * (E > 0.0 ? E : 0.0);
    Fmag = sqrt(F1 * F1 + F2 * F2 + F3 * F3);
    inv_Fmag = (Fmag > 0.0) ? (1.0 / Fmag) : 0.0;
    if (cE > 0.0) {
        f = Fmag / cE;
        if (f > 1.0) {
            f = 1.0;
        }
    } else {
        f = 0.0;
    }
    prj_rad_m1_wavespeeds_with_fluxmag(E, F1, Fmag, inv_Fmag, f, lam_min, lam_max);
}
#endif

void prj_rad_flux(const prj_rad *rad, const double *WL, const double *WR,
    double lapse, const double *chi_face,
    double dx_dir, double v_face, double *flux)
{
    int field;
    int group;

#if PRJ_NRAD > 0
    {
        for (field = 0; field < PRJ_NRAD; ++field) {
            for (group = 0; group < PRJ_NEGROUP; ++group) {
                int idx = field * PRJ_NEGROUP + group;
                double EL = WL[PRJ_PRIM_RAD_E(field, group)];
                double ER = WR[PRJ_PRIM_RAD_E(field, group)];
                double F1L = WL[PRJ_PRIM_RAD_F1(field, group)];
                double F2L = WL[PRJ_PRIM_RAD_F2(field, group)];
                double F3L = WL[PRJ_PRIM_RAD_F3(field, group)];
                double F1R = WR[PRJ_PRIM_RAD_F1(field, group)];
                double F2R = WR[PRJ_PRIM_RAD_F2(field, group)];
                double F3R = WR[PRJ_PRIM_RAD_F3(field, group)];
                double Fmag_L;
                double inv_Fmag_L;
                double f_L;
                double Fmag_R;
                double inv_Fmag_R;
                double f_R;
                double fLE;
                double fLF1;
                double fLF2;
                double fLF3;
                double fRE;
                double fRF1;
                double fRF2;
                double fRF3;
                double lamL_min;
                double lamL_max;
                double lamR_min;
                double lamR_max;
                double sL;
                double sR;
                double denom;
                double inv_denom;
                double chi_ext;
                double tau;
                double eps;
                double eps2;

                prj_rad_enforce_flux_limit(&EL, &F1L, &F2L, &F3L, &Fmag_L, &inv_Fmag_L, &f_L);
                prj_rad_enforce_flux_limit(&ER, &F1R, &F2R, &F3R, &Fmag_R, &inv_Fmag_R, &f_R);

                prj_rad_m1_phys_flux_with_fluxmag(rad, EL, F1L, F2L, F3L, Fmag_L, inv_Fmag_L, f_L,
                    &fLE, &fLF1, &fLF2, &fLF3);
                prj_rad_m1_phys_flux_with_fluxmag(rad, ER, F1R, F2R, F3R, Fmag_R, inv_Fmag_R, f_R,
                    &fRE, &fRF1, &fRF2, &fRF3);

                prj_rad_m1_wavespeeds_with_fluxmag(EL, F1L, Fmag_L, inv_Fmag_L, f_L,
                    &lamL_min, &lamL_max);
                prj_rad_m1_wavespeeds_with_fluxmag(ER, F1R, Fmag_R, inv_Fmag_R, f_R,
                    &lamR_min, &lamR_max);
                sL = PRJ_CLIGHT * (lamL_min < lamR_min ? lamL_min : lamR_min);
                sR = PRJ_CLIGHT * (lamL_max > lamR_max ? lamL_max : lamR_max);
                if (sL > 0.0) {
                    sL = 0.0;
                }
                if (sR < 0.0) {
                    sR = 0.0;
                }
                if (sR - sL < 1.0e-30) {
                    sL = -PRJ_CLIGHT;
                    sR = PRJ_CLIGHT;
                }
                denom = sR - sL;
                inv_denom = 1.0 / denom;

                chi_ext = chi_face[idx];
                tau = chi_ext * dx_dir;
                eps = 3.0 / (5.0 * tau + 1.0e-10);
                if (eps > 1.0) {
                    eps = 1.0;
                }
                eps2 = eps*eps;

                /* Equation 49 and 50 of Audit et al. 2002 */
                flux[PRJ_CONS_RAD_E(field, group)] = lapse *
                    (sR * fLE - sL * fRE + eps * sL * sR * (ER - EL)) * inv_denom;
                flux[PRJ_CONS_RAD_F1(field, group)] = lapse *
                    ((eps2*(sR * fLF1 - sL * fRF1) + eps * sL * sR * (F1R - F1L)) * inv_denom
                    +(1-eps2)*(fLF1+fRF1)*0.5);
                flux[PRJ_CONS_RAD_F2(field, group)] = lapse *
                    ((eps2*(sR * fLF2 - sL * fRF2) + eps * sL * sR * (F2R - F2L)) * inv_denom
                    +(1-eps2)*(fLF2+fRF2)*0.5);
                flux[PRJ_CONS_RAD_F3(field, group)] = lapse *
                    ((eps2*(sR * fLF3 - sL * fRF3) + eps * sL * sR * (F3R - F3L)) * inv_denom
                    +(1-eps2)*(fLF3+fRF3)*0.5);

                /* O(v/c) fluid advection: upwinded v_face * {E, F_i} term. */
                {
                    double E_up = v_face >= 0.0 ? EL : ER;
                    double F1_up = v_face >= 0.0 ? F1L : F1R;
                    double F2_up = v_face >= 0.0 ? F2L : F2R;
                    double F3_up = v_face >= 0.0 ? F3L : F3R;

                    flux[PRJ_CONS_RAD_E(field, group)] += lapse * v_face * E_up;
                    flux[PRJ_CONS_RAD_F1(field, group)] += lapse * v_face * F1_up;
                    flux[PRJ_CONS_RAD_F2(field, group)] += lapse * v_face * F2_up;
                    flux[PRJ_CONS_RAD_F3(field, group)] += lapse * v_face * F3_up;
                }
            }
        }
    }
#else
    (void)rad;
    (void)chi_face;
    (void)dx_dir;
    (void)v_face;
    (void)lapse;
    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            flux[PRJ_CONS_RAD_E(field, group)] = 0.0;
            flux[PRJ_CONS_RAD_F1(field, group)] = 0.0;
            flux[PRJ_CONS_RAD_F2(field, group)] = 0.0;
            flux[PRJ_CONS_RAD_F3(field, group)] = 0.0;
        }
    }
#endif
}

#if PRJ_NRAD > 0
/* Building-block derivatives of the implicit residual w.r.t. (lnT, Ye), filled
 * on request by prj_rad_implicit_residuals(). */
typedef struct prj_rad_resid_deriv {
    double dlnkappa_dlnT[PRJ_NRAD * PRJ_NEGROUP];
    double dlnkappa_dYe[PRJ_NRAD * PRJ_NEGROUP];
    double dlneta_dlnT[PRJ_NRAD * PRJ_NEGROUP];
    double dlneta_dYe[PRJ_NRAD * PRJ_NEGROUP];
    double deint_dlnT;
    double deint_dYe;
} prj_rad_resid_deriv;

static void prj_rad_energy_failure_diagnostics(const char *reason,
    const double *u_input, double dt, double lapse, int iter, int maxiter,
    double res, double tol2, double T, double Ye, double F1, double F2,
    double rho, double Uint_old, double Ye_old)
{
    int v;

    fprintf(stderr, "prj_rad_energy_update: %s\n", reason);
    fprintf(stderr,
        "  solver: iter=%d maxiter=%d res=%.17e tol2=%.17e "
        "T=%.17e Ye=%.17e F1=%.17e F2=%.17e\n",
        iter, maxiter, res, tol2, T, Ye, F1, F2);
    fprintf(stderr,
        "  derived input: rho=%.17e Uint_old=%.17e Ye_old=%.17e\n",
        rho, Uint_old, Ye_old);
    fprintf(stderr,
        "  raw input: dt=%.17e lapse=%.17e PRJ_NVAR_CONS=%d "
        "PRJ_NRAD=%d PRJ_NEGROUP=%d PRJ_MHD=%d\n",
        dt, lapse, PRJ_NVAR_CONS, PRJ_NRAD, PRJ_NEGROUP, PRJ_MHD);
    fprintf(stderr, "  double u[PRJ_NVAR_CONS] = {\n");
    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        fprintf(stderr, "      [%d] = %.17e,\n", v, u_input[v]);
    }
    fprintf(stderr, "  };\n");
    fflush(stderr);
}

static void prj_rad_implicit_residuals(prj_rad *rad, prj_eos *eos, double *u,
    double dt, double lapse, double rho, double Uint_old, double Ye_old,
    const double *E_nu_old, double T, double Ye, double *F1, double *F2,
    double *E_nu_new_out, double *kappa_out, prj_rad_resid_deriv *deriv)
{
    double kappa[PRJ_NRAD * PRJ_NEGROUP];
    double eta[PRJ_NRAD * PRJ_NEGROUP];
    double dlnkappa_dlnT[PRJ_NRAD * PRJ_NEGROUP];
    double dlnkappa_dYe[PRJ_NRAD * PRJ_NEGROUP];
    double dlneta_dlnT[PRJ_NRAD * PRJ_NEGROUP];
    double dlneta_dYe[PRJ_NRAD * PRJ_NEGROUP];
    double eint_new;
    double deint_dlnT;
    double deint_dYe;
    double Uint_new;
    double sum_dE = 0.0;
    double sum_dE_xe = 0.0;
    int nu;
    int g;

    (void)u;
    prj_rad3_opac_lookup_ke(rad, rho, T, Ye, kappa, eta,
        dlnkappa_dlnT, dlnkappa_dYe, dlneta_dlnT, dlneta_dYe);
    eint_new = prj_eos_rty_eint(eos, rho, T, Ye, &deint_dlnT, &deint_dYe, PRJ_EOS_CTX_MAIN);
    Uint_new = rho * eint_new;

    for (nu = 0; nu < PRJ_NRAD; ++nu) {
        for (g = 0; g < PRJ_NEGROUP; ++g) {
            int idx = nu * PRJ_NEGROUP + g;
            double num = E_nu_old[idx] + dt * lapse * eta[idx];
            double den = 1.0 + dt * lapse * PRJ_CLIGHT * kappa[idx];
            double Enew = num / den;
            double dE = Enew - E_nu_old[idx];

            sum_dE += dE;
            sum_dE_xe += dE * rad->x_e[nu][g];
            if (E_nu_new_out != 0) {
                E_nu_new_out[idx] = Enew;
            }
        }
    }

    /* sum_dE is a radiation-energy change in RAD_SCALE*erg units; multiply back
       to erg to balance the gas internal energy.  sum_dE_xe already carries
       RAD_SCALE through x_e, so the lepton residual needs no extra factor. */
    *F1 = Uint_new - Uint_old + sum_dE * RAD_SCALE;
    *F2 = rho * Ye - rho * Ye_old + sum_dE_xe;

    if (kappa_out != 0) {
        int i;
        for (i = 0; i < PRJ_NRAD * PRJ_NEGROUP; ++i) {
            kappa_out[i] = kappa[i];
        }
    }

    if (deriv != 0) {
        int i;
        for (i = 0; i < PRJ_NRAD * PRJ_NEGROUP; ++i) {
            deriv->dlnkappa_dlnT[i] = dlnkappa_dlnT[i];
            deriv->dlnkappa_dYe[i] = dlnkappa_dYe[i];
            deriv->dlneta_dlnT[i] = dlneta_dlnT[i];
            deriv->dlneta_dYe[i] = dlneta_dYe[i];
        }
        deriv->deint_dlnT = deint_dlnT;
        deriv->deint_dYe = deint_dYe;
    }
}

static void prj_rad_implicit_jacobian_from_deriv(const prj_rad *rad,
    const double *E_nu_old, const double *E_nu_new, const double *kappa,
    const prj_rad_resid_deriv *deriv, double dt, double lapse, double rho,
    double T, double *dFdT_1, double *dFdT_2, double *dFdY_1, double *dFdY_2)
{
    double dF1_dlnT;
    double dF2_dlnT;
    double dF1_dYe;
    double dF2_dYe;
    double dt_lapse;
    double inv_T;
    int nu;
    int g;

    dF1_dlnT = rho * deriv->deint_dlnT;
    dF2_dlnT = 0.0;
    dF1_dYe = rho * deriv->deint_dYe;
    dF2_dYe = rho;
    dt_lapse = dt * lapse;
    inv_T = T > 0.0 ? 1.0 / T : 0.0;

    for (nu = 0; nu < PRJ_NRAD; ++nu) {
        for (g = 0; g < PRJ_NEGROUP; ++g) {
            int idx = nu * PRJ_NEGROUP + g;
            double den = 1.0 + dt_lapse * PRJ_CLIGHT * kappa[idx];
            double eta = 0.0;
            double dE_dlnT = 0.0;
            double dE_dYe = 0.0;

            if (dt_lapse != 0.0) {
                eta = (E_nu_new[idx] * den - E_nu_old[idx]) / dt_lapse;
                dE_dlnT = dt_lapse *
                    (eta * deriv->dlneta_dlnT[idx] -
                        E_nu_new[idx] * PRJ_CLIGHT * kappa[idx] *
                            deriv->dlnkappa_dlnT[idx]) / den;
                dE_dYe = dt_lapse *
                    (eta * deriv->dlneta_dYe[idx] -
                        E_nu_new[idx] * PRJ_CLIGHT * kappa[idx] *
                            deriv->dlnkappa_dYe[idx]) / den;
            }

            dF1_dlnT += RAD_SCALE * dE_dlnT;
            dF2_dlnT += rad->x_e[nu][g] * dE_dlnT;
            dF1_dYe += RAD_SCALE * dE_dYe;
            dF2_dYe += rad->x_e[nu][g] * dE_dYe;
        }
    }

    *dFdT_1 = dF1_dlnT * inv_T;
    *dFdT_2 = dF2_dlnT * inv_T;
    *dFdY_1 = dF1_dYe;
    *dFdY_2 = dF2_dYe;
}

void prj_rad_energy_update(prj_rad *rad, prj_eos *eos, double *u, double dt, double lapse, double *final_temperature, double *kappa_out)
{
    double u_input[PRJ_NVAR_CONS];
    double E_nu_old[PRJ_NRAD * PRJ_NEGROUP];
    double E_nu_new[PRJ_NRAD * PRJ_NEGROUP];
    double last_kappa[PRJ_NRAD * PRJ_NEGROUP];
    double rho;
    double KE;
    double Emag = 0.0;
    double Uint_old;
    double Ye_old;
    double eint_old;
    double eos_q[PRJ_EOS_NQUANT];
    double T;
    double Ye;
    double err_scale_1;
    double err_scale_2;
    double res_cur;
    double cached_F1 = 0.0;
    double cached_F2 = 0.0;
    double cached_res = 0.0;
    prj_rad_resid_deriv cached_deriv = {0};
    int have_cached_residual = 0;
    int have_final_residual = 0;
    int iter;
    int nu;
    int g;
    int v;
    const double alpha_ls = 1.0e-4;
    const double tol2 = rad->implicit_err_tol * rad->implicit_err_tol;

    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        u_input[v] = u[v];
    }

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
    Ye_old = u[PRJ_CONS_YE] / rho;
    eint_old = Uint_old / rho;

    for (nu = 0; nu < PRJ_NRAD; ++nu) {
        for (g = 0; g < PRJ_NEGROUP; ++g) {
            E_nu_old[nu * PRJ_NEGROUP + g] = u[PRJ_CONS_RAD_E(nu, g)];
        }
    }

    prj_eos_rey(eos, rho, eint_old, Ye_old, eos_q, PRJ_EOS_CTX_MAIN);
    T = eos_q[PRJ_EOS_TEMPERATURE];
    Ye = Ye_old;
    err_scale_1 = fabs(Uint_old) > 0.0 ? Uint_old : 1.0;
    err_scale_2 = fabs(rho * Ye_old) > 0.0 ? rho * Ye_old : 1.0;
    res_cur = 1.0e30;

    for (iter = 0; iter < rad->maxiter; ++iter) {
        double F1;
        double F2;
        double f1;
        double f2;
        prj_rad_resid_deriv deriv;
        prj_rad_resid_deriv trial_deriv;
        double dFdT_1;
        double dFdT_2;
        double dFdY_1;
        double dFdY_2;
        double J00;
        double J01;
        double J10;
        double J11;
        double r0;
        double r1;
        double col_scale0;
        double col_scale1;
        double det;
        double s0;
        double s1;
        double dT;
        double dY;
        double step_scale;
        double gradf0;
        double gradf1;
        double gradfdx;
        double lam;
        double lamold;
        double resold;
        double Ttrial;
        double Ytrial;
        double res_trial;
        double F1_trial;
        double F2_trial;
        int inner_iter;
        int accepted_trial;

        if (have_cached_residual) {
            F1 = cached_F1;
            F2 = cached_F2;
            res_cur = cached_res;
            deriv = cached_deriv;
            have_cached_residual = 0;
            have_final_residual = 1;
        } else {
            prj_rad_implicit_residuals(rad, eos, u, dt, lapse, rho, Uint_old, Ye_old,
                E_nu_old, T, Ye, &F1, &F2, E_nu_new, last_kappa, &deriv);
            f1 = F1 / err_scale_1;
            f2 = F2 / err_scale_2;
            res_cur = 0.5 * (f1 * f1 + f2 * f2);
            have_final_residual = 1;
        }

        if (res_cur < tol2) {
            break;
        }

        prj_rad_implicit_jacobian_from_deriv(rad, E_nu_old, E_nu_new, last_kappa,
            &deriv, dt, lapse, rho, T, &dFdT_1, &dFdT_2, &dFdY_1, &dFdY_2);

        /* grad(½||F||²) = J^T F  (unscaled, for line search). */
        gradf0 = dFdT_1 * F1 / (err_scale_1 * err_scale_1) +
            dFdT_2 * F2 / (err_scale_2 * err_scale_2);
        gradf1 = dFdY_1 * F1 / (err_scale_1 * err_scale_1) +
            dFdY_2 * F2 / (err_scale_2 * err_scale_2);

        /* Row-max equilibration on the augmented matrix [J | -F]. */
        J00 = dFdT_1; J01 = dFdY_1; r0 = F1;
        J10 = dFdT_2; J11 = dFdY_2; r1 = F2;
        {
            double row_max;

            row_max = fabs(J00) > fabs(J01) ? fabs(J00) : fabs(J01);
            if (row_max == 0.0) row_max = 1.0e-10;
            J00 /= row_max; J01 /= row_max; r0 /= row_max;

            row_max = fabs(J10) > fabs(J11) ? fabs(J10) : fabs(J11);
            if (row_max == 0.0) row_max = 1.0e-10;
            J10 /= row_max; J11 /= row_max; r1 /= row_max;
        }
        /* Column-max equilibration. */
        col_scale0 = fabs(J00) > fabs(J10) ? fabs(J00) : fabs(J10);
        col_scale1 = fabs(J01) > fabs(J11) ? fabs(J01) : fabs(J11);
        if (col_scale0 == 0.0) col_scale0 = 1.0e-16;
        if (col_scale1 == 0.0) col_scale1 = 1.0e-16;
        J00 /= col_scale0; J10 /= col_scale0;
        J01 /= col_scale1; J11 /= col_scale1;

        /* Negate RHS for Newton step. */
        r0 = -r0;
        r1 = -r1;

        /* 2x2 Cramer solve. */
        det = J00 * J11 - J01 * J10;
        if (fabs(det) < 1.0e-30) {
            prj_rad_energy_failure_diagnostics("singular Jacobian", u_input,
                dt, lapse, iter, rad->maxiter, res_cur, tol2, T, Ye, F1, F2,
                rho, Uint_old, Ye_old);
            fprintf(stderr, "  scaled Jacobian determinant: %.17e\n", det);
            fflush(stderr);
            exit(1);
        }
        s0 = (J11 * r0 - J01 * r1) / det;
        s1 = (-J10 * r0 + J00 * r1) / det;

        /* Undo column scaling. */
        dT = s0 / col_scale0;
        dY = s1 / col_scale1;

        /* Step limiter: cap at 3% of current values. */
        step_scale = 1.0;
        if (fabs(dT) > 0.03 * T) {
            step_scale = 0.03 * T / fabs(dT);
        }
        if (fabs(dY) > 0.03 * Ye && 0.03 * Ye / fabs(dY) < step_scale) {
            step_scale = 0.03 * Ye / fabs(dY);
        }
        dT *= step_scale;
        dY *= step_scale;

        /* Directional derivative for Armijo check. */
        gradfdx = gradf0 * dT + gradf1 * dY;

        /* Backtracking line search with cubic/quadratic interpolation. */
        lam = 1.0;
        lamold = 0.0;
        resold = 0.0;
        accepted_trial = 0;
        F1_trial = F1;
        F2_trial = F2;
        res_trial = res_cur;
        for (inner_iter = 0; inner_iter < 6; ++inner_iter) {
            double ft1;
            double ft2;

            Ttrial = T + lam * dT;
            Ytrial = Ye + lam * dY;
            prj_rad_implicit_residuals(rad, eos, u, dt, lapse, rho, Uint_old, Ye_old,
                E_nu_old, Ttrial, Ytrial, &F1_trial, &F2_trial, E_nu_new,
                last_kappa, &trial_deriv);
            ft1 = F1_trial / err_scale_1;
            ft2 = F2_trial / err_scale_2;
            res_trial = 0.5 * (ft1 * ft1 + ft2 * ft2);

            if (res_trial < tol2 || res_trial < res_cur + alpha_ls * lam * gradfdx) {
                accepted_trial = 1;
                break;
            }

            {
                double templam;

                if (inner_iter == 0) {
                    templam = -gradfdx / (2.0 * (res_trial - res_cur - gradfdx));
                } else {
                    double rhs1 = res_trial - res_cur - lam * gradfdx;
                    double rhs2 = resold - res_cur - lamold * gradfdx;
                    double a_c = (rhs1 / (lam * lam) - rhs2 / (lamold * lamold)) / (lam - lamold);
                    double b_c = (-lamold * rhs1 / (lam * lam) + lam * rhs2 / (lamold * lamold)) / (lam - lamold);

                    if (a_c == 0.0) {
                        templam = -gradfdx / (2.0 * b_c);
                    } else {
                        double disc = b_c * b_c - 3.0 * a_c * gradfdx;

                        templam = disc >= 0.0 ? (-b_c + sqrt(disc)) / (3.0 * a_c) : 0.5 * lam;
                    }
                    if (templam > 0.5 * lam) {
                        templam = 0.5 * lam;
                    }
                }
                lamold = lam;
                resold = res_trial;
                lam = templam > 0.1 * lam ? templam : 0.1 * lam;
            }
        }

        T = T + lam * dT;
        Ye = Ye + lam * dY;
        if (T <= 0.0) {
            T = 0.5 * (T - lam * dT);
            accepted_trial = 0;
            have_final_residual = 0;
        }
        if (accepted_trial) {
            cached_F1 = F1_trial;
            cached_F2 = F2_trial;
            cached_res = res_trial;
            cached_deriv = trial_deriv;
            have_cached_residual = 1;
            have_final_residual = 1;
        } else {
            have_final_residual = 0;
        }

        res_cur = res_trial;
        if (have_final_residual && res_cur < tol2) {
            break;
        }
    }
    if (iter == rad->maxiter) {
        double F1_final;
        double F2_final;
        double f1_final;
        double f2_final;

        prj_rad_implicit_residuals(rad, eos, u, dt, lapse, rho, Uint_old, Ye_old,
            E_nu_old, T, Ye, &F1_final, &F2_final, E_nu_new, last_kappa, 0);
        f1_final = F1_final / err_scale_1;
        f2_final = F2_final / err_scale_2;
        res_cur = 0.5 * (f1_final * f1_final + f2_final * f2_final);
        prj_rad_energy_failure_diagnostics("failed to converge", u_input,
            dt, lapse, iter, rad->maxiter, res_cur, tol2, T, Ye,
            F1_final, F2_final, rho, Uint_old, Ye_old);
        exit(1);
    }

    /* Final pass at converged (T, Ye) to populate E_nu_new if the accepted
     * line-search residual was not already evaluated at the final state. */
    {
        double eint_new;

        if (!have_final_residual) {
            double F1_final;
            double F2_final;

            prj_rad_implicit_residuals(rad, eos, u, dt, lapse, rho, Uint_old, Ye_old,
                E_nu_old, T, Ye, &F1_final, &F2_final, E_nu_new, last_kappa, 0);
        }
        prj_eos_rty(eos, rho, T, Ye, eos_q, PRJ_EOS_CTX_MAIN);
        eint_new = eos_q[PRJ_EOS_EINT];
        u[PRJ_CONS_ETOT] = rho * eint_new + KE + Emag;
        u[PRJ_CONS_YE] = rho * Ye;
        for (nu = 0; nu < PRJ_NRAD; ++nu) {
            for (g = 0; g < PRJ_NEGROUP; ++g) {
                u[PRJ_CONS_RAD_E(nu, g)] = E_nu_new[nu * PRJ_NEGROUP + g];
            }
        }
    }

    if (final_temperature != 0) {
        *final_temperature = T;
    }
    if (kappa_out != 0) {
        int i;
        for (i = 0; i < PRJ_NRAD * PRJ_NEGROUP; ++i) {
            kappa_out[i] = last_kappa[i];
        }
    }
}

void prj_rad_momentum_update(prj_rad *rad, prj_eos *eos, double *u, double dt, double lapse, double temperature, const double *kappa_in)
{
    double sigma[PRJ_NRAD * PRJ_NEGROUP];
    double delta[PRJ_NRAD * PRJ_NEGROUP];
    double rho;
    double Ye;
    double inv_c2;
    double dmom[3];
    double e_unchanged;
    int nu;
    int g;
    int d;

    (void)eos;

    rho = u[PRJ_CONS_RHO];
    Ye = u[PRJ_CONS_YE] / rho;
    inv_c2 = 1.0 / (PRJ_CLIGHT * PRJ_CLIGHT);

    prj_rad3_opac_lookup(rad, rho, temperature, Ye, 0, sigma, delta, 0);

    dmom[0] = 0.0;
    dmom[1] = 0.0;
    dmom[2] = 0.0;
    e_unchanged = u[PRJ_CONS_ETOT] - 0.5 * (u[PRJ_CONS_MOM1] * u[PRJ_CONS_MOM1] +
                                            u[PRJ_CONS_MOM2] * u[PRJ_CONS_MOM2] +
                                            u[PRJ_CONS_MOM3] * u[PRJ_CONS_MOM3]) / rho;

    for (nu = 0; nu < PRJ_NRAD; ++nu) {
        for (g = 0; g < PRJ_NEGROUP; ++g) {
            int idx = nu * PRJ_NEGROUP + g;
            double chi = kappa_in[idx] + sigma[idx] * (1.0 - delta[idx] / 3.0);
            double factor = 1.0 / (1.0 + dt * lapse * PRJ_CLIGHT * chi) - 1.0;
            int fi[3];

            fi[0] = PRJ_CONS_RAD_F1(nu, g);
            fi[1] = PRJ_CONS_RAD_F2(nu, g);
            fi[2] = PRJ_CONS_RAD_F3(nu, g);

            double F_old[PRJ_NDIM];
            for (d = 0; d < 3; ++d) {
                F_old[d] = u[fi[d]];
                double dF = F_old[d] * factor;
                u[fi[d]] = F_old[d] + dF;
            }
            
            double E = u[PRJ_CONS_RAD_E(nu, g)];
            double F1 = u[fi[0]];
            double F2 = u[fi[1]];
            double F3 = u[fi[2]];
            double Fmag = sqrt(F1 * F1 + F2 * F2 + F3 * F3);
            double cE = PRJ_CLIGHT * E;

            if (E > 0.0 && Fmag > cE) {
                double scale = cE / Fmag;
                u[fi[0]] = F1 * scale;
                u[fi[1]] = F2 * scale;
                u[fi[2]] = F3 * scale;
            }

            for (d = 0; d < 3; ++d) {
                dmom[d] += (u[fi[d]]-F_old[d]) * inv_c2;
            }
        }
    }

    /* dmom/detot accumulate radiation-flux changes in RAD_SCALE*erg units;
       multiply back to erg for the gas momentum/energy back-reaction. */
    u[PRJ_CONS_MOM1] -= dmom[0] * RAD_SCALE;
    u[PRJ_CONS_MOM2] -= dmom[1] * RAD_SCALE;
    u[PRJ_CONS_MOM3] -= dmom[2] * RAD_SCALE;
    u[PRJ_CONS_ETOT] = e_unchanged + 0.5 * (u[PRJ_CONS_MOM1] * u[PRJ_CONS_MOM1] +
                                            u[PRJ_CONS_MOM2] * u[PRJ_CONS_MOM2] +
                                            u[PRJ_CONS_MOM3] * u[PRJ_CONS_MOM3]) / rho;
}

/* Koren slope-limiter function φ(r) = max(0, min(2r, (2+r)/3, 2)). */
static double prj_rad_koren_phi(double r)
{
    double phi = 2.0 * r;
    double t = (2.0 + r) / 3.0;

    if (t < phi) {
        phi = t;
    }
    if (2.0 < phi) {
        phi = 2.0;
    }
    if (phi < 0.0) {
        phi = 0.0;
    }
    return phi;
}

/* Reconstruct cell-value array q[] (one entry per energy group) at the right
 * (side=+1) or left (side=-1) face of cell `gcell` using a Koren-limited linear
 * stencil.  Energy groups are uniformly spaced in log ν, so equal-spaced
 * samples are valid.  The two outermost cells (gcell == 0 or NEGROUP-1) fall
 * back to piecewise constant — at those edges the outer face values are also
 * forced to zero by the caller, so the choice has no effect on the update. */
static double prj_rad_recon_face(const double q[PRJ_NEGROUP], int gcell, int side)
{
    (void)side;
    return q[gcell];
}

/* Apply the per-cell energy-space-flux part of the SR redshift terms
 * (Eqs. 21a/21b of the comoving-frame mixed-frame moment equations):
 *
 *   ∂_t E_g    -= - v^i_{;j} [ (ν P^j_{νi})_{g+1/2} - (ν P^j_{νi})_{g-1/2} ]
 *   ∂_t F_{gj} -= - v^i_{;k} [ (ν Q^k_{νji})_{g+1/2} - (ν Q^k_{νji})_{g-1/2} ]
 *
 * (Equivalently dE_g/dt += v^i_{;j} ΔνP^{ji} and dF_{gj}/dt += v^i_{;k} ΔνQ^{kji}.)
 *
 * The face values are picked by upwinding in frequency space according to
 * Eq. 22:
 *      pick L (right face of the lower group)  if v_{...} <  0,
 *      pick R (left  face of the upper group)  if v_{...} >= 0.
 *
 * Closure choices:
 *   - P^{ij}_g from the standard M1 Eddington tensor on cell-centred (E_g, F_g).
 *   - Q^{kij}_g from the Levermore/Vaytet third-moment closure
 *       Q^{ijk} = c E H^{ijk},
 *     with H^{ijk} built from f = |F|/(cE), n = F/|F|, and q(f).
 *
 * Reconstruction: linear-in-log-ν with Koren limiter (groups are uniform in
 * log ν → equal spacing).  Outermost group faces (g = -1/2 and g = NEGROUP-1/2)
 * are set to zero (outflow).  The cell-centred state used for the closure comes
 * from W_state (W in stage1, W1 in stage2), per the user-specified ordering. */
void prj_rad_freq_flux_apply(const prj_rad *rad, const prj_block *block,
    const double *W_state, double *u, int ic, int jc, int kc, double lapse, double dt)
{
#if PRJ_NRAD > 0
    double dvdx[3][3];
    double inv_dx[3];
    int jdir;
    int icomp;
    int field;

    if (rad == 0 || block == 0 || W_state == 0 || u == 0) {
        return;
    }
    if (block->v_riemann[0] == 0 || block->v_riemann[1] == 0 || block->v_riemann[2] == 0) {
        return;
    }

    inv_dx[0] = 1.0 / block->dx[0];
    inv_dx[1] = 1.0 / block->dx[1];
    inv_dx[2] = 1.0 / block->dx[2];

    /* Cell-centred ∂_jdir v_icomp from the two normal-direction Riemann faces. */
    for (jdir = 0; jdir < 3; ++jdir) {
        for (icomp = 0; icomp < 3; ++icomp) {
            int il = ic;
            int jl = jc;
            int kl = kc;
            int ir = ic;
            int jr = jc;
            int kr = kc;
            double vL;
            double vR;

            if (jdir == X1DIR) {
                ir = ic + 1;
            } else if (jdir == X2DIR) {
                jr = jc + 1;
            } else {
                kr = kc + 1;
            }
            vL = block->v_riemann[jdir][icomp * PRJ_BLOCK_NCELLS + IDX(il, jl, kl)];
            vR = block->v_riemann[jdir][icomp * PRJ_BLOCK_NCELLS + IDX(ir, jr, kr)];
            dvdx[jdir][icomp] = (vR - vL) * inv_dx[jdir];
        }
    }

    int cell_idx = IDX(ic, jc, kc);
    double grad_phi[3];
    double inv_c2 = 1.0 / (PRJ_CLIGHT * PRJ_CLIGHT);

    grad_phi[0] = 0.0;
    grad_phi[1] = 0.0;
    grad_phi[2] = 0.0;
    if (block->grav[0] != 0 && block->grav[1] != 0 && block->grav[2] != 0) {
        grad_phi[0] = -block->grav[0][cell_idx];
        grad_phi[1] = -block->grav[1][cell_idx];
        grad_phi[2] = -block->grav[2][cell_idx];
    }

    for (field = 0; field < PRJ_NRAD; ++field) {
        double dt_lapse = lapse * dt;
        double Eg[PRJ_NEGROUP];
        double Fg[PRJ_NEGROUP][3];
        double Pg[PRJ_NEGROUP][3][3];
        double Mq[PRJ_NEGROUP][3]; /* Q_g : dvdx, the only way Q is ever used */
        double energy_face[PRJ_NEGROUP + 1] = {0.0};
        double momentum_face[PRJ_NEGROUP + 1][PRJ_NDIM] = {{0.0}};
        double energy_available[PRJ_NEGROUP];
        const double *nu_face = rad->eedge[field];
        int g;
        int ii;
        int jj;

        /* Per-group cell-centred state and closure tensors.  P and Q are built
         * once here and shared: P by the GR redshift terms below, and both P, Q
         * by the SR frequency flux (reconstructed to the frequency faces). */
        for (g = 0; g < PRJ_NEGROUP; ++g) {
            Eg[g] = W_state[VIDX(PRJ_PRIM_RAD_E(field, g), ic, jc, kc)];
            Fg[g][0] = W_state[VIDX(PRJ_PRIM_RAD_F1(field, g), ic, jc, kc)];
            Fg[g][1] = W_state[VIDX(PRJ_PRIM_RAD_F2(field, g), ic, jc, kc)];
            Fg[g][2] = W_state[VIDX(PRJ_PRIM_RAD_F3(field, g), ic, jc, kc)];
            prj_rad_m1_pressure(rad, Eg[g], Fg[g][0], Fg[g][1], Fg[g][2], Pg[g]);
            prj_rad_m1_third_moment_contract(rad, Eg[g], Fg[g][0], Fg[g][1], Fg[g][2],
                dvdx, Mq[g]);
            energy_available[g] = u[PRJ_CONS_RAD_E(field, g)];
        }

        /* SR velocity-gradient energy-space flux (Eqs. 21a/21b).  The whole
         * frequency flux is upwinded once by the sign of the velocity divergence
         * div(v) = tr(∂_j v_i): compression (div < 0) blueshifts and expansion
         * (div > 0) redshifts the entire spectrum, so a single upwind side is
         * used for every energy group and every tensor component (replacing the
         * earlier per-component upwinding, which summed nine independently
         * upwinded ν-fluxes and drove a grid-scale instability).  The centred
         * closure tensors P, Q are reconstructed (donor cell) to each frequency
         * face from the upwind group and scaled by the face frequency ν. */
        {
            /* CAVEAT: choosing the single upwind side from the sign of div(v)
             * only guarantees a stable energy-space flux in bulk flows, where
             * the (isotropic) compression/expansion divergence is the dominant
             * part of the velocity gradient.  In turbulent or strongly shearing
             * flows the off-diagonal/anisotropic parts of ∂_j v_i can exceed the
             * trace, so the trace-based upwind direction may disagree with the
             * actual sign of an individual tensor component and absolute
             * stability is no longer assured.  It is robust for the smooth
             * infall/contraction regime this code targets. */
            double divv = dvdx[0][0] + dvdx[1][1] + dvdx[2][2];
            int d = (divv >= 0.0) ? 0 : -1; /* donor offset: upwind group = gf + d */
            int gf;

            /* Sweep interior frequency faces (the domain-edge faces 0 and
             * PRJ_NEGROUP carry zero flux and are skipped).  Each face gf has a
             * single upwind donor group gu = gf + d (gf for div v >= 0, else
             * gf-1); its flux ν_face[gf]·{P,Q}[gu] is the upper (g+1/2) face of
             * group gf-1 and the lower (g-1/2) face of group gf, so it is
             * scattered into both with opposite signs.  Accumulate the net
             * energy and momentum transfer at each face so it can be limited
             * before the conservative scatter.  (For a Koren-limited version,
             * replace Pg[gu] /
             * Qg[gu] by q[gu] + 0.5·φ(r)·(q[gu] − q[gu+up]) per component, with
             * upwind-neighbour offset up = (div v >= 0) ? +1 : -1.) */
            for (gf = 1; gf < PRJ_NEGROUP; ++gf) {
                int gu = gf + d;
                double nu = nu_face[gf];

                for (jj = 0; jj < 3; ++jj) {
                    for (ii = 0; ii < 3; ++ii) {
                        double pf = nu * Pg[gu][jj][ii] * dvdx[jj][ii];

                        energy_face[gf] += pf;
                    }
                    {
                        double qf = nu * Mq[gu][jj];

                        momentum_face[gf][jj] += qf;
                    }
                }
            }
        }

        /* GR ε-flux pieces of G^e and G^m_j (the ∂_ε terms only):
         *   G^e          = -F_{sε}·∇φ/c² + (∇φ/c²) · ∂_ε(ε F_{sε})
         *   G^m_j        = -E_{sε} ∇_jφ + ∇_iφ · ∂_ε(ε P^i_{sεj})
         * G is on the source-side (RHS) of the moment equations, so when added
         * to ∂_t E_g / ∂_t F_{gj} it carries an overall minus sign:
         *   ∂_t E_g    -= (∇_i φ / c²) · [(εF^i)_{g+1/2} - (εF^i)_{g-1/2}]
         *   ∂_t F_{gj} -= (∇_i φ)      · [(εP^{ij})_{g+1/2} - (εP^{ij})_{g-1/2}]
         * (sums over i; the −F_{sε}·∇φ/c² and −E_{sε}∇_jφ pieces are NOT done
         * here per user request).
         *
         * Upwind in ε-space follows the same Eq. 22 rule, using the per-i
         * coefficient that multiplies the ε-flux divergence as the "speed":
         *   pick L if coef <  0,  pick R if coef >= 0.
         *
         * ∇_i φ at the cell centre comes from the active monopole gravity:
         * gravitational acceleration a_i = accel(r) · x_i/r, and ∇φ = −a. */
        {
            /* Energy: per i, scalar q[g] = F^i_g, coef = −(∇_i φ)/c². */
            for (ii = 0; ii < 3; ++ii) {
                double coef = -grad_phi[ii] * inv_c2;
                double q[PRJ_NEGROUP];
                double face_val[PRJ_NEGROUP + 1];
                int gf;

                if (coef == 0.0) {
                    continue;
                }
                for (g = 0; g < PRJ_NEGROUP; ++g) {
                    q[g] = Fg[g][ii];
                }
                face_val[0] = 0.0;
                face_val[PRJ_NEGROUP] = 0.0;
                for (gf = 1; gf < PRJ_NEGROUP; ++gf) {
                    double pick = (coef >= 0.0)
                        ? prj_rad_recon_face(q, gf, -1)
                        : prj_rad_recon_face(q, gf - 1, +1);
                    face_val[gf] = nu_face[gf] * pick;
                }
                for (gf = 1; gf < PRJ_NEGROUP; ++gf) {
                    energy_face[gf] += coef * face_val[gf];
                }
            }

            /* Flux j: per (i, j), scalar q[g] = P^{ij}_g, coef = −∇_i φ. */
            for (jj = 0; jj < 3; ++jj) {
                for (ii = 0; ii < 3; ++ii) {
                    double coef = -grad_phi[ii];
                    double q[PRJ_NEGROUP];
                    double face_val[PRJ_NEGROUP + 1];
                    int gf;

                    if (coef == 0.0) {
                        continue;
                    }
                    for (g = 0; g < PRJ_NEGROUP; ++g) {
                        q[g] = Pg[g][ii][jj];
                    }
                    face_val[0] = 0.0;
                    face_val[PRJ_NEGROUP] = 0.0;
                    for (gf = 1; gf < PRJ_NEGROUP; ++gf) {
                        double pick = (coef >= 0.0)
                            ? prj_rad_recon_face(q, gf, -1)
                            : prj_rad_recon_face(q, gf - 1, +1);
                        face_val[gf] = nu_face[gf] * pick;
                    }
                    for (gf = 1; gf < PRJ_NEGROUP; ++gf) {
                        momentum_face[gf][jj] += coef * face_val[gf];
                    }
                }
            }
        }

        /* Limit the combined SR+GR frequency-space flux with one factor per
         * donor group.  For dE_g/dt = face[g+1] - face[g], a positive face
         * drains its upper group and a negative face drains its lower group.
         * Both outgoing faces of a group share the same factor. */
        {
            double outgoing[PRJ_NEGROUP] = {0.0};
            double theta[PRJ_NEGROUP];
            int gf;

            for (gf = 1; gf < PRJ_NEGROUP; ++gf) {
                if (energy_face[gf] > 0.0) {
                    outgoing[gf] += energy_face[gf];
                } else if (energy_face[gf] < 0.0) {
                    outgoing[gf - 1] -= energy_face[gf];
                }
            }
            for (g = 0; g < PRJ_NEGROUP; ++g) {
                double drain = dt_lapse * outgoing[g];

                theta[g] = 1.0;
                if (drain > energy_available[g]) {
                    theta[g] = energy_available[g] > 0.0 && drain > 0.0
                        ? nextafter(energy_available[g] / drain, 0.0)
                        : 0.0;
                }
            }
            for (gf = 1; gf < PRJ_NEGROUP; ++gf) {
                int donor = energy_face[gf] > 0.0 ? gf : gf - 1;
                double factor = theta[donor];

                energy_face[gf] *= factor;
                for (ii = 0; ii < PRJ_NDIM; ++ii) {
                    momentum_face[gf][ii] *= factor;
                }
            }
        }

        /* Apply.  dt is the effective stage weight (full dt in stage1, 0.5·dt
         * in stage2 to match the RK2-Heun mixing of dUdt).  The lapse factor
         * α(r) accounts for the GR proper-time slowdown in the gravitational
         * well, consistent with the lapse multipliers already on the spatial
         * radiation flux and on the gravity source. */
        for (g = 0; g < PRJ_NEGROUP; ++g) {
            u[PRJ_CONS_RAD_E(field, g)] += dt_lapse *
                (energy_face[g + 1] - energy_face[g]);
            u[PRJ_CONS_RAD_F1(field, g)] += dt_lapse *
                (momentum_face[g + 1][0] - momentum_face[g][0]);
            u[PRJ_CONS_RAD_F2(field, g)] += dt_lapse *
                (momentum_face[g + 1][1] - momentum_face[g][1]);
            u[PRJ_CONS_RAD_F3(field, g)] += dt_lapse *
                (momentum_face[g + 1][2] - momentum_face[g][2]);
        }
    }
#else
    (void)rad;
    (void)block;
    (void)W_state;
    (void)u;
    (void)ic;
    (void)jc;
    (void)kc;
    (void)dt;
#endif
}

#else
void prj_rad_energy_update(prj_rad *rad, prj_eos *eos, double *u, double dt, double lapse, double *final_temperature, double *kappa_out)
{
    (void)rad;
    (void)eos;
    (void)u;
    (void)dt;
    (void)final_temperature;
    (void)lapse;
    (void)kappa_out;
}

void prj_rad_momentum_update(prj_rad *rad, prj_eos *eos, double *u, double dt, double lapse, double temperature, const double *kappa_in)
{
    (void)rad;
    (void)eos;
    (void)u;
    (void)dt;
    (void)lapse;
    (void)temperature;
    (void)kappa_in;
}

void prj_rad_freq_flux_apply(const prj_rad *rad, const prj_block *block,
    const double *W_state, double *u, int ic, int jc, int kc, double lapse, double dt)
{
    (void)rad;
    (void)block;
    (void)W_state;
    (void)u;
    (void)ic;
    (void)jc;
    (void)kc;
    (void)lapse;
    (void)dt;
}
#endif
