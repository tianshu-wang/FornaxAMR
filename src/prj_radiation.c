#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "prj.h"

void prj_rad_init(prj_rad *rad)
{
#if PRJ_NRAD > 0
    prj_rad3_opac_init(rad);
#else
    (void)rad;
#endif
}

void prj_rad_prim2cons(const double *W, double *U)
{
    int field;
    int group;

    if (W == 0 || U == 0) {
        return;
    }

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

    if (U == 0 || W == 0) {
        return;
    }

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
void prj_rad_m1_pressure(double E, double F1, double F2, double F3, double P[3][3])
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
    chi = (3.0 + 4.0 * f * f) / (5.0 + 2.0 * sqrt(4.0 - 3.0 * f * f));
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
/* Levermore/Vaytet third-moment scalar q(f), evaluated through the equivalent
 * boost parameter β = 3f / (2 + √(4 - 3f²)).  This form is numerically stable
 * at moderate and large f; for very small β we switch to the series expansion
 *
 *   q = 4β/5 + 4β^3/21 - 4β^5/315 + O(β^7),
 *
 * which matches the exact boosted-isotropic closure and avoids catastrophic
 * cancellation in the raw closed form. */
static double prj_rad_levermore_q_factor(double f)
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

/* Levermore third moment Q^{ijk} = c E H^{ijk}, with
 *
 *   H^{ijk} = (5q - 3f)/2 n^i n^j n^k
 *           + (f - q)/2 (n^i δ^{jk} + n^j δ^{ki} + n^k δ^{ij}),
 *
 * where f = |F|/(cE), n = F/|F|, and q = q(f) above.  The tensor is fully
 * symmetric in its three indices. */
static void prj_rad_m1_third_moment(double E, double F1, double F2, double F3,
    double Q[3][3][3])
{
    double E_pos;
    double Fmag;
    double cE;
    double f;
    double q_fac;
    double n[3];
    double coef_nnn;
    double coef_mix;
    int a;
    int b;
    int c;

    E_pos = E > 0.0 ? E : 0.0;
    Fmag = sqrt(F1 * F1 + F2 * F2 + F3 * F3);
    cE = PRJ_CLIGHT * E_pos;

    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            for (c = 0; c < 3; ++c) {
                Q[a][b][c] = 0.0;
            }
        }
    }
    if (cE <= 0.0 || Fmag <= 0.0) {
        return;
    }

    f = Fmag / cE;
    if (f > 1.0) {
        f = 1.0;
    }
    q_fac = prj_rad_levermore_q_factor(f);
    n[0] = F1 / Fmag;
    n[1] = F2 / Fmag;
    n[2] = F3 / Fmag;
    coef_nnn = 0.5 * cE * (5.0 * q_fac - 3.0 * f);
    coef_mix = 0.5 * cE * (f - q_fac);

    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            for (c = 0; c < 3; ++c) {
                Q[a][b][c] = coef_nnn * n[a] * n[b] * n[c] +
                    coef_mix * (n[a] * (b == c ? 1.0 : 0.0) +
                        n[b] * (c == a ? 1.0 : 0.0) +
                        n[c] * (a == b ? 1.0 : 0.0));
            }
        }
    }
}

static void prj_rad_m1_phys_flux(double E, double F1, double F2, double F3,
    double *fE, double *fF1, double *fF2, double *fF3)
{
    double Fmag;
    double f;
    double chi;
    double n1;
    double n2;
    double n3;
    double D11;
    double D12;
    double D13;
    double c2;
    double cE;

    c2 = PRJ_CLIGHT * PRJ_CLIGHT;
    Fmag = sqrt(F1 * F1 + F2 * F2 + F3 * F3);
    cE = PRJ_CLIGHT * (E > 0.0 ? E : 0.0);

    if (cE <= 0.0 || Fmag <= 0.0) {
        /* isotropic: P = (E/3) I */
        *fE = F1;
        *fF1 = c2 * E / 3.0;
        *fF2 = 0.0;
        *fF3 = 0.0;
        return;
    }

    f = Fmag / cE;
    if (f > 1.0) {
        f = 1.0;
    }
    chi = (3.0 + 4.0 * f * f) / (5.0 + 2.0 * sqrt(4.0 - 3.0 * f * f));

    n1 = F1 / Fmag;
    n2 = F2 / Fmag;
    n3 = F3 / Fmag;

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

static void prj_rad_enforce_flux_limit(double *E, double *F1, double *F2, double *F3)
{
    double Fmag;
    double cE;
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
    }
}

static void prj_rad_m1_wavespeeds(double E, double F1, double F2, double F3,
    double *lam_min, double *lam_max)
{
    double Fmag;
    double cE;
    double f;
    double mu;
    double fsq;
    double ffac;
    double lterm;

    cE = PRJ_CLIGHT * (E > 0.0 ? E : 0.0);
    Fmag = sqrt(F1 * F1 + F2 * F2 + F3 * F3);

    if (cE <= 0.0 || Fmag <= 0.0) {
        *lam_min = -1.0 / sqrt(3.0);
        *lam_max = 1.0 / sqrt(3.0);
        return;
    }

    f = Fmag / cE;
    if (f > 1.0) {
        f = 1.0;
    }
    mu = F1 / Fmag;
    if (mu > 1.0) {
        mu = 1.0;
    } else if (mu < -1.0) {
        mu = -1.0;
    }

    fsq = f * f;
    ffac = sqrt(4.0 - 3.0 * fsq);
    lterm = sqrt(fabs((2.0 / 3.0) * (4.0 - 3.0 * fsq - ffac) + 2.0 * mu * mu * (2.0 - fsq - ffac)));
    *lam_min = (mu * f - lterm) / ffac;
    *lam_max = (mu * f + lterm) / ffac;
    if (*lam_min < -1.0) {
        *lam_min = -1.0;
    }
    if (*lam_max > 1.0) {
        *lam_max = 1.0;
    }
}
#endif

void prj_rad_flux(const double *WL, const double *WR,
    const prj_grav_mono *grav_mono, const double *x_face, const double *chi_face,
    double dx_dir, double v_face, double *flux)
{
    int field;
    int group;
    double lapse;

    if (WL == 0 || WR == 0 || flux == 0) {
        return;
    }

    if (grav_mono != 0 && x_face != 0) {
        double r = sqrt(x_face[0] * x_face[0] + x_face[1] * x_face[1] + x_face[2] * x_face[2]);
        lapse = prj_gravity_interp_lapse(grav_mono, r);
    } else {
        lapse = 1.0;
    }

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

                prj_rad_enforce_flux_limit(&EL, &F1L, &F2L, &F3L);
                prj_rad_enforce_flux_limit(&ER, &F1R, &F2R, &F3R);
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
                double chi_ext;
                double tau;
                double eps;

                prj_rad_m1_phys_flux(EL, F1L, F2L, F3L, &fLE, &fLF1, &fLF2, &fLF3);
                prj_rad_m1_phys_flux(ER, F1R, F2R, F3R, &fRE, &fRF1, &fRF2, &fRF3);

                prj_rad_m1_wavespeeds(EL, F1L, F2L, F3L, &lamL_min, &lamL_max);
                prj_rad_m1_wavespeeds(ER, F1R, F2R, F3R, &lamR_min, &lamR_max);
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

                chi_ext = chi_face != 0 ? chi_face[idx] : 0.0;
                tau = chi_ext * dx_dir;
                eps = 3.0 / (5.0 * tau + 1.0e-10);
                if (eps > 1.0) {
                    eps = 1.0;
                }

                flux[PRJ_CONS_RAD_E(field, group)] = lapse *
                    (sR * fLE - sL * fRE + eps * sL * sR * (ER - EL)) / denom;
                flux[PRJ_CONS_RAD_F1(field, group)] = lapse *
                    (sR * fLF1 - sL * fRF1 + sL * sR * (F1R - F1L)) / denom;
                flux[PRJ_CONS_RAD_F2(field, group)] = lapse *
                    (sR * fLF2 - sL * fRF2 + sL * sR * (F2R - F2L)) / denom;
                flux[PRJ_CONS_RAD_F3(field, group)] = lapse *
                    (sR * fLF3 - sL * fRF3 + sL * sR * (F3R - F3L)) / denom;

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
static void prj_rad_implicit_residuals(prj_rad *rad, prj_eos *eos, double *u,
    double dt, double lapse, double rho, double Uint_old, double Ye_old,
    const double *E_nu_old, double T, double Ye, double *F1, double *F2,
    double *E_nu_new_out)
{
    double kappa[PRJ_NRAD * PRJ_NEGROUP];
    double eta[PRJ_NRAD * PRJ_NEGROUP];
    double eos_q[PRJ_EOS_NQUANT];
    double eint_new;
    double Uint_new;
    double sum_dE = 0.0;
    double sum_dE_xe = 0.0;
    int nu;
    int g;

    (void)u;
    prj_rad3_opac_lookup(rad, rho, T, Ye, kappa, 0, 0, eta);
    prj_eos_rty(eos, rho, T, Ye, eos_q);
    eint_new = eos_q[PRJ_EOS_EINT];
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

    *F1 = Uint_new - Uint_old + sum_dE;
    *F2 = rho * Ye - rho * Ye_old + sum_dE_xe;
}

void prj_rad_energy_update(prj_rad *rad, prj_eos *eos, double *u, double dt, double lapse, double *final_temperature)
{
    double E_nu_old[PRJ_NRAD * PRJ_NEGROUP];
    double E_nu_new[PRJ_NRAD * PRJ_NEGROUP];
    double rho;
    double KE;
    double Uint_old;
    double Ye_old;
    double eint_old;
    double eos_q[PRJ_EOS_NQUANT];
    double T;
    double Ye;
    double err_scale_1;
    double err_scale_2;
    double res_cur;
    int iter;
    int nu;
    int g;
    const double alpha_ls = 1.0e-4;

    rho = u[PRJ_CONS_RHO];
    KE = 0.5 * (u[PRJ_CONS_MOM1] * u[PRJ_CONS_MOM1] +
        u[PRJ_CONS_MOM2] * u[PRJ_CONS_MOM2] +
        u[PRJ_CONS_MOM3] * u[PRJ_CONS_MOM3]) / rho;
    Uint_old = u[PRJ_CONS_ETOT] - KE;
    Ye_old = u[PRJ_CONS_YE] / rho;
    eint_old = Uint_old / rho;

    for (nu = 0; nu < PRJ_NRAD; ++nu) {
        for (g = 0; g < PRJ_NEGROUP; ++g) {
            E_nu_old[nu * PRJ_NEGROUP + g] = u[PRJ_CONS_RAD_E(nu, g)];
        }
    }

    prj_eos_rey(eos, rho, eint_old, Ye_old, eos_q);
    T = eos_q[PRJ_EOS_TEMPERATURE];
    Ye = Ye_old;
    err_scale_1 = fabs(Uint_old) > 0.0 ? Uint_old : 1.0;
    err_scale_2 = fabs(rho * Ye_old) > 0.0 ? rho * Ye_old : 1.0;
    res_cur = 1.0e30;

    for (iter = 0; iter < rad->maxiter; ++iter) {
        double F1;
        double F2;
        double F1_p;
        double F2_p;
        double f1;
        double f2;
        double hT;
        double hY;
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
        int inner_iter;

        prj_rad_implicit_residuals(rad, eos, u, dt, lapse, rho, Uint_old, Ye_old,
            E_nu_old, T, Ye, &F1, &F2, E_nu_new);
        f1 = F1 / err_scale_1;
        f2 = F2 / err_scale_2;
        res_cur = 0.5 * (f1 * f1 + f2 * f2);

        /* Jacobian via finite differences with power-of-2 step. */
        hT = 3.0e-8 * T;
        hT = pow(2.0, round(log(hT) / log(2.0)));
        hY = 3.0e-8 * Ye;
        hY = pow(2.0, round(log(hY) / log(2.0)));

        prj_rad_implicit_residuals(rad, eos, u, dt, lapse, rho, Uint_old, Ye_old,
            E_nu_old, T + hT, Ye, &F1_p, &F2_p, 0);
        dFdT_1 = (F1_p - F1) / hT;
        dFdT_2 = (F2_p - F2) / hT;

        prj_rad_implicit_residuals(rad, eos, u, dt, lapse, rho, Uint_old, Ye_old,
            E_nu_old, T, Ye + hY, &F1_p, &F2_p, 0);
        dFdY_1 = (F1_p - F1) / hY;
        dFdY_2 = (F2_p - F2) / hY;

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
            fprintf(stderr, "prj_rad_energy_update: singular Jacobian at iter=%d\n", iter);
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
        for (inner_iter = 0; inner_iter < 6; ++inner_iter) {
            double F1t;
            double F2t;
            double ft1;
            double ft2;

            Ttrial = T + lam * dT;
            Ytrial = Ye + lam * dY;
            prj_rad_implicit_residuals(rad, eos, u, dt, lapse, rho, Uint_old, Ye_old,
                E_nu_old, Ttrial, Ytrial, &F1t, &F2t, E_nu_new);
            ft1 = F1t / err_scale_1;
            ft2 = F2t / err_scale_2;
            res_trial = 0.5 * (ft1 * ft1 + ft2 * ft2);

            if (res_trial < rad->implicit_err_tol * rad->implicit_err_tol ||
                res_trial < res_cur + alpha_ls * lam * gradfdx) {
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
        }

        /* Convergence: residual norm AND step-size test. */
        if (res_trial < rad->implicit_err_tol * rad->implicit_err_tol &&
            fabs(lam * dT / T) < rad->implicit_err_tol &&
            fabs(lam * dY / Ye) < rad->implicit_err_tol) {
            res_cur = res_trial;
            break;
        }
        res_cur = res_trial;
    }
    if (iter == rad->maxiter) {
        fprintf(stderr, "prj_rad_energy_update: failed to converge (res=%e)\n", res_cur);
        exit(1);
    }

    /* Final pass at converged (T, Ye) to populate E_nu_new. */
    {
        double F1_final;
        double F2_final;
        double eint_new;

        prj_rad_implicit_residuals(rad, eos, u, dt, lapse, rho, Uint_old, Ye_old,
            E_nu_old, T, Ye, &F1_final, &F2_final, E_nu_new);
        prj_eos_rty(eos, rho, T, Ye, eos_q);
        eint_new = eos_q[PRJ_EOS_EINT];
        u[PRJ_CONS_ETOT] = rho * eint_new + KE;
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
}

void prj_rad_momentum_update(prj_rad *rad, prj_eos *eos, double *u, double dt, double lapse, double temperature)
{
    double kappa[PRJ_NRAD * PRJ_NEGROUP];
    double sigma[PRJ_NRAD * PRJ_NEGROUP];
    double delta[PRJ_NRAD * PRJ_NEGROUP];
    double rho;
    double Ye;
    double v[3];
    double inv_c2;
    double dmom[3];
    double detot;
    int nu;
    int g;
    int d;

    (void)eos;

    rho = u[PRJ_CONS_RHO];
    Ye = u[PRJ_CONS_YE] / rho;
    v[0] = u[PRJ_CONS_MOM1] / rho;
    v[1] = u[PRJ_CONS_MOM2] / rho;
    v[2] = u[PRJ_CONS_MOM3] / rho;
    inv_c2 = 1.0 / (PRJ_CLIGHT * PRJ_CLIGHT);

    prj_rad3_opac_lookup(rad, rho, temperature, Ye, kappa, sigma, delta, 0);

    dmom[0] = 0.0;
    dmom[1] = 0.0;
    dmom[2] = 0.0;
    detot = 0.0;

    for (nu = 0; nu < PRJ_NRAD; ++nu) {
        for (g = 0; g < PRJ_NEGROUP; ++g) {
            int idx = nu * PRJ_NEGROUP + g;
            double chi = kappa[idx] + sigma[idx] * (1.0 - delta[idx] / 3.0);
            double factor = 1.0 / (1.0 + dt * lapse * PRJ_CLIGHT * chi) - 1.0;
            int fi[3];

            fi[0] = PRJ_CONS_RAD_F1(nu, g);
            fi[1] = PRJ_CONS_RAD_F2(nu, g);
            fi[2] = PRJ_CONS_RAD_F3(nu, g);
            for (d = 0; d < 3; ++d) {
                double F_old = u[fi[d]];
                double dF = F_old * factor;
                u[fi[d]] = F_old + dF;
                dmom[d] += dF * inv_c2;
                detot += v[d] * dF * inv_c2;
            }
        }
    }

    u[PRJ_CONS_MOM1] -= dmom[0];
    u[PRJ_CONS_MOM2] -= dmom[1];
    u[PRJ_CONS_MOM3] -= dmom[2];
    u[PRJ_CONS_ETOT] -= detot;

    for (nu = 0; nu < PRJ_NRAD; ++nu) {
        for (g = 0; g < PRJ_NEGROUP; ++g) {
            double E = u[PRJ_CONS_RAD_E(nu, g)];
            double F1 = u[PRJ_CONS_RAD_F1(nu, g)];
            double F2 = u[PRJ_CONS_RAD_F2(nu, g)];
            double F3 = u[PRJ_CONS_RAD_F3(nu, g)];
            double Fmag = sqrt(F1 * F1 + F2 * F2 + F3 * F3);
            double cE = PRJ_CLIGHT * E;

            if (E > 0.0 && Fmag > cE) {
                double scale = 0.99999 * cE / Fmag;
                u[PRJ_CONS_RAD_F1(nu, g)] = F1 * scale;
                u[PRJ_CONS_RAD_F2(nu, g)] = F2 * scale;
                u[PRJ_CONS_RAD_F3(nu, g)] = F3 * scale;
            }
        }
    }
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
    double dqp;
    double dqm;
    double r;
    double slope;

    if (gcell <= 0 || gcell >= PRJ_NEGROUP - 1) {
        return q[gcell];
    }
    dqp = q[gcell + 1] - q[gcell];
    dqm = q[gcell] - q[gcell - 1];
    if (dqp == 0.0) {
        slope = 0.0;
    } else {
        r = dqm / dqp;
        slope = prj_rad_koren_phi(r) * dqp;
    }
    return q[gcell] + 0.5 * (double)side * slope;
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

    for (field = 0; field < PRJ_NRAD; ++field) {
        double Eg[PRJ_NEGROUP];
        double Fg[PRJ_NEGROUP][3];
        double Pg[PRJ_NEGROUP][3][3];
        double Qg[PRJ_NEGROUP][3][3][3];
        double dE_acc[PRJ_NEGROUP];
        double dF_acc[PRJ_NEGROUP][3];
        const double *nu_face = rad->eedge[field];
        int g;
        int ii;
        int jj;
        int kk;

        /* Per-group cell-centred state and M1 pressure tensor. */
        for (g = 0; g < PRJ_NEGROUP; ++g) {
            Eg[g] = W_state[VIDX(PRJ_PRIM_RAD_E(field, g), ic, jc, kc)];
            Fg[g][0] = W_state[VIDX(PRJ_PRIM_RAD_F1(field, g), ic, jc, kc)];
            Fg[g][1] = W_state[VIDX(PRJ_PRIM_RAD_F2(field, g), ic, jc, kc)];
            Fg[g][2] = W_state[VIDX(PRJ_PRIM_RAD_F3(field, g), ic, jc, kc)];
            prj_rad_m1_pressure(Eg[g], Fg[g][0], Fg[g][1], Fg[g][2], Pg[g]);
            prj_rad_m1_third_moment(Eg[g], Fg[g][0], Fg[g][1], Fg[g][2], Qg[g]);
            dE_acc[g] = 0.0;
            dF_acc[g][0] = 0.0;
            dF_acc[g][1] = 0.0;
            dF_acc[g][2] = 0.0;
        }

        /* Energy: dE_g/dt += Σ_{i,j} v^i_{;j} · [(ν P^{ji})_{g+1/2} - (ν P^{ji})_{g-1/2}]. */
        for (jj = 0; jj < 3; ++jj) {
            for (ii = 0; ii < 3; ++ii) {
                double vij = dvdx[jj][ii];
                double q[PRJ_NEGROUP];
                double face_val[PRJ_NEGROUP + 1];
                int gf;

                if (vij == 0.0) {
                    continue;
                }
                for (g = 0; g < PRJ_NEGROUP; ++g) {
                    q[g] = Pg[g][jj][ii];
                }

                face_val[0] = 0.0;
                face_val[PRJ_NEGROUP] = 0.0;
                for (gf = 1; gf < PRJ_NEGROUP; ++gf) {
                    /* Face gf sits between cells gf-1 (lower ν) and gf (upper ν). */
                    double pick = (vij >= 0.0)
                        ? prj_rad_recon_face(q, gf, -1)      /* R: left edge of cell gf  */
                        : prj_rad_recon_face(q, gf - 1, +1); /* L: right edge of cell gf-1 */
                    face_val[gf] = nu_face[gf] * pick;
                }

                for (g = 0; g < PRJ_NEGROUP; ++g) {
                    dE_acc[g] += vij * (face_val[g + 1] - face_val[g]);
                }
            }
        }

        /* Flux: dF_{gj}/dt += Σ_{i,k} v^i_{;k} · [(ν Q^{kji})_{g+1/2} - (ν Q^{kji})_{g-1/2}].
         * Closure: Q^{kji}_g from the Levermore/Vaytet third moment. */
        for (jj = 0; jj < 3; ++jj) {
            for (kk = 0; kk < 3; ++kk) {
                for (ii = 0; ii < 3; ++ii) {
                    double vik = dvdx[kk][ii];
                    double q[PRJ_NEGROUP];
                    double face_val[PRJ_NEGROUP + 1];
                    int gf;

                    if (vik == 0.0) {
                        continue;
                    }
                    for (g = 0; g < PRJ_NEGROUP; ++g) {
                        q[g] = Qg[g][kk][jj][ii];
                    }

                    face_val[0] = 0.0;
                    face_val[PRJ_NEGROUP] = 0.0;
                    for (gf = 1; gf < PRJ_NEGROUP; ++gf) {
                        double pick = (vik >= 0.0)
                            ? prj_rad_recon_face(q, gf, -1)
                            : prj_rad_recon_face(q, gf - 1, +1);
                        face_val[gf] = nu_face[gf] * pick;
                    }

                    for (g = 0; g < PRJ_NEGROUP; ++g) {
                        dF_acc[g][jj] += vik * (face_val[g + 1] - face_val[g]);
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
            const prj_grav_mono *gm = prj_gravity_active_monopole();
            double xc1 = block->xmin[0] + ((double)ic + 0.5) * block->dx[0];
            double xc2 = block->xmin[1] + ((double)jc + 0.5) * block->dx[1];
            double xc3 = block->xmin[2] + ((double)kc + 0.5) * block->dx[2];
            double r_cell = sqrt(xc1 * xc1 + xc2 * xc2 + xc3 * xc3);
            double grad_phi[3];
            double inv_c2 = 1.0 / (PRJ_CLIGHT * PRJ_CLIGHT);

            grad_phi[0] = 0.0;
            grad_phi[1] = 0.0;
            grad_phi[2] = 0.0;
            if (gm != 0 && r_cell > 0.0) {
                double accel = prj_gravity_interp_accel(gm, r_cell);

                grad_phi[0] = -accel * xc1 / r_cell;
                grad_phi[1] = -accel * xc2 / r_cell;
                grad_phi[2] = -accel * xc3 / r_cell;
            }

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
                for (g = 0; g < PRJ_NEGROUP; ++g) {
                    dE_acc[g] += coef * (face_val[g + 1] - face_val[g]);
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
                    for (g = 0; g < PRJ_NEGROUP; ++g) {
                        dF_acc[g][jj] += coef * (face_val[g + 1] - face_val[g]);
                    }
                }
            }
        }

        /* Apply.  dt is the effective stage weight (full dt in stage1, 0.5·dt
         * in stage2 to match the RK2-Heun mixing of dUdt).  The lapse factor
         * α(r) accounts for the GR proper-time slowdown in the gravitational
         * well, consistent with the lapse multipliers already on the spatial
         * radiation flux and on the gravity source. */
        for (g = 0; g < PRJ_NEGROUP; ++g) {
            u[PRJ_CONS_RAD_E(field, g)] += lapse * dt * dE_acc[g];
            u[PRJ_CONS_RAD_F1(field, g)] += lapse * dt * dF_acc[g][0];
            u[PRJ_CONS_RAD_F2(field, g)] += lapse * dt * dF_acc[g][1];
            u[PRJ_CONS_RAD_F3(field, g)] += lapse * dt * dF_acc[g][2];
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
void prj_rad_energy_update(prj_rad *rad, prj_eos *eos, double *u, double dt, double lapse, double *final_temperature)
{
    (void)rad;
    (void)eos;
    (void)u;
    (void)dt;
    (void)final_temperature;
    (void)lapse;
}

void prj_rad_momentum_update(prj_rad *rad, prj_eos *eos, double *u, double dt, double lapse, double temperature)
{
    (void)rad;
    (void)eos;
    (void)u;
    (void)dt;
    (void)lapse;
    (void)temperature;
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
