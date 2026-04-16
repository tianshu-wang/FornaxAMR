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

void prj_rad_flux(const double *WL, const double *WR, const prj_eos *eos, const prj_rad *rad,
    const prj_grav_mono *grav_mono, const double *x_face, double dx_dir, double *flux)
{
    int field;
    int group;
    double lapse;

    (void)eos;

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
        double kappa_avg[PRJ_NRAD * PRJ_NEGROUP];
        double sigma_avg[PRJ_NRAD * PRJ_NEGROUP];
        double delta_avg[PRJ_NRAD * PRJ_NEGROUP];
        double eta_avg[PRJ_NRAD * PRJ_NEGROUP];
        double rho_avg;
        double eint_avg;
        double ye_avg;
        double T_avg;
        double eos_q[PRJ_EOS_NQUANT];

        rho_avg = 0.5 * (WL[PRJ_PRIM_RHO] + WR[PRJ_PRIM_RHO]);
        eint_avg = 0.5 * (WL[PRJ_PRIM_EINT] + WR[PRJ_PRIM_EINT]);
        ye_avg = 0.5 * (WL[PRJ_PRIM_YE] + WR[PRJ_PRIM_YE]);
        prj_eos_rey((prj_eos *)eos, rho_avg, eint_avg, ye_avg, eos_q);
        T_avg = eos_q[PRJ_EOS_TEMPERATURE];
        prj_rad3_opac_lookup(rad, rho_avg, T_avg, ye_avg, kappa_avg, sigma_avg, delta_avg, eta_avg);

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

                chi_ext = kappa_avg[idx] + sigma_avg[idx];
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
            }
        }
    }
#else
    (void)rad;
    (void)dx_dir;
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
    double sigma[PRJ_NRAD * PRJ_NEGROUP];
    double delta[PRJ_NRAD * PRJ_NEGROUP];
    double eta[PRJ_NRAD * PRJ_NEGROUP];
    double eos_q[PRJ_EOS_NQUANT];
    double eint_new;
    double Uint_new;
    double sum_dE = 0.0;
    double sum_dE_xe = 0.0;
    int nu;
    int g;

    (void)u;
    prj_rad3_opac_lookup(rad, rho, T, Ye, kappa, sigma, delta, eta);
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
    double eta[PRJ_NRAD * PRJ_NEGROUP];
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

    prj_rad3_opac_lookup(rad, rho, temperature, Ye, kappa, sigma, delta, eta);

    dmom[0] = 0.0;
    dmom[1] = 0.0;
    dmom[2] = 0.0;
    detot = 0.0;

    for (nu = 0; nu < PRJ_NRAD; ++nu) {
        for (g = 0; g < PRJ_NEGROUP; ++g) {
            int idx = nu * PRJ_NEGROUP + g;
            double chi = kappa[idx] + sigma[idx] * (1.0 - delta[idx] / 3.0);
            double factor = 1.0 / (1.0 + dt * lapse * chi) - 1.0;
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
#endif
