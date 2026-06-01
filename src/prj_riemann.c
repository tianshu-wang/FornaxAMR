#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "prj.h"

static double prj_riemann_theta_limiter(double cfmax,
    double deltau, double deltav, double deltaw)
{
    double theta_r;

    theta_r = (-PRJ_MIN(deltau, 0.0) + cfmax) /
        (-PRJ_MIN(PRJ_MIN(deltav, deltaw), 0.0) + cfmax);
    theta_r *= theta_r;
    return PRJ_MIN(1.0, theta_r * theta_r);
}

#if PRJ_MHD
#define PRJ_HLLD_SMALL_NUMBER 1.0e-8
#ifndef PRJ_RIEMANN_DEBUG
#define PRJ_RIEMANN_DEBUG 0
#endif

typedef struct prj_hlld_state {
    double rho;
    double vx;
    double vy;
    double vz;
    double bx;
    double by;
    double bz;
    double p;
    double gamma;
    double ye;
    double e;
    double pt;
    double U[PRJ_NHYDRO];
    double F[PRJ_NHYDRO];
} prj_hlld_state;

typedef struct prj_hlld_star {
    double rho;
    double vx;
    double vy;
    double vz;
    double bx;
    double by;
    double bz;
    double ye;
    double e;
    double U[PRJ_NHYDRO];
} prj_hlld_star;

static void prj_riemann_hlld_fail(const char *message)
{
    fprintf(stderr, "prj_riemann_hlld: %s\n", message);
    exit(1);
}

#if PRJ_RIEMANN_DEBUG
static void prj_hlld_require_finite(double x, const char *name)
{
    if (!isfinite(x)) {
        fprintf(stderr, "prj_riemann_hlld: non-finite %s\n", name);
        exit(1);
    }
}
#else
#define prj_hlld_require_finite(x, name) ((void)(x))
#endif

static void prj_hlld_fill_conserved(double rho, double vx, double vy, double vz,
    double bx, double by, double bz, double ye, double e, double *U)
{
    U[PRJ_CONS_RHO] = rho;
    U[PRJ_CONS_MOM1] = rho * vx;
    U[PRJ_CONS_MOM2] = rho * vy;
    U[PRJ_CONS_MOM3] = rho * vz;
    U[PRJ_CONS_ETOT] = e;
    U[PRJ_CONS_YE] = rho * ye;
    U[PRJ_CONS_B1] = bx;
    U[PRJ_CONS_B2] = by;
    U[PRJ_CONS_B3] = bz;
}

static void prj_hlld_fill_flux(double rho, double vx, double vy, double vz,
    double bx, double by, double bz, double p, double ye, double e, double *F)
{
    double b2;
    double pt;
    double vdotb;

    b2 = bx * bx + by * by + bz * bz;
    pt = p + 0.5 * b2;
    vdotb = vx * bx + vy * by + vz * bz;

    F[PRJ_CONS_RHO] = rho * vx;
    F[PRJ_CONS_MOM1] = rho * vx * vx + pt - bx * bx;
    F[PRJ_CONS_MOM2] = rho * vx * vy - bx * by;
    F[PRJ_CONS_MOM3] = rho * vx * vz - bx * bz;
    F[PRJ_CONS_ETOT] = (e + pt) * vx - bx * vdotb;
    F[PRJ_CONS_YE] = rho * ye * vx;
    F[PRJ_CONS_B1] = 0.0;
    F[PRJ_CONS_B2] = by * vx - bx * vy;
    F[PRJ_CONS_B3] = bz * vx - bx * vz;
}

static double prj_hlld_fast_speed(const prj_hlld_state *s)
{
    double a2;
    double bx2_rho;
    double b2_rho;
    double disc;
    double cf2;

    a2 = s->gamma * s->p / s->rho;
    bx2_rho = s->bx * s->bx / s->rho;
    b2_rho = (s->bx * s->bx + s->by * s->by + s->bz * s->bz) / s->rho;
    disc = (a2 + b2_rho) * (a2 + b2_rho) - 4.0 * a2 * bx2_rho;
    if (disc < 0.0) {
        if (disc < -1.0e-12 * (a2 + b2_rho) * (a2 + b2_rho)) {
            prj_riemann_hlld_fail("negative fast-speed discriminant");
        }
        disc = 0.0;
    }
    cf2 = 0.5 * (a2 + b2_rho + sqrt(disc));
    if (cf2 < 0.0 || !isfinite(cf2)) {
        prj_riemann_hlld_fail("invalid fast speed");
    }
    return sqrt(cf2);
}

static void prj_hlld_state_from_prim(const double *W, double p, double gamma,
    double bn, prj_hlld_state *s)
{
    double kinetic;
    double magnetic;

    s->rho = W[PRJ_PRIM_RHO];
    s->vx = W[PRJ_PRIM_V1];
    s->vy = W[PRJ_PRIM_V2];
    s->vz = W[PRJ_PRIM_V3];
    s->bx = bn;
    s->by = W[PRJ_PRIM_B2];
    s->bz = W[PRJ_PRIM_B3];
    s->p = p;
    s->gamma = gamma;
    s->ye = W[PRJ_PRIM_YE];

    if (s->rho <= 0.0) {
        prj_riemann_hlld_fail("non-positive density");
    }
    if (s->p <= 0.0) {
        prj_riemann_hlld_fail("non-positive gas pressure");
    }
    if (s->gamma <= 0.0) {
        prj_riemann_hlld_fail("non-positive adiabatic index");
    }

    kinetic = 0.5 * s->rho * (s->vx * s->vx + s->vy * s->vy + s->vz * s->vz);
    magnetic = 0.5 * (s->bx * s->bx + s->by * s->by + s->bz * s->bz);
    s->e = s->rho * W[PRJ_PRIM_EINT] + kinetic + magnetic;
    s->pt = s->p + magnetic;

    prj_hlld_require_finite(s->e, "total energy");
    prj_hlld_fill_conserved(s->rho, s->vx, s->vy, s->vz, s->bx, s->by, s->bz,
        s->ye, s->e, s->U);
    prj_hlld_fill_flux(s->rho, s->vx, s->vy, s->vz, s->bx, s->by, s->bz,
        s->p, s->ye, s->e, s->F);
}

static void prj_hlld_fill_star_conserved(prj_hlld_star *s)
{
    prj_hlld_fill_conserved(s->rho, s->vx, s->vy, s->vz, s->bx, s->by, s->bz,
        s->ye, s->e, s->U);
}

static void prj_hlld_outer_star(const prj_hlld_state *s, double S, double SM,
    double pt_star, prj_hlld_star *star)
{
    double alpha;
    double smdiff;
    double denom;
    double vdotb;
    double vstar_dotbstar;

    alpha = S - s->vx;
    smdiff = S - SM;
    if (fabs(smdiff) <= 1.0e-14 * (fabs(S) + fabs(SM) + 1.0)) {
        prj_riemann_hlld_fail("outer wave coincides with contact wave");
    }

    star->rho = s->rho * alpha / smdiff;
    if (star->rho <= 0.0 || !isfinite(star->rho)) {
        prj_riemann_hlld_fail("invalid star density");
    }
    star->vx = SM;
    star->bx = s->bx;
    star->ye = s->ye;

    denom = s->rho * alpha * smdiff - s->bx * s->bx;
    if (fabs(denom) < PRJ_HLLD_SMALL_NUMBER * pt_star) {
        star->vy = s->vy;
        star->vz = s->vz;
        star->by = s->by;
        star->bz = s->bz;
    } else {
        double sm_diff_vx_over_denom = (SM - s->vx) / denom;
        double mag_factor = (s->rho * alpha * alpha - s->bx * s->bx) / denom;
        star->vy = s->vy - s->bx * s->by * sm_diff_vx_over_denom;
        star->vz = s->vz - s->bx * s->bz * sm_diff_vx_over_denom;
        star->by = s->by * mag_factor;
        star->bz = s->bz * mag_factor;
    }

    vdotb = s->vx * s->bx + s->vy * s->by + s->vz * s->bz;
    vstar_dotbstar = star->vx * star->bx + star->vy * star->by + star->vz * star->bz;
    star->e = (alpha * s->e - s->pt * s->vx + pt_star * SM +
        s->bx * (vdotb - vstar_dotbstar)) / smdiff;
    prj_hlld_require_finite(star->e, "star total energy");

    prj_hlld_fill_star_conserved(star);
}

static void prj_hlld_double_star_one(const prj_hlld_star *left, const prj_hlld_star *right,
    double bx, int right_side, prj_hlld_star *ss)
{
    double sqrt_l;
    double sqrt_r;
    double inv_sum;
    double sign_bx;
    double vy_ss;
    double vz_ss;
    double by_ss;
    double bz_ss;
    double vbdot_ss;

    sqrt_l = sqrt(left->rho);
    sqrt_r = sqrt(right->rho);
    if (sqrt_l + sqrt_r <= 0.0) {
        prj_riemann_hlld_fail("invalid double-star density weights");
    }
    inv_sum = 1.0 / (sqrt_l + sqrt_r);
    sign_bx = bx < 0.0 ? -1.0 : 1.0;

    vy_ss = (sqrt_l * left->vy + sqrt_r * right->vy +
        (right->by - left->by) * sign_bx) * inv_sum;
    vz_ss = (sqrt_l * left->vz + sqrt_r * right->vz +
        (right->bz - left->bz) * sign_bx) * inv_sum;
    by_ss = (sqrt_l * right->by + sqrt_r * left->by +
        sqrt_l * sqrt_r * (right->vy - left->vy) * sign_bx) * inv_sum;
    bz_ss = (sqrt_l * right->bz + sqrt_r * left->bz +
        sqrt_l * sqrt_r * (right->vz - left->vz) * sign_bx) * inv_sum;

    vbdot_ss = (right_side ? right->vx : left->vx) * bx + vy_ss * by_ss + vz_ss * bz_ss;

    if (!right_side) {
        double vbdot_l = left->vx * bx + left->vy * left->by + left->vz * left->bz;
        *ss = *left;
        ss->e = left->e - sqrt_l * sign_bx * (vbdot_l - vbdot_ss);
        prj_hlld_require_finite(ss->e, "left double-star total energy");
    } else {
        double vbdot_r = right->vx * bx + right->vy * right->by + right->vz * right->bz;
        *ss = *right;
        ss->e = right->e + sqrt_r * sign_bx * (vbdot_r - vbdot_ss);
        prj_hlld_require_finite(ss->e, "right double-star total energy");
    }
    ss->vy = vy_ss;
    ss->vz = vz_ss;
    ss->by = by_ss;
    ss->bz = bz_ss;
    prj_hlld_fill_star_conserved(ss);
}

static void prj_hlld_flux_from_jump(const double *F0, double S,
    const double *U1, const double *U0, double *F)
{
    int v;

    for (v = 0; v < PRJ_NHYDRO; ++v) {
        F[v] = F0[v] + S * (U1[v] - U0[v]);
    }
}

static void prj_hlld_flux_from_double_jump(const double *F0, double S0,
    const double *Ustar, const double *U0, double Sstar,
    const double *Uss, double *F)
{
    int v;

    for (v = 0; v < PRJ_NHYDRO; ++v) {
        F[v] = F0[v] + S0 * (Ustar[v] - U0[v]) + Sstar * (Uss[v] - Ustar[v]);
    }
}

static void prj_hlld_face_outputs(double vx, double vy, double vz,
    double bx, double by, double bz, double v_face[3], double *Bv1, double *Bv2)
{
    if (v_face != 0) {
        v_face[0] = vx;
        v_face[1] = vy;
        v_face[2] = vz;
    }
    if (Bv1 != 0) {
        *Bv1 = bz * vx - bx * vz;
    }
    if (Bv2 != 0) {
        *Bv2 = bx * vy - by * vx;
    }
}

void prj_riemann_hlld(const double *WL, const double *WR,
    double pL, double pR, double gL, double gR,
    const prj_eos *eos, double bn, double *flux, double v_face[3],
    double *Bv1, double *Bv2, double deltau, double deltav, double deltaw)
{
    prj_hlld_state L;
    prj_hlld_state R;
    prj_hlld_star Ls;
    prj_hlld_star Rs;
    prj_hlld_star ss;
    double cfL;
    double cfR;
    double SL;
    double SR;
    double SM;
    double SLs;
    double SRs;
    double denom;
    double pt_star_l;
    double pt_star_r;
    double pt_star;
    int v;
    int bn_small;

    /* HLLD owns only hydro/MHD fluxes; radiation fluxes are built separately. */
    (void)eos;

    if (WL == 0 || WR == 0 || flux == 0) {
        prj_riemann_hlld_fail("null input");
    }

    prj_hlld_state_from_prim(WL, pL, gL, bn, &L);
    prj_hlld_state_from_prim(WR, pR, gR, bn, &R);
    cfL = prj_hlld_fast_speed(&L);
    cfR = prj_hlld_fast_speed(&R);
    SL = PRJ_MIN(L.vx - cfL, R.vx - cfR);
    SR = PRJ_MAX(L.vx + cfL, R.vx + cfR);
    if (!(SL < SR)) {
        prj_riemann_hlld_fail("invalid fast-wave ordering");
    }

    if (0.0 <= SL) {
        for (v = 0; v < PRJ_NHYDRO; ++v) {
            flux[v] = L.F[v];
        }
        prj_hlld_face_outputs(L.vx, L.vy, L.vz, L.bx, L.by, L.bz, v_face, Bv1, Bv2);
        return;
    }
    if (SR <= 0.0) {
        for (v = 0; v < PRJ_NHYDRO; ++v) {
            flux[v] = R.F[v];
        }
        prj_hlld_face_outputs(R.vx, R.vy, R.vz, R.bx, R.by, R.bz, v_face, Bv1, Bv2);
        return;
    }

    denom = R.rho * (SR - R.vx) - L.rho * (SL - L.vx);
    if (fabs(denom) <= 1.0e-14 *
        (fabs(R.rho * (SR - R.vx)) + fabs(L.rho * (SL - L.vx)) + 1.0)) {
        prj_riemann_hlld_fail("degenerate contact-speed denominator");
    }
    double cfmax = PRJ_MAX(cfL,cfR);
    double theta = prj_riemann_theta_limiter(cfmax, deltau, deltav, deltaw);
    SM = (R.rho * R.vx * (SR - R.vx) -
        L.rho * L.vx * (SL - L.vx) + theta*(L.pt - R.pt)) / denom;
    prj_hlld_require_finite(SM, "contact speed");

    pt_star_l = L.pt + L.rho * (SL - L.vx) * (SM - L.vx);
    pt_star_r = R.pt + R.rho * (SR - R.vx) * (SM - R.vx);
    pt_star = 0.5 * (pt_star_l + pt_star_r);
    if (pt_star <= 0.0 || !isfinite(pt_star)) {
        prj_riemann_hlld_fail("invalid star total pressure");
    }
    bn_small = (0.5 * bn * bn < PRJ_HLLD_SMALL_NUMBER * pt_star);

    if (0.0 <= SM) {
        prj_hlld_outer_star(&L, SL, SM, pt_star, &Ls);
        SLs = SM - fabs(bn) / sqrt(Ls.rho);
        prj_hlld_require_finite(SLs, "left Alfven speed");
        if (0.0 <= SLs || bn_small) {
            prj_hlld_flux_from_jump(L.F, SL, Ls.U, L.U, flux);
            prj_hlld_face_outputs(Ls.vx, Ls.vy, Ls.vz, Ls.bx, Ls.by, Ls.bz, v_face, Bv1, Bv2);
        } else {
            prj_hlld_outer_star(&R, SR, SM, pt_star, &Rs);
            prj_hlld_double_star_one(&Ls, &Rs, bn, 0, &ss);
            prj_hlld_flux_from_double_jump(L.F, SL, Ls.U, L.U, SLs, ss.U, flux);
            prj_hlld_face_outputs(ss.vx, ss.vy, ss.vz, ss.bx, ss.by, ss.bz, v_face, Bv1, Bv2);
        }
    } else {
        prj_hlld_outer_star(&R, SR, SM, pt_star, &Rs);
        SRs = SM + fabs(bn) / sqrt(Rs.rho);
        prj_hlld_require_finite(SRs, "right Alfven speed");
        if (SRs <= 0.0 || bn_small) {
            prj_hlld_flux_from_jump(R.F, SR, Rs.U, R.U, flux);
            prj_hlld_face_outputs(Rs.vx, Rs.vy, Rs.vz, Rs.bx, Rs.by, Rs.bz, v_face, Bv1, Bv2);
        } else {
            prj_hlld_outer_star(&L, SL, SM, pt_star, &Ls);
            prj_hlld_double_star_one(&Ls, &Rs, bn, 1, &ss);
            prj_hlld_flux_from_double_jump(R.F, SR, Rs.U, R.U, SRs, ss.U, flux);
            prj_hlld_face_outputs(ss.vx, ss.vy, ss.vz, ss.bx, ss.by, ss.bz, v_face, Bv1, Bv2);
        }
    }
}
#endif

static void prj_riemann_state(const double *W, double pressure, double gamma,
    const prj_eos *eos, double *U, double *F, double *cs)
{
    double rho;
    double v1;
    double v2;
    double v3;
    double eint;
    double etot;
#if PRJ_MHD
    double b1;
    double b2;
    double b3;
#endif
    int v;

    for (v = 0; v < PRJ_NHYDRO; ++v) {
        U[v] = 0.0;
        F[v] = 0.0;
    }

    /* HLLC owns only hydro/MHD fluxes; radiation fluxes are built separately. */
    (void)eos;
    rho = W[PRJ_PRIM_RHO];
    v1 = W[PRJ_PRIM_V1];
    v2 = W[PRJ_PRIM_V2];
    v3 = W[PRJ_PRIM_V3];
    eint = W[PRJ_PRIM_EINT];
#if PRJ_MHD
    b1 = W[PRJ_PRIM_B1];
    b2 = W[PRJ_PRIM_B2];
    b3 = W[PRJ_PRIM_B3];
#endif
    etot = rho * eint + 0.5 * rho * (v1 * v1 + v2 * v2 + v3 * v3)
#if PRJ_MHD
        + 0.5 * (b1 * b1 + b2 * b2 + b3 * b3)
#endif
        ;

    U[PRJ_CONS_RHO] = rho;
    U[PRJ_CONS_MOM1] = rho * v1;
    U[PRJ_CONS_MOM2] = rho * v2;
    U[PRJ_CONS_MOM3] = rho * v3;
    U[PRJ_CONS_ETOT] = etot;
    U[PRJ_CONS_YE] = rho * W[PRJ_PRIM_YE];
#if PRJ_MHD
    U[PRJ_CONS_B1] = b1;
    U[PRJ_CONS_B2] = b2;
    U[PRJ_CONS_B3] = b3;
#endif

    *cs = sqrt(gamma * pressure / rho);

    F[PRJ_CONS_RHO] = rho * v1;
    F[PRJ_CONS_MOM1] = rho * v1 * v1 + pressure;
    F[PRJ_CONS_MOM2] = rho * v1 * v2;
    F[PRJ_CONS_MOM3] = rho * v1 * v3;
    F[PRJ_CONS_ETOT] = (etot + pressure) * v1;
    F[PRJ_CONS_YE] = rho * W[PRJ_PRIM_YE] * v1;
}

void prj_riemann_hllc(const double *WL, const double *WR,
    double pL, double pR, double gL, double gR,
    const prj_eos *eos, double *flux, double v_face[3],
    double deltau, double deltav, double deltaw)
{
    double UL[PRJ_NHYDRO];
    double UR[PRJ_NHYDRO];
    double FL[PRJ_NHYDRO];
    double FR[PRJ_NHYDRO];
    double Ustar[PRJ_NHYDRO];
    double csL;
    double csR;
    double SL;
    double SR;
    double SM;
    double cfmax;
    double theta;
    double rho_star;
    double e_over_rho;
    int v;

    prj_riemann_state(WL, pL, gL, eos, UL, FL, &csL);
    prj_riemann_state(WR, pR, gR, eos, UR, FR, &csR);

    SL = PRJ_MIN(WL[PRJ_PRIM_V1] - csL, WR[PRJ_PRIM_V1] - csR);
    SR = PRJ_MAX(WL[PRJ_PRIM_V1] + csL, WR[PRJ_PRIM_V1] + csR);
    cfmax = PRJ_MAX(csL, csR);
    theta = prj_riemann_theta_limiter(cfmax, deltau, deltav, deltaw);

    /* Toro (2009), Eq. 10.37 with the HLLD theta pressure-jump limiter. */
    SM = (theta * (pR - pL) +
        UL[PRJ_CONS_MOM1] * (SL - WL[PRJ_PRIM_V1]) -
        UR[PRJ_CONS_MOM1] * (SR - WR[PRJ_PRIM_V1])) /
        (UL[PRJ_CONS_RHO] * (SL - WL[PRJ_PRIM_V1]) -
            UR[PRJ_CONS_RHO] * (SR - WR[PRJ_PRIM_V1]));

    if (0.0 <= SL) {
        for (v = 0; v < PRJ_NHYDRO; ++v) {
            flux[v] = FL[v];
        }
        if (v_face != 0) {
            v_face[0] = WL[PRJ_PRIM_V1];
            v_face[1] = WL[PRJ_PRIM_V2];
            v_face[2] = WL[PRJ_PRIM_V3];
        }
        return;
    }
    if (SR <= 0.0) {
        for (v = 0; v < PRJ_NHYDRO; ++v) {
            flux[v] = FR[v];
        }
        if (v_face != 0) {
            v_face[0] = WR[PRJ_PRIM_V1];
            v_face[1] = WR[PRJ_PRIM_V2];
            v_face[2] = WR[PRJ_PRIM_V3];
        }
        return;
    }

    if (v_face != 0) {
        /* HLLC star state: contact normal speed SM, transverse velocities preserved from upwind side. */
        v_face[0] = SM;
        if (0.0 <= SM) {
            v_face[1] = WL[PRJ_PRIM_V2];
            v_face[2] = WL[PRJ_PRIM_V3];
        } else {
            v_face[1] = WR[PRJ_PRIM_V2];
            v_face[2] = WR[PRJ_PRIM_V3];
        }
    }
    if (0.0 <= SM) {
        /* Toro (2009), Eqs. 10.33-10.36: left star state. */
        for (v = 0; v < PRJ_NHYDRO; ++v) {
            Ustar[v] = UL[v];
        }
        rho_star = UL[PRJ_CONS_RHO] * (SL - WL[PRJ_PRIM_V1]) / (SL - SM);
        Ustar[PRJ_CONS_RHO] = rho_star;
        Ustar[PRJ_CONS_MOM1] = rho_star * SM;
        Ustar[PRJ_CONS_MOM2] = rho_star * WL[PRJ_PRIM_V2];
        Ustar[PRJ_CONS_MOM3] = rho_star * WL[PRJ_PRIM_V3];
        e_over_rho = UL[PRJ_CONS_ETOT] / UL[PRJ_CONS_RHO] +
            (SM - WL[PRJ_PRIM_V1]) * (SM + pL / (UL[PRJ_CONS_RHO] * (SL - WL[PRJ_PRIM_V1])));
        Ustar[PRJ_CONS_ETOT] = rho_star * e_over_rho;
        Ustar[PRJ_CONS_YE] = rho_star * WL[PRJ_PRIM_YE];
        for (v = 0; v < PRJ_NHYDRO; ++v) {
            flux[v] = FL[v] + SL * (Ustar[v] - UL[v]);
        }
    } else {
        /* Toro (2009), Eqs. 10.33-10.36 with right star state. */
        for (v = 0; v < PRJ_NHYDRO; ++v) {
            Ustar[v] = UR[v];
        }
        rho_star = UR[PRJ_CONS_RHO] * (SR - WR[PRJ_PRIM_V1]) / (SR - SM);
        Ustar[PRJ_CONS_RHO] = rho_star;
        Ustar[PRJ_CONS_MOM1] = rho_star * SM;
        Ustar[PRJ_CONS_MOM2] = rho_star * WR[PRJ_PRIM_V2];
        Ustar[PRJ_CONS_MOM3] = rho_star * WR[PRJ_PRIM_V3];
        e_over_rho = UR[PRJ_CONS_ETOT] / UR[PRJ_CONS_RHO] +
            (SM - WR[PRJ_PRIM_V1]) * (SM + pR / (UR[PRJ_CONS_RHO] * (SR - WR[PRJ_PRIM_V1])));
        Ustar[PRJ_CONS_ETOT] = rho_star * e_over_rho;
        Ustar[PRJ_CONS_YE] = rho_star * WR[PRJ_PRIM_YE];
        for (v = 0; v < PRJ_NHYDRO; ++v) {
            flux[v] = FR[v] + SR * (Ustar[v] - UR[v]);
        }
    }
}

int prj_blocks_overlap_open(double amin, double amax, double bmin, double bmax)
{
    const double tol = 1.0e-12;
    return PRJ_MIN(amax, bmax) - PRJ_MAX(amin, bmin) > tol;
}

double prj_overlap_length(double amin, double amax, double bmin, double bmax)
{
    double overlap = PRJ_MIN(amax, bmax) - PRJ_MAX(amin, bmin);

    return overlap > 0.0 ? overlap : 0.0;
}

/* prj_riemann_set_mesh kept for API compatibility; no longer needed internally */
void prj_riemann_set_mesh(prj_mesh *mesh)
{
    (void)mesh;
}

/*
 * Helper: restrict fine-block face fluxes onto one coarse face cell and write them
 * into dst[v] for v in [0, PRJ_NVAR_CONS).  Returns 1 if weight_sum > 0, else 0.
 */
int prj_riemann_restrict_one(
    const prj_block *fine, int axis, int side,
    double cmin0, double cmax0, double cmin1, double cmax1,
    int tan0, int tan1,
    double *dst)
{
    int fi0s = (int)((cmin0 - fine->xmin[tan0]) / fine->dx[tan0]);
    int fi0e = (int)((cmax0 - fine->xmin[tan0]) / fine->dx[tan0]);
    int fi1s = (int)((cmin1 - fine->xmin[tan1]) / fine->dx[tan1]);
    int fi1e = (int)((cmax1 - fine->xmin[tan1]) / fine->dx[tan1]);
    int fine_ax_idx = side == 1 ? PRJ_BLOCK_SIZE : 0;
    int r0, r1, v;

    if (fi0s < 0) fi0s = 0;
    if (fi1s < 0) fi1s = 0;
    if (fi0e >= PRJ_BLOCK_SIZE) fi0e = PRJ_BLOCK_SIZE - 1;
    if (fi1e >= PRJ_BLOCK_SIZE) fi1e = PRJ_BLOCK_SIZE - 1;

    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        double restricted = 0.0;
        double weight_sum = 0.0;
        int fine_face[3] = {0, 0, 0};

        fine_face[axis] = fine_ax_idx;
        for (r0 = fi0s; r0 <= fi0e; ++r0) {
            for (r1 = fi1s; r1 <= fi1e; ++r1) {
                double fm0 = fine->xmin[tan0] + (double)r0 * fine->dx[tan0];
                double fm1 = fine->xmin[tan1] + (double)r1 * fine->dx[tan1];
                double ov0 = prj_overlap_length(cmin0, cmax0, fm0, fm0 + fine->dx[tan0]);
                double ov1 = prj_overlap_length(cmin1, cmax1, fm1, fm1 + fine->dx[tan1]);
                double w = ov0 * ov1;

                if (w <= 0.0) continue;
                fine_face[tan0] = r0;
                fine_face[tan1] = r1;
                restricted += fine->flux[axis][VIDX(v, fine_face[0], fine_face[1], fine_face[2])] * w;
                weight_sum += w;
            }
        }
        if (weight_sum <= 0.0) return 0;
        dst[v] = restricted / weight_sum;
    }
    return 1;
}

/*
 * Helper: find the face axis and side for block→slot, with tangential overlap check.
 * Returns axis in [0,2] (side set via *side_out), or -1 if not a valid face neighbor.
 */
int prj_riemann_face_axis(const prj_block *block, const double slot_xmin[3],
    const double slot_xmax[3], int *side_out)
{
    int axis = -1;
    int side = -1;
    int d;

    for (d = 0; d < 3; ++d) {
        if (fabs(block->xmax[d] - slot_xmin[d]) < 1.0e-12*block->dx[d]) {
            if (axis >= 0) return -1;
            axis = d; side = 1;
        } else if (fabs(slot_xmax[d] - block->xmin[d]) < 1.0e-12*block->dx[d]) {
            if (axis >= 0) return -1;
            axis = d; side = 0;
        }
    }
    if (axis < 0) return -1;
    for (d = 0; d < 3; ++d) {
        if (d == axis) continue;
        if (!prj_blocks_overlap_open(block->xmin[d], block->xmax[d],
                slot_xmin[d], slot_xmax[d]))
            return -1;
    }
    *side_out = side;
    return axis;
}

void prj_riemann_flux_send(prj_mesh *mesh, const prj_mpi *mpi)
{
    int my_rank;
    int bidx;

    if (mesh == 0) return;

    my_rank = (mpi != 0) ? mpi->rank : 0;

    /* ---- Phase 1: local (same-rank) coarse-fine correction ---- */
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        int n;

        if (block->id < 0 || block->active != 1 || block->rank != my_rank) continue;

        for (n = 0; n < 56; ++n) {
            const prj_neighbor *slot = &block->slot[n];
            if (slot->type!=PRJ_NEIGHBOR_FACE||slot->rel_level>=0||slot->rank!=my_rank) continue;
            int nid = slot->id;
            if (nid < 0 || nid >= mesh->nblocks) continue;
            prj_block *neighbor = &mesh->blocks[nid];
            
            int axis, side, tan0, tan1, it0, it1;
            axis = prj_riemann_face_axis(block, neighbor->xmin, neighbor->xmax, &side);
            if (axis < 0) continue;

            tan0 = (axis + 1) % 3;
            tan1 = (axis + 2) % 3;

            for (it0 = 0; it0 < slot->recv_loc_end[tan0]-slot->recv_loc_start[tan0]; ++it0) {
                for (it1 = 0; it1 < slot->recv_loc_end[tan1]-slot->recv_loc_start[tan1]; ++it1) {
                    int it_recv[3] = {0, 0, 0};
                    int it_send0[3] = {0, 0, 0};
                    int it_send1[3] = {0, 0, 0};
                    int it_send2[3] = {0, 0, 0};
                    int it_send3[3] = {0, 0, 0};
                    it_recv[axis] = side == 1 ? 0 : PRJ_BLOCK_SIZE;
                    it_recv[tan0] = it0+slot->recv_loc_start[tan0];
                    it_recv[tan1] = it1+slot->recv_loc_start[tan1];
                    
                    it_send0[axis] = side == 1 ? PRJ_BLOCK_SIZE : 0;
                    it_send0[tan0] = 2*it0+slot->send_loc_start[tan0];
                    it_send0[tan1] = 2*it1+slot->send_loc_start[tan1];
                    
                    it_send1[axis] = it_send0[axis];
                    it_send1[tan0] = 2*it0+slot->send_loc_start[tan0]+1;
                    it_send1[tan1] = 2*it1+slot->send_loc_start[tan1];
                    
                    it_send2[axis] = it_send0[axis];
                    it_send2[tan0] = 2*it0+slot->send_loc_start[tan0];
                    it_send2[tan1] = 2*it1+slot->send_loc_start[tan1]+1;

                    it_send3[axis] = it_send0[axis];
                    it_send3[tan0] = 2*it0+slot->send_loc_start[tan0]+1;
                    it_send3[tan1] = 2*it1+slot->send_loc_start[tan1]+1;
                    for (int v = 0; v < PRJ_NVAR_CONS; ++v) {
                        neighbor->flux[axis][VIDX(v, it_recv[0], it_recv[1], it_recv[2])] = 
                            0.25 * (block->flux[axis][VIDX(v, it_send0[0], it_send0[1], it_send0[2])]
                                   +block->flux[axis][VIDX(v, it_send1[0], it_send1[1], it_send1[2])]
                                   +block->flux[axis][VIDX(v, it_send2[0], it_send2[1], it_send2[2])]
                                   +block->flux[axis][VIDX(v, it_send3[0], it_send3[1], it_send3[2])]
                                   );
                    }
                }
            }
        }
    }
}
