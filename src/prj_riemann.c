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

/* Reciprocals of c as compile-time constants: the GR HLLD hot path divides
 * primitives by c / c^2 on every face. Strict-IEEE builds (no -ffast-math)
 * keep `x / c` as a division; multiplying by the folded reciprocal turns those
 * into multiplies. Changes results by <=0.5 ULP (bit-identity intentionally
 * not preserved here). */
#define PRJ_GR_INV_CLIGHT (1.0 / PRJ_CLIGHT)
#define PRJ_GR_INV_CLIGHT2 (1.0 / (PRJ_CLIGHT * PRJ_CLIGHT))

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
    (void)deltau;
    (void)deltav;
    (void)deltaw;

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
    SM = (R.rho * R.vx * (SR - R.vx) -
        L.rho * L.vx * (SL - L.vx) + (L.pt - R.pt)) / denom;
    prj_hlld_require_finite(SM, "contact speed");

    pt_star_l = L.pt + L.rho * (SL - L.vx) * (SM - L.vx);
    pt_star_r = R.pt + R.rho * (SR - R.vx) * (SM - R.vx);
    pt_star = 0.5 * (pt_star_l + pt_star_r);
    if (pt_star <= 0.0 || !isfinite(pt_star)) {
        /* Unphysical star pressure (strong rarefaction): fall back to the
         * positivity-robust HLL flux over SL..SR instead of aborting.  Only
         * reached with SL < 0 < SR, so SR - SL > 0 and the HLL density > 0;
         * face velocity/EMF come from the HLL average state. */
        double inv = 1.0 / (SR - SL);
        double Uhll[PRJ_NHYDRO];
        double rho;

        for (v = 0; v < PRJ_NHYDRO; ++v) {
            flux[v] = (SR * L.F[v] - SL * R.F[v] + SL * SR * (R.U[v] - L.U[v])) * inv;
            Uhll[v] = (SR * R.U[v] - SL * L.U[v] + L.F[v] - R.F[v]) * inv;
        }
        rho = Uhll[PRJ_CONS_RHO];
        prj_hlld_face_outputs(Uhll[PRJ_CONS_MOM1] / rho, Uhll[PRJ_CONS_MOM2] / rho,
            Uhll[PRJ_CONS_MOM3] / rho, Uhll[PRJ_CONS_B1], Uhll[PRJ_CONS_B2],
            Uhll[PRJ_CONS_B3], v_face, Bv1, Bv2);
        return;
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

#if PRJ_DYNAMIC_GR
#define PRJ_GR_HLLD_MAX_ITER 32
#define PRJ_GR_HLLD_EPS 1.0e-12
/* Weak-normal-field degeneracy: relativistic HLLD's rotational (Alfven) fan
 * collapses when the normal field is dynamically negligible, leaving no
 * constructible pressure root. Detect it by magnetization B_n^2/ptot rather
 * than raw |B_n|, and route such faces to HLLE. */
#define PRJ_GR_HLLD_BN2_DEGEN 1.0e-14

enum prj_gr_hlld_q {
    PRJ_GRQ_D = 0,
    PRJ_GRQ_J1 = 1,
    PRJ_GRQ_J2 = 2,
    PRJ_GRQ_J3 = 3,
    PRJ_GRQ_H = 4,
    PRJ_GRQ_YE = 5,
    PRJ_GRQ_B2 = 6,
    PRJ_GRQ_B3 = 7,
    PRJ_GRQ_N = 8
};

typedef struct prj_gr_hlld_tetrad {
    double lower[3][3];
    double upper[3][3];
    double etx;
    double eup_t[3];
    double exx;
} prj_gr_hlld_tetrad;

typedef struct prj_gr_hlld_state {
    double rho;
    double eint_m;
    double p;
    double gamma;
    double ye;
    double v[3];
    double B[3];
    double v2;
    double wlor;
    double wlor2;
    double B2;
    double Bv;
    double b2;
    double rhoh;
    double ptot;
    double D;
    double H;
    double J[3];
    double q[PRJ_GRQ_N];
    double f[PRJ_GRQ_N];
    double lam_m;
    double lam_p;
} prj_gr_hlld_state;

typedef struct prj_gr_hlld_side {
    const prj_gr_hlld_state *base;
    double lambda;
    double RD;
    double RJ[3];
    double RH;
    double RB[3];
} prj_gr_hlld_side;

typedef struct prj_gr_hlld_fan_state {
    double v[3];
    double B[3];
    double D;
    double H;
    double J[3];
    double rhohtot;
    double eta;
    double K[3];
    double lambda_a;
    double ye;
    double q[PRJ_GRQ_N];
} prj_gr_hlld_fan_state;

typedef struct prj_gr_hlld_fan {
    double lambda_l;
    double lambda_r;
    double lambda_al;
    double lambda_ar;
    double lambda_c;
    double ptot;
    prj_gr_hlld_fan_state aL;
    prj_gr_hlld_fan_state aR;
    prj_gr_hlld_fan_state cL;
    prj_gr_hlld_fan_state cR;
    double f_aL[PRJ_GRQ_N];
    double f_aR[PRJ_GRQ_N];
    double f_cL[PRJ_GRQ_N];
    double f_cR[PRJ_GRQ_N];
} prj_gr_hlld_fan;

static double prj_gr_hlld_dot3(const double a[3], const double b[3])
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static int prj_gr_hlld_finite_array(const double *a, int n)
{
    int i;

    for (i = 0; i < n; ++i) {
        if (!isfinite(a[i])) {
            return 0;
        }
    }
    return 1;
}

static void prj_gr_hlld_apply_lower(const prj_gr_hlld_tetrad *tet,
    const double vcoord[3], double vhat[3])
{
    vhat[0] = tet->lower[0][0] * vcoord[0];
    vhat[1] = tet->lower[1][0] * vcoord[0] +
        tet->lower[1][1] * vcoord[1];
    vhat[2] = tet->lower[2][0] * vcoord[0] +
        tet->lower[2][1] * vcoord[1] +
        tet->lower[2][2] * vcoord[2];
}

static void prj_gr_hlld_apply_upper(const prj_gr_hlld_tetrad *tet,
    const double vhat[3], double vcoord[3])
{
    vcoord[0] = tet->upper[0][0] * vhat[0];
    vcoord[1] = tet->upper[0][1] * vhat[0] +
        tet->upper[1][1] * vhat[1];
    vcoord[2] = tet->upper[0][2] * vhat[0] +
        tet->upper[1][2] * vhat[1] +
        tet->upper[2][2] * vhat[2];
}

static int prj_gr_hlld_build_tetrad(const double gamma[3][3],
    double alpha, const double beta[3], prj_gr_hlld_tetrad *tet)
{
    double gu[3][3];
    double det;
    double minor_yz;
    double bhat;
    double chat;
    double dhat;
    double inv_det;
#if PRJ_RIEMANN_DEBUG
    int a;
    int b;
#endif

    if (gamma == 0 || beta == 0 || tet == 0 || alpha <= 0.0 ||
        !isfinite(alpha)) {
        return 0;
    }
#if PRJ_RIEMANN_DEBUG
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            if (!isfinite(gamma[a][b])) {
                return 0;
            }
        }
    }
#endif
    det = gamma[0][0] * (gamma[1][1] * gamma[2][2] -
            gamma[1][2] * gamma[2][1]) -
        gamma[0][1] * (gamma[1][0] * gamma[2][2] -
            gamma[1][2] * gamma[2][0]) +
        gamma[0][2] * (gamma[1][0] * gamma[2][1] -
            gamma[1][1] * gamma[2][0]);
    if (!isfinite(det) || det <= 0.0) {
        return 0;
    }
    inv_det = 1.0 / det;
    gu[0][0] = (gamma[1][1] * gamma[2][2] -
        gamma[1][2] * gamma[2][1]) * inv_det;
    gu[0][1] = (gamma[0][2] * gamma[2][1] -
        gamma[0][1] * gamma[2][2]) * inv_det;
    gu[0][2] = (gamma[0][1] * gamma[1][2] -
        gamma[0][2] * gamma[1][1]) * inv_det;
    gu[1][0] = (gamma[1][2] * gamma[2][0] -
        gamma[1][0] * gamma[2][2]) * inv_det;
    gu[1][1] = (gamma[0][0] * gamma[2][2] -
        gamma[0][2] * gamma[2][0]) * inv_det;
    gu[1][2] = (gamma[0][2] * gamma[1][0] -
        gamma[0][0] * gamma[1][2]) * inv_det;
    gu[2][0] = (gamma[1][0] * gamma[2][1] -
        gamma[1][1] * gamma[2][0]) * inv_det;
    gu[2][1] = (gamma[0][1] * gamma[2][0] -
        gamma[0][0] * gamma[2][1]) * inv_det;
    gu[2][2] = (gamma[0][0] * gamma[1][1] -
        gamma[0][1] * gamma[1][0]) * inv_det;

    minor_yz = gamma[1][1] * gamma[2][2] - gamma[1][2] * gamma[1][2];
    if (gu[0][0] <= 0.0 || gamma[2][2] <= 0.0 || minor_yz <= 0.0) {
        return 0;
    }
    bhat = 1.0 / sqrt(gu[0][0]);
    chat = 1.0 / sqrt(gamma[2][2]);
    dhat = 1.0 / sqrt(gamma[2][2] * minor_yz);
    memset(tet, 0, sizeof(*tet));

    tet->upper[0][0] = bhat * gu[0][0];
    tet->upper[0][1] = bhat * gu[0][1];
    tet->upper[0][2] = bhat * gu[0][2];
    tet->upper[1][1] = dhat * gamma[2][2];
    tet->upper[1][2] = -dhat * gamma[1][2];
    tet->upper[2][2] = chat;

    tet->lower[0][0] = bhat;
    tet->lower[1][0] = dhat *
        (gamma[0][1] * gamma[2][2] - gamma[0][2] * gamma[1][2]);
    tet->lower[1][1] = dhat * minor_yz;
    tet->lower[2][0] = chat * gamma[0][2];
    tet->lower[2][1] = chat * gamma[1][2];
    tet->lower[2][2] = chat * gamma[2][2];

    tet->eup_t[0] = -beta[0] / alpha;
    tet->eup_t[1] = -beta[1] / alpha;
    tet->eup_t[2] = -beta[2] / alpha;
    tet->etx = tet->eup_t[0];
    tet->exx = tet->upper[0][0];
#if PRJ_RIEMANN_DEBUG
    return prj_gr_hlld_finite_array(&tet->lower[0][0], 9) &&
        prj_gr_hlld_finite_array(&tet->upper[0][0], 9) &&
        prj_gr_hlld_finite_array(tet->eup_t, 3) &&
        isfinite(tet->etx) && isfinite(tet->exx);
#else
    return 1;
#endif
}

static void prj_gr_hlld_fill_qf(prj_gr_hlld_state *s)
{
    /* 1/W^2 == 1 - v^2 exactly (W = 1/sqrt(1-v^2)); use the identity to avoid a
     * division. */
    double inv_w2 = 1.0 - s->v2;
    double vB = s->Bv;

    s->q[PRJ_GRQ_D] = s->D;
    s->q[PRJ_GRQ_J1] = s->J[0];
    s->q[PRJ_GRQ_J2] = s->J[1];
    s->q[PRJ_GRQ_J3] = s->J[2];
    s->q[PRJ_GRQ_H] = s->H;
    s->q[PRJ_GRQ_YE] = s->D * s->ye;
    s->q[PRJ_GRQ_B2] = s->B[1];
    s->q[PRJ_GRQ_B3] = s->B[2];
    s->f[PRJ_GRQ_D] = s->D * s->v[0];
    s->f[PRJ_GRQ_J1] = s->J[0] * s->v[0] + s->ptot -
        s->B[0] * (s->B[0] * inv_w2 + vB * s->v[0]);
    s->f[PRJ_GRQ_J2] = s->J[1] * s->v[0] -
        s->B[0] * (s->B[1] * inv_w2 + vB * s->v[1]);
    s->f[PRJ_GRQ_J3] = s->J[2] * s->v[0] -
        s->B[0] * (s->B[2] * inv_w2 + vB * s->v[2]);
    s->f[PRJ_GRQ_H] = (s->H + s->ptot) * s->v[0] - vB * s->B[0];
    s->f[PRJ_GRQ_YE] = s->f[PRJ_GRQ_D] * s->ye;
    s->f[PRJ_GRQ_B2] = s->v[0] * s->B[1] - s->v[1] * s->B[0];
    s->f[PRJ_GRQ_B3] = s->v[0] * s->B[2] - s->v[2] * s->B[0];
}

static int prj_gr_hlld_state_from_prim(const double *W, double pressure,
    double gas_gamma, const prj_gr_hlld_tetrad *tet, double bn_phys,
    prj_gr_hlld_state *s)
{
    double beta_coord[3];
    double B_coord_m[3];
    double sqrt_one_minus_v2;
    double cs2;
    double va2;
    double zeta;
    double disc;
    double den;
    double inv_den;
    int d;

    if (W == 0 || tet == 0 || s == 0 || pressure <= 0.0 ||
        !isfinite(pressure)) {
        return 0;
    }
    s->rho = W[PRJ_PRIM_RHO];
    s->eint_m = W[PRJ_PRIM_EINT] * PRJ_GR_INV_CLIGHT2;
    s->p = pressure * PRJ_GR_INV_CLIGHT2;
    s->gamma = isfinite(gas_gamma) && gas_gamma > 0.0 ? gas_gamma : 5.0 / 3.0;
    s->ye = W[PRJ_PRIM_YE];
    if (s->rho <= 0.0 || s->eint_m < 0.0 || s->p <= 0.0 ||
        !isfinite(s->rho) || !isfinite(s->eint_m) || !isfinite(s->ye)) {
        return 0;
    }
    beta_coord[0] = W[PRJ_PRIM_V1] * PRJ_GR_INV_CLIGHT;
    beta_coord[1] = W[PRJ_PRIM_V2] * PRJ_GR_INV_CLIGHT;
    beta_coord[2] = W[PRJ_PRIM_V3] * PRJ_GR_INV_CLIGHT;
    B_coord_m[0] = bn_phys * PRJ_GR_INV_CLIGHT;
    B_coord_m[1] = W[PRJ_PRIM_B2] * PRJ_GR_INV_CLIGHT;
    B_coord_m[2] = W[PRJ_PRIM_B3] * PRJ_GR_INV_CLIGHT;
    prj_gr_hlld_apply_lower(tet, beta_coord, s->v);
    prj_gr_hlld_apply_lower(tet, B_coord_m, s->B);
    s->v2 = prj_gr_hlld_dot3(s->v, s->v);
    if (!isfinite(s->v2) || s->v2 < 0.0 || s->v2 >= 1.0) {
        return 0;
    }
    sqrt_one_minus_v2 = sqrt(1.0 - s->v2);
    s->wlor = 1.0 / sqrt_one_minus_v2;
    s->wlor2 = s->wlor * s->wlor;
    s->B2 = prj_gr_hlld_dot3(s->B, s->B);
    s->Bv = prj_gr_hlld_dot3(s->B, s->v);
    s->b2 = s->B2 * (1.0 - s->v2) + s->Bv * s->Bv;
    s->rhoh = s->rho * (1.0 + s->eint_m) + s->p;
    s->ptot = s->p + 0.5 * s->b2;
    s->D = s->rho * s->wlor;
    s->H = s->rhoh * s->wlor2 - s->p + s->B2 - 0.5 * s->b2;
    for (d = 0; d < 3; ++d) {
        s->J[d] = (s->rhoh + s->B2) * s->wlor2 * s->v[d] -
            s->Bv * s->B[d];
    }
    if (s->rhoh <= 0.0 || s->ptot <= 0.0 || s->D <= 0.0 ||
        !isfinite(s->b2) || !isfinite(s->H)) {
        return 0;
    }
    prj_gr_hlld_fill_qf(s);

    cs2 = s->gamma * s->p / s->rhoh;
    if (!isfinite(cs2) || cs2 < 0.0) {
        cs2 = 0.0;
    }
    if (cs2 > 1.0) {
        cs2 = 1.0;
    }
    va2 = s->b2 / (s->rhoh + s->b2);
    if (!isfinite(va2) || va2 < 0.0) {
        va2 = 0.0;
    }
    if (va2 > 1.0) {
        va2 = 1.0;
    }
    zeta = va2 + cs2 - va2 * cs2;
    if (!isfinite(zeta) || zeta < 0.0) {
        zeta = 0.0;
    }
    if (zeta > 1.0) {
        zeta = 1.0;
    }
    den = 1.0 - s->v2 * zeta;
    disc = (1.0 - s->v2) *
        (1.0 - s->v2 * zeta - (1.0 - zeta) * s->v[0] * s->v[0]);
    if (disc < 0.0 && disc > -1.0e-13) {
        disc = 0.0;
    }
    if (den <= 0.0 || disc < 0.0 || !isfinite(disc)) {
        return 0;
    }
    disc = sqrt(zeta) * sqrt(disc);
    inv_den = 1.0 / den;
    s->lam_m = (s->v[0] * (1.0 - zeta) - disc) * inv_den;
    s->lam_p = (s->v[0] * (1.0 - zeta) + disc) * inv_den;
    if (!isfinite(s->lam_m) || !isfinite(s->lam_p)) {
        return 0;
    }
    if (s->lam_m < -1.0) s->lam_m = -1.0;
    if (s->lam_p > 1.0) s->lam_p = 1.0;
    return 1;
}

static void prj_gr_hlld_side_invariants(const prj_gr_hlld_state *s,
    double lambda, prj_gr_hlld_side *side)
{
    int d;

    side->base = s;
    side->lambda = lambda;
    side->RD = lambda * s->q[PRJ_GRQ_D] - s->f[PRJ_GRQ_D];
    side->RH = lambda * s->q[PRJ_GRQ_H] - s->f[PRJ_GRQ_H];
    for (d = 0; d < 3; ++d) {
        side->RJ[d] = lambda * s->q[PRJ_GRQ_J1 + d] -
            s->f[PRJ_GRQ_J1 + d];
    }
    side->RB[0] = lambda * s->B[0];
    side->RB[1] = lambda * s->q[PRJ_GRQ_B2] - s->f[PRJ_GRQ_B2];
    side->RB[2] = lambda * s->q[PRJ_GRQ_B3] - s->f[PRJ_GRQ_B3];
}

static void prj_gr_hlld_fill_fan_q(const prj_gr_hlld_fan_state *s,
    double q[PRJ_GRQ_N])
{
    q[PRJ_GRQ_D] = s->D;
    q[PRJ_GRQ_J1] = s->J[0];
    q[PRJ_GRQ_J2] = s->J[1];
    q[PRJ_GRQ_J3] = s->J[2];
    q[PRJ_GRQ_H] = s->H;
    q[PRJ_GRQ_YE] = s->D * s->ye;
    q[PRJ_GRQ_B2] = s->B[1];
    q[PRJ_GRQ_B3] = s->B[2];
}

static int prj_gr_hlld_build_outer_state(const prj_gr_hlld_side *side,
    double Bx, double ptot, prj_gr_hlld_fan_state *a)
{
    double lambda = side->lambda;
    double A;
    double G;
    double C;
    double Q;
    double X;
    double den;
    double vB;
    int d;

    a->B[0] = Bx;
    a->ye = side->base->ye;
    A = side->RJ[0] - lambda * side->RH + ptot * (1.0 - lambda * lambda);
    G = side->RB[1] * side->RB[1] + side->RB[2] * side->RB[2];
    C = side->RJ[1] * side->RB[1] + side->RJ[2] * side->RB[2];
    Q = -A - G + Bx * Bx * (1.0 - lambda * lambda);
    X = Bx * (A * lambda * Bx + C) -
        (A + G) * (lambda * ptot + side->RH);
    if (fabs(X) <= PRJ_GR_HLLD_EPS *
        (fabs(Bx * (A * lambda * Bx + C)) +
            fabs((A + G) * (lambda * ptot + side->RH)) + 1.0)) {
        return 0;
    }
    a->v[0] = (Bx * (A * Bx + lambda * C) -
        (A + G) * (ptot + side->RJ[0])) / X;
    a->v[1] = (Q * side->RJ[1] +
        side->RB[1] * (C + Bx * (lambda * side->RJ[0] - side->RH))) / X;
    a->v[2] = (Q * side->RJ[2] +
        side->RB[2] * (C + Bx * (lambda * side->RJ[0] - side->RH))) / X;
    if (!prj_gr_hlld_finite_array(a->v, 3) ||
        prj_gr_hlld_dot3(a->v, a->v) >= 1.0) {
        return 0;
    }
    den = lambda - a->v[0];
    if (fabs(den) <= PRJ_GR_HLLD_EPS * (fabs(lambda) + fabs(a->v[0]) + 1.0)) {
        return 0;
    }
    a->B[1] = (side->RB[1] - Bx * a->v[1]) / den;
    a->B[2] = (side->RB[2] - Bx * a->v[2]) / den;
    vB = prj_gr_hlld_dot3(a->v, a->B);
    a->rhohtot = ptot + (side->RH -
        prj_gr_hlld_dot3(a->v, side->RJ)) / den;
    a->D = side->RD / den;
    a->H = (side->RH + ptot * a->v[0] - vB * Bx) / den;
    for (d = 0; d < 3; ++d) {
        a->J[d] = (a->H + ptot) * a->v[d] - vB * a->B[d];
    }
    if (a->rhohtot <= 0.0 || a->D <= 0.0 || !isfinite(a->H) ||
        !prj_gr_hlld_finite_array(a->B, 3) ||
        !prj_gr_hlld_finite_array(a->J, 3)) {
        return 0;
    }
    prj_gr_hlld_fill_fan_q(a, a->q);
    return prj_gr_hlld_finite_array(a->q, PRJ_GRQ_N);
}

static int prj_gr_hlld_contact_one(const prj_gr_hlld_side *side,
    const prj_gr_hlld_fan_state *a, const double Bc[3],
    double ptot, int right_side, prj_gr_hlld_fan_state *c)
{
    double sign_bx = side->base->B[0] < 0.0 ? -1.0 : 1.0;
    double den;
    double K2;
    double KB;
    double fac;
    double vB;
    double wave_den;
    int d;

    c->eta = (right_side ? 1.0 : -1.0) * sign_bx * sqrt(a->rhohtot);
    den = side->lambda * ptot + side->RH + side->base->B[0] * c->eta;
    if (fabs(den) <= PRJ_GR_HLLD_EPS *
        (fabs(side->lambda * ptot) + fabs(side->RH) +
            fabs(side->base->B[0] * c->eta) + 1.0)) {
        return 0;
    }
    for (d = 0; d < 3; ++d) {
        c->K[d] = (side->RJ[d] + (d == 0 ? ptot : 0.0) +
            side->RB[d] * c->eta) / den;
        c->B[d] = Bc[d];
    }
    c->lambda_a = c->K[0];
    K2 = prj_gr_hlld_dot3(c->K, c->K);
    KB = prj_gr_hlld_dot3(c->K, c->B);
    den = c->eta - KB;
    if (fabs(den) <= PRJ_GR_HLLD_EPS * (fabs(c->eta) + fabs(KB) + 1.0)) {
        return 0;
    }
    fac = (1.0 - K2) / den;
    for (d = 0; d < 3; ++d) {
        c->v[d] = c->K[d] - c->B[d] * fac;
    }
    if (!prj_gr_hlld_finite_array(c->v, 3) ||
        prj_gr_hlld_dot3(c->v, c->v) >= 1.0) {
        return 0;
    }
    c->D = a->D * (c->lambda_a - a->v[0]) /
        (c->lambda_a - c->v[0]);
    vB = prj_gr_hlld_dot3(c->v, c->B);
    wave_den = c->lambda_a - c->v[0];
    if (fabs(wave_den) <= PRJ_GR_HLLD_EPS *
        (fabs(c->lambda_a) + fabs(c->v[0]) + 1.0)) {
        return 0;
    }
    c->H = (c->lambda_a * a->H - a->J[0] +
        ptot * c->v[0] - vB * c->B[0]) / wave_den;
    c->rhohtot = a->rhohtot;
    c->ye = a->ye;
    for (d = 0; d < 3; ++d) {
        c->J[d] = (c->H + ptot) * c->v[d] - vB * c->B[d];
    }
    if (c->D <= 0.0 || c->H <= 0.0 ||
        !prj_gr_hlld_finite_array(c->J, 3)) {
        return 0;
    }
    prj_gr_hlld_fill_fan_q(c, c->q);
    return prj_gr_hlld_finite_array(c->q, PRJ_GRQ_N);
}

static int prj_gr_hlld_build_contact_states(const prj_gr_hlld_side *L,
    const prj_gr_hlld_side *R, const prj_gr_hlld_fan_state *aL,
    const prj_gr_hlld_fan_state *aR, double ptot,
    prj_gr_hlld_fan_state *cL, prj_gr_hlld_fan_state *cR,
    double *residual)
{
    double Bc[3];
    double lambda_al;
    double lambda_ar;
    double den;
    int k;

    if (!prj_gr_hlld_contact_one(L, aL, aL->B, ptot, 0, cL) ||
        !prj_gr_hlld_contact_one(R, aR, aR->B, ptot, 1, cR)) {
        return 0;
    }
    lambda_al = cL->lambda_a;
    lambda_ar = cR->lambda_a;
    den = lambda_ar - lambda_al;
    if (fabs(den) <= PRJ_GR_HLLD_EPS *
        (fabs(lambda_ar) + fabs(lambda_al) + 1.0)) {
        return 0;
    }
    Bc[0] = L->base->B[0];
    for (k = 1; k < 3; ++k) {
        Bc[k] = (aR->B[k] * (R->lambda - aR->v[0]) +
            Bc[0] * aR->v[k] -
            aL->B[k] * (L->lambda - aL->v[0]) -
            Bc[0] * aL->v[k]) / den;
    }
    if (!prj_gr_hlld_contact_one(L, aL, Bc, ptot, 0, cL) ||
        !prj_gr_hlld_contact_one(R, aR, Bc, ptot, 1, cR)) {
        return 0;
    }
    if (residual != 0) {
        *residual = cL->v[0] - cR->v[0];
    }
    return 1;
}

static void prj_gr_hlld_flux_from_jump(const double f0[PRJ_GRQ_N],
    double lambda, const double q1[PRJ_GRQ_N], const double q0[PRJ_GRQ_N],
    double f[PRJ_GRQ_N])
{
    int q;

    for (q = 0; q < PRJ_GRQ_N; ++q) {
        f[q] = f0[q] + lambda * (q1[q] - q0[q]);
    }
}

static int prj_gr_hlld_build_fan(const prj_gr_hlld_side *L,
    const prj_gr_hlld_side *R,
    double ptot, prj_gr_hlld_fan *fan, double *residual)
{
    int q;

    if (ptot <= 0.0 || !isfinite(ptot)) {
        return 0;
    }
    if (!prj_gr_hlld_build_outer_state(L, L->base->B[0], ptot, &fan->aL) ||
        !prj_gr_hlld_build_outer_state(R, L->base->B[0], ptot, &fan->aR) ||
        !prj_gr_hlld_build_contact_states(L, R, &fan->aL, &fan->aR,
            ptot, &fan->cL, &fan->cR, residual)) {
        return 0;
    }
    fan->lambda_l = L->lambda;
    fan->lambda_r = R->lambda;
    fan->lambda_al = fan->cL.lambda_a;
    fan->lambda_ar = fan->cR.lambda_a;
    fan->lambda_c = 0.5 * (fan->cL.v[0] + fan->cR.v[0]);
    fan->ptot = ptot;
    if (!(fan->lambda_l < fan->lambda_r) ||
        !(fan->lambda_l <= fan->lambda_al) ||
        !(fan->lambda_al <= fan->lambda_ar) ||
        !(fan->lambda_ar <= fan->lambda_r)) {
        return 0;
    }
    prj_gr_hlld_flux_from_jump(L->base->f, L->lambda, fan->aL.q, L->base->q,
        fan->f_aL);
    prj_gr_hlld_flux_from_jump(R->base->f, R->lambda, fan->aR.q, R->base->q,
        fan->f_aR);
    for (q = 0; q < PRJ_GRQ_N; ++q) {
        fan->f_cL[q] = fan->f_aL[q] +
            fan->lambda_al * (fan->cL.q[q] - fan->aL.q[q]);
        fan->f_cR[q] = fan->f_aR[q] +
            fan->lambda_ar * (fan->cR.q[q] - fan->aR.q[q]);
    }
    return prj_gr_hlld_finite_array(fan->f_aL, PRJ_GRQ_N) &&
        prj_gr_hlld_finite_array(fan->f_aR, PRJ_GRQ_N) &&
        prj_gr_hlld_finite_array(fan->f_cL, PRJ_GRQ_N) &&
        prj_gr_hlld_finite_array(fan->f_cR, PRJ_GRQ_N);
}

static double prj_gr_hlld_pressure_guess(const prj_gr_hlld_state *L,
    const prj_gr_hlld_state *R)
{
    double guess = 0.5 * (L->ptot + R->ptot);

    if (!(guess > 0.0) || !isfinite(guess)) {
        guess = PRJ_MAX(L->ptot, R->ptot);
    }
    return guess;
}

static int prj_gr_hlld_same_primitive_state(const double *WL, const double *WR,
    double pL, double pR, double gL, double gR)
{
    const int vars[] = {
        PRJ_PRIM_RHO,
        PRJ_PRIM_V1,
        PRJ_PRIM_V2,
        PRJ_PRIM_V3,
        PRJ_PRIM_EINT,
        PRJ_PRIM_YE,
        PRJ_PRIM_B2,
        PRJ_PRIM_B3
    };
    int n;

    if (pL != pR || gL != gR) {
        return 0;
    }
    for (n = 0; n < (int)(sizeof(vars) / sizeof(vars[0])); ++n) {
        int v = vars[n];
        if (WL[v] != WR[v]) {
            return 0;
        }
    }
    return 1;
}

static int prj_gr_hlld_fan_converged(double residual, const prj_gr_hlld_fan *fan)
{
    return fabs(residual) <= 1.0e-10 *
        (fabs(fan->lambda_ar - fan->lambda_al) + 1.0);
}

static int prj_gr_hlld_solve_pressure(const prj_gr_hlld_state *L,
    const prj_gr_hlld_state *R, double lambda_l, double lambda_r,
    prj_gr_hlld_fan *fan)
{
    const double scale = PRJ_MAX(PRJ_MAX(L->ptot, R->ptot), 1.0);
    prj_gr_hlld_side side_l;
    prj_gr_hlld_side side_r;
    double guesses[5];
    int ig;

    prj_gr_hlld_side_invariants(L, lambda_l, &side_l);
    prj_gr_hlld_side_invariants(R, lambda_r, &side_r);
    guesses[0] = prj_gr_hlld_pressure_guess(L, R);
    guesses[1] = L->ptot;
    guesses[2] = R->ptot;
    guesses[3] = 2.0 * guesses[0];
    guesses[4] = 0.5 * guesses[0];
    for (ig = 0; ig < 5; ++ig) {
        double p0 = guesses[ig];
        double p1;
        double f0;
        double f1;
        double dp;
        prj_gr_hlld_fan trial;
        int iter;

        if (!(p0 > 0.0) || !isfinite(p0)) {
            continue;
        }
        /* Bootstrap the secant with two points built entirely from this face's
         * own L/R states -- p0 is a face-local seed and p1 is a fixed relative
         * offset from it. No dependence on any other face's result, so the flux
         * loop stays free of loop-carried state. */
        if (!prj_gr_hlld_build_fan(&side_l, &side_r, p0, &trial, &f0)) {
            continue;
        }
        if (prj_gr_hlld_fan_converged(f0, &trial)) {
            *fan = trial;
            return 1;
        }
        dp = 1.0e-6 * PRJ_MAX(fabs(p0), scale);
        p1 = p0 + dp;
        if (!prj_gr_hlld_build_fan(&side_l, &side_r, p1, &trial, &f1)) {
            p1 = p0 - dp;
            if (!(p1 > 0.0) ||
                !prj_gr_hlld_build_fan(&side_l, &side_r, p1, &trial, &f1)) {
                continue;
            }
        }
        if (prj_gr_hlld_fan_converged(f1, &trial)) {
            *fan = trial;
            return 1;
        }
        /* Secant iteration: one new residual evaluation per step, reusing the
         * (p0,f0),(p1,f1) pair -- half the build_fan calls of the previous
         * finite-difference Newton, which spent a second evaluation per step
         * purely to estimate df/dp. */
        for (iter = 0; iter < PRJ_GR_HLLD_MAX_ITER; ++iter) {
            double denom = f1 - f0;
            double pnew;
            double fnew;
            double step;
            int ok;
            int bt;

            if (fabs(denom) <= 1.0e-300 || !isfinite(denom)) {
                break;
            }
            pnew = p1 - f1 * (p1 - p0) / denom;
            if (!(pnew > 0.0) || !isfinite(pnew)) {
                pnew = 0.5 * p1;
            }
            if (pnew > 10.0 * p1) {
                pnew = 10.0 * p1;
            }
            if (pnew < 0.1 * p1) {
                pnew = 0.1 * p1;
            }
            ok = prj_gr_hlld_build_fan(&side_l, &side_r, pnew, &trial, &fnew);
            /* If the trial pressure yields an invalid fan, backtrack toward the
             * last good point p1 (bounded) rather than abandoning the seed. */
            step = pnew - p1;
            for (bt = 0; !ok && bt < 8; ++bt) {
                step *= 0.5;
                pnew = p1 + step;
                if (!(pnew > 0.0)) {
                    break;
                }
                ok = prj_gr_hlld_build_fan(&side_l, &side_r, pnew, &trial, &fnew);
            }
            if (!ok) {
                break;
            }
            p0 = p1;
            f0 = f1;
            p1 = pnew;
            f1 = fnew;
            if (prj_gr_hlld_fan_converged(f1, &trial)) {
                *fan = trial;
                return 1;
            }
        }
    }
    return 0;
}

static void prj_gr_hlld_emit_flux(const prj_gr_hlld_tetrad *tet,
    double sqrt_gamma, double alpha, const double beta[3],
    const double q[PRJ_GRQ_N], const double f[PRJ_GRQ_N],
    const double Bhat[3], const double vhat[3],
    double *flux, double v_face[3], double *Bv1, double *Bv2)
{
    double c = PRJ_CLIGHT;
    double c2 = c * c;
    double c3 = c2 * c;
    double flux_b2;
    double flux_b3;
    double mom0;
    double mom1;
    double mom2;

    flux[PRJ_CONS_RHO] = sqrt_gamma * c * alpha *
        (tet->etx * q[PRJ_GRQ_D] + tet->exx * f[PRJ_GRQ_D]);
    flux[PRJ_CONS_YE] = sqrt_gamma * c * alpha *
        (tet->etx * q[PRJ_GRQ_YE] + tet->exx * f[PRJ_GRQ_YE]);
    mom0 = tet->etx * q[PRJ_GRQ_J1] + tet->exx * f[PRJ_GRQ_J1];
    mom1 = tet->etx * q[PRJ_GRQ_J2] + tet->exx * f[PRJ_GRQ_J2];
    mom2 = tet->etx * q[PRJ_GRQ_J3] + tet->exx * f[PRJ_GRQ_J3];
    flux[PRJ_CONS_MOM1] = sqrt_gamma * c2 * alpha *
        (tet->upper[0][0] * mom0);
    flux[PRJ_CONS_MOM2] = sqrt_gamma * c2 * alpha *
        (tet->upper[0][1] * mom0 + tet->upper[1][1] * mom1);
    flux[PRJ_CONS_MOM3] = sqrt_gamma * c2 * alpha *
        (tet->upper[0][2] * mom0 + tet->upper[1][2] * mom1 +
            tet->upper[2][2] * mom2);
    flux[PRJ_CONS_ETOT] = sqrt_gamma * c3 * alpha *
        (tet->etx * (q[PRJ_GRQ_H] - q[PRJ_GRQ_D]) +
            tet->exx * (f[PRJ_GRQ_H] - f[PRJ_GRQ_D]));

    flux[PRJ_CONS_B1] = 0.0;
    flux_b2 = tet->upper[0][1] * tet->etx * Bhat[0] +
        tet->upper[1][1] * tet->etx * q[PRJ_GRQ_B2] -
        tet->eup_t[1] * tet->exx * Bhat[0] +
        tet->upper[1][1] * tet->exx * f[PRJ_GRQ_B2];
    flux_b3 = tet->upper[0][2] * tet->etx * Bhat[0] +
        tet->upper[1][2] * tet->etx * q[PRJ_GRQ_B2] +
        tet->upper[2][2] * tet->etx * q[PRJ_GRQ_B3] -
        tet->eup_t[2] * tet->exx * Bhat[0] +
        tet->upper[1][2] * tet->exx * f[PRJ_GRQ_B2] +
        tet->upper[2][2] * tet->exx * f[PRJ_GRQ_B3];
    flux[PRJ_CONS_B2] = sqrt_gamma * c2 * alpha * flux_b2;
    flux[PRJ_CONS_B3] = sqrt_gamma * c2 * alpha * flux_b3;
    if (v_face != 0) {
        double beta_coord[3] = {0.0, 0.0, 0.0};
        int d;

        if (vhat != 0) {
            prj_gr_hlld_apply_upper(tet, vhat, beta_coord);
        } else if (fabs(q[PRJ_GRQ_D]) > 0.0) {
            double approx[3] = {0.0, 0.0, 0.0};

            approx[0] = f[PRJ_GRQ_D] / q[PRJ_GRQ_D];
            if (approx[0] > 1.0) approx[0] = 1.0;
            if (approx[0] < -1.0) approx[0] = -1.0;
            prj_gr_hlld_apply_upper(tet, approx, beta_coord);
        }
        for (d = 0; d < 3; ++d) {
            v_face[d] = c * (alpha * beta_coord[d] - beta[d]);
        }
    }
    if (Bv1 != 0) {
        *Bv1 = flux[PRJ_CONS_B3];
    }
    if (Bv2 != 0) {
        *Bv2 = -flux[PRJ_CONS_B2];
    }
}

static void prj_gr_hlld_hlle_fallback(const prj_gr_hlld_state *L,
    const prj_gr_hlld_state *R, const prj_gr_hlld_tetrad *tet,
    double sqrt_gamma, double alpha, const double beta[3],
    double lambda_l, double lambda_r, double s_interface,
    double *flux, double v_face[3], double *Bv1, double *Bv2)
{
    double qhll[PRJ_GRQ_N];
    double fhll[PRJ_GRQ_N];
    double Bhll[3];
    double inv;
    int q;

    if (s_interface <= lambda_l) {
        prj_gr_hlld_emit_flux(tet, sqrt_gamma, alpha, beta, L->q, L->f,
            L->B, L->v, flux, v_face, Bv1, Bv2);
        return;
    }
    if (s_interface >= lambda_r) {
        prj_gr_hlld_emit_flux(tet, sqrt_gamma, alpha, beta, R->q, R->f,
            R->B, R->v, flux, v_face, Bv1, Bv2);
        return;
    }
    inv = 1.0 / (lambda_r - lambda_l);
    for (q = 0; q < PRJ_GRQ_N; ++q) {
        qhll[q] = (lambda_r * R->q[q] - lambda_l * L->q[q] +
            L->f[q] - R->f[q]) * inv;
        fhll[q] = (lambda_r * L->f[q] - lambda_l * R->f[q] +
            lambda_l * lambda_r * (R->q[q] - L->q[q])) * inv;
    }
    Bhll[0] = L->B[0];
    Bhll[1] = qhll[PRJ_GRQ_B2];
    Bhll[2] = qhll[PRJ_GRQ_B3];
    prj_gr_hlld_emit_flux(tet, sqrt_gamma, alpha, beta, qhll, fhll,
        Bhll, 0, flux, v_face, Bv1, Bv2);
}

static void prj_riemann_gr_hlld_impl(const double *WL, const double *WR,
    double pL, double pR, double gL, double gR,
    prj_eos *eos, const double gamma[3][3], double sqrt_gamma,
    double alpha, const double beta[3], double bn_tilde, double *flux,
    double v_face[3], double *Bv1, double *Bv2)
{
    prj_gr_hlld_tetrad tet;
    prj_gr_hlld_state L;
    prj_gr_hlld_state R;
    prj_gr_hlld_fan fan;
    const double *qsel = 0;
    const double *fsel = 0;
    const double *Bsel = 0;
    const double *vsel = 0;
    double bn_phys;
    double lambda_l;
    double lambda_r;
    double s_interface;

    (void)eos;
    if (WL == 0 || WR == 0 || gamma == 0 || beta == 0 || flux == 0 ||
        sqrt_gamma <= 0.0 || alpha <= 0.0 || !isfinite(sqrt_gamma) ||
        !isfinite(alpha)) {
        prj_riemann_hlld_fail("invalid GR HLLD input");
    }
    if (!prj_gr_hlld_build_tetrad(gamma, alpha, beta, &tet)) {
        prj_riemann_hlld_fail("invalid GR HLLD tetrad");
    }
    bn_phys = bn_tilde / sqrt_gamma;
    if (!prj_gr_hlld_state_from_prim(WL, pL, gL, &tet, bn_phys, &L) ||
        !prj_gr_hlld_state_from_prim(WR, pR, gR, &tet, bn_phys, &R)) {
        prj_riemann_hlld_fail("invalid GR HLLD primitive state");
    }
    if (prj_gr_hlld_same_primitive_state(WL, WR, pL, pR, gL, gR)) {
        prj_gr_hlld_emit_flux(&tet, sqrt_gamma, alpha, beta, L.q, L.f,
            L.B, L.v, flux, v_face, Bv1, Bv2);
        return;
    }
    lambda_l = PRJ_MIN(L.lam_m, R.lam_m);
    lambda_r = PRJ_MAX(L.lam_p, R.lam_p);
    if (!(lambda_l < lambda_r) || !isfinite(lambda_l) || !isfinite(lambda_r)) {
        prj_riemann_hlld_fail("invalid GR HLLD wave ordering");
    }
    /* Interface speed in the tetrad (orthonormal normal) frame: the shift
     * projects onto the normal one-form with norm sqrt(gamma^xx), which the
     * tetrad already carries as tet.exx = sqrt(gu[0][0]). 1/sqrt(gamma[0][0])
     * only equals sqrt(gamma^xx) for a diagonal metric. */
    s_interface = beta[0] * tet.exx / alpha;
    if (!isfinite(s_interface)) {
        s_interface = 0.0;
    }
    if (s_interface <= lambda_l) {
        prj_gr_hlld_emit_flux(&tet, sqrt_gamma, alpha, beta, L.q, L.f,
            L.B, L.v, flux, v_face, Bv1, Bv2);
        return;
    }
    if (s_interface >= lambda_r) {
        prj_gr_hlld_emit_flux(&tet, sqrt_gamma, alpha, beta, R.q, R.f,
            R.B, R.v, flux, v_face, Bv1, Bv2);
        return;
    }
    if (L.B[0] * L.B[0] < PRJ_GR_HLLD_BN2_DEGEN * PRJ_MAX(L.ptot, R.ptot)) {
        prj_gr_hlld_hlle_fallback(&L, &R, &tet, sqrt_gamma, alpha, beta,
            lambda_l, lambda_r, s_interface, flux, v_face, Bv1, Bv2);
        return;
    }
    if (!prj_gr_hlld_solve_pressure(&L, &R, lambda_l, lambda_r, &fan)) {
        /* The B_normal~=0 degeneracy above is handled by HLLE (relativistic
         * HLLD is undefined there). Reaching here means HLLD *should* apply but
         * the pressure root find did not converge -- surface it rather than
         * silently degrading to HLLE. */
        prj_riemann_hlld_fail("GR HLLD pressure solve failed to converge");
    }
    if (s_interface <= fan.lambda_l) {
        qsel = L.q;
        fsel = L.f;
        Bsel = L.B;
        vsel = L.v;
    } else if (s_interface <= fan.lambda_al) {
        qsel = fan.aL.q;
        fsel = fan.f_aL;
        Bsel = fan.aL.B;
        vsel = fan.aL.v;
    } else if (s_interface <= fan.lambda_c) {
        qsel = fan.cL.q;
        fsel = fan.f_cL;
        Bsel = fan.cL.B;
        vsel = fan.cL.v;
    } else if (s_interface <= fan.lambda_ar) {
        qsel = fan.cR.q;
        fsel = fan.f_cR;
        Bsel = fan.cR.B;
        vsel = fan.cR.v;
    } else if (s_interface <= fan.lambda_r) {
        qsel = fan.aR.q;
        fsel = fan.f_aR;
        Bsel = fan.aR.B;
        vsel = fan.aR.v;
    } else {
        qsel = R.q;
        fsel = R.f;
        Bsel = R.B;
        vsel = R.v;
    }
    prj_gr_hlld_emit_flux(&tet, sqrt_gamma, alpha, beta, qsel, fsel,
        Bsel, vsel, flux, v_face, Bv1, Bv2);
}

void prj_riemann_gr_hlld(const double *WL, const double *WR,
    double pL, double pR, double gL, double gR,
    prj_eos *eos, const double gamma[3][3], double sqrt_gamma,
    double alpha, const double beta[3], double bn_tilde, double *flux,
    double v_face[3], double *Bv1, double *Bv2,
    double deltau, double deltav, double deltaw)
{
    (void)deltau;
    (void)deltav;
    (void)deltaw;
    prj_riemann_gr_hlld_impl(WL, WR, pL, pR, gL, gR, eos, gamma,
        sqrt_gamma, alpha, beta, bn_tilde, flux, v_face, Bv1, Bv2);
}
#endif
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
