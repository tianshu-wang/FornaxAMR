#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "prj.h"

static prj_mesh *prj_riemann_flux_mesh = 0;

double prj_riemann_min_double(double a, double b)
{
    return a < b ? a : b;
}

double prj_riemann_max_double(double a, double b)
{
    return a > b ? a : b;
}

static double prj_abs_double(double x)
{
    return x < 0.0 ? -x : x;
}

#if PRJ_MHD
#define PRJ_HLLD_SMALL_NUMBER 1.0e-8

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
    double U[PRJ_NVAR_CONS];
    double F[PRJ_NVAR_CONS];
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
    double U[PRJ_NVAR_CONS];
} prj_hlld_star;

static void prj_riemann_hlld_fail(const char *message)
{
    fprintf(stderr, "prj_riemann_hlld: %s\n", message);
    exit(1);
}

static void prj_hlld_require_finite(double x, const char *name)
{
    if (!isfinite(x)) {
        fprintf(stderr, "prj_riemann_hlld: non-finite %s\n", name);
        exit(1);
    }
}

static void prj_hlld_fill_conserved(double rho, double vx, double vy, double vz,
    double bx, double by, double bz, double ye, double e, double *U)
{
    int v;

    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        U[v] = 0.0;
    }
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
    int v;

    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        F[v] = 0.0;
    }

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
    if (prj_abs_double(smdiff) <= 1.0e-14 * (prj_abs_double(S) + prj_abs_double(SM) + 1.0)) {
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
    if (prj_abs_double(denom) < PRJ_HLLD_SMALL_NUMBER * pt_star) {
        star->vy = s->vy;
        star->vz = s->vz;
        star->by = s->by;
        star->bz = s->bz;
    } else {
        star->vy = s->vy - s->bx * s->by * (SM - s->vx) / denom;
        star->vz = s->vz - s->bx * s->bz * (SM - s->vx) / denom;
        star->by = s->by * (s->rho * alpha * alpha - s->bx * s->bx) / denom;
        star->bz = s->bz * (s->rho * alpha * alpha - s->bx * s->bx) / denom;
    }

    vdotb = s->vx * s->bx + s->vy * s->by + s->vz * s->bz;
    vstar_dotbstar = star->vx * star->bx + star->vy * star->by + star->vz * star->bz;
    star->e = (alpha * s->e - s->pt * s->vx + pt_star * SM +
        s->bx * (vdotb - vstar_dotbstar)) / smdiff;
    prj_hlld_require_finite(star->e, "star total energy");

    prj_hlld_fill_star_conserved(star);
}

static void prj_hlld_double_star(const prj_hlld_star *left, const prj_hlld_star *right,
    double bx, prj_hlld_star *left_ss, prj_hlld_star *right_ss)
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
    double vbdot_l;
    double vbdot_r;

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

    *left_ss = *left;
    *right_ss = *right;
    left_ss->vy = vy_ss;
    left_ss->vz = vz_ss;
    left_ss->by = by_ss;
    left_ss->bz = bz_ss;
    right_ss->vy = vy_ss;
    right_ss->vz = vz_ss;
    right_ss->by = by_ss;
    right_ss->bz = bz_ss;

    vbdot_ss = left->vx * bx + vy_ss * by_ss + vz_ss * bz_ss;
    vbdot_l = left->vx * bx + left->vy * left->by + left->vz * left->bz;
    vbdot_r = right->vx * bx + right->vy * right->by + right->vz * right->bz;
    left_ss->e = left->e - sqrt_l * sign_bx * (vbdot_l - vbdot_ss);
    right_ss->e = right->e + sqrt_r * sign_bx * (vbdot_r - vbdot_ss);
    prj_hlld_require_finite(left_ss->e, "left double-star total energy");
    prj_hlld_require_finite(right_ss->e, "right double-star total energy");

    prj_hlld_fill_star_conserved(left_ss);
    prj_hlld_fill_star_conserved(right_ss);
}

static void prj_hlld_flux_from_jump(const double *F0, double S,
    const double *U1, const double *U0, double *F)
{
    int v;

    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        F[v] = F0[v] + S * (U1[v] - U0[v]);
    }
}

static void prj_hlld_flux_from_double_jump(const double *F0, double S0,
    const double *Ustar, const double *U0, double Sstar,
    const double *Uss, double *F)
{
    double Fstar[PRJ_NVAR_CONS];
    int v;

    prj_hlld_flux_from_jump(F0, S0, Ustar, U0, Fstar);
    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        F[v] = Fstar[v] + Sstar * (Uss[v] - Ustar[v]);
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
    double *Bv1, double *Bv2)
{
    prj_hlld_state L;
    prj_hlld_state R;
    prj_hlld_star Ls;
    prj_hlld_star Rs;
    prj_hlld_star Lss;
    prj_hlld_star Rss;
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

    (void)eos;

    if (WL == 0 || WR == 0 || flux == 0) {
        prj_riemann_hlld_fail("null input");
    }

    prj_hlld_state_from_prim(WL, pL, gL, bn, &L);
    prj_hlld_state_from_prim(WR, pR, gR, bn, &R);
    cfL = prj_hlld_fast_speed(&L);
    cfR = prj_hlld_fast_speed(&R);
    SL = prj_riemann_min_double(L.vx - cfL, R.vx - cfR);
    SR = prj_riemann_max_double(L.vx + cfL, R.vx + cfR);
    if (!(SL < SR)) {
        prj_riemann_hlld_fail("invalid fast-wave ordering");
    }

    denom = R.rho * (SR - R.vx) - L.rho * (SL - L.vx);
    if (prj_abs_double(denom) <= 1.0e-14 *
        (prj_abs_double(R.rho * (SR - R.vx)) + prj_abs_double(L.rho * (SL - L.vx)) + 1.0)) {
        prj_riemann_hlld_fail("degenerate contact-speed denominator");
    }
    SM = (R.rho * R.vx * (SR - R.vx) -
        L.rho * L.vx * (SL - L.vx) + L.pt - R.pt) / denom;
    prj_hlld_require_finite(SM, "contact speed");

    pt_star_l = L.pt + L.rho * (SL - L.vx) * (SM - L.vx);
    pt_star_r = R.pt + R.rho * (SR - R.vx) * (SM - R.vx);
    pt_star = 0.5 * (pt_star_l + pt_star_r);
    if (pt_star <= 0.0 || !isfinite(pt_star)) {
        prj_riemann_hlld_fail("invalid star total pressure");
    }

    prj_hlld_outer_star(&L, SL, SM, pt_star, &Ls);
    prj_hlld_outer_star(&R, SR, SM, pt_star, &Rs);
    SLs = SM - prj_abs_double(bn) / sqrt(Ls.rho);
    SRs = SM + prj_abs_double(bn) / sqrt(Rs.rho);
    prj_hlld_require_finite(SLs, "left Alfven speed");
    prj_hlld_require_finite(SRs, "right Alfven speed");
    if (0.5 * bn * bn < PRJ_HLLD_SMALL_NUMBER * pt_star) {
        Lss = Ls;
        Rss = Rs;
    } else {
        prj_hlld_double_star(&Ls, &Rs, bn, &Lss, &Rss);
    }

    if (0.0 <= SL) {
        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
            flux[v] = L.F[v];
        }
        prj_hlld_face_outputs(L.vx, L.vy, L.vz, L.bx, L.by, L.bz, v_face, Bv1, Bv2);
    } else if (0.0 <= SLs) {
        prj_hlld_flux_from_jump(L.F, SL, Ls.U, L.U, flux);
        prj_hlld_face_outputs(Ls.vx, Ls.vy, Ls.vz, Ls.bx, Ls.by, Ls.bz, v_face, Bv1, Bv2);
    } else if (0.0 <= SM) {
        prj_hlld_flux_from_double_jump(L.F, SL, Ls.U, L.U, SLs, Lss.U, flux);
        prj_hlld_face_outputs(Lss.vx, Lss.vy, Lss.vz, Lss.bx, Lss.by, Lss.bz, v_face, Bv1, Bv2);
    } else if (0.0 <= SRs) {
        prj_hlld_flux_from_double_jump(R.F, SR, Rs.U, R.U, SRs, Rss.U, flux);
        prj_hlld_face_outputs(Rss.vx, Rss.vy, Rss.vz, Rss.bx, Rss.by, Rss.bz, v_face, Bv1, Bv2);
    } else if (0.0 <= SR) {
        prj_hlld_flux_from_jump(R.F, SR, Rs.U, R.U, flux);
        prj_hlld_face_outputs(Rs.vx, Rs.vy, Rs.vz, Rs.bx, Rs.by, Rs.bz, v_face, Bv1, Bv2);
    } else {
        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
            flux[v] = R.F[v];
        }
        prj_hlld_face_outputs(R.vx, R.vy, R.vz, R.bx, R.by, R.bz, v_face, Bv1, Bv2);
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
    double etot;
    int v;

    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        U[v] = 0.0;
        F[v] = 0.0;
    }

    prj_eos_prim2cons((prj_eos *)eos, (double *)W, U);

    rho = W[PRJ_PRIM_RHO];
    v1 = W[PRJ_PRIM_V1];
    v2 = W[PRJ_PRIM_V2];
    v3 = W[PRJ_PRIM_V3];
    etot = U[PRJ_CONS_ETOT];

    *cs = sqrt(gamma * pressure / rho);

    F[PRJ_CONS_RHO] = rho * v1;
    F[PRJ_CONS_MOM1] = rho * v1 * v1 + pressure;
    F[PRJ_CONS_MOM2] = rho * v1 * v2;
    F[PRJ_CONS_MOM3] = rho * v1 * v3;
    F[PRJ_CONS_ETOT] = (etot + pressure) * v1;
    F[PRJ_CONS_YE] = rho * W[PRJ_PRIM_YE] * v1;
}

void prj_riemann_hlle(const double *WL, const double *WR,
    double pL, double pR, double gL, double gR,
    const prj_eos *eos, double *flux, double v_face[3])
{
    double UL[PRJ_NVAR_CONS];
    double UR[PRJ_NVAR_CONS];
    double FL[PRJ_NVAR_CONS];
    double FR[PRJ_NVAR_CONS];
    double csL;
    double csR;
    double SL;
    double SR;
    int v;

    prj_riemann_state(WL, pL, gL, eos, UL, FL, &csL);
    prj_riemann_state(WR, pR, gR, eos, UR, FR, &csR);

    SL = prj_riemann_min_double(WL[PRJ_PRIM_V1] - csL, WR[PRJ_PRIM_V1] - csR);
    SR = prj_riemann_max_double(WL[PRJ_PRIM_V1] + csL, WR[PRJ_PRIM_V1] + csR);

    if (0.0 <= SL) {
        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
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
        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
            flux[v] = FR[v];
        }
        if (v_face != 0) {
            v_face[0] = WR[PRJ_PRIM_V1];
            v_face[1] = WR[PRJ_PRIM_V2];
            v_face[2] = WR[PRJ_PRIM_V3];
        }
        return;
    }

    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        flux[v] = (SR * FL[v] - SL * FR[v] + SL * SR * (UR[v] - UL[v])) / (SR - SL);
    }
    if (v_face != 0) {
        double rho_hll = (SR * UR[PRJ_CONS_RHO] - SL * UL[PRJ_CONS_RHO] -
            (FR[PRJ_CONS_RHO] - FL[PRJ_CONS_RHO])) / (SR - SL);
        double mom1_hll = (SR * UR[PRJ_CONS_MOM1] - SL * UL[PRJ_CONS_MOM1] -
            (FR[PRJ_CONS_MOM1] - FL[PRJ_CONS_MOM1])) / (SR - SL);
        double mom2_hll = (SR * UR[PRJ_CONS_MOM2] - SL * UL[PRJ_CONS_MOM2] -
            (FR[PRJ_CONS_MOM2] - FL[PRJ_CONS_MOM2])) / (SR - SL);
        double mom3_hll = (SR * UR[PRJ_CONS_MOM3] - SL * UL[PRJ_CONS_MOM3] -
            (FR[PRJ_CONS_MOM3] - FL[PRJ_CONS_MOM3])) / (SR - SL);
        if (rho_hll != 0.0) {
            v_face[0] = mom1_hll / rho_hll;
            v_face[1] = mom2_hll / rho_hll;
            v_face[2] = mom3_hll / rho_hll;
        } else {
            v_face[0] = 0.0;
            v_face[1] = 0.0;
            v_face[2] = 0.0;
        }
    }
}

void prj_riemann_hllc(const double *WL, const double *WR,
    double pL, double pR, double gL, double gR,
    const prj_eos *eos, double *flux, double v_face[3])
{
    double UL[PRJ_NVAR_CONS];
    double UR[PRJ_NVAR_CONS];
    double FL[PRJ_NVAR_CONS];
    double FR[PRJ_NVAR_CONS];
    double Ustar[PRJ_NVAR_CONS];
    double csL;
    double csR;
    double SL;
    double SR;
    double SM;
    double rho_star;
    double e_over_rho;
    int v;

    prj_riemann_state(WL, pL, gL, eos, UL, FL, &csL);
    prj_riemann_state(WR, pR, gR, eos, UR, FR, &csR);

    SL = prj_riemann_min_double(WL[PRJ_PRIM_V1] - csL, WR[PRJ_PRIM_V1] - csR);
    SR = prj_riemann_max_double(WL[PRJ_PRIM_V1] + csL, WR[PRJ_PRIM_V1] + csR);

    /* Toro (2009), Eq. 10.37: contact wave speed from Rankine-Hugoniot conditions. */
    SM = (pR - pL +
        UL[PRJ_CONS_MOM1] * (SL - WL[PRJ_PRIM_V1]) -
        UR[PRJ_CONS_MOM1] * (SR - WR[PRJ_PRIM_V1])) /
        (UL[PRJ_CONS_RHO] * (SL - WL[PRJ_PRIM_V1]) -
            UR[PRJ_CONS_RHO] * (SR - WR[PRJ_PRIM_V1]));

    if (0.0 <= SL) {
        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
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
        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
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
        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
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
        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
            flux[v] = FL[v] + SL * (Ustar[v] - UL[v]);
        }
    } else {
        /* Toro (2009), Eqs. 10.33-10.36 with right star state. */
        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
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
        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
            flux[v] = FR[v] + SR * (Ustar[v] - UR[v]);
        }
    }
}

static int prj_blocks_overlap_open(double amin, double amax, double bmin, double bmax)
{
    const double tol = 1.0e-12;
    return prj_riemann_min_double(amax, bmax) - prj_riemann_max_double(amin, bmin) > tol;
}

static double prj_overlap_length(double amin, double amax, double bmin, double bmax)
{
    double overlap = prj_riemann_min_double(amax, bmax) - prj_riemann_max_double(amin, bmin);

    return overlap > 0.0 ? overlap : 0.0;
}

void prj_riemann_set_mesh(prj_mesh *mesh)
{
    prj_riemann_flux_mesh = mesh;
}

int prj_riemann_detect_shock(const double *WL, const double *WR, double pL, double pR)
{
    double pressure_ratio;
    double compression[3];
    int best_dir;

    if (WL[PRJ_PRIM_EINT] <= 0.0 || WR[PRJ_PRIM_EINT] <= 0.0 ||
        WL[PRJ_PRIM_RHO] <= 0.0 || WR[PRJ_PRIM_RHO] <= 0.0) {
        return -1;
    }

    pressure_ratio = prj_riemann_max_double(pL, pR) /
        prj_riemann_max_double(prj_riemann_min_double(pL, pR), 1.0e-12);
    if (pressure_ratio < 1.5 || WR[PRJ_PRIM_V1] - WL[PRJ_PRIM_V1] >= 0.0) {
        return -1;
    }

    compression[0] = prj_abs_double(WR[PRJ_PRIM_V1] - WL[PRJ_PRIM_V1]);
    compression[1] = prj_abs_double(WR[PRJ_PRIM_V2] - WL[PRJ_PRIM_V2]);
    compression[2] = prj_abs_double(WR[PRJ_PRIM_V3] - WL[PRJ_PRIM_V3]);
    best_dir = X1DIR;
    if (compression[X2DIR] > compression[best_dir]) {
        best_dir = X2DIR;
    }
    if (compression[X3DIR] > compression[best_dir]) {
        best_dir = X3DIR;
    }
    if (compression[best_dir] < 2.0 * prj_riemann_max_double(
            compression[(best_dir + 1) % 3],
            compression[(best_dir + 2) % 3])) {
        return -1;
    }
    return best_dir;
}

void prj_riemann_flux_send(prj_block *block)
{
    int n;

    if (prj_riemann_flux_mesh == 0 || block == 0) {
        return;
    }

    for (n = 0; n < 56; ++n) {
        int nid = block->slot[n].id;

        if (nid >= 0 && nid < prj_riemann_flux_mesh->nblocks) {
            prj_block *neighbor = &prj_riemann_flux_mesh->blocks[nid];
            int axis = -1;
            int side = -1;
            int d;

            if (neighbor->id < 0 || neighbor->active != 1) {
                continue;
            }
            for (d = 0; d < 3; ++d) {
                if (prj_abs_double(block->xmax[d] - neighbor->xmin[d]) < 1.0e-12) {
                    if (axis >= 0) {
                        axis = -2;
                        break;
                    }
                    axis = d;
                    side = 1;
                } else if (prj_abs_double(neighbor->xmax[d] - block->xmin[d]) < 1.0e-12) {
                    if (axis >= 0) {
                        axis = -2;
                        break;
                    }
                    axis = d;
                    side = 0;
                }
            }
            if (axis < 0) {
                continue;
            }
            for (d = 0; d < 3; ++d) {
                if (d == axis) {
                    continue;
                }
                if (!prj_blocks_overlap_open(block->xmin[d], block->xmax[d], neighbor->xmin[d], neighbor->xmax[d])) {
                    axis = -1;
                    break;
                }
            }
            if (axis < 0) {
                continue;
            }
            if (neighbor->dx[axis] <= block->dx[axis]) {
                continue;
            }
            if (neighbor->rank == block->rank) {
                int ratio = (int)(neighbor->dx[axis] / block->dx[axis] + 0.5);
                int it0;
                int it1;

                if (ratio < 1) {
                    ratio = 1;
                }

                for (it0 = 0; it0 < PRJ_BLOCK_SIZE; ++it0) {
                    for (it1 = 0; it1 < PRJ_BLOCK_SIZE; ++it1) {
                        int coarse_face[3] = {0, 0, 0};
                        int tan0 = (axis + 1) % 3;
                        int tan1 = (axis + 2) % 3;
                        double coarse_min0;
                        double coarse_max0;
                        double coarse_min1;
                        double coarse_max1;
                        int fine_i0_start;
                        int fine_i0_end;
                        int fine_i1_start;
                        int fine_i1_end;
                        int r0;
                        int r1;

                        coarse_face[axis] = side == 1 ? 0 : PRJ_BLOCK_SIZE;
                        coarse_face[tan0] = it0;
                        coarse_face[tan1] = it1;
                        coarse_min0 = neighbor->xmin[tan0] + (double)it0 * neighbor->dx[tan0];
                        coarse_max0 = coarse_min0 + neighbor->dx[tan0];
                        coarse_min1 = neighbor->xmin[tan1] + (double)it1 * neighbor->dx[tan1];
                        coarse_max1 = coarse_min1 + neighbor->dx[tan1];

                        if (!prj_blocks_overlap_open(coarse_min0, coarse_max0, block->xmin[tan0], block->xmax[tan0]) ||
                            !prj_blocks_overlap_open(coarse_min1, coarse_max1, block->xmin[tan1], block->xmax[tan1])) {
                            continue;
                        }

                        fine_i0_start = (int)((coarse_min0 - block->xmin[tan0]) / block->dx[tan0]);
                        fine_i0_end = (int)((coarse_max0 - block->xmin[tan0]) / block->dx[tan0]);
                        fine_i1_start = (int)((coarse_min1 - block->xmin[tan1]) / block->dx[tan1]);
                        fine_i1_end = (int)((coarse_max1 - block->xmin[tan1]) / block->dx[tan1]);
                        if (fine_i0_start < 0) {
                            fine_i0_start = 0;
                        }
                        if (fine_i1_start < 0) {
                            fine_i1_start = 0;
                        }
                        if (fine_i0_end >= PRJ_BLOCK_SIZE) {
                            fine_i0_end = PRJ_BLOCK_SIZE - 1;
                        }
                        if (fine_i1_end >= PRJ_BLOCK_SIZE) {
                            fine_i1_end = PRJ_BLOCK_SIZE - 1;
                        }
                        {
                            int v;

                            for (v = 0; v < PRJ_NVAR_CONS; ++v) {
                                double restricted_flux = 0.0;
                                double weight_sum = 0.0;

                                for (r0 = fine_i0_start; r0 <= fine_i0_end; ++r0) {
                                    for (r1 = fine_i1_start; r1 <= fine_i1_end; ++r1) {
                                        int fine_face[3] = {0, 0, 0};
                                        double fine_min0 = block->xmin[tan0] + (double)r0 * block->dx[tan0];
                                        double fine_max0 = fine_min0 + block->dx[tan0];
                                        double fine_min1 = block->xmin[tan1] + (double)r1 * block->dx[tan1];
                                        double fine_max1 = fine_min1 + block->dx[tan1];
                                        double overlap0 = prj_overlap_length(coarse_min0, coarse_max0, fine_min0, fine_max0);
                                        double overlap1 = prj_overlap_length(coarse_min1, coarse_max1, fine_min1, fine_max1);
                                        double weight = overlap0 * overlap1;

                                        if (weight <= 0.0) {
                                            continue;
                                        }
                                        fine_face[axis] = side == 1 ? PRJ_BLOCK_SIZE : 0;
                                        fine_face[tan0] = r0;
                                        fine_face[tan1] = r1;
                                        restricted_flux +=
                                            block->flux[axis][VIDX(v, fine_face[0], fine_face[1], fine_face[2])] * weight;
                                        weight_sum += weight;
                                    }
                                }
                                if (weight_sum > 0.0) {
                                    neighbor->flux[axis][VIDX(v, coarse_face[0], coarse_face[1], coarse_face[2])] =
                                        restricted_flux / weight_sum;
                                }
                            }
                        }
                    }
                }
            } else {
                /* MPI flux exchange buffer path not implemented yet. */
            }
        }
    }
}
