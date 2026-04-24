#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#if defined(PRJ_ENABLE_MPI)
#include <mpi.h>
#endif

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

static void prj_riemann_state_hydro(const double *W, double pressure, double gamma,
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

#if PRJ_MHD
static void prj_riemann_abort(const char *message)
{
    prj_mpi *mpi = prj_mpi_current();

    fprintf(stderr, "%s\n", message);
#if defined(PRJ_ENABLE_MPI)
    if (mpi != 0 && mpi->totrank > 1) {
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
#endif
    exit(1);
}

static double prj_riemann_fast_speed_mhd(const double *W, double pressure, double gamma)
{
    double rho;
    double bx;
    double by;
    double bz;
    double a2;
    double b2;
    double disc;
    double cf2;

    rho = W[PRJ_PRIM_RHO];
    bx = W[PRJ_PRIM_B1];
    by = W[PRJ_PRIM_B2];
    bz = W[PRJ_PRIM_B3];
    if (rho <= 0.0 || pressure <= 0.0 || gamma <= 0.0) {
        prj_riemann_abort("prj_riemann_hlld: invalid input state");
    }

    a2 = gamma * pressure / rho;
    b2 = (bx * bx + by * by + bz * bz) / rho;
    disc = (a2 + b2) * (a2 + b2) - 4.0 * a2 * bx * bx / rho;
    if (disc < 0.0 && disc > -1.0e-12 * (a2 + b2) * (a2 + b2)) {
        disc = 0.0;
    }
    if (disc < 0.0) {
        prj_riemann_abort("prj_riemann_hlld: negative fast-wave discriminant");
    }
    cf2 = 0.5 * (a2 + b2 + sqrt(disc));
    if (cf2 < 0.0) {
        prj_riemann_abort("prj_riemann_hlld: negative fast-wave speed");
    }
    return sqrt(cf2);
}

static void prj_riemann_state_mhd(const double *W, double pressure, double gamma,
    const prj_eos *eos, double *U, double *F, double *cf)
{
    double rho;
    double v1;
    double v2;
    double v3;
    double b1;
    double b2;
    double b3;
    double etot;
    double ptot;
    double vdotb;
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
    b1 = W[PRJ_PRIM_B1];
    b2 = W[PRJ_PRIM_B2];
    b3 = W[PRJ_PRIM_B3];
    etot = U[PRJ_CONS_ETOT];
    ptot = pressure + 0.5 * (b1 * b1 + b2 * b2 + b3 * b3);
    vdotb = v1 * b1 + v2 * b2 + v3 * b3;

    *cf = prj_riemann_fast_speed_mhd(W, pressure, gamma);

    F[PRJ_CONS_RHO] = rho * v1;
    F[PRJ_CONS_MOM1] = rho * v1 * v1 + ptot - b1 * b1;
    F[PRJ_CONS_MOM2] = rho * v1 * v2 - b1 * b2;
    F[PRJ_CONS_MOM3] = rho * v1 * v3 - b1 * b3;
    F[PRJ_CONS_ETOT] = (etot + ptot) * v1 - b1 * vdotb;
    F[PRJ_CONS_YE] = rho * W[PRJ_PRIM_YE] * v1;
    F[PRJ_CONS_B1] = 0.0;
    F[PRJ_CONS_B2] = b2 * v1 - v2 * b1;
    F[PRJ_CONS_B3] = b3 * v1 - v3 * b1;
}

static void prj_riemann_store_face_state(double v1, double v2, double v3,
    double b1, double b2, double b3, double v_face[3], double emf_face[2])
{
    if (v_face != 0) {
        v_face[0] = v1;
        v_face[1] = v2;
        v_face[2] = v3;
    }
    if (emf_face != 0) {
        emf_face[0] = v3 * b1 - v1 * b3;
        emf_face[1] = v1 * b2 - v2 * b1;
    }
}
#endif

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

    prj_riemann_state_hydro(WL, pL, gL, eos, UL, FL, &csL);
    prj_riemann_state_hydro(WR, pR, gR, eos, UR, FR, &csR);

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

    prj_riemann_state_hydro(WL, pL, gL, eos, UL, FL, &csL);
    prj_riemann_state_hydro(WR, pR, gR, eos, UR, FR, &csR);

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

#if PRJ_MHD
void prj_riemann_hlld(const double *WL, const double *WR,
    double pL, double pR, double gL, double gR,
    const prj_eos *eos, double *flux, double v_face[3], double emf_face[2])
{
    const double tiny = 1.0e-12;
    double UL[PRJ_NVAR_CONS];
    double UR[PRJ_NVAR_CONS];
    double FL[PRJ_NVAR_CONS];
    double FR[PRJ_NVAR_CONS];
    double UstarL[PRJ_NVAR_CONS];
    double UstarR[PRJ_NVAR_CONS];
    double FStarL[PRJ_NVAR_CONS];
    double FStarR[PRJ_NVAR_CONS];
    double UssL[PRJ_NVAR_CONS];
    double UssR[PRJ_NVAR_CONS];
    double FssL[PRJ_NVAR_CONS];
    double FssR[PRJ_NVAR_CONS];
    double rhoL;
    double rhoR;
    double vxL;
    double vxR;
    double vyL;
    double vyR;
    double vzL;
    double vzR;
    double yeL;
    double yeR;
    double bx;
    double byL;
    double byR;
    double bzL;
    double bzR;
    double ptotL;
    double ptotR;
    double cfL;
    double cfR;
    double SL;
    double SR;
    double denom_sm;
    double SM;
    double ptotStarL;
    double ptotStarR;
    double ptotStar;
    double rhoStarL;
    double rhoStarR;
    double denL;
    double denR;
    double vyStarL;
    double vyStarR;
    double vzStarL;
    double vzStarR;
    double byStarL;
    double byStarR;
    double bzStarL;
    double bzStarR;
    double vdotbL;
    double vdotbR;
    double vdotbStarL;
    double vdotbStarR;
    double EstarL;
    double EstarR;
    double bxabs;
    int use_rotational;
    int v;

    prj_riemann_state_mhd(WL, pL, gL, eos, UL, FL, &cfL);
    prj_riemann_state_mhd(WR, pR, gR, eos, UR, FR, &cfR);

    rhoL = WL[PRJ_PRIM_RHO];
    rhoR = WR[PRJ_PRIM_RHO];
    vxL = WL[PRJ_PRIM_V1];
    vxR = WR[PRJ_PRIM_V1];
    vyL = WL[PRJ_PRIM_V2];
    vyR = WR[PRJ_PRIM_V2];
    vzL = WL[PRJ_PRIM_V3];
    vzR = WR[PRJ_PRIM_V3];
    yeL = WL[PRJ_PRIM_YE];
    yeR = WR[PRJ_PRIM_YE];
    bx = WL[PRJ_PRIM_B1];
    if (prj_abs_double(WR[PRJ_PRIM_B1] - bx) > 1.0e-10 *
            prj_riemann_max_double(1.0, prj_abs_double(bx))) {
        prj_riemann_abort("prj_riemann_hlld: left/right normal magnetic field mismatch");
    }
    byL = WL[PRJ_PRIM_B2];
    byR = WR[PRJ_PRIM_B2];
    bzL = WL[PRJ_PRIM_B3];
    bzR = WR[PRJ_PRIM_B3];
    ptotL = pL + 0.5 * (bx * bx + byL * byL + bzL * bzL);
    ptotR = pR + 0.5 * (bx * bx + byR * byR + bzR * bzR);
    SL = prj_riemann_min_double(vxL - cfL, vxR - cfR);
    SR = prj_riemann_max_double(vxL + cfL, vxR + cfR);

    if (0.0 <= SL) {
        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
            flux[v] = FL[v];
        }
        prj_riemann_store_face_state(vxL, vyL, vzL, bx, byL, bzL, v_face, emf_face);
        return;
    }
    if (SR <= 0.0) {
        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
            flux[v] = FR[v];
        }
        prj_riemann_store_face_state(vxR, vyR, vzR, bx, byR, bzR, v_face, emf_face);
        return;
    }

    denom_sm = rhoR * (SR - vxR) - rhoL * (SL - vxL);
    if (prj_abs_double(denom_sm) < tiny) {
        prj_riemann_abort("prj_riemann_hlld: singular contact-wave denominator");
    }
    SM = (rhoR * vxR * (SR - vxR) - rhoL * vxL * (SL - vxL) - ptotR + ptotL) / denom_sm;
    ptotStarL = ptotL + rhoL * (SL - vxL) * (SM - vxL);
    ptotStarR = ptotR + rhoR * (SR - vxR) * (SM - vxR);
    ptotStar = 0.5 * (ptotStarL + ptotStarR);
    rhoStarL = rhoL * (SL - vxL) / (SL - SM);
    rhoStarR = rhoR * (SR - vxR) / (SR - SM);
    if (rhoStarL <= 0.0 || rhoStarR <= 0.0) {
        prj_riemann_abort("prj_riemann_hlld: non-positive star density");
    }

    denL = rhoL * (SL - vxL) * (SL - SM) - bx * bx;
    denR = rhoR * (SR - vxR) * (SR - SM) - bx * bx;
    if (prj_abs_double(denL) > tiny) {
        vyStarL = vyL - bx * byL * (SM - vxL) / denL;
        vzStarL = vzL - bx * bzL * (SM - vxL) / denL;
        byStarL = byL * (rhoL * (SL - vxL) * (SL - vxL) - bx * bx) / denL;
        bzStarL = bzL * (rhoL * (SL - vxL) * (SL - vxL) - bx * bx) / denL;
    } else {
        vyStarL = vyL;
        vzStarL = vzL;
        byStarL = byL;
        bzStarL = bzL;
    }
    if (prj_abs_double(denR) > tiny) {
        vyStarR = vyR - bx * byR * (SM - vxR) / denR;
        vzStarR = vzR - bx * bzR * (SM - vxR) / denR;
        byStarR = byR * (rhoR * (SR - vxR) * (SR - vxR) - bx * bx) / denR;
        bzStarR = bzR * (rhoR * (SR - vxR) * (SR - vxR) - bx * bx) / denR;
    } else {
        vyStarR = vyR;
        vzStarR = vzR;
        byStarR = byR;
        bzStarR = bzR;
    }

    vdotbL = vxL * bx + vyL * byL + vzL * bzL;
    vdotbR = vxR * bx + vyR * byR + vzR * bzR;
    vdotbStarL = SM * bx + vyStarL * byStarL + vzStarL * bzStarL;
    vdotbStarR = SM * bx + vyStarR * byStarR + vzStarR * bzStarR;
    EstarL = ((SL - vxL) * UL[PRJ_CONS_ETOT] - ptotL * vxL +
        ptotStar * SM + bx * (vdotbL - vdotbStarL)) / (SL - SM);
    EstarR = ((SR - vxR) * UR[PRJ_CONS_ETOT] - ptotR * vxR +
        ptotStar * SM + bx * (vdotbR - vdotbStarR)) / (SR - SM);

    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        UstarL[v] = 0.0;
        UstarR[v] = 0.0;
    }
    UstarL[PRJ_CONS_RHO] = rhoStarL;
    UstarL[PRJ_CONS_MOM1] = rhoStarL * SM;
    UstarL[PRJ_CONS_MOM2] = rhoStarL * vyStarL;
    UstarL[PRJ_CONS_MOM3] = rhoStarL * vzStarL;
    UstarL[PRJ_CONS_ETOT] = EstarL;
    UstarL[PRJ_CONS_YE] = rhoStarL * yeL;
    UstarL[PRJ_CONS_B1] = bx;
    UstarL[PRJ_CONS_B2] = byStarL;
    UstarL[PRJ_CONS_B3] = bzStarL;

    UstarR[PRJ_CONS_RHO] = rhoStarR;
    UstarR[PRJ_CONS_MOM1] = rhoStarR * SM;
    UstarR[PRJ_CONS_MOM2] = rhoStarR * vyStarR;
    UstarR[PRJ_CONS_MOM3] = rhoStarR * vzStarR;
    UstarR[PRJ_CONS_ETOT] = EstarR;
    UstarR[PRJ_CONS_YE] = rhoStarR * yeR;
    UstarR[PRJ_CONS_B1] = bx;
    UstarR[PRJ_CONS_B2] = byStarR;
    UstarR[PRJ_CONS_B3] = bzStarR;

    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        FStarL[v] = FL[v] + SL * (UstarL[v] - UL[v]);
        FStarR[v] = FR[v] + SR * (UstarR[v] - UR[v]);
    }

    bxabs = prj_abs_double(bx);
    use_rotational = bxabs > tiny;
    if (!use_rotational) {
        if (SM >= 0.0) {
            for (v = 0; v < PRJ_NVAR_CONS; ++v) {
                flux[v] = FStarL[v];
            }
            prj_riemann_store_face_state(SM, vyStarL, vzStarL, bx, byStarL, bzStarL, v_face, emf_face);
        } else {
            for (v = 0; v < PRJ_NVAR_CONS; ++v) {
                flux[v] = FStarR[v];
            }
            prj_riemann_store_face_state(SM, vyStarR, vzStarR, bx, byStarR, bzStarR, v_face, emf_face);
        }
        return;
    }

    {
        double sqrtRhoL = sqrt(rhoStarL);
        double sqrtRhoR = sqrt(rhoStarR);
        double signBx = bx >= 0.0 ? 1.0 : -1.0;
        double vySS;
        double vzSS;
        double bySS;
        double bzSS;
        double vdotbSS;
        double EssL;
        double EssR;
        double SAL;
        double SAR;
        double inv_sum = 1.0 / (sqrtRhoL + sqrtRhoR);

        vySS = (sqrtRhoL * vyStarL + sqrtRhoR * vyStarR +
            (byStarR - byStarL) * signBx) * inv_sum;
        vzSS = (sqrtRhoL * vzStarL + sqrtRhoR * vzStarR +
            (bzStarR - bzStarL) * signBx) * inv_sum;
        bySS = (sqrtRhoL * byStarR + sqrtRhoR * byStarL +
            sqrtRhoL * sqrtRhoR * (vyStarR - vyStarL) * signBx) * inv_sum;
        bzSS = (sqrtRhoL * bzStarR + sqrtRhoR * bzStarL +
            sqrtRhoL * sqrtRhoR * (vzStarR - vzStarL) * signBx) * inv_sum;
        vdotbSS = SM * bx + vySS * bySS + vzSS * bzSS;
        EssL = EstarL - sqrtRhoL * (vdotbStarL - vdotbSS) * signBx;
        EssR = EstarR + sqrtRhoR * (vdotbStarR - vdotbSS) * signBx;
        SAL = SM - bxabs / sqrtRhoL;
        SAR = SM + bxabs / sqrtRhoR;

        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
            UssL[v] = 0.0;
            UssR[v] = 0.0;
        }
        UssL[PRJ_CONS_RHO] = rhoStarL;
        UssL[PRJ_CONS_MOM1] = rhoStarL * SM;
        UssL[PRJ_CONS_MOM2] = rhoStarL * vySS;
        UssL[PRJ_CONS_MOM3] = rhoStarL * vzSS;
        UssL[PRJ_CONS_ETOT] = EssL;
        UssL[PRJ_CONS_YE] = rhoStarL * yeL;
        UssL[PRJ_CONS_B1] = bx;
        UssL[PRJ_CONS_B2] = bySS;
        UssL[PRJ_CONS_B3] = bzSS;

        UssR[PRJ_CONS_RHO] = rhoStarR;
        UssR[PRJ_CONS_MOM1] = rhoStarR * SM;
        UssR[PRJ_CONS_MOM2] = rhoStarR * vySS;
        UssR[PRJ_CONS_MOM3] = rhoStarR * vzSS;
        UssR[PRJ_CONS_ETOT] = EssR;
        UssR[PRJ_CONS_YE] = rhoStarR * yeR;
        UssR[PRJ_CONS_B1] = bx;
        UssR[PRJ_CONS_B2] = bySS;
        UssR[PRJ_CONS_B3] = bzSS;

        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
            FssL[v] = FStarL[v] + SAL * (UssL[v] - UstarL[v]);
            FssR[v] = FStarR[v] + SAR * (UssR[v] - UstarR[v]);
        }

        if (SAL >= 0.0) {
            for (v = 0; v < PRJ_NVAR_CONS; ++v) {
                flux[v] = FStarL[v];
            }
            prj_riemann_store_face_state(SM, vyStarL, vzStarL, bx, byStarL, bzStarL, v_face, emf_face);
        } else if (SM >= 0.0) {
            for (v = 0; v < PRJ_NVAR_CONS; ++v) {
                flux[v] = FssL[v];
            }
            prj_riemann_store_face_state(SM, vySS, vzSS, bx, bySS, bzSS, v_face, emf_face);
        } else if (SAR > 0.0) {
            for (v = 0; v < PRJ_NVAR_CONS; ++v) {
                flux[v] = FssR[v];
            }
            prj_riemann_store_face_state(SM, vySS, vzSS, bx, bySS, bzSS, v_face, emf_face);
        } else {
            for (v = 0; v < PRJ_NVAR_CONS; ++v) {
                flux[v] = FStarR[v];
            }
            prj_riemann_store_face_state(SM, vyStarR, vzStarR, bx, byStarR, bzStarR, v_face, emf_face);
        }
    }
}
#endif

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
