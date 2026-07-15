/* Round-trip check of the GR prim<->cons transforms against the TABULAR SFHo EOS
 * at the exact state of the failing cell (block 14364, cell 2,1,7), plus a sweep
 * in eint to find where the round trip breaks. Standalone diagnostic. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "prj.h"

#define EOS_FILE "/Users/tianshuw/vibe-coding/test-codex/eos_tmp/" \
    "SFHoEOS__ye__0.035_0.56_50__logT_-4.793_2.176_500__logrho_-8.699_15.5_500_extend.dat"

static const char *const NM[6] = {"RHO", "MOM1", "MOM2", "MOM3", "ETOT", "YE"};

static void set_geom(prj_eos_gr_geom *geom, double chi)
{
    /* CFC: gamma_ij = delta_ij / chi  (matches prj_z4c_load_hydro_geom). */
    memset(geom, 0, sizeof(*geom));
    geom->gamma[0][0] = 1.0 / chi;
    geom->gamma[1][1] = 1.0 / chi;
    geom->gamma[2][2] = 1.0 / chi;
}

static void set_prim(double *W, double rho, double v1, double v2, double v3,
    double eint, double ye)
{
    int v;
    for (v = 0; v < PRJ_NVAR_PRIM; ++v) W[v] = 0.0;
    W[PRJ_PRIM_RHO] = rho; W[PRJ_PRIM_V1] = v1; W[PRJ_PRIM_V2] = v2;
    W[PRJ_PRIM_V3] = v3; W[PRJ_PRIM_EINT] = eint; W[PRJ_PRIM_YE] = ye;
}

static int roundtrip(prj_eos *eos, const prj_eos_gr_geom *geom, const double *W,
    double *Wout, int verbose)
{
    double U[PRJ_NVAR_CONS];
    int sp = prj_eos_gr_prim2cons(eos, geom, W, U, PRJ_EOS_CTX_MAIN);
    int sc;
    if (sp != PRJ_EOS_GR_OK) {
        if (verbose) fprintf(stderr, "  prim2cons status=%d\n", sp);
        return sp;
    }
    if (verbose) {
        int v;
        fprintf(stderr, "  U (densitized):\n");
        for (v = 0; v < 6; ++v)
            fprintf(stderr, "    U[%d] %-5s = %.17e\n", v, NM[v], U[v]);
    }
    sc = prj_eos_gr_cons2prim(eos, geom, U, Wout, PRJ_EOS_CTX_MAIN);
    if (sc != PRJ_EOS_GR_OK) {
        if (verbose) fprintf(stderr, "  cons2prim status=%d\n", sc);
        return 100 + sc;
    }
    return 0;
}

int main(void)
{
    prj_eos eos;
    prj_eos_gr_geom geom;
    double W[PRJ_NVAR_PRIM], Wout[PRJ_NVAR_PRIM];
    double chi = 0.980646;              /* failing cell metric (dump) */
    double rho = 1.596616e10, ye = 0.430379;
    double v1 = -3.548e6, v2 = -2.129e6, v3 = 7.096e5;
    double eint0 = 2.909635e18;
    int st, k;

    memset(&eos, 0, sizeof(eos));
    eos.kind = PRJ_EOS_KIND_TABLE;
    strncpy(eos.filename, EOS_FILE, sizeof(eos.filename) - 1);
    prj_eos_init(&eos, 0);
    if (eos.table_loaded != 1) { fprintf(stderr, "EOS table failed to load\n"); return 1; }
    set_geom(&geom, chi);

    fprintf(stderr, "=== round trip at the exact INITIAL failing-cell state ===\n");
    fprintf(stderr, "rho=%.6e eint=%.6e ye=%.6f v=(%.3e,%.3e,%.3e) chi=%.6f\n",
        rho, eint0, ye, v1, v2, v3, chi);
    set_prim(W, rho, v1, v2, v3, eint0, ye);
    st = roundtrip(&eos, &geom, W, Wout, 1);
    if (st == 0) {
        fprintf(stderr, "  ROUND TRIP OK. deltas:\n");
        fprintf(stderr, "    d_rho =%.3e (rel %.3e)\n",
            Wout[PRJ_PRIM_RHO]-rho, (Wout[PRJ_PRIM_RHO]-rho)/rho);
        fprintf(stderr, "    d_eint=%.3e (rel %.3e)\n",
            Wout[PRJ_PRIM_EINT]-eint0, (Wout[PRJ_PRIM_EINT]-eint0)/eint0);
        fprintf(stderr, "    d_v1  =%.3e  d_v2=%.3e  d_v3=%.3e\n",
            Wout[PRJ_PRIM_V1]-v1, Wout[PRJ_PRIM_V2]-v2, Wout[PRJ_PRIM_V3]-v3);
        fprintf(stderr, "    d_ye  =%.3e\n", Wout[PRJ_PRIM_YE]-ye);
    } else {
        fprintf(stderr, "  ROUND TRIP FAILED status=%d\n", st);
    }

    fprintf(stderr, "\n=== eint sweep (does round trip survive as eint drops?) ===\n");
    fprintf(stderr, " eint          rel_to_init   status  eint_out(rel err)\n");
    for (k = 0; k <= 20; ++k) {
        double frac = 1.0 - 0.005 * k;      /* 100% down to 90% of eint0 */
        double e = eint0 * frac;
        set_prim(W, rho, v1, v2, v3, e, ye);
        st = roundtrip(&eos, &geom, W, Wout, 0);
        if (st == 0)
            fprintf(stderr, " %.6e  %6.2f%%       OK      %.6e (%.2e)\n",
                e, (frac-1)*100, Wout[PRJ_PRIM_EINT], (Wout[PRJ_PRIM_EINT]-e)/e);
        else
            fprintf(stderr, " %.6e  %6.2f%%       FAIL(%d)\n", e, (frac-1)*100, st);
    }

    /* The exact EVOLVED conserved state from the crash -> cons2prim, swept over
     * the metric chi (the one unknown at step 1). eint=U[ETOT]/U[RHO] is
     * metric-independent, but rho=U[RHO]/sqrt_det (hence the cold floor) is not. */
    fprintf(stderr, "\n=== crash's evolved U -> cons2prim vs metric chi ===\n");
    fprintf(stderr, " (production used the step-1 metric; dump step-0 chi=0.980646)\n");
    fprintf(stderr, "   chi      sqrt_det     rho=D/w        status  rho_out      eint_out\n");
    {
        double U[PRJ_NVAR_CONS];
        int v, kc;
        for (v = 0; v < PRJ_NVAR_CONS; ++v) U[v] = 0.0;
        U[PRJ_CONS_RHO]  = 1.80182020941170197e10;
        U[PRJ_CONS_MOM1] = -1.03404871635359888e17;
        U[PRJ_CONS_MOM2] = -5.75520776456552080e16;
        U[PRJ_CONS_MOM3] = 4.51447346641542640e16;
        U[PRJ_CONS_ETOT] = 4.82337588310982163e28;
        U[PRJ_CONS_YE]   = 7.75469187605319881e9;
        for (kc = 0; kc <= 24; ++kc) {
            double c = 0.90 + 0.005 * kc;    /* chi from 0.90 to 1.02 */
            double sd, rho_approx;
            set_geom(&geom, c);
            sd = pow(c, -1.5);
            rho_approx = U[PRJ_CONS_RHO] / sd;
            st = prj_eos_gr_cons2prim(&eos, &geom, U, Wout, PRJ_EOS_CTX_MAIN);
            if (st == PRJ_EOS_GR_OK)
                fprintf(stderr, "  %.4f   %.6f   %.6e   OK      %.6e %.6e\n",
                    c, sd, rho_approx, Wout[PRJ_PRIM_RHO], Wout[PRJ_PRIM_EINT]);
            else
                fprintf(stderr, "  %.4f   %.6f   %.6e   FAIL(%d)\n",
                    c, sd, rho_approx, st);
        }
    }
    return 0;
}
