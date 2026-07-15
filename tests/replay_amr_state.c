#include <stdio.h>
#include <string.h>

#include "prj.h"

int main(void)
{
    prj_eos eos;
    prj_eos_gr_geom geom;
    double U[PRJ_NVAR_CONS] = {0.0};
    double W[PRJ_NVAR_PRIM] = {0.0};
    double Urt[PRJ_NVAR_CONS] = {0.0};
    int status;
    int v;

    memset(&eos, 0, sizeof(eos));
    eos.kind = PRJ_EOS_KIND_TABLE;
    strncpy(eos.filename,
        "../eos_tmp/SFHoEOS__ye__0.035_0.56_50__logT_-4.793_2.176_500__logrho_-8.699_15.5_500_extend.dat",
        sizeof(eos.filename) - 1);
    prj_eos_init(&eos, 0);

    memset(&geom, 0, sizeof(geom));
    geom.gamma[0][0] = 1.01270104755041346;
    geom.gamma[1][1] = 1.01270104755041346;
    geom.gamma[2][2] = 1.01270104755041346;

    U[PRJ_CONS_RHO] = 1.93689471268773985e9;
    U[PRJ_CONS_MOM1] = -4.23851859448649920e16;
    U[PRJ_CONS_MOM2] = -4.23851859448649920e16;
    U[PRJ_CONS_MOM3] = -4.23851859448649920e16;
    U[PRJ_CONS_ETOT] = 2.80445979841485171e27;
    U[PRJ_CONS_YE] = 8.84555906614651322e8;

    status = prj_eos_gr_cons2prim(&eos, &geom, U, W, PRJ_EOS_CTX_AMR);
    printf("cons2prim status=%d\n", status);
    if (status != PRJ_EOS_GR_OK) {
        return 1;
    }
    printf("W rho=%.17e v=(%.17e %.17e %.17e) eint=%.17e Ye=%.17e\n",
        W[PRJ_PRIM_RHO], W[PRJ_PRIM_V1], W[PRJ_PRIM_V2],
        W[PRJ_PRIM_V3], W[PRJ_PRIM_EINT], W[PRJ_PRIM_YE]);
    status = prj_eos_gr_prim2cons(&eos, &geom, W, Urt, PRJ_EOS_CTX_AMR);
    printf("prim2cons status=%d\n", status);
    for (v = 0; v < PRJ_NHYDRO; ++v) {
        double rel = U[v] != 0.0 ? (Urt[v] - U[v]) / U[v] : Urt[v];

        printf("U[%d] in=%.17e out=%.17e rel=%.9e\n", v, U[v], Urt[v], rel);
    }
    return 0;
}
