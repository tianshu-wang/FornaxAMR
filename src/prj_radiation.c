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

void prj_rad_flux(const double *WL, const double *WR, const prj_eos *eos, const prj_rad *rad, double *flux)
{
    int field;
    int group;

    (void)eos;
    (void)rad;

    if (WL == 0 || WR == 0 || flux == 0) {
        return;
    }

    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            double el = WL[PRJ_PRIM_RAD_E(field, group)];
            double er = WR[PRJ_PRIM_RAD_E(field, group)];
            double f1l = WL[PRJ_PRIM_RAD_F1(field, group)];
            double f1r = WR[PRJ_PRIM_RAD_F1(field, group)];
            double ul_e = el;
            double ur_e = er;
            double ul_f1 = f1l;
            double ur_f1 = f1r;
            double fl_e = f1l;
            double fr_e = f1r;
            double fl_f1 = PRJ_CLIGHT * PRJ_CLIGHT * el;
            double fr_f1 = PRJ_CLIGHT * PRJ_CLIGHT * er;
            double sl = -PRJ_CLIGHT;
            double sr = PRJ_CLIGHT;

            flux[PRJ_CONS_RAD_E(field, group)] =
                (sr * fl_e - sl * fr_e + sl * sr * (ur_e - ul_e)) / (sr - sl);
            flux[PRJ_CONS_RAD_F1(field, group)] =
                (sr * fl_f1 - sl * fr_f1 + sl * sr * (ur_f1 - ul_f1)) / (sr - sl);
            flux[PRJ_CONS_RAD_F2(field, group)] = 0.0;
            flux[PRJ_CONS_RAD_F3(field, group)] = 0.0;
        }
    }
}
