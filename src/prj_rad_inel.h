#ifndef PRJ_RAD_INEL_H
#define PRJ_RAD_INEL_H

#include "prj_defs.h"
#include "prj_types.h"

#if PRJ_NRAD > 0

void prj_rad_eleinel_init(prj_rad *rad);
void prj_rad_eleinel_free(prj_rad *rad);
void prj_rad_eleinel_lookup(const prj_rad *rad,
    double rho, double T, double Ye,
    double etael,
    const double *je, const double *he,
    double *source, double *sink, double *scatt);
void prj_rad_eleinel_step(prj_rad *rad, prj_eos *eos, double *u, double dt, double T_cell);

int prj_rad_nucinel_step(prj_rad *rad, prj_eos *eos, double *u, double dt, double T_cell);

#if PRJ_USE_RADIATION_FSA
void prj_rad_eleinel_fsa(prj_rad *rad, prj_eos *eos, double *u, double dt, double T_cell);
int prj_rad_nucinel_fsa(prj_rad *rad, prj_eos *eos, double *u, double dt, double T_cell);
#endif

#endif

#endif
