#ifndef PRJ_RAD3_OPAC_H
#define PRJ_RAD3_OPAC_H

#include "prj_defs.h"
#include "prj_types.h"

#if PRJ_NRAD > 0

#define PRJ_MEV_TO_ERG 1.602176634e-6
#define PRJ_AVOGADRO 6.02214076e23

void prj_rad3_opac_init(prj_rad *rad);
void prj_rad3_opac_free(prj_rad *rad);
void prj_rad3_opac_lookup(const prj_rad *rad, double rho, double temp, double ye,
    double *kappa, double *sigma, double *delta, double *eta);

#endif

#endif
