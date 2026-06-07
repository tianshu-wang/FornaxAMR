#ifndef PRJ_RAD3_OPAC_H
#define PRJ_RAD3_OPAC_H

#include "prj_defs.h"
#include "prj_types.h"

#if PRJ_NRAD > 0

#define PRJ_MEV_TO_ERG 1.602176634e-6
#define PRJ_AVOGADRO 6.02214076e23

void prj_rad3_opac_init(prj_rad *rad);
void prj_rad3_opac_free(prj_rad *rad);
/* Any output pointer may be NULL; unrequested fields are skipped. */
void prj_rad3_opac_lookup(const prj_rad *rad, double rho, double temp, double ye,
    double *kappa, double *sigma, double *delta, double *eta);
/* Implicit-solver hot path: always writes kappa, eta and their log-derivatives
 * d(ln kappa)/d(lnT), d(ln kappa)/d(Ye), d(ln eta)/d(lnT), d(ln eta)/d(Ye)
 * (no NULL checks).  The table temperature axis is natural log(T), so the lnT
 * derivatives need no log10 conversion. */
void prj_rad3_opac_lookup_ke(const prj_rad *rad, double rho, double temp, double ye,
    double *kappa, double *eta,
    double *dlnkappa_dlnT, double *dlnkappa_dYe,
    double *dlneta_dlnT, double *dlneta_dYe);

#endif

#endif
