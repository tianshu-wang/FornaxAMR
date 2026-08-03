#ifndef PRJ_RAD3_OPAC_H
#define PRJ_RAD3_OPAC_H

#include "prj_defs.h"
#include "prj_types.h"

#if PRJ_NRAD > 0

#define PRJ_MEV_TO_ERG 1.602176634e-6
#define PRJ_AVOGADRO 6.02214076e23
#define PRJ_RAD3_OPAC_NGROUPS (PRJ_NRAD * PRJ_NEGROUP)

typedef struct prj_rad3_opac_interp_result {
    double kappa[PRJ_RAD3_OPAC_NGROUPS];
    double sigma[PRJ_RAD3_OPAC_NGROUPS];
    double delta[PRJ_RAD3_OPAC_NGROUPS];
    double eta[PRJ_RAD3_OPAC_NGROUPS];
    double kappa_raw_slope[3][PRJ_RAD3_OPAC_NGROUPS];
    double sigma_raw_slope[3][PRJ_RAD3_OPAC_NGROUPS];
    double delta_raw_slope[3][PRJ_RAD3_OPAC_NGROUPS];
    double eta_raw_slope[3][PRJ_RAD3_OPAC_NGROUPS];
    double coord_scale[3];
    double inv_rho;
    double inv_temp;
    int physical_derivatives;
} prj_rad3_opac_interp_result;

#if PRJ_OPAC_PROVIDER == PRJ_PROVIDER_USER
typedef void (*prj_user_opac_fn)(const prj_rad *rad, int species, int group,
    double rho, double temperature, double ye, double *value,
    double derivatives[3]);
void opac(const prj_rad *rad, int species, int group, double rho,
    double temperature, double ye, double *value, double derivatives[3]);
void emis(const prj_rad *rad, int species, int group, double rho,
    double temperature, double ye, double *value, double derivatives[3]);
void scat(const prj_rad *rad, int species, int group, double rho,
    double temperature, double ye, double *value, double derivatives[3]);
void delta(const prj_rad *rad, int species, int group, double rho,
    double temperature, double ye, double *value, double derivatives[3]);
#endif

void prj_rad3_opac_init(prj_rad *rad);
void prj_rad3_opac_free(prj_rad *rad);
/* Any output pointer may be NULL; unrequested fields are skipped. */
void prj_rad3_opac_lookup(const prj_rad *rad, double rho, double temp, double ye,
    double *kappa, double *sigma, double *delta, double *eta);
void prj_rad3_opac_lookup_interp(const prj_rad *rad, double rho, double temp,
    double ye, prj_rad3_opac_interp_result *result);
int prj_rad3_opac_interp_group_derivs(
    const prj_rad3_opac_interp_result *result, int group,
    double dkappa[3], double dsigma[3], double ddelta[3], double deta[3]);
void prj_rad3_opac_lookup_derivs(const prj_rad *rad, double rho, double temp,
    double ye, double *kappa, double *sigma, double *delta, double *eta,
    double *dkappa_drho, double *dkappa_dT, double *dkappa_dYe,
    double *dsigma_drho, double *dsigma_dT, double *dsigma_dYe,
    double *ddelta_drho, double *ddelta_dT, double *ddelta_dYe,
    double *deta_drho, double *deta_dT, double *deta_dYe);
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
