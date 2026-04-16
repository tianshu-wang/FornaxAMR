#ifndef PRJ_FLUX_H
#define PRJ_FLUX_H

void prj_flux_update(prj_eos *eos, prj_rad *rad, const prj_block *block, double *W, double *eosvar, double *flux[3]);
void prj_flux_div(double *flux[3], double area[3], double vol, double *dUdt);

#endif
