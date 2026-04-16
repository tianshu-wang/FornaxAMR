#ifndef PRJ_RADIATION_H
#define PRJ_RADIATION_H

#define PRJ_CLIGHT 2.99792458e10

void prj_rad_init(prj_rad *rad);
void prj_rad_prim2cons(const double *W, double *U);
void prj_rad_cons2prim(const double *U, double *W);
void prj_rad_flux(const double *WL, const double *WR, const prj_eos *eos, const prj_rad *rad, double *flux);
void prj_rad_energy_update(prj_rad *rad, prj_eos *eos, double *u, double dt, double lapse, double *final_temperature);
void prj_rad_momentum_update(prj_rad *rad, prj_eos *eos, double *u, double dt, double lapse, double temperature);

#endif
