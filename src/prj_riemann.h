#ifndef PRJ_RIEMANN_H
#define PRJ_RIEMANN_H

#include "prj_types.h"

double prj_riemann_min_double(double a, double b);
double prj_riemann_max_double(double a, double b);
void prj_riemann_set_mesh(prj_mesh *mesh);
void prj_riemann_hlle(const double *WL, const double *WR, const prj_eos *eos, double *flux, double v_face[3]);
void prj_riemann_hllc(const double *WL, const double *WR, const prj_eos *eos, double *flux, double v_face[3]);
int prj_riemann_detect_shock(const double *WL, const double *WR, double pL, double pR);
void prj_riemann_flux_send(prj_block *block);

#endif
