#ifndef PRJ_RIEMANN_H
#define PRJ_RIEMANN_H

#include "prj_types.h"

void prj_riemann_set_mesh(prj_mesh *mesh);
void prj_riemann_hlle(const double *WL, const double *WR,
    double pL, double pR, double gL, double gR,
    const prj_eos *eos, double *flux, double v_face[3]);
void prj_riemann_hllc(const double *WL, const double *WR,
    double pL, double pR, double gL, double gR,
    const prj_eos *eos, double *flux, double v_face[3]);
#if PRJ_MHD
void prj_riemann_hlld(const double *WL, const double *WR,
    double pL, double pR, double gL, double gR,
    const prj_eos *eos, double bn, double *flux, double v_face[3],
    double *Bv1, double *Bv2, double deltau, double deltav, double deltaw);
#endif
int prj_riemann_detect_shock(const double *WL, const double *WR, double pL, double pR);
void prj_riemann_flux_send(prj_mesh *mesh);

#endif
