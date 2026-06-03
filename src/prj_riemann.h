#ifndef PRJ_RIEMANN_H
#define PRJ_RIEMANN_H

#include "prj_types.h"

void prj_riemann_set_mesh(prj_mesh *mesh);
void prj_riemann_hllc(const double *WL, const double *WR,
    double pL, double pR, double gL, double gR,
    const prj_eos *eos, double *flux, double v_face[3],
    double deltau, double deltav, double deltaw);
#if PRJ_MHD
void prj_riemann_hlld(const double *WL, const double *WR,
    double pL, double pR, double gL, double gR,
    const prj_eos *eos, double bn, double *flux, double v_face[3],
    double *Bv1, double *Bv2, double deltau, double deltav, double deltaw);
void prj_riemann_lhlld(const double *WL, const double *WR,
    double pL, double pR, double gL, double gR,
    const prj_eos *eos, double bn, double *flux, double v_face[3],
    double *Bv1, double *Bv2, double deltau, double deltav, double deltaw);
#endif
void prj_riemann_flux_send(prj_mesh *mesh, const prj_mpi *mpi);
int prj_blocks_overlap_open(double amin, double amax, double bmin, double bmax);
double prj_overlap_length(double amin, double amax, double bmin, double bmax);
int prj_riemann_face_axis(const prj_block *block, const double slot_xmin[3],
    const double slot_xmax[3], int *side_out);
int prj_riemann_restrict_one(
    const prj_block *fine, int axis, int side,
    double cmin0, double cmax0, double cmin1, double cmax1,
    int tan0, int tan1,
    double *dst);

#endif
