#ifndef PRJ_SRC_H
#define PRJ_SRC_H

void prj_src_geom(prj_eos *eos, double *W, double *dUdt);
void prj_src_user(prj_eos *eos, double *W, double *dUdt);
void prj_src_monopole_gravity(prj_mesh *mesh, const prj_grav_mono *grav_mono,
    double *restrict W, double *restrict U, double *restrict dUdt);
void prj_src_update(prj_mesh *mesh, prj_eos *eos, double *restrict W,
    double *restrict U, double *restrict dUdt);

#endif
