#ifndef PRJ_SRC_H
#define PRJ_SRC_H

void prj_src_geom(prj_eos *eos, double *W, double *dUdt);
void prj_src_user(prj_eos *eos, double *W, double *dUdt);
void prj_src_monopole_gravity(const prj_block *block, const prj_grav_mono *grav_mono,
    double *restrict W, double *restrict dUdt);
void prj_src_radiation_vel_grad(const prj_block *block,
    double *restrict W, double *restrict dUdt);
void prj_src_update(prj_eos *eos, const prj_block *block, double *restrict W,
    double *restrict dUdt);

#endif
