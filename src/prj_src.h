#ifndef PRJ_SRC_H
#define PRJ_SRC_H

void prj_src_geom(prj_eos *eos, double *W_mhd, double *W_rad,
    double *mhd_rhs, double *rad_rhs);
void prj_src_user(prj_eos *eos, double *W_mhd, double *W_rad,
    double *mhd_rhs, double *rad_rhs);
void prj_src_monopole_gravity(const prj_rad *rad, const prj_block *block,
    const prj_grav *grav, double *restrict W_mhd, double *restrict W_rad,
    double *restrict mhd_rhs, double *restrict rad_rhs);
void prj_src_radiation_vel_grad(const prj_rad *rad, const prj_block *block,
    double *restrict W_rad, double *restrict rad_rhs);
void prj_src_update(prj_eos *eos, const prj_rad *rad, const prj_grav *grav,
    const prj_mesh *mesh, const prj_block *block, int z4c_stage,
    double *restrict W_mhd, double *restrict W_rad,
    double *restrict mhd_rhs, double *restrict rad_rhs);

#endif
