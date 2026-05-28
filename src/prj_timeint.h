#ifndef PRJ_TIMEINT_H
#define PRJ_TIMEINT_H

double prj_timeint_calc_dt(const prj_mesh *mesh, prj_eos *eos, const prj_mpi *mpi, double cfl);
void prj_timeint_stage1(prj_mesh *mesh, const prj_coord *coord, const prj_bc *bc, prj_eos *eos,
    prj_rad *rad, prj_grav *grav, prj_mpi *mpi, double dt, double *dt_src, prj_timer *timer);
void prj_timeint_stage2(prj_mesh *mesh, const prj_coord *coord, const prj_bc *bc, prj_eos *eos,
    prj_rad *rad, prj_grav *grav, prj_mpi *mpi, double dt, double *dt_src, prj_timer *timer);
void prj_timeint_step(prj_mesh *mesh, const prj_coord *coord, const prj_bc *bc, prj_eos *eos,
    prj_rad *rad, prj_grav *grav, prj_mpi *mpi, double dt, double *dt_src, prj_timer *timer);

#endif
