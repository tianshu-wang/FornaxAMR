#ifndef PRJ_TIMEINT_H
#define PRJ_TIMEINT_H

typedef struct prj_timeint_imex_tableau {
    int nstages;
    const double *a_ex;
    const double *a_im;
    const double *b_ex;
    const double *b_im;
    double r;
} prj_timeint_imex_tableau;

extern const prj_timeint_imex_tableau PRJ_TIMEINT_TABLEAU_NAME;

void prj_timeint_init(const prj_timeint_imex_tableau *tableau);
double prj_timeint_calc_dt(const prj_mesh *mesh, prj_eos *eos, const prj_mpi *mpi,
    double cfl, const prj_timeint_imex_tableau *tableau);
void prj_timeint_stage1(prj_mesh *mesh, const prj_coord *coord, const prj_bc *bc, prj_eos *eos,
    prj_rad *rad, prj_grav *grav, prj_mpi *mpi, double dt, double *dt_src, prj_timer *timer);
void prj_timeint_stage2(prj_mesh *mesh, const prj_coord *coord, const prj_bc *bc, prj_eos *eos,
    prj_rad *rad, prj_grav *grav, prj_mpi *mpi, double dt, double *dt_src, prj_timer *timer);
void prj_timeint_eSSPRK_step(prj_mesh *mesh, const prj_coord *coord, const prj_bc *bc, prj_eos *eos,
    prj_rad *rad, prj_grav *grav, prj_mpi *mpi, double dt, double *dt_src, prj_timer *timer);
void prj_timeint_eSSPRK_final(prj_mesh *mesh, const prj_coord *coord, const prj_bc *bc, prj_eos *eos,
    prj_rad *rad, prj_grav *grav, prj_mpi *mpi, prj_timer *timer);
void prj_timeint_step_im(prj_mesh *mesh, const prj_coord *coord, const prj_bc *bc, prj_eos *eos,
    prj_rad *rad, prj_grav *grav, prj_mpi *mpi, const prj_timeint_imex_tableau *tableau,
    int stage, double dt, double dt_implicit, double *dt_src, prj_timer *timer);
void prj_timeint_step_ex(prj_mesh *mesh, const prj_coord *coord, const prj_bc *bc, prj_eos *eos,
    prj_rad *rad, prj_grav *grav, prj_mpi *mpi, const prj_timeint_imex_tableau *tableau,
    int stage, double dt, double *dt_src, prj_timer *timer);
void prj_timeint_step(prj_mesh *mesh, const prj_coord *coord, const prj_bc *bc, prj_eos *eos,
    prj_rad *rad, prj_grav *grav, prj_mpi *mpi, const prj_timeint_imex_tableau *tableau,
    double dt, double *dt_src, prj_timer *timer);

#endif
