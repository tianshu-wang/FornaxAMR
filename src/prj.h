#ifndef PRJ_H
#define PRJ_H

#include "prj_defs.h"
#include "prj_types.h"
#include "prj_amr.h"
#include "prj_boundary.h"
#include "prj_diagnostics.h"
#include "prj_eos.h"
#include "prj_flux.h"
#include "prj_gravity.h"
#include "prj_io.h"
#include "prj_mesh.h"
#include "prj_mhd.h"
#include "prj_mpi.h"
#include "prj_perturbation.h"
#include "prj_radiation.h"
#include "prj_rad3_opac.h"
#include "prj_rad_inel.h"
#include "prj_reconstruct.h"
#include "prj_riemann.h"
#include "prj_src.h"
#include "prj_timer.h"
#include "prj_timeint.h"
#include "prj_utils.h"
#include "prj_z4c.h"

void prj_problem_initial_condition(double x1, double x2, double x3, double *data);
void prj_problem_general(prj_sim *sim, prj_mpi *mpi);
void prj_problem_cc(prj_sim *sim, prj_mpi *mpi);
void prj_problem_ccsn(prj_sim *sim, prj_mpi *mpi);
void prj_problem_sedov(prj_sim *sim, prj_mpi *mpi);
void prj_problem_sedov_offcenter(prj_sim *sim, prj_mpi *mpi);
void prj_problem_shock1d(prj_sim *sim, prj_mpi *mpi);
void prj_problem_kh(prj_sim *sim, prj_mpi *mpi);
void prj_problem_shocktube(prj_sim *sim, prj_mpi *mpi);
void prj_problem_rad_free_streaming(prj_sim *sim, prj_mpi *mpi);
void prj_problem_rad_diffusive_source(prj_sim *sim, prj_mpi *mpi);
void prj_problem_rad_picket_fence(prj_sim *sim, prj_mpi *mpi);
void prj_problem_user_boundary(const prj_mesh *mesh, const prj_block *block,
    double *W, int axis, int side, int i, int j, int k,
    const double position[3], double time_seconds);
void prj_problem_user_source(const prj_mesh *mesh, const prj_block *block,
    double *W_mhd, double *W_rad, double *mhd_rhs, double *rad_rhs);
void prj_problem_z4c_one_puncture(prj_sim *sim, prj_mpi *mpi);
void prj_problem_z4c_two_puncture(prj_sim *sim, prj_mpi *mpi);

#endif
