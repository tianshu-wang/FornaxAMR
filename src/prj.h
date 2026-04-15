#ifndef PRJ_H
#define PRJ_H

#include "prj_defs.h"
#include "prj_types.h"
#include "prj_amr.h"
#include "prj_boundary.h"
#include "prj_eos.h"
#include "prj_flux.h"
#include "prj_gravity.h"
#include "prj_io.h"
#include "prj_mesh.h"
#include "prj_mpi.h"
#include "prj_radiation.h"
#include "prj_rad3_opac.h"
#include "prj_reconstruct.h"
#include "prj_riemann.h"
#include "prj_src.h"
#include "prj_timeint.h"
#include "prj_utils.h"

void prj_problem_initial_condition(double x1, double x2, double x3, double *data);
void prj_problem_general(prj_sim *sim);
void prj_problem_cc(prj_sim *sim);
void prj_problem_sedov(prj_sim *sim);
void prj_problem_sedov_offcenter(prj_sim *sim);
void prj_problem_shock1d(prj_sim *sim);

#endif
