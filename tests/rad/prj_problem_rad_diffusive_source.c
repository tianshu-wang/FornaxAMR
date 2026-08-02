#include "prj_problem_rad_common.h"

void prj_problem_rad_diffusive_source(prj_sim *sim, prj_mpi *mpi)
{
    (void)mpi;
    prj_rad_test_fill_problem(sim, 2);
}
