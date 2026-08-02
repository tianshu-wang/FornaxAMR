#include "prj_problem_rad_common.h"

void prj_problem_rad_free_streaming(prj_sim *sim, prj_mpi *mpi)
{
    (void)mpi;
    prj_rad_test_fill_problem(sim, 1);
}
