#ifndef PRJ_PERTURBATION_H
#define PRJ_PERTURBATION_H

#include "prj_types.h"

void prj_set_perturbation(prj_mesh *mesh, prj_eos *eos, const prj_mpi *mpi,
    double gaussian_norm, unsigned long long seed);

#endif
