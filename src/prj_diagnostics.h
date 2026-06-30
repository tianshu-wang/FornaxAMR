#ifndef PRJ_DIAGNOSTICS_H
#define PRJ_DIAGNOSTICS_H

#include "prj_types.h"

void prj_diagnostics_truncate_dqdt(const prj_mpi *mpi);
void prj_diagnostics_write_dqdt(const prj_mesh *mesh, const prj_mpi *mpi, double time);

#endif
