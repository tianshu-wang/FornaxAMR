#ifndef PRJ_IO_H
#define PRJ_IO_H

#include <stddef.h>

void prj_io_write_restart(const prj_mesh *mesh, const prj_mpi *mpi, double time, int step, int dump_count,
    double last_output_time, double last_restart_time, double dt);
void prj_io_read_restart(prj_mesh *mesh, const prj_eos *eos, prj_mpi *mpi, const char *filename,
    double *time, int *step, int *dump_count, double *last_output_time, double *last_restart_time,
    double *dt);
void prj_io_write_dump(const prj_mesh *mesh, const prj_grav *grav, const prj_mpi *mpi,
    int dump_index, int step, double time);
void prj_io_parser(prj_sim *sim, char *filename);
int prj_io_find_latest_restart(const char *dir, char *out_filename, size_t out_size, int *out_id);

#endif
