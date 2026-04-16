#ifndef PRJ_IO_H
#define PRJ_IO_H

void prj_io_write_restart(const prj_mesh *mesh, double time, int step, int dump_count,
    double next_output_time, double next_restart_time);
void prj_io_read_restart(prj_mesh *mesh, const prj_eos *eos, const char *filename,
    double *time, int *step, int *dump_count, double *next_output_time, double *next_restart_time);
void prj_io_write_dump(const prj_mesh *mesh, const char *basename, int dump_index, int step, double time);
void prj_io_parser(prj_sim *sim, char *filename);

#endif
