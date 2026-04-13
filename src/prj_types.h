#ifndef PRJ_TYPES_H
#define PRJ_TYPES_H

#include <stddef.h>

#include "prj_defs.h"

typedef struct prj_block prj_block;
typedef struct prj_mesh prj_mesh;
typedef struct prj_sim prj_sim;
typedef struct prj_eos prj_eos;
typedef struct prj_coord prj_coord;
typedef struct prj_neighbor prj_neighbor;
typedef struct prj_bc prj_bc;
typedef struct prj_grav_mono prj_grav_mono;
typedef struct prj_rad prj_rad;
typedef struct prj_mpi_buffer prj_mpi_buffer;
typedef struct prj_mpi prj_mpi;
typedef void (*prj_problem_init_fn)(prj_sim *sim);

struct prj_coord {
    double x1min;
    double x1max;
    double x2min;
    double x2max;
    double x3min;
    double x3max;
};

struct prj_neighbor {
    int id;
    int rank;
    double xmin[3];
    double xmax[3];
    double dx[3];
};

struct prj_bc {
    int bc_x1_inner;
    int bc_x1_outer;
    int bc_x2_inner;
    int bc_x2_outer;
    int bc_x3_inner;
    int bc_x3_outer;
};

struct prj_block {
    int id;
    int rank;
    int level;
    int active;
    int refine_flag;
    int base_block;
    double xmin[3];
    double xmax[3];
    double dx[3];
    double *W;
    double *W1;
    double *eosvar;
    double *U;
    double *dUdt;
    double *flux[3];
    double vol;
    double area[3];
    int parent;
    int children[8];
    prj_neighbor slot[56];
};

struct prj_mesh {
    int nblocks;
    int nblocks_max;
    int max_level;
    int root_nx[3];
    prj_block *blocks;
    prj_coord coord;
    double amr_refine_thresh;
    double amr_derefine_thresh;
    double amr_eps;
    int amr_estimator;
    double E_floor;
};

struct prj_eos {
    char filename[100];
    int kind;
    int table_loaded;
    int table_is_mmap;
    int nt;
    int nr;
    int ny;
    double r1;
    double r2;
    double t1;
    double t2;
    double y1c;
    double y2c;
    double dlogrho;
    double dlogT;
    double dYe;
    size_t table_bytes;
    double *table;
};

struct prj_grav_mono {
    int nbins;
    double dr_min;
    double *rf;
    double *ms;
    double *phi;
    double *accel;
    double *lapse;
    double *vol;
    double *rho_avg;
    double *vr_avg;
    double *pgas_avg;
    double *eint_avg;
    double *erad_avg;
    double *prad_avg;
    double *vdotF_avg;
};

struct prj_rad {
    char filename[100];
};

struct prj_sim {
    prj_mesh mesh;
    prj_eos eos;
    prj_rad rad;
    prj_grav_mono monopole_grav;
    prj_coord coord;
    prj_bc bc;
    double time;
    double dt;
    int step;
    double cfl;
    double t_end;
    int max_steps;
    int output_interval;
    int restart_interval;
    char output_dir[256];
    char progenitor_file[256];
    int amr_interval;
};

struct prj_mpi_buffer {
    int receiver_rank;
    int number;
    int *receiver_blocks;
    int *cell_data_size_send;
    int *cell_data_idx_send[3];
    double *cell_buffer_send;
    int *face_data_size_send;
    int *face_data_idx_send[3];
    double *face_buffer_send;
    int *cell_data_size_recv;
    int *cell_data_idx_recv[3];
    double *cell_buffer_recv;
    int *face_data_size_recv;
    int *face_data_idx_recv[3];
    double *face_buffer_recv;
};

struct prj_mpi {
    int totrank;
    int rank;
    int neighbor_number;
    prj_mpi_buffer *neighbor_buffer;
};

#endif
