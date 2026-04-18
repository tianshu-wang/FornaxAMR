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
    int *eos_done;
    double *U;
    double *dUdt;
    double *flux[3];
    double *v_riemann[3];
    int *ridx;
    double *fr;
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
    double min_dx;
    int max_active_level;
    int root_nx[3];
    prj_block *blocks;
    prj_coord coord;
    double amr_refine_thresh[PRJ_AMR_N];
    double amr_derefine_thresh[PRJ_AMR_N];
    double amr_lohner_eps[PRJ_AMR_N];
    int amr_estimator[PRJ_AMR_N];
    int amr_lohner_var[PRJ_AMR_N];
    int amr_criterion_set[PRJ_AMR_N];
    int use_amr_angle_resolution;
    double amr_angle_resolution_limit;
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
    double rho_min;
    double rho_max;
    double temp_min;
    double temp_max;
    double inv_dlogrho;
    double inv_dlogT;
    double inv_dYe;
    double ln10_t1;
    double ln10_dlogT;
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
    double *uavg_int;
    double *erad_avg;
    double *prad_avg;
    double *vdotF_avg;
};

struct prj_rad {
    char filename[100];
    int maxiter;
    double implicit_err_tol;
#if PRJ_NRAD > 0
    double emin[PRJ_NRAD];
    double emax[PRJ_NRAD];
    double *egroup[PRJ_NRAD];
    double *eedge[PRJ_NRAD];
    double *egroup_erg[PRJ_NRAD];
    double *degroup_erg[PRJ_NRAD];
    double *x_e[PRJ_NRAD];

    char table_param_file[100];
    char table_file[100];
    int nromax;
    int ntmax;
    int nyemax;
    int ngmax;
    double romin;
    double romax;
    double tmin;
    double tmax;
    double yemin;
    double yemax;
    double *enuk;
    double *absopac[PRJ_NRAD];
    double *scaopac[PRJ_NRAD];
    double *emis[PRJ_NRAD];
    double *sdelta[PRJ_NRAD];
#endif
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
    double dt_factor;
    double t_end;
    double output_dt;
    double restart_dt;
    int max_steps;
    int output_interval;
    int restart_interval;
    char output_dir[256];
    char progenitor_file[256];
    char problem_name[64];
    int restart_from_file;
    char restart_file_name[256];
    int amr_interval;
    int dump_count;
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
