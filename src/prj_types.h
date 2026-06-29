#ifndef PRJ_TYPES_H
#define PRJ_TYPES_H

#include <stddef.h>

#include "prj_defs.h"

#ifndef LMAX
#define LMAX 12
#endif
#if LMAX < 1
#error "LMAX must be at least 1"
#endif

typedef struct prj_block prj_block;
typedef struct prj_mesh prj_mesh;
typedef struct prj_sim prj_sim;
typedef struct prj_eos prj_eos;
typedef struct prj_coord prj_coord;
typedef struct prj_neighbor prj_neighbor;
typedef struct prj_morton_lookup_entry prj_morton_lookup_entry;
typedef struct prj_bc prj_bc;
typedef struct prj_grav prj_grav;
typedef struct prj_rad prj_rad;
typedef struct prj_mpi_buffer prj_mpi_buffer;
typedef struct prj_mpi prj_mpi;
typedef void (*prj_problem_init_fn)(prj_sim *sim, prj_mpi *mpi);
typedef int (*prj_amr_init_refine_fn)(const prj_block *block, void *userdata);

struct prj_coord {
    double x1min;
    double x1max;
    double x2min;
    double x2max;
    double x3min;
    double x3max;
};

enum prj_neighbor_type {
    PRJ_NEIGHBOR_NONE = 0,
    PRJ_NEIGHBOR_FACE = 1,
    PRJ_NEIGHBOR_EDGE = 2,
    PRJ_NEIGHBOR_CORNER = 3
};

struct prj_neighbor {
    int id;
    int rank;
    int rel_level;
    int type;
    int send_loc_start[3];
    int send_loc_end[3];
    int recv_loc_start[3];
    int recv_loc_end[3];
    int send_loc_start_rad[3];
    int send_loc_end_rad[3];
    int recv_loc_start_rad[3];
    int recv_loc_end_rad[3];
    double xmin[3];
    double xmax[3];
    double dx[3];
};

struct prj_morton_lookup_entry {
    int occupied;
    int level;
    int coord[3];
    int id;
    unsigned long long key;
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
    int can_refine;
    double xmin[3];
    double xmax[3];
    double dx[3];
    double *W;
    double *eosvar;
    int *cell_derived_done;
    double *U;
    double *dUdt;
    double *flux[3];
    double *v_riemann[3];
    double *kappa_cell;
    double *sigma_cell;
    double *lapse;
    double *grav[3];
    double *r_com;
    /* Cell-major: Ylm[cell * (LMAX * LMAX) + yidx]. */
    double *Ylm;
#if PRJ_MHD
    int *face_fidelity[3];
    int *edge_fidelity[3];
    double *Bf[3];
    double *Bv1[3];
    double *Bv2[3];
    double *emf[3];
#endif
    int *ridx;
    double *fr;
    double vol;
    double area[3];
    int parent;
    int children[8];
    prj_neighbor slot[56];
};

static inline double *prj_block_prim_stage(prj_block *block, int stage)
{
    return block != 0 && stage >= 0 && stage < PRJ_BLOCK_NSTAGES && block->W != 0 ?
        PRJ_BLOCK_STAGE_W(block->W, stage) : 0;
}

static inline const double *prj_block_prim_stage_const(const prj_block *block, int stage)
{
    return block != 0 && stage >= 0 && stage < PRJ_BLOCK_NSTAGES && block->W != 0 ?
        PRJ_BLOCK_STAGE_W(block->W, stage) : 0;
}

#if PRJ_MHD
static inline double *prj_block_bf_stage(prj_block *block, int dir, int stage)
{
    return block != 0 && dir >= 0 && dir < 3 &&
        stage >= 0 && stage < PRJ_BLOCK_NSTAGES && block->Bf[dir] != 0 ?
        PRJ_BLOCK_STAGE_BF(block->Bf[dir], stage) : 0;
}

static inline const double *prj_block_bf_stage_const(const prj_block *block, int dir, int stage)
{
    return block != 0 && dir >= 0 && dir < 3 &&
        stage >= 0 && stage < PRJ_BLOCK_NSTAGES && block->Bf[dir] != 0 ?
        PRJ_BLOCK_STAGE_BF(block->Bf[dir], stage) : 0;
}
#endif

struct prj_mesh {
    int nblocks;
    int nblocks_max;
    int max_blocks;
    int max_level;
    double min_dx;
    double min_allowable_cell_size;
    int max_active_level;
    int root_nx[3];
    double x_com[3];
    prj_block *blocks;
    prj_morton_lookup_entry *morton_lookup;
    int morton_lookup_capacity;
    int morton_lookup_count;
    prj_coord coord;
    double amr_refine_thresh[PRJ_AMR_N];
    double amr_derefine_thresh[PRJ_AMR_N];
    double amr_lohner_eps[PRJ_AMR_N];
    int amr_estimator[PRJ_AMR_N];
    int amr_lohner_var[PRJ_AMR_N];
    int amr_fractional_jump_var[PRJ_AMR_N];
    int amr_criterion_set[PRJ_AMR_N];
    int use_amr_angular_resolution_limit;
    int use_BJ;
    double amr_init_scale_factor;
    double E_floor;
    prj_amr_init_refine_fn amr_init_refine_fn;
    void *amr_init_refine_userdata;
};

struct prj_eos {
    char filename[PRJ_PATH_MAX];
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
    double E_injected;
};

struct prj_grav {
    int use_multipole_gravity;
    int nbins;
    double dr_min;
    double rmax;
    double min_cell;
    double *rf;
    double *ms;
    double *phi;
    double *lapse;
    /* Radial-bin-major: {C,D}lm[bin * (LMAX * LMAX) + yidx]. */
    double *Clm;
    double *Dlm;
    double *vol;
    double *rho_avg;
    double *vr_avg;
    double *pgas_avg;
    double *uavg_int;
    double *erad_avg;
    double *prad_avg;
    double *vdotF_avg;
    double *reduce_avg_buf;
    double *reduce_lm_buf;
};

struct prj_rad {
    char filename[100];
    int maxiter;
    double implicit_err_tol;
#if PRJ_NRAD > 0
    double chi[NCLOSURE + 1];
    double q[NCLOSURE + 1];
    double emin[PRJ_NRAD];
    double emax[PRJ_NRAD];
    double *egroup[PRJ_NRAD];
    double *eedge[PRJ_NRAD];
    double *egroup_erg[PRJ_NRAD];
    double *degroup_erg[PRJ_NRAD];
    double *x_e[PRJ_NRAD];

    char table_param_file[PRJ_PATH_MAX];
    char table_file[PRJ_PATH_MAX];
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
    double *log_enuk;
    double *log_egroup[PRJ_NRAD];
    double log_romin;
    double log_romax;
    double log_tmin;
    double log_tmax;
    double inv_logrho_span;
    double inv_logtemp_span;
    double inv_ye_span;
    double *absopac[PRJ_NRAD];
    double *scaopac[PRJ_NRAD];
    double *emis[PRJ_NRAD];
    double *sdelta[PRJ_NRAD];
    char eleinel_table_dir[PRJ_PATH_MAX];
    int eleinel_table_loaded;
    double *eleinel_phi_ee[PRJ_NRAD];

    double *spec_factor[PRJ_NRAD];
    double kom_epsilon;
    double kom_delta;
    double kom_dtmin;
    double kom_rhocut;
    double kom_nucinel_const;
    double kom_nucinel_prot;
    double kom_nucinel_neut;
    double kom_Ecut[PRJ_NRAD];
    double min_inel_density;
    double eleinel_factf;
    double eleinel_constin;
    double eleinel_constout;
    double eleinel_freqe3[PRJ_NRAD * PRJ_NEGROUP];
    double eleinel_factf_over_freqe3[PRJ_NRAD * PRJ_NEGROUP];
    double eleinel_freqe2_dnue[PRJ_NRAD * PRJ_NEGROUP];
    double expe[PRJ_NRAD][PRJ_NEGROUP * PRJ_NEGROUP * PRJ_EXPE_NT];
#endif
};

struct prj_sim {
    prj_mesh mesh;
    prj_eos eos;
    prj_rad rad;
    prj_grav grav;
    prj_coord coord;
    prj_bc bc;
    double time;
    double dt;
    int step;
    double cfl;
    double dt_factor;
    double x_com_err_tol;
    double t_end;
    double output_dt;
    double restart_dt;
    int max_steps;
    int output_interval;
    int restart_interval;
    int timer_interval;
    char progenitor_file[256];
    char problem_name[64];
    int restart_from_file;
    int restart_from_latest;
    char restart_file_name[256];
    int amr_interval;
    int dump_count;
    double perturbation_gaussian_norm;
    unsigned long long perturbation_seed;
#if PRJ_MHD
    int mhd_init_type;
    double mhd_B_norm;
    double mhd_B_scale;
#endif
};

struct prj_mpi_buffer {
    int receiver_rank;
    int number;
    int *receiver_blocks;
    int *cell_data_size_send;
    int *cell_data_idx_send[3];
    double *cell_buffer_send;
    int cell_send_count_by_kind[6];
    int cell_send_count_rad_by_kind[6];
    double *cell_buffer_send_by_kind[6];
#if PRJ_NRAD > 0 && PRJ_MIXED_PRECISION
    float *cell_buffer_send_rad_by_kind[6];
#endif
    int *face_data_size_send;
    int *face_data_idx_send[3];
    double *face_buffer_send;
    int *cell_data_size_recv;
    int *cell_data_idx_recv[3];
    double *cell_buffer_recv;
    int cell_recv_count_by_kind[6];
    int cell_recv_count_rad_by_kind[6];
    double *cell_buffer_recv_by_kind[6];
#if PRJ_NRAD > 0 && PRJ_MIXED_PRECISION
    float *cell_buffer_recv_rad_by_kind[6];
#endif
    int *face_data_size_recv;
    int *face_data_idx_recv[3];
    double *face_buffer_recv;
    /* Radiation stream covers only cells in the narrower radiation ghost band. */
    int *cell_data_size_send_rad;
    int *cell_data_idx_send_rad[3];
    int *cell_data_size_recv_rad;
    int *cell_data_idx_recv_rad[3];
    int flux_send_count;
    int flux_recv_count;
    int *flux_idx_send;
    double *flux_value_send;
    int *flux_idx_recv;
    double *flux_value_recv;
#if PRJ_MHD
    int bf_send_record_count[6];
    int bf_send_value_count[6];
    int bf_recv_record_count[6];
    int bf_recv_value_count[6];
    int bf_send_record_capacity;
    int bf_send_value_capacity;
    int bf_recv_record_capacity;
    int bf_recv_value_capacity;
    int *bf_headers_send;
    double *bf_values_send;
    int *bf_headers_recv;
    double *bf_values_recv;
    int emf_send_count;
    int emf_recv_count;
    int *emf_idx_send[3];
    int *emf_idx_recv[3];
    double *emf_value_send;
    double *emf_value_recv;
    int *emf_src_block;
    int *emf_src_dir;
    int *emf_src_idx[3];
    int amr_bf_send_value_count;
    int amr_bf_recv_record_count;
    int amr_bf_recv_value_count;
    int amr_bf_send_value_capacity;
    int amr_bf_recv_record_capacity;
    int amr_bf_recv_value_capacity;
    double *amr_bf_values_send;
    int *amr_bf_headers_recv;
    double *amr_bf_values_recv;
#endif
};

struct prj_mpi {
    int totrank;
    int rank;
    int neighbor_number;
    prj_mpi_buffer *neighbor_buffer;
    void *request_buffer;
    int request_capacity;
    int ghost_request_count;
#if PRJ_MHD
    int *amr_bf_headers;
    double *amr_bf_values;
    int amr_bf_record_count;
    int amr_bf_value_count;
    int amr_bf_record_capacity;
    int amr_bf_value_capacity;
#endif
};

#endif
