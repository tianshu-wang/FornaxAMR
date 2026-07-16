#ifndef PRJ_Z4C_H
#define PRJ_Z4C_H

#include "prj_types.h"

const char *prj_z4c_var_name(int var);
const char *prj_z4c_tmunu_var_name(int var);
int prj_z4c_runtime_enabled(const prj_mesh *mesh);
void prj_z4c_init_params(prj_z4c_params *params);

typedef struct prj_z4c_hydro_geom {
    double gamma[3][3];
    double gamma_inv[3][3];
    double dgamma[3][3][3];
    double alpha;
    double beta[3];
    double dalpha[3];
    double dbeta[3][3];
    double sqrt_gamma;
    double K_dd[3][3];
} prj_z4c_hydro_geom;

int prj_z4c_load_hydro_geom(const prj_mesh *mesh, const prj_block *block,
    int stage, int i, int j, int k, prj_z4c_hydro_geom *geom);
int prj_z4c_load_hydro_metric_geom(const prj_mesh *mesh, const prj_block *block,
    int stage, int i, int j, int k, prj_z4c_hydro_geom *geom);
/* Lightweight sqrt(gamma) only: same value load_hydro_geom returns in
 * geom.sqrt_gamma, but without the metric inverse, lapse/shift, or extrinsic
 * curvature. Returns 1 on success. */
int prj_z4c_cell_sqrt_gamma(const prj_mesh *mesh, const prj_block *block,
    int stage, int i, int j, int k, double *sqrt_gamma_out);

double prj_z4c_calc_dt_seconds(const prj_mesh *mesh, const prj_mpi *mpi, double cfl);
void prj_z4c_init_mesh_flat(prj_mesh *mesh, const prj_mpi *mpi);
void prj_z4c_init_mesh_spherical(prj_mesh *mesh, const prj_mpi *mpi,
    const double x_com[3], const double *r_grid, const double *chi_grid,
    const double *alpha_grid, int npts);
void prj_z4c_fill_ghosts(prj_mesh *mesh, prj_mpi *mpi, const prj_bc *bc, int stage);
void prj_z4c_init_punctures(prj_mesh *mesh, const prj_mpi *mpi, int npunctures,
    const double centers_cm[][3], const double masses_cm[],
    const double momenta_cm[][3], double floor_radius_cm);

void prj_z4c_compute_rhs(prj_mesh *mesh, const prj_mpi *mpi,
    const prj_rad *rad, int state_stage, int rhs_stage, double tau_cm);
void prj_z4c_apply_sommerfeld_rhs(prj_mesh *mesh, const prj_mpi *mpi,
    const prj_bc *bc, int state_stage, int rhs_stage);
void prj_z4c_update_linear(prj_mesh *mesh, const prj_mpi *mpi,
    int dst_stage, int a_stage, double a_w, int b_stage, double b_w,
    int rhs_stage, double dtau_cm);
/* Pointwise single-cell analog of prj_z4c_update_linear: advances one cell's
 * Z4c state in place (used to fuse the stage geometry advance into the per-cell
 * hydro update, between densitization and cons2prim recovery). Callers guard on
 * the dynamic-GR predicate. */
void prj_z4c_update_linear_cell(prj_block *block,
    int dst_stage, int a_stage, double a_w, int b_stage, double b_w,
    int rhs_stage, double dtau_cm, int i, int j, int k);
void prj_z4c_finalize_stage(prj_mesh *mesh, prj_mpi *mpi, const prj_bc *bc, int stage);

void prj_z4c_save_stage(prj_mesh *mesh, const prj_mpi *mpi, int dst_stage, int src_stage);
void prj_z4c_blend_with_saved(prj_mesh *mesh, prj_mpi *mpi, const prj_bc *bc,
    int saved_stage, double saved_weight);

void prj_z4c_amr_prolongate_child(const prj_block *parent, prj_block *child, int child_oct);
void prj_z4c_amr_restrict_parent(const prj_block *children[8], prj_block *parent);

#endif
