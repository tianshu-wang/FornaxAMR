#ifndef PRJ_MHD_H
#define PRJ_MHD_H

#include "prj_types.h"

typedef struct prj_mhd_face_patch {
    double value[3][3][2][2];
    unsigned char fixed[3][3][2][2];
} prj_mhd_face_patch;

void prj_mhd_face_patch_clear(prj_mhd_face_patch *patch);
void prj_mhd_bf_prolong_direct(const prj_block *coarse, const double *Bface[3],
    int i, int j, int k, prj_mhd_face_patch *patch);
double *prj_mhd_bface_stage(prj_block *block, int stage, int dir);
const double *prj_mhd_bface_stage_const(const prj_block *block, int stage, int dir);
int prj_mhd_face_is_interior(int dir, int i, int j, int k);
int prj_mhd_face_is_owned(int dir, int i, int j, int k);
void prj_mhd_face_position(const prj_block *block, int dir, int i, int j, int k, double x[3]);
int prj_mhd_face_index_from_position(const prj_block *block, int dir, const double x[3], int *i, int *j, int *k);
int prj_mhd_block_owns_face_position(const prj_block *block, int dir, const double x[3], int *i, int *j, int *k);
int prj_mhd_edge_is_interior(int dir, int i, int j, int k);
int prj_mhd_edge_is_owned(int dir, int i, int j, int k);
void prj_mhd_edge_position(const prj_block *block, int dir, int i, int j, int k, double x[3]);
int prj_mhd_edge_index_from_position(const prj_block *block, int dir, const double x[3], int *i, int *j, int *k);
int prj_mhd_block_owns_edge_position(const prj_block *block, int dir, const double x[3], int *i, int *j, int *k);
void prj_mhd_apply_bf_patch(prj_mesh *mesh, prj_block *fine, int stage,
    const prj_block *coarse, const double *Bface[3], int i, int j, int k);
void prj_mhd_init(prj_sim *sim);
void prj_mhd_bf2bc(prj_block *block, int stage);
double prj_mhd_emf_upwind(const double emf_face[4], const double emf_cell[4], const double face_velocity[4]);
void prj_mhd_emf_send(prj_mesh *mesh);
void prj_mhd_debug_check_emf(const prj_mesh *mesh);
void prj_mhd_debug_check_divergence(const prj_mesh *mesh, int stage);

#endif
