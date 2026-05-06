#ifndef PRJ_MHD_H
#define PRJ_MHD_H

#include "prj_types.h"

#if PRJ_MHD
void prj_mhd_init(prj_sim *sim);
void prj_mhd_bf2bc(prj_eos *eos, prj_block *block, int use_bf1);
void prj_mhd_prolong_bf_from_buffer(const double *buf[3],
    const int buf_lo[3][3], const int buf_n[3][3], const double coarse_dx[3],
    prj_block *fine, int ci, int cj, int ck, int fi, int fj, int fk,
    int use_bf1);
double prj_mhd_emf_upwind(prj_block *block, int dir, int i, int j, int k,
    const double emf_face[4], const double emf_cell[4], const double v_norm[4]);
void prj_mhd_emf_send(prj_mesh *mesh);
void prj_mhd_debug_check_emf(const prj_mesh *mesh);
void prj_mhd_debug_check_divb(const prj_mesh *mesh, int use_bf1);
#else
static inline void prj_mhd_init(prj_sim *sim)
{
    (void)sim;
}

static inline void prj_mhd_bf2bc(prj_eos *eos, prj_block *block, int use_bf1)
{
    (void)eos;
    (void)block;
    (void)use_bf1;
}

static inline void prj_mhd_prolong_bf_from_buffer(const double *buf[3],
    const int buf_lo[3][3], const int buf_n[3][3], const double coarse_dx[3],
    prj_block *fine, int ci, int cj, int ck, int fi, int fj, int fk,
    int use_bf1)
{
    (void)buf;
    (void)buf_lo;
    (void)buf_n;
    (void)coarse_dx;
    (void)fine;
    (void)ci;
    (void)cj;
    (void)ck;
    (void)fi;
    (void)fj;
    (void)fk;
    (void)use_bf1;
}

static inline double prj_mhd_emf_upwind(prj_block *block, int dir, int i, int j, int k,
    const double emf_face[4], const double emf_cell[4], const double v_norm[4])
{
    (void)block;
    (void)dir;
    (void)i;
    (void)j;
    (void)k;
    (void)emf_face;
    (void)emf_cell;
    (void)v_norm;
    return 0.0;
}

static inline void prj_mhd_emf_send(prj_mesh *mesh)
{
    (void)mesh;
}

static inline void prj_mhd_debug_check_emf(const prj_mesh *mesh)
{
    (void)mesh;
}

static inline void prj_mhd_debug_check_divb(const prj_mesh *mesh, int use_bf1)
{
    (void)mesh;
    (void)use_bf1;
}
#endif

#endif
