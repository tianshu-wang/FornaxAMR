#ifndef PRJ_MHD_H
#define PRJ_MHD_H

#include "prj_types.h"

#if PRJ_MHD
void prj_mhd_init(prj_sim *sim);
void prj_mhd_bf2bc(prj_eos *eos, prj_block *block, int use_bf1);
void prj_mhd_bf_prolongate(const prj_block *coarse, prj_block *fine,
    int ci, int cj, int ck, int fi, int fj, int fk, int use_bf1);
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

static inline void prj_mhd_bf_prolongate(const prj_block *coarse, prj_block *fine,
    int ci, int cj, int ck, int fi, int fj, int fk, int use_bf1)
{
    (void)coarse;
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
