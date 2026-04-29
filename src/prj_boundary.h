#ifndef PRJ_BOUNDARY_H
#define PRJ_BOUNDARY_H

void prj_boundary_get_prim(const prj_block *block, int stage, double x1, double x2, double x3, double *w);
void prj_boundary_get_eosvar(const prj_block *block, double x1, double x2, double x3, double *eosv);
void prj_boundary_send(prj_block *block, int stage, int fill_kind);
void prj_boundary_physical(const prj_bc *bc, prj_block *block, int stage, int mode);
void prj_boundary_mpi_recv(prj_mesh *mesh, int stage, int fill_kind);
void prj_boundary_fill_ghosts(prj_mesh *mesh, const prj_bc *bc, int stage);
#if PRJ_MHD
void prj_boundary_send_bf(prj_block *block, int use_bf1, int fill_kind);
void prj_boundary_fill_bf(prj_mesh *mesh, const prj_bc *bc, int use_bf1);
#endif

#endif
