#ifndef PRJ_AMR_H
#define PRJ_AMR_H

void prj_amr_tag(prj_mesh *mesh, prj_eos *eos, const prj_mpi *mpi);
int prj_amr_criteria_need_eosvar(const prj_mesh *mesh);
void prj_amr_enforce_two_to_one(prj_mesh *mesh, const prj_mpi *mpi);
int prj_amr_refine_marked_blocks(prj_mesh *mesh, const prj_mpi *mpi, const prj_grav *grav);
void prj_amr_refine_block(prj_mesh *mesh, const prj_mpi *mpi, int block_id, const prj_grav *grav);
int prj_amr_coarsen_block(prj_mesh *mesh, const prj_mpi *mpi, int parent_id);
void prj_amr_prolongate(const prj_mesh *mesh, const prj_mpi *mpi, const prj_block *parent,
    prj_block *child, int child_oct);
void prj_amr_restrict(const prj_block *children[8], prj_block *parent);
int prj_amr_adapt(prj_mesh *mesh, prj_eos *eos, prj_mpi *mpi, const prj_grav *grav);
void prj_amr_init_neighbors(prj_mesh *mesh);
double prj_amr_angular_resolution_limit(double x1, double x2, double x3, double r_com);

#endif
