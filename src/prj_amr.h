#ifndef PRJ_AMR_H
#define PRJ_AMR_H

void prj_amr_tag(prj_mesh *mesh, prj_eos *eos);
int prj_amr_criteria_need_eosvar(const prj_mesh *mesh);
int prj_amr_criteria_need_gravity(const prj_mesh *mesh);
void prj_amr_enforce_two_to_one(prj_mesh *mesh);
int prj_amr_refine_marked_blocks(prj_mesh *mesh);
void prj_amr_refine_block(prj_mesh *mesh, int block_id);
int prj_amr_coarsen_block(prj_mesh *mesh, int parent_id);
void prj_amr_prolongate(const prj_block *parent, prj_block *child, int child_oct, double E_floor);
void prj_amr_restrict(const prj_block *children[8], prj_block *parent);
void prj_amr_adapt(prj_mesh *mesh, prj_eos *eos);
void prj_amr_init_neighbors(prj_mesh *mesh);

#endif
