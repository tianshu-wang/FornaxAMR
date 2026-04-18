#ifndef PRJ_MESH_H
#define PRJ_MESH_H

#include "prj_types.h"

int prj_mesh_init(prj_mesh *mesh, int root_nx1, int root_nx2, int root_nx3, int max_level, const prj_coord *coord);
int prj_block_alloc_data(prj_block *b);
void prj_block_free_data(prj_block *b);
void prj_block_setup_geometry(prj_block *b, const prj_coord *coord);
int prj_block_cache_index(int i, int j, int k);
prj_block *prj_mesh_get_block(prj_mesh *mesh, int id);
int prj_mesh_count_active(const prj_mesh *mesh);
double prj_mesh_min_cell_size(const prj_mesh *mesh);
void prj_mesh_update_max_active_level(prj_mesh *mesh);
void prj_mesh_mark_base_blocks(prj_mesh *mesh);
void prj_mesh_destroy(prj_mesh *mesh);

#endif
