#include <stddef.h>
#include <stdlib.h>

#if defined(PRJ_ENABLE_MPI)
#include <mpi.h>
#endif

#include "prj.h"

static int prj_neighbor_slot_index(int ox, int oy, int oz)
{
    int index;

    if (ox < -1 || ox > 1 || oy < -1 || oy > 1 || oz < -1 || oz > 1) {
        return -1;
    }
    if (ox == 0 && oy == 0 && oz == 0) {
        return -1;
    }

    index = (ox + 1) * 9 + (oy + 1) * 3 + (oz + 1);
    if (index > 13) {
        index -= 1;
    }
    return index;
}

static void prj_block_init_empty(prj_block *b)
{
    int n;

    b->id = -1;
    b->rank = 0;
    b->level = 0;
    b->active = 1;
    b->refine_flag = 0;
    b->base_block = 0;
    b->W = 0;
    b->W1 = 0;
    b->eosvar = 0;
    b->eos_done = 0;
    b->U = 0;
    b->dUdt = 0;
    b->flux[0] = 0;
    b->flux[1] = 0;
    b->flux[2] = 0;
    b->v_riemann[0] = 0;
    b->v_riemann[1] = 0;
    b->v_riemann[2] = 0;
    b->ridx = 0;
    b->fr = 0;
    b->vol = 0.0;
    b->area[0] = 0.0;
    b->area[1] = 0.0;
    b->area[2] = 0.0;
    b->parent = -1;

    for (n = 0; n < 8; ++n) {
        b->children[n] = -1;
    }
    for (n = 0; n < 56; ++n) {
        b->slot[n].id = -1;
        b->slot[n].rank = 0;
        b->slot[n].xmin[0] = 0.0;
        b->slot[n].xmin[1] = 0.0;
        b->slot[n].xmin[2] = 0.0;
        b->slot[n].xmax[0] = 0.0;
        b->slot[n].xmax[1] = 0.0;
        b->slot[n].xmax[2] = 0.0;
        b->slot[n].dx[0] = 0.0;
        b->slot[n].dx[1] = 0.0;
        b->slot[n].dx[2] = 0.0;
    }
}

int prj_block_alloc_data(prj_block *b)
{
    size_t prim_count;
    size_t eosvar_count;
    size_t cons_count;
    size_t total_count;
    double *base;
    int *eos_done;
    int *ridx;
    double *fr;

    if (b == 0) {
        return 1;
    }

    prj_block_free_data(b);

    prim_count = (size_t)PRJ_NVAR_PRIM * (size_t)PRJ_BLOCK_NCELLS;
    eosvar_count = (size_t)PRJ_NVAR_EOSVAR * (size_t)PRJ_BLOCK_NCELLS;
    cons_count = (size_t)PRJ_NVAR_CONS * (size_t)PRJ_BLOCK_NCELLS;
    total_count = 2U * prim_count + eosvar_count + 5U * cons_count + 9U * (size_t)PRJ_BLOCK_NCELLS;

    base = (double *)malloc(total_count * sizeof(*base));
    if (base == 0) {
        return 2;
    }
    eos_done = (int *)calloc((size_t)PRJ_BLOCK_NCELLS, sizeof(*eos_done));
    ridx = (int *)malloc((size_t)PRJ_BLOCK_NCELLS * sizeof(*ridx));
    fr = (double *)malloc((size_t)PRJ_BLOCK_NCELLS * sizeof(*fr));
    if (eos_done == 0 || ridx == 0 || fr == 0) {
        free(fr);
        free(ridx);
        free(eos_done);
        free(base);
        return 2;
    }

    b->W = base;
    base += prim_count;
    b->W1 = base;
    base += prim_count;
    b->eosvar = base;
    base += eosvar_count;
    b->eos_done = eos_done;
    b->U = base;
    base += cons_count;
    b->dUdt = base;
    base += cons_count;
    b->flux[0] = base;
    base += cons_count;
    b->flux[1] = base;
    base += cons_count;
    b->flux[2] = base;
    base += cons_count;
    b->v_riemann[0] = base;
    base += 3U * (size_t)PRJ_BLOCK_NCELLS;
    b->v_riemann[1] = base;
    base += 3U * (size_t)PRJ_BLOCK_NCELLS;
    b->v_riemann[2] = base;
    b->ridx = ridx;
    b->fr = fr;
    prj_gravity_cache_block(b);
    return 0;
}

void prj_block_free_data(prj_block *b)
{
    if (b == 0) {
        return;
    }

    free(b->W);
    free(b->eos_done);
    free(b->ridx);
    free(b->fr);
    b->W = 0;
    b->W1 = 0;
    b->eosvar = 0;
    b->eos_done = 0;
    b->U = 0;
    b->dUdt = 0;
    b->flux[0] = 0;
    b->flux[1] = 0;
    b->flux[2] = 0;
    b->v_riemann[0] = 0;
    b->v_riemann[1] = 0;
    b->v_riemann[2] = 0;
    b->ridx = 0;
    b->fr = 0;
}

void prj_block_setup_geometry(prj_block *b, const prj_coord *coord)
{
    (void)coord;

    if (b == 0) {
        return;
    }

    b->vol = b->dx[0] * b->dx[1] * b->dx[2];
    b->area[0] = b->dx[1] * b->dx[2];
    b->area[1] = b->dx[0] * b->dx[2];
    b->area[2] = b->dx[0] * b->dx[1];
    prj_gravity_cache_block(b);
}

int prj_block_cache_index(int i, int j, int k)
{
    return IDX(i, j, k);
}

prj_block *prj_mesh_get_block(prj_mesh *mesh, int id)
{
    if (mesh == 0 || id < 0 || id >= mesh->nblocks) {
        return 0;
    }

    return &mesh->blocks[id];
}

int prj_mesh_count_active(const prj_mesh *mesh)
{
    int i;
    int count;

    if (mesh == 0) {
        return 0;
    }

    count = 0;
    for (i = 0; i < mesh->nblocks; ++i) {
        if (mesh->blocks[i].active == 1) {
            count += 1;
        }
    }
    return count;
}

void prj_mesh_update_max_active_level(prj_mesh *mesh)
{
    int i;
    int local_max;

    if (mesh == 0) {
        return;
    }

    local_max = -1;
    for (i = 0; i < mesh->nblocks; ++i) {
        prj_block *b = &mesh->blocks[i];

        if (b->active == 0) {
            continue;
        }
        if (b->level > local_max) {
            local_max = b->level;
        }
    }

#if defined(PRJ_ENABLE_MPI)
    MPI_Allreduce(&local_max, &mesh->max_active_level, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#else
    mesh->max_active_level = local_max;
#endif
}

void prj_mesh_mark_base_blocks(prj_mesh *mesh)
{
    int i;

    if (mesh == 0) {
        return;
    }
    for (i = 0; i < mesh->nblocks; ++i) {
        prj_block *block = &mesh->blocks[i];

        block->base_block = (block->id >= 0 && block->active == 1) ? 1 : 0;
    }
}

int prj_mesh_init(prj_mesh *mesh, int root_nx1, int root_nx2, int root_nx3, int max_level, const prj_coord *coord)
{
    int i;
    int j;
    int k;
    int id;
    int nroot;
    int capacity;
    double block_dx[3];
    double saved_amr_refine_thresh[PRJ_AMR_N];
    double saved_amr_derefine_thresh[PRJ_AMR_N];
    int saved_amr_estimator[PRJ_AMR_N];
    int saved_amr_criterion_set[PRJ_AMR_N];
    double saved_amr_eps;
    int saved_use_amr_angle_resolution;
    double saved_amr_angle_resolution_limit;
    double saved_E_floor;
    int amr_idx;

    if (mesh == 0 || coord == 0 || root_nx1 <= 0 || root_nx2 <= 0 || root_nx3 <= 0 || max_level < 0) {
        return 1;
    }

    for (amr_idx = 0; amr_idx < PRJ_AMR_N; ++amr_idx) {
        saved_amr_refine_thresh[amr_idx] = mesh->amr_refine_thresh[amr_idx];
        saved_amr_derefine_thresh[amr_idx] = mesh->amr_derefine_thresh[amr_idx];
        saved_amr_estimator[amr_idx] = mesh->amr_estimator[amr_idx];
        saved_amr_criterion_set[amr_idx] = mesh->amr_criterion_set[amr_idx];
    }
    saved_amr_eps = mesh->amr_eps;
    saved_use_amr_angle_resolution = mesh->use_amr_angle_resolution;
    saved_amr_angle_resolution_limit = mesh->amr_angle_resolution_limit;
    saved_E_floor = mesh->E_floor;

    mesh->nblocks = 0;
    mesh->nblocks_max = 0;
    mesh->max_level = max_level;
    mesh->max_active_level = -1;
    mesh->root_nx[0] = root_nx1;
    mesh->root_nx[1] = root_nx2;
    mesh->root_nx[2] = root_nx3;
    mesh->coord = *coord;
    for (amr_idx = 0; amr_idx < PRJ_AMR_N; ++amr_idx) {
        mesh->amr_refine_thresh[amr_idx] = saved_amr_refine_thresh[amr_idx];
        mesh->amr_derefine_thresh[amr_idx] = saved_amr_derefine_thresh[amr_idx];
        mesh->amr_estimator[amr_idx] = saved_amr_estimator[amr_idx];
        mesh->amr_criterion_set[amr_idx] = saved_amr_criterion_set[amr_idx];
    }
    mesh->amr_eps = saved_amr_eps;
    mesh->use_amr_angle_resolution = saved_use_amr_angle_resolution;
    mesh->amr_angle_resolution_limit = saved_amr_angle_resolution_limit;
    mesh->E_floor = saved_E_floor;
    mesh->blocks = 0;

    nroot = root_nx1 * root_nx2 * root_nx3;
    capacity = 65536;

    mesh->blocks = (prj_block *)malloc((size_t)capacity * sizeof(*mesh->blocks));
    if (mesh->blocks == 0) {
        return 2;
    }
    mesh->nblocks_max = capacity;

    for (i = 0; i < capacity; ++i) {
        prj_block_init_empty(&mesh->blocks[i]);
    }

    block_dx[0] = (coord->x1max - coord->x1min) / (double)root_nx1;
    block_dx[1] = (coord->x2max - coord->x2min) / (double)root_nx2;
    block_dx[2] = (coord->x3max - coord->x3min) / (double)root_nx3;

    id = 0;
    for (i = 0; i < root_nx1; ++i) {
        for (j = 0; j < root_nx2; ++j) {
            for (k = 0; k < root_nx3; ++k) {
                prj_block *b = &mesh->blocks[id];

                b->id = id;
                b->rank = 0;
                b->level = 0;
                b->active = 1;
                b->refine_flag = 0;
                b->base_block = 0;
                b->xmin[0] = coord->x1min + (double)i * block_dx[0];
                b->xmax[0] = b->xmin[0] + block_dx[0];
                b->xmin[1] = coord->x2min + (double)j * block_dx[1];
                b->xmax[1] = b->xmin[1] + block_dx[1];
                b->xmin[2] = coord->x3min + (double)k * block_dx[2];
                b->xmax[2] = b->xmin[2] + block_dx[2];
                b->dx[0] = block_dx[0] / (double)PRJ_BLOCK_SIZE;
                b->dx[1] = block_dx[1] / (double)PRJ_BLOCK_SIZE;
                b->dx[2] = block_dx[2] / (double)PRJ_BLOCK_SIZE;
                prj_block_setup_geometry(b, coord);
                if (prj_block_alloc_data(b) != 0) {
                    prj_mesh_destroy(mesh);
                    return 3;
                }
                id += 1;
            }
        }
    }
    mesh->nblocks = nroot;
    prj_mesh_update_max_active_level(mesh);

    for (i = 0; i < root_nx1; ++i) {
        for (j = 0; j < root_nx2; ++j) {
            for (k = 0; k < root_nx3; ++k) {
                prj_block *b = &mesh->blocks[(i * root_nx2 + j) * root_nx3 + k];
                int ox;
                int oy;
                int oz;

                for (ox = -1; ox <= 1; ++ox) {
                    for (oy = -1; oy <= 1; ++oy) {
                        for (oz = -1; oz <= 1; ++oz) {
                            int ni = i + ox;
                            int nj = j + oy;
                            int nk = k + oz;
                            int slot_index = prj_neighbor_slot_index(ox, oy, oz);

                            if (slot_index < 0) {
                                continue;
                            }
                            if (ni < 0 || ni >= root_nx1 || nj < 0 || nj >= root_nx2 || nk < 0 || nk >= root_nx3) {
                                continue;
                            }

                            id = (ni * root_nx2 + nj) * root_nx3 + nk;
                            b->slot[slot_index].id = id;
                            b->slot[slot_index].rank = 0;
                            b->slot[slot_index].xmin[0] = mesh->blocks[id].xmin[0];
                            b->slot[slot_index].xmin[1] = mesh->blocks[id].xmin[1];
                            b->slot[slot_index].xmin[2] = mesh->blocks[id].xmin[2];
                            b->slot[slot_index].xmax[0] = mesh->blocks[id].xmax[0];
                            b->slot[slot_index].xmax[1] = mesh->blocks[id].xmax[1];
                            b->slot[slot_index].xmax[2] = mesh->blocks[id].xmax[2];
                            b->slot[slot_index].dx[0] = mesh->blocks[id].dx[0];
                            b->slot[slot_index].dx[1] = mesh->blocks[id].dx[1];
                            b->slot[slot_index].dx[2] = mesh->blocks[id].dx[2];
                        }
                    }
                }
            }
        }
    }

    return 0;
}

void prj_mesh_destroy(prj_mesh *mesh)
{
    int i;

    if (mesh == 0) {
        return;
    }

    for (i = 0; i < mesh->nblocks; ++i) {
        prj_block_free_data(&mesh->blocks[i]);
    }
    free(mesh->blocks);
    mesh->blocks = 0;
    mesh->nblocks = 0;
    mesh->nblocks_max = 0;
    mesh->max_active_level = -1;
}
