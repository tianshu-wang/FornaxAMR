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

static double prj_neighbor_abs(double x)
{
    return x < 0.0 ? -x : x;
}

void prj_neighbor_compute_geometry(const prj_block *a, const prj_block *b, prj_neighbor *slot)
{
    double tol=0;
    int axisrel[3];
    int touching;
    int d;

    if (a == 0 || b == 0 || slot == 0) {
        return;
    }

    touching = 0;
    for (d = 0; d < 3; ++d) {
        tol = 1.0e-2*PRJ_MIN(a->dx[d],b->dx[d]);
        if (prj_neighbor_abs(a->xmax[d] - b->xmin[d]) < tol) {
            axisrel[d] = 1;
            touching += 1;
        } else if (prj_neighbor_abs(b->xmax[d] - a->xmin[d]) < tol) {
            axisrel[d] = -1;
            touching += 1;
        } else {
            axisrel[d] = 0;
        }
    }

    if (touching == 1) {
        slot->type = PRJ_NEIGHBOR_FACE;
    } else if (touching == 2) {
        slot->type = PRJ_NEIGHBOR_EDGE;
    } else if (touching == 3) {
        slot->type = PRJ_NEIGHBOR_CORNER;
    } else {
        slot->type = PRJ_NEIGHBOR_NONE;
    }

    slot->rel_level = b->level - a->level;

    for (d = 0; d < 3; ++d) {
        int ri_lo = -PRJ_NGHOST;
        int ri_hi = PRJ_BLOCK_SIZE + PRJ_NGHOST;
        int si_lo = -PRJ_NGHOST;
        int si_hi = PRJ_BLOCK_SIZE + PRJ_NGHOST;
        tol = 1.0e-2*PRJ_MIN(a->dx[d],b->dx[d]);

        while (ri_lo < ri_hi &&
               b->xmin[d] + ((double)ri_lo + 0.5) * b->dx[d] < a->xmin[d] - tol)
            ri_lo++;
        while (ri_hi > ri_lo &&
               b->xmin[d] + ((double)ri_hi - 0.5) * b->dx[d] >= a->xmax[d] + tol)
            ri_hi--;

        slot->recv_loc_start[d] = ri_lo;
        slot->recv_loc_end[d] = ri_hi;

        if (ri_hi > ri_lo) {
            double rpos_lo = b->xmin[d] + ((double)ri_lo + 0.5) * b->dx[d];
            double rpos_hi = b->xmin[d] + ((double)ri_hi - 0.5) * b->dx[d];
            while (si_lo < si_hi &&
                   a->xmin[d] + ((double)si_lo + 0.5) * a->dx[d] < rpos_lo - tol)
                si_lo++;
            while (si_hi > si_lo &&
                   a->xmin[d] + ((double)si_hi - 0.5) * a->dx[d] > rpos_hi + tol)
                si_hi--;
        } else {
            si_lo = 0;
            si_hi = 0;
        }
        slot->send_loc_start[d] = si_lo;
        slot->send_loc_end[d] = si_hi;
        
    }
}

static void prj_neighbor_clear_derived(prj_neighbor *slot)
{
    int d;

    if (slot == 0) {
        return;
    }
    slot->rel_level = 0;
    slot->type = PRJ_NEIGHBOR_NONE;
    for (d = 0; d < 3; ++d) {
        slot->send_loc_start[d] = 0;
        slot->send_loc_end[d] = 0;
        slot->recv_loc_start[d] = 0;
        slot->recv_loc_end[d] = 0;
    }
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
#if PRJ_MHD
    for (n = 0; n < 3; ++n) {
        b->face_fidelity[n] = 0;
        b->edge_fidelity[n] = 0;
        b->Bf[n] = 0;
        b->Bf1[n] = 0;
        b->Bv1[n] = 0;
        b->Bv2[n] = 0;
        b->emf[n] = 0;
    }
#endif
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
        prj_neighbor_clear_derived(&b->slot[n]);
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
#if PRJ_MHD
    int *face_fidelity[3] = {0, 0, 0};
    int *edge_fidelity[3] = {0, 0, 0};
#endif
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
#if PRJ_MHD
    total_count += 15U * (size_t)PRJ_BLOCK_NCELLS;
#endif

    base = (double *)malloc(total_count * sizeof(*base));
    if (base == 0) {
        return 2;
    }
    eos_done = (int *)calloc((size_t)PRJ_BLOCK_NCELLS, sizeof(*eos_done));
#if PRJ_MHD
    for (int d = 0; d < 3; ++d) {
        face_fidelity[d] = (int *)calloc((size_t)PRJ_BLOCK_NCELLS, sizeof(*face_fidelity[d]));
        edge_fidelity[d] = (int *)calloc((size_t)PRJ_BLOCK_NCELLS, sizeof(*edge_fidelity[d]));
    }
#endif
    ridx = (int *)malloc((size_t)PRJ_BLOCK_NCELLS * sizeof(*ridx));
    fr = (double *)malloc((size_t)PRJ_BLOCK_NCELLS * sizeof(*fr));
    if (eos_done == 0 ||
#if PRJ_MHD
        face_fidelity[0] == 0 || face_fidelity[1] == 0 || face_fidelity[2] == 0 ||
        edge_fidelity[0] == 0 || edge_fidelity[1] == 0 || edge_fidelity[2] == 0 ||
#endif
        ridx == 0 || fr == 0) {
        free(fr);
        free(ridx);
#if PRJ_MHD
        for (int d = 0; d < 3; ++d) {
            free(edge_fidelity[d]);
            free(face_fidelity[d]);
        }
#endif
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
#if PRJ_MHD
    base += 3U * (size_t)PRJ_BLOCK_NCELLS;
    for (int d = 0; d < 3; ++d) {
        b->face_fidelity[d] = face_fidelity[d];
        b->edge_fidelity[d] = edge_fidelity[d];
        b->Bf[d] = base;
        base += (size_t)PRJ_BLOCK_NCELLS;
    }
    for (int d = 0; d < 3; ++d) {
        b->Bf1[d] = base;
        base += (size_t)PRJ_BLOCK_NCELLS;
    }
    for (int d = 0; d < 3; ++d) {
        b->Bv1[d] = base;
        base += (size_t)PRJ_BLOCK_NCELLS;
    }
    for (int d = 0; d < 3; ++d) {
        b->Bv2[d] = base;
        base += (size_t)PRJ_BLOCK_NCELLS;
    }
    for (int d = 0; d < 3; ++d) {
        b->emf[d] = base;
        base += (size_t)PRJ_BLOCK_NCELLS;
    }
#endif
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
#if PRJ_MHD
    for (int d = 0; d < 3; ++d) {
        free(b->face_fidelity[d]);
        free(b->edge_fidelity[d]);
    }
#endif
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
#if PRJ_MHD
    for (int d = 0; d < 3; ++d) {
        b->face_fidelity[d] = 0;
        b->edge_fidelity[d] = 0;
        b->Bf[d] = 0;
        b->Bf1[d] = 0;
        b->Bv1[d] = 0;
        b->Bv2[d] = 0;
        b->emf[d] = 0;
    }
#endif
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

static double prj_mesh_block_cell_size(const prj_block *block)
{
    double cell_size;

    if (block == 0) {
        return 0.0;
    }
    cell_size = block->dx[0];
    if (block->dx[1] > cell_size) {
        cell_size = block->dx[1];
    }
    if (block->dx[2] > cell_size) {
        cell_size = block->dx[2];
    }
    return cell_size;
}

double prj_mesh_min_cell_size(const prj_mesh *mesh)
{
    double min_cell_size = 1.0e99;
    int i;
    int found = 0;

    if (mesh == 0) {
        return 0.0;
    }
    for (i = 0; i < mesh->nblocks; ++i) {
        if (mesh->blocks[i].active == 1) {
            double cell_size = prj_mesh_block_cell_size(&mesh->blocks[i]);

            if (cell_size < min_cell_size) {
                min_cell_size = cell_size;
            }
            found = 1;
        }
    }
    return found != 0 ? min_cell_size : 0.0;
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
    int saved_amr_lohner_var[PRJ_AMR_N];
    int saved_amr_criterion_set[PRJ_AMR_N];
    double saved_amr_lohner_eps[PRJ_AMR_N];
    int saved_use_amr_angle_resolution;
    double saved_amr_angle_resolution_limit;
    double saved_E_floor;
    double saved_min_dx;
    int amr_idx;

    if (mesh == 0 || coord == 0 || root_nx1 <= 0 || root_nx2 <= 0 || root_nx3 <= 0) {
        return 1;
    }

    for (amr_idx = 0; amr_idx < PRJ_AMR_N; ++amr_idx) {
        saved_amr_refine_thresh[amr_idx] = mesh->amr_refine_thresh[amr_idx];
        saved_amr_derefine_thresh[amr_idx] = mesh->amr_derefine_thresh[amr_idx];
        saved_amr_estimator[amr_idx] = mesh->amr_estimator[amr_idx];
        saved_amr_lohner_var[amr_idx] = mesh->amr_lohner_var[amr_idx];
        saved_amr_lohner_eps[amr_idx] = mesh->amr_lohner_eps[amr_idx];
        saved_amr_criterion_set[amr_idx] = mesh->amr_criterion_set[amr_idx];
    }
    saved_use_amr_angle_resolution = mesh->use_amr_angle_resolution;
    saved_amr_angle_resolution_limit = mesh->amr_angle_resolution_limit;
    saved_E_floor = mesh->E_floor;
    saved_min_dx = mesh->min_dx;

    mesh->nblocks = 0;
    mesh->nblocks_max = 0;
    mesh->max_level = max_level;
    mesh->min_dx = saved_min_dx;
    mesh->max_active_level = -1;
    mesh->root_nx[0] = root_nx1;
    mesh->root_nx[1] = root_nx2;
    mesh->root_nx[2] = root_nx3;
    mesh->coord = *coord;
    for (amr_idx = 0; amr_idx < PRJ_AMR_N; ++amr_idx) {
        mesh->amr_refine_thresh[amr_idx] = saved_amr_refine_thresh[amr_idx];
        mesh->amr_derefine_thresh[amr_idx] = saved_amr_derefine_thresh[amr_idx];
        mesh->amr_estimator[amr_idx] = saved_amr_estimator[amr_idx];
        mesh->amr_lohner_var[amr_idx] = saved_amr_lohner_var[amr_idx];
        mesh->amr_lohner_eps[amr_idx] = saved_amr_lohner_eps[amr_idx];
        mesh->amr_criterion_set[amr_idx] = saved_amr_criterion_set[amr_idx];
    }
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
                            prj_neighbor_compute_geometry(b, &mesh->blocks[id], &b->slot[slot_index]);
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
