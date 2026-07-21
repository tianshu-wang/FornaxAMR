#include <math.h>
#include <limits.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include "prj.h"

#if defined(PRJ_ENABLE_MPI)
#include <mpi.h>
#endif

double prj_amr_angular_resolution_limit(double x1, double x2, double x3, double r_com)
{
    (void)x1;
    (void)x2;
    (void)x3;
    (void)r_com;
    return 0.024;
}

static double prj_abs_double(double x)
{
    return x < 0.0 ? -x : x;
}

static double prj_min_double(double a, double b)
{
    return a < b ? a : b;
}

static double prj_max_double(double a, double b)
{
    return a > b ? a : b;
}

static double prj_sqrt_double(double x)
{
    return sqrt(x);
}

static double prj_block_cell_size(const prj_block *block)
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

#if defined(PRJ_ENABLE_MPI)
#define PRJ_AMR_CHILD_TRANSFER_TAG 300
#endif

static void prj_amr_fatal(const char *message)
{
    if (message != 0) {
        fprintf(stderr, "%s\n", message);
    }
#if defined(PRJ_ENABLE_MPI)
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#endif
    exit(EXIT_FAILURE);
}

static int prj_clamp_storage_index(int idx)
{
    if (idx < -PRJ_NGHOST) {
        return -PRJ_NGHOST;
    }
    if (idx >= PRJ_BLOCK_SIZE + PRJ_NGHOST) {
        return PRJ_BLOCK_SIZE + PRJ_NGHOST - 1;
    }
    return idx;
}

static int prj_is_active_block(const prj_block *b)
{
    return b != 0 && b->id >= 0 && b->active == 1;
}

static int prj_is_local_active_block(const prj_mpi *mpi, const prj_block *b)
{
    return prj_is_active_block(b) && (mpi == 0 || b->rank == mpi->rank);
}

static int prj_is_local_block_owner(const prj_mpi *mpi, const prj_block *b)
{
    return b != 0 && b->id >= 0 && (mpi == 0 || b->rank == mpi->rank);
}

enum {
    PRJ_AMR_BOUNDARY_XMIN = 1 << 0,
    PRJ_AMR_BOUNDARY_XMAX = 1 << 1,
    PRJ_AMR_BOUNDARY_YMIN = 1 << 2,
    PRJ_AMR_BOUNDARY_YMAX = 1 << 3,
    PRJ_AMR_BOUNDARY_ZMIN = 1 << 4,
    PRJ_AMR_BOUNDARY_ZMAX = 1 << 5
};

static unsigned int prj_amr_refine_boundary_mask_for_cell(int i, int j, int k)
{
    const int boundary_buffer = PRJ_AMR_BUFFER_ZONE;
    unsigned int mask = 0U;

    if (i < boundary_buffer) {
        mask |= PRJ_AMR_BOUNDARY_XMIN;
    }
    if (i >= PRJ_BLOCK_SIZE - boundary_buffer) {
        mask |= PRJ_AMR_BOUNDARY_XMAX;
    }
    if (j < boundary_buffer) {
        mask |= PRJ_AMR_BOUNDARY_YMIN;
    }
    if (j >= PRJ_BLOCK_SIZE - boundary_buffer) {
        mask |= PRJ_AMR_BOUNDARY_YMAX;
    }
    if (k < boundary_buffer) {
        mask |= PRJ_AMR_BOUNDARY_ZMIN;
    }
    if (k >= PRJ_BLOCK_SIZE - boundary_buffer) {
        mask |= PRJ_AMR_BOUNDARY_ZMAX;
    }
    return mask;
}

static int prj_amr_boundary_mask_has_side(unsigned int mask, int axis, int sign)
{
    if (sign == 0) {
        return 1;
    }
    if (axis == 0) {
        return (mask & (sign < 0 ? PRJ_AMR_BOUNDARY_XMIN : PRJ_AMR_BOUNDARY_XMAX)) != 0U;
    }
    if (axis == 1) {
        return (mask & (sign < 0 ? PRJ_AMR_BOUNDARY_YMIN : PRJ_AMR_BOUNDARY_YMAX)) != 0U;
    }
    return (mask & (sign < 0 ? PRJ_AMR_BOUNDARY_ZMIN : PRJ_AMR_BOUNDARY_ZMAX)) != 0U;
}

static void prj_amr_sync_refine_flags(prj_mesh *mesh, const prj_mpi *mpi)
{
#if defined(PRJ_ENABLE_MPI)
    int *local_pos;
    int *global_pos;
    int *local_neg;
    int *global_neg;
    int i;

    if (mesh == 0 || mpi == 0 || mpi->totrank <= 1) {
        return;
    }
    local_pos = (int *)prj_calloc((size_t)mesh->nblocks, sizeof(*local_pos));
    global_pos = (int *)prj_calloc((size_t)mesh->nblocks, sizeof(*global_pos));
    local_neg = (int *)prj_calloc((size_t)mesh->nblocks, sizeof(*local_neg));
    global_neg = (int *)prj_calloc((size_t)mesh->nblocks, sizeof(*global_neg));
    if (local_pos == 0 || global_pos == 0 || local_neg == 0 || global_neg == 0) {
        free(local_pos);
        free(global_pos);
        free(local_neg);
        free(global_neg);
        return;
    }
    for (i = 0; i < mesh->nblocks; ++i) {
        if (mesh->blocks[i].refine_flag > 0 && mesh->blocks[i].can_refine != 0) {
            if (prj_is_active_block(&mesh->blocks[i])) {
                local_pos[i] = 1;
            }
        } else if (mesh->blocks[i].refine_flag < 0 &&
            prj_is_local_active_block(mpi, &mesh->blocks[i])) {
            local_neg[i] = -1;
        }
    }
    /*
     * Synchronize three-state AMR tags across ranks with clear precedence:
     * refine (+1) wins over derefine (-1), and derefine wins over neutral (0).
     * Each active block has a single owner, so a derefine request is produced by
     * that owner only, while MPI_MAX/MPI_MIN broadcast the agreed state to every rank.
     */
    MPI_Allreduce(local_pos, global_pos, mesh->nblocks, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(local_neg, global_neg, mesh->nblocks, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    for (i = 0; i < mesh->nblocks; ++i) {
        if (global_pos[i] > 0 && mesh->blocks[i].can_refine != 0) {
            mesh->blocks[i].refine_flag = 1;
        } else if (global_neg[i] < 0) {
            mesh->blocks[i].refine_flag = -1;
        } else {
            mesh->blocks[i].refine_flag = 0;
        }
    }
    free(local_pos);
    free(global_pos);
    free(local_neg);
    free(global_neg);
#else
    (void)mesh;
    (void)mpi;
#endif
}

static void prj_zero_block_arrays(prj_block *b)
{
    size_t n;
    size_t total;

    if (b == 0 || b->W_mhd == 0) {
        return;
    }

    total = prj_block_data_count();
    for (n = 0; n < total; ++n) {
        b->W_mhd[n] = 0.0;
    }
}

static void prj_amr_move_children_to_parent_rank(prj_mesh *mesh, const prj_mpi *mpi, prj_block *parent)
{
#if defined(PRJ_ENABLE_MPI)
    size_t data_count;
    int parent_rank;
    int oct;

    if (mesh == 0 || parent == 0) {
        return;
    }
    if (mpi == 0 || mpi->totrank <= 1) {
        return;
    }

    data_count = prj_block_data_count();
    parent_rank = parent->rank;
    for (oct = 0; oct < 8; ++oct) {
        int child_id;
        prj_block *child;
        int source_rank;

        child_id = parent->children[oct];
        if (child_id < 0 || child_id >= mesh->nblocks) {
            continue;
        }
        child = &mesh->blocks[child_id];
        source_rank = child->rank;
        if (source_rank == parent_rank) {
            continue;
        }

        if (mpi->rank == parent_rank) {
            if (child->W_mhd == 0 && prj_block_alloc_data(child) != 0) {
                prj_amr_fatal("prj_amr_move_children_to_parent_rank: failed to allocate receiving child data");
            }
            MPI_Recv(child->W_mhd, (int)data_count, MPI_DOUBLE, source_rank,
                PRJ_AMR_CHILD_TRANSFER_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else if (mpi->rank == source_rank) {
            if (child->W_mhd == 0) {
                prj_amr_fatal("prj_amr_move_children_to_parent_rank: source child is missing block data");
            }
            MPI_Send(child->W_mhd, (int)data_count, MPI_DOUBLE, parent_rank,
                PRJ_AMR_CHILD_TRANSFER_TAG, MPI_COMM_WORLD);
            prj_block_free_data(child);
        }
        child->rank = parent_rank;
    }
#else
    (void)mesh;
    (void)mpi;
    (void)parent;
#endif
}

static void prj_clear_neighbors(prj_block *b)
{
    int n;
    int d;

    for (n = 0; n < 56; ++n) {
        b->slot[n].id = -1;
        b->slot[n].rank = 0;
        b->slot[n].rel_level = 0;
        b->slot[n].type = PRJ_NEIGHBOR_NONE;
        for (d = 0; d < 3; ++d) {
            b->slot[n].xmin[d] = 0.0;
            b->slot[n].xmax[d] = 0.0;
            b->slot[n].dx[d] = 0.0;
            b->slot[n].send_loc_start[d] = 0;
            b->slot[n].send_loc_end[d] = 0;
            b->slot[n].recv_loc_start[d] = 0;
            b->slot[n].recv_loc_end[d] = 0;
        }
    }
}

static void prj_amr_mark_dirty_block_and_neighbors(const prj_mesh *mesh,
    int block_id, int *dirty)
{
    const prj_block *block;
    int n;

    if (mesh == 0 || dirty == 0 || block_id < 0 || block_id >= mesh->nblocks) {
        return;
    }
    dirty[block_id] = 1;
    block = &mesh->blocks[block_id];
    for (n = 0; n < 56; ++n) {
        int nid = block->slot[n].id;

        if (nid >= 0 && nid < mesh->nblocks) {
            dirty[nid] = 1;
        }
    }
}

static void prj_amr_mark_coarsen_dirty(const prj_mesh *mesh, int parent_id,
    int *dirty)
{
    const prj_block *parent;
    int oct;

    if (mesh == 0 || dirty == 0 || parent_id < 0 || parent_id >= mesh->nblocks) {
        return;
    }
    dirty[parent_id] = 1;
    parent = &mesh->blocks[parent_id];
    for (oct = 0; oct < 8; ++oct) {
        prj_amr_mark_dirty_block_and_neighbors(mesh, parent->children[oct],
            dirty);
    }
}

static void prj_amr_clear_dirty_mask(int *dirty, int count)
{
    int i;

    if (dirty == 0) {
        return;
    }
    for (i = 0; i < count; ++i) {
        dirty[i] = 0;
    }
}

static void prj_reset_children(prj_block *b)
{
    int n;

    for (n = 0; n < 8; ++n) {
        b->children[n] = -1;
    }
}

static int prj_find_free_block_slot(prj_mesh *mesh)
{
    int i;

    for (i = 0; i < mesh->nblocks; ++i) {
        if (mesh->blocks[i].id < 0) {
            return i;
        }
    }
    if (mesh->nblocks < mesh->nblocks_max) {
        i = mesh->nblocks;
        mesh->nblocks += 1;
        return i;
    }
    fprintf(stderr,
        "prj_find_free_block_slot: block capacity exceeded "
        "(nblocks=%d >= max_blocks=%d). Increase max_blocks in the param file.\n",
        mesh->nblocks, mesh->nblocks_max);
#if defined(PRJ_ENABLE_MPI)
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#endif
    exit(EXIT_FAILURE);
}

static int prj_boxes_overlap_or_touch(double amin, double amax, double bmin, double bmax, double tol)
{
    return prj_min_double(amax, bmax) >= prj_max_double(amin, bmin) - tol;
}

static int prj_boxes_touch_face_edge_corner(const prj_block *a, const prj_block *b)
{
    const double tol = 1.0e-12;
    int axis;
    int touches = 0;

    for (axis = 0; axis < 3; ++axis) {
        if (!prj_boxes_overlap_or_touch(a->xmin[axis], a->xmax[axis], b->xmin[axis], b->xmax[axis], tol)) {
            return 0;
        }
        if (prj_abs_double(a->xmax[axis] - b->xmin[axis]) < tol ||
            prj_abs_double(b->xmax[axis] - a->xmin[axis]) < tol) {
            touches += 1;
        }
    }
    if (touches == 0) {
        return 0;
    }
    if (prj_abs_double(a->xmin[0] - b->xmin[0]) < tol &&
        prj_abs_double(a->xmax[0] - b->xmax[0]) < tol &&
        prj_abs_double(a->xmin[1] - b->xmin[1]) < tol &&
        prj_abs_double(a->xmax[1] - b->xmax[1]) < tol &&
        prj_abs_double(a->xmin[2] - b->xmin[2]) < tol &&
        prj_abs_double(a->xmax[2] - b->xmax[2]) < tol) {
        return 0;
    }
    return 1;
}

static int prj_blocks_overlap_on_axis(const prj_block *a, const prj_block *b, int axis)
{
    const double tol = 1.0e-12;
    double overlap = prj_min_double(a->xmax[axis], b->xmax[axis]) -
        prj_max_double(a->xmin[axis], b->xmin[axis]);

    return overlap > tol;
}

static int prj_blocks_are_face_neighbors(const prj_block *a, const prj_block *b)
{
    const double tol = 1.0e-12;
    int axis;
    int touching_axes = 0;

    for (axis = 0; axis < 3; ++axis) {
        if (prj_abs_double(a->xmax[axis] - b->xmin[axis]) < tol ||
            prj_abs_double(b->xmax[axis] - a->xmin[axis]) < tol) {
            touching_axes += 1;
        } else if (!prj_blocks_overlap_on_axis(a, b, axis)) {
            return 0;
        }
    }
    return touching_axes == 1;
}

static int prj_add_neighbor(prj_block *a, const prj_block *b)
{
    int n;

    for (n = 0; n < 56; ++n) {
        if (a->slot[n].id == b->id) {
            return 0;
        }
    }

    for (n = 0; n < 56; ++n) {
        if (a->slot[n].id < 0) {
            int d;

            a->slot[n].id = b->id;
            a->slot[n].rank = b->rank;
            for (d = 0; d < 3; ++d) {
                a->slot[n].xmin[d] = b->xmin[d];
                a->slot[n].xmax[d] = b->xmax[d];
                a->slot[n].dx[d] = b->dx[d];
            }
            prj_neighbor_compute_geometry(a, b, &a->slot[n]);
            return 0;
        }
    }
    return 1;
}

static int prj_amr_neighbor_lookup_range(const prj_mesh *mesh,
    const prj_block *block, int level, int axis, int *start, int *end)
{
    double coord_min[3];
    double coord_max[3];
    double extent;
    double scale;
    double lo;
    double hi;
    double level_count;

    if (mesh == 0 || block == 0 || level < 0 || axis < 0 || axis >= 3 ||
        start == 0 || end == 0 || mesh->root_nx[axis] <= 0) {
        return 0;
    }
    coord_min[0] = mesh->coord.x1min;
    coord_min[1] = mesh->coord.x2min;
    coord_min[2] = mesh->coord.x3min;
    coord_max[0] = mesh->coord.x1max;
    coord_max[1] = mesh->coord.x2max;
    coord_max[2] = mesh->coord.x3max;

    extent = coord_max[axis] - coord_min[axis];
    level_count = ldexp((double)mesh->root_nx[axis], level);
    if (extent <= 0.0 || !isfinite(level_count) || level_count <= 0.0 ||
        level_count > (double)INT_MAX) {
        return 0;
    }
    scale = level_count / extent;
    lo = (block->xmin[axis] - coord_min[axis]) * scale;
    hi = (block->xmax[axis] - coord_min[axis]) * scale;
    if (!isfinite(lo) || !isfinite(hi) || lo > hi ||
        lo < (double)INT_MIN || hi > (double)INT_MAX) {
        return 0;
    }
    *start = (int)floor(lo) - 1;
    *end = (int)ceil(hi);
    return *start <= *end;
}

static void prj_amr_try_neighbor_candidate(prj_mesh *mesh, prj_block *a, int nid)
{
    prj_block *b;

    if (mesh == 0 || a == 0 || nid < 0 || nid >= mesh->nblocks || nid == a->id) {
        return;
    }
    b = &mesh->blocks[nid];
    if (!prj_is_active_block(b)) {
        return;
    }
    if (prj_boxes_touch_face_edge_corner(a, b)) {
        prj_add_neighbor(a, b);
        prj_add_neighbor(b, a);
    }
}

static void prj_amr_probe_neighbor_level(prj_mesh *mesh, prj_block *block, int level)
{
    int is;
    int ie;
    int js;
    int je;
    int ks;
    int ke;
    int ix;
    int iy;
    int iz;

    if (!prj_amr_neighbor_lookup_range(mesh, block, level, 0, &is, &ie) ||
        !prj_amr_neighbor_lookup_range(mesh, block, level, 1, &js, &je) ||
        !prj_amr_neighbor_lookup_range(mesh, block, level, 2, &ks, &ke)) {
        return;
    }

    for (ix = is; ix <= ie; ++ix) {
        for (iy = js; iy <= je; ++iy) {
            for (iz = ks; iz <= ke; ++iz) {
                int nid = prj_mesh_morton_lookup_block(mesh, level, ix, iy, iz);

                prj_amr_try_neighbor_candidate(mesh, block, nid);
            }
        }
    }
}

static void prj_amr_init_neighbors_pairwise(prj_mesh *mesh,
    const int *rebuild_mask)
{
    int i;
    int j;

    for (i = 0; i < mesh->nblocks; ++i) {
        prj_block *a = &mesh->blocks[i];

        if (!prj_is_active_block(a)) {
            continue;
        }
        for (j = i + 1; j < mesh->nblocks; ++j) {
            prj_block *b = &mesh->blocks[j];

            if (!prj_is_active_block(b)) {
                continue;
            }
            if (rebuild_mask != 0 && rebuild_mask[i] == 0 &&
                rebuild_mask[j] == 0) {
                continue;
            }
            if (prj_boxes_touch_face_edge_corner(a, b)) {
                prj_add_neighbor(a, b);
                prj_add_neighbor(b, a);
            }
        }
    }
}

static void prj_amr_init_neighbors_morton(prj_mesh *mesh,
    const int *rebuild_mask)
{
    int i;

    for (i = 0; i < mesh->nblocks; ++i) {
        prj_block *block = &mesh->blocks[i];
        int level;

        if (!prj_is_active_block(block)) {
            continue;
        }
        if (rebuild_mask != 0 && rebuild_mask[i] == 0) {
            continue;
        }
        for (level = block->level - 1; level <= block->level + 1; ++level) {
            if (level >= 0) {
                prj_amr_probe_neighbor_level(mesh, block, level);
            }
        }
    }
}

static int prj_amr_clamp_storage_index(int idx)
{
    if (idx < -PRJ_NGHOST) {
        return -PRJ_NGHOST;
    }
    if (idx >= PRJ_BLOCK_SIZE + PRJ_NGHOST) {
        return PRJ_BLOCK_SIZE + PRJ_NGHOST - 1;
    }
    return idx;
}

static double prj_block_primitive_at(const prj_block *b, int v, int i, int j, int k)
{
    i = prj_amr_clamp_storage_index(i);
    j = prj_amr_clamp_storage_index(j);
    k = prj_amr_clamp_storage_index(k);
    return prj_block_prim_value_const(b, 0, v, i, j, k);
}

static double prj_amr_metric_dot(const double metric[3][3], const double a[3], const double b[3])
{
    double dot = 0.0;
    int m;
    int n;

    for (m = 0; m < 3; ++m) {
        for (n = 0; n < 3; ++n) {
            dot += metric[m][n] * a[m] * b[n];
        }
    }
    return dot;
}

static double prj_amr_gr_magnetic_pressure_at(
    const prj_z4c_hydro_geom *geom, const prj_block *b, int i, int j, int k)
{
#if PRJ_MHD
    double Bcon[3];
    double beta_con[3];
    double Bsq;
    double beta2;
    double Bbeta;
    double wlor2;
    double pmag;
    int d;

    if (geom == 0 || b == 0 || b->W_mhd == 0) {
        prj_amr_fatal("AMR error: full-GR magnetic pressure requires Z4c geometry and MHD primitives");
    }

    for (d = 0; d < 3; ++d) {
        Bcon[d] = prj_block_primitive_at(b, PRJ_PRIM_B1 + d, i, j, k);
        beta_con[d] = prj_block_primitive_at(b, PRJ_PRIM_V1 + d, i, j, k) / PRJ_CLIGHT;
    }

    Bsq = prj_amr_metric_dot(geom->gamma, Bcon, Bcon);
    beta2 = prj_amr_metric_dot(geom->gamma, beta_con, beta_con);
    Bbeta = prj_amr_metric_dot(geom->gamma, Bcon, beta_con);
    if (!isfinite(Bsq) || !isfinite(beta2) || !isfinite(Bbeta) ||
        Bsq < 0.0 || beta2 < 0.0 || beta2 >= 1.0) {
        prj_amr_fatal("AMR error: invalid full-GR metric contraction in magnetic pressure");
    }

    wlor2 = 1.0 / (1.0 - beta2);
    pmag = 0.5 * (Bsq / wlor2 + Bbeta * Bbeta);
    if (!isfinite(pmag) || pmag < 0.0) {
        prj_amr_fatal("AMR error: invalid full-GR magnetic pressure");
    }
    return pmag;
#else
    (void)geom;
    (void)b;
    (void)i;
    (void)j;
    (void)k;
    return 0.0;
#endif
}

static double prj_amr_total_pressure_with_gr_geom_at(
    const prj_block *b, const prj_z4c_hydro_geom *geom, int i, int j, int k)
{
    double pressure;

    if (b == 0 || geom == 0 || b->eosvar == 0) {
        prj_amr_fatal("AMR error: full-GR pressure requires EOS variables and Z4c geometry");
    }
    pressure = b->eosvar[EIDX(PRJ_EOSVAR_PRESSURE, i, j, k)];
    pressure += prj_amr_gr_magnetic_pressure_at(geom, b, i, j, k);
    return pressure;
}

static double prj_amr_total_pressure_non_gr_at(const prj_block *b, int i, int j, int k)
{
    double pressure = b->eosvar[EIDX(PRJ_EOSVAR_PRESSURE, i, j, k)];
#if PRJ_MHD
    double b1 = b->W_mhd[WIDX(PRJ_PRIM_B1, i, j, k)];
    double b2 = b->W_mhd[WIDX(PRJ_PRIM_B2, i, j, k)];
    double b3 = b->W_mhd[WIDX(PRJ_PRIM_B3, i, j, k)];

    pressure += 0.5 * (b1 * b1 + b2 * b2 + b3 * b3);
#endif
    return pressure;
}

/* Total pressure used by the AMR estimators.  Non-GR runs keep the historical
 * gas + 0.5*B^2 pressure.  Full dynamic GR uses the Z4c spatial metric and the
 * same GRMHD magnetic-pressure convention as prj_eos_grmhd_state_from_prim. */
static double prj_amr_total_pressure_at(
    const prj_mesh *mesh, const prj_block *b, int stage, int i, int j, int k)
{
#if PRJ_DYNAMIC_GR
    if (prj_eos_full_dynamic_gr_enabled(mesh)) {
        prj_z4c_hydro_geom geom;

        if (!prj_z4c_load_hydro_geom(mesh, b, stage, i, j, k, &geom)) {
            prj_amr_fatal("AMR error: failed to load full-GR geometry for pressure estimator");
        }
        return prj_amr_total_pressure_with_gr_geom_at(b, &geom, i, j, k);
    }
#else
    (void)mesh;
    (void)stage;
#endif
    return prj_amr_total_pressure_non_gr_at(b, i, j, k);
}

static double prj_block_sound_speed_at(
    const prj_mesh *mesh, const prj_block *b, prj_eos *eos, int i, int j, int k)
{
    double rho;
    double eint;
    double pressure;
    double gamma;

    rho = prj_block_primitive_at(b, PRJ_PRIM_RHO, i, j, k);
    eint = prj_block_primitive_at(b, PRJ_PRIM_EINT, i, j, k);
    (void)eos;
    if (rho <= 0.0 || eint < 0.0) {
        return 0.0;
    }
    /* Use total (gas + magnetic) pressure so the velocity estimator normalises
     * by a fast-magnetosonic-like speed in MHD runs. */
    pressure = prj_amr_total_pressure_at(mesh, b, 0, i, j, k);
    gamma = b->eosvar[EIDX(PRJ_EOSVAR_GAMMA, i, j, k)];
    if (pressure <= 0.0 || gamma <= 0.0) {
        return 0.0;
    }
    return prj_sqrt_double(gamma * pressure / rho);
}

static double prj_block_conserved_at(const prj_block *b, int v, int i, int j, int k)
{
    i = prj_clamp_storage_index(i);
    j = prj_clamp_storage_index(j);
    k = prj_clamp_storage_index(k);
    return prj_block_cons_value_const(b, v, i, j, k);
}

static void prj_apply_eint_floor(prj_eos *eos, double E_floor, double cell_vol,
    double *U, double *W, double *e_injected)
{
    double rho;
    double eint_floor;
    double v1;
    double v2;
    double v3;
    double kinetic;
    double etot_old;
    double etot_new;
#if PRJ_MHD
    double magnetic;
#endif

    if (E_floor <= 0.0 || U == 0 || W == 0) {
        return;
    }

    rho = W[PRJ_PRIM_RHO];
    if (rho <= 0.0) {
        return;
    }

    /* For the tabulated EOS the floor is measured relative to the table's
     * low-temperature boundary internal energy (a function of rho and Ye);
     * the helper returns 0 for the ideal-gas EOS so the floor is unchanged. */
    eint_floor = E_floor + prj_eos_low_temp_eint(eos, rho, W[PRJ_PRIM_YE], PRJ_EOS_CTX_AMR);

    if (W[PRJ_PRIM_EINT] >= eint_floor) {
        return;
    }

    v1 = W[PRJ_PRIM_V1];
    v2 = W[PRJ_PRIM_V2];
    v3 = W[PRJ_PRIM_V3];
    kinetic = 0.5 * (v1 * v1 + v2 * v2 + v3 * v3);
#if PRJ_MHD
    magnetic = 0.5 * (W[PRJ_PRIM_B1] * W[PRJ_PRIM_B1] +
        W[PRJ_PRIM_B2] * W[PRJ_PRIM_B2] +
        W[PRJ_PRIM_B3] * W[PRJ_PRIM_B3]);
#endif
    etot_old = U[PRJ_CONS_ETOT];
    etot_new = rho * (eint_floor + kinetic)
#if PRJ_MHD
        + magnetic
#endif
        ;
    W[PRJ_PRIM_EINT] = eint_floor;
    U[PRJ_CONS_ETOT] = etot_new;

    /* U holds conserved energy density; weight the increase by cell volume so
     * the accumulator tracks the actual extra energy added by the floor. */
    if (e_injected != 0) {
        *e_injected += (etot_new - etot_old) * cell_vol;
    }
}

static double prj_loehner_cell_value(const prj_mesh *mesh, const prj_block *b, int lohner_var, int i, int j, int k)
{
    if (b == 0) {
        return 0.0;
    }
    if (lohner_var == PRJ_LOHNER_VAR_DENSITY) {
        return prj_block_primitive_at(b, PRJ_PRIM_RHO, i, j, k);
    }
    if (lohner_var == PRJ_LOHNER_VAR_TEMPERATURE) {
        return b->eosvar[EIDX(PRJ_EOSVAR_TEMPERATURE, i, j, k)];
    }
    return prj_amr_total_pressure_at(mesh, b, 0, i, j, k);
}

static double prj_loehner_cell_indicator(
    const prj_mesh *mesh, const prj_block *b, prj_eos *eos, int lohner_var, double lohner_eps, int i, int j, int k)
{
    const double small = 1.0e-14;
    double alpha = lohner_eps;
    double u0;
    double uxm;
    double uxp;
    double uym;
    double uyp;
    double uzm;
    double uzp;
    double num_x;
    double num_y;
    double num_z;
    double num_xy;
    double num_xz;
    double num_yz;
    double den_x;
    double den_y;
    double den_z;
    double numerator;
    double denominator;

    (void)eos;
    u0 = prj_loehner_cell_value(mesh, b, lohner_var, i, j, k);
    uxm = prj_loehner_cell_value(mesh, b, lohner_var, i - 1, j, k);
    uxp = prj_loehner_cell_value(mesh, b, lohner_var, i + 1, j, k);
    uym = prj_loehner_cell_value(mesh, b, lohner_var, i, j - 1, k);
    uyp = prj_loehner_cell_value(mesh, b, lohner_var, i, j + 1, k);
    uzm = prj_loehner_cell_value(mesh, b, lohner_var, i, j, k - 1);
    uzp = prj_loehner_cell_value(mesh, b, lohner_var, i, j, k + 1);
    num_x = uxp - 2.0 * u0 + uxm;
    num_y = uyp - 2.0 * u0 + uym;
    num_z = uzp - 2.0 * u0 + uzm;
    num_xy = 0.25 * (
        prj_loehner_cell_value(mesh, b, lohner_var, i + 1, j + 1, k) -
        prj_loehner_cell_value(mesh, b, lohner_var, i + 1, j - 1, k) -
        prj_loehner_cell_value(mesh, b, lohner_var, i - 1, j + 1, k) +
        prj_loehner_cell_value(mesh, b, lohner_var, i - 1, j - 1, k));
    num_xz = 0.25 * (
        prj_loehner_cell_value(mesh, b, lohner_var, i + 1, j, k + 1) -
        prj_loehner_cell_value(mesh, b, lohner_var, i + 1, j, k - 1) -
        prj_loehner_cell_value(mesh, b, lohner_var, i - 1, j, k + 1) +
        prj_loehner_cell_value(mesh, b, lohner_var, i - 1, j, k - 1));
    num_yz = 0.25 * (
        prj_loehner_cell_value(mesh, b, lohner_var, i, j + 1, k + 1) -
        prj_loehner_cell_value(mesh, b, lohner_var, i, j + 1, k - 1) -
        prj_loehner_cell_value(mesh, b, lohner_var, i, j - 1, k + 1) +
        prj_loehner_cell_value(mesh, b, lohner_var, i, j - 1, k - 1));
    den_x = prj_abs_double(uxp - u0) + prj_abs_double(u0 - uxm) +
        alpha * (prj_abs_double(uxp) + 2.0 * prj_abs_double(u0) + prj_abs_double(uxm));
    den_y = prj_abs_double(uyp - u0) + prj_abs_double(u0 - uym) +
        alpha * (prj_abs_double(uyp) + 2.0 * prj_abs_double(u0) + prj_abs_double(uym));
    den_z = prj_abs_double(uzp - u0) + prj_abs_double(u0 - uzm) +
        alpha * (prj_abs_double(uzp) + 2.0 * prj_abs_double(u0) + prj_abs_double(uzm));
    numerator = num_x * num_x + num_y * num_y + num_z * num_z +
        2.0 * (num_xy * num_xy + num_xz * num_xz + num_yz * num_yz);
    denominator = den_x * den_x + den_y * den_y + den_z * den_z + small;

    return prj_sqrt_double(numerator / denominator);
}

static double prj_velocity_cell_indicator(
    const prj_mesh *mesh, const prj_block *b, prj_eos *eos, int i, int j, int k)
{
    const double small = 1.0e-14;
    double v1;
    double v2;
    double v3;
    double cs;
    double max_indicator = 0.0;
    int axis;

    v1 = prj_block_primitive_at(b, PRJ_PRIM_V1, i, j, k);
    v2 = prj_block_primitive_at(b, PRJ_PRIM_V2, i, j, k);
    v3 = prj_block_primitive_at(b, PRJ_PRIM_V3, i, j, k);
    cs = prj_block_sound_speed_at(mesh, b, eos, i, j, k);
    for (axis = 0; axis < 3; ++axis) {
        int side;

        for (side = -1; side <= 1; side += 2) {
            int di = axis == 0 ? side : 0;
            int dj = axis == 1 ? side : 0;
            int dk = axis == 2 ? side : 0;
            double dv1 = v1 - prj_block_primitive_at(b, PRJ_PRIM_V1, i + di, j + dj, k + dk);
            double dv2 = v2 - prj_block_primitive_at(b, PRJ_PRIM_V2, i + di, j + dj, k + dk);
            double dv3 = v3 - prj_block_primitive_at(b, PRJ_PRIM_V3, i + di, j + dj, k + dk);
            double dv = prj_sqrt_double(dv1 * dv1 + dv2 * dv2 + dv3 * dv3) / (cs + small);

            max_indicator = prj_max_double(max_indicator, dv);
        }
    }
    return max_indicator;
}

static double prj_pressure_scale_height_cell_indicator(
    const prj_mesh *mesh, const prj_block *b, int i, int j, int k)
{
    int cache_idx;
    double rho;
    double eint;
    double pressure_gas;
    double pressure;
    double accel;
    double g1;
    double g2;
    double g3;
    double cell_size;
    double Hp;

    if (prj_eos_full_dynamic_gr_enabled(mesh)) {
        prj_z4c_hydro_geom geom;
        double dalpha2;
        double c2 = PRJ_CLIGHT * PRJ_CLIGHT;
        double rho_eff;

        if (b == 0 || b->eosvar == 0) {
            return 0.0;
        }
        if (!prj_z4c_load_hydro_geom(mesh, b, 0, i, j, k, &geom)) {
            prj_amr_fatal("AMR error: failed to load full-GR geometry for pressure-scale-height estimator");
        }

        rho = prj_block_primitive_at(b, PRJ_PRIM_RHO, i, j, k);
        eint = prj_block_primitive_at(b, PRJ_PRIM_EINT, i, j, k);
        pressure_gas = b->eosvar[EIDX(PRJ_EOSVAR_PRESSURE, i, j, k)];
        pressure = prj_amr_total_pressure_with_gr_geom_at(b, &geom, i, j, k);
        dalpha2 = prj_amr_metric_dot((const double (*)[3])geom.gamma_inv,
            geom.dalpha, geom.dalpha);
        if (!isfinite(dalpha2) || dalpha2 < 0.0 || !isfinite(geom.alpha) || geom.alpha <= 0.0) {
            prj_amr_fatal("AMR error: invalid full-GR lapse-gradient contraction in pressure-scale-height estimator");
        }
        accel = c2 * prj_sqrt_double(dalpha2) / geom.alpha;
        rho_eff = rho + (rho * eint + pressure_gas) / c2;
        cell_size = prj_block_cell_size(b);

        if (rho_eff <= 0.0 || pressure <= 0.0 || accel <= 0.0 || cell_size <= 0.0) {
            return 0.0;
        }
        Hp = pressure / (rho_eff * accel);
        if (!isfinite(Hp) || Hp <= 0.0) {
            return 0.0;
        }
        return cell_size / Hp;
    }

    if (b == 0 || b->grav[0] == 0 || b->grav[1] == 0 || b->grav[2] == 0) {
        return 0.0;
    }

    cache_idx = prj_block_cache_index(i, j, k);
    rho = prj_block_primitive_at(b, PRJ_PRIM_RHO, i, j, k);
    pressure = prj_amr_total_pressure_non_gr_at(b, i, j, k);
    g1 = b->grav[0][cache_idx];
    g2 = b->grav[1][cache_idx];
    g3 = b->grav[2][cache_idx];
    accel = prj_sqrt_double(g1 * g1 + g2 * g2 + g3 * g3);
    cell_size = b->dx[0];
    if (b->dx[1] > cell_size) {
        cell_size = b->dx[1];
    }
    if (b->dx[2] > cell_size) {
        cell_size = b->dx[2];
    }

    if (rho <= 0.0 || pressure <= 0.0 || accel <= 0.0) {
        return 0.0;
    }

    Hp = pressure / (rho * accel);
    if (Hp <= 0.0) {
        return 0.0;
    }
    return cell_size / Hp;
}

static double prj_fractional_jump_cell_value(
    const prj_mesh *mesh, const prj_block *b, int jump_var, int i, int j, int k)
{
    if (jump_var == PRJ_FRACTIONAL_JUMP_VAR_PRESSURE) {
        return prj_amr_total_pressure_at(mesh, b, 0, i, j, k);
    }
    return b->W_mhd[WIDX(PRJ_PRIM_RHO, i, j, k)];
}

static double prj_fractional_jump_cell_indicator(
    const prj_mesh *mesh, const prj_block *b, int jump_var, int i, int j, int k)
{
    const double small = 1.0e-14;
    double q0;
    double max_indicator = 0.0;
    int di;
    int dj;
    int dk;

    if (b == 0) {
        return 0.0;
    }
    if (jump_var == PRJ_FRACTIONAL_JUMP_VAR_PRESSURE && b->eosvar == 0) {
        return 0.0;
    }
    if (jump_var != PRJ_FRACTIONAL_JUMP_VAR_PRESSURE && b->W_mhd == 0) {
        return 0.0;
    }

    q0 = prj_fractional_jump_cell_value(mesh, b, jump_var, i, j, k);
    if (q0 <= 0.0) {
        return 0.0;
    }

    for (di = -1; di <= 1; ++di) {
        for (dj = -1; dj <= 1; ++dj) {
            for (dk = -1; dk <= 1; ++dk) {
                double qnei;
                double denom;
                double jump;

                if (di == 0 && dj == 0 && dk == 0) {
                    continue;
                }
                qnei = prj_fractional_jump_cell_value(mesh, b, jump_var, i + di, j + dj, k + dk);
                if (qnei <= 0.0) {
                    continue;
                }
                denom = q0 < qnei ? q0 : qnei;
                jump = prj_abs_double(q0 - qnei) / (denom + small);
                max_indicator = prj_max_double(max_indicator, jump);
            }
        }
    }

    return max_indicator;
}

static double prj_amr_cell_indicator_for_estimator(
    const prj_mesh *mesh, const prj_block *b, prj_eos *eos, int amr_idx, int estimator, int i, int j, int k)
{
    if (mesh != 0 && estimator == PRJ_AMR_ESTIMATOR_FRACTIONAL_JUMP) {
        return prj_fractional_jump_cell_indicator(mesh, b, mesh->amr_fractional_jump_var[amr_idx], i, j, k);
    }
    if (mesh != 0 && estimator == PRJ_AMR_ESTIMATOR_PRESSURE_SCALE_HEIGHT) {
        return prj_pressure_scale_height_cell_indicator(mesh, b, i, j, k);
    }
    if (mesh != 0 && estimator == PRJ_AMR_ESTIMATOR_VELOCITY) {
        return prj_velocity_cell_indicator(mesh, b, eos, i, j, k);
    }
    return prj_loehner_cell_indicator(
        mesh, b, eos, mesh->amr_lohner_var[amr_idx], mesh->amr_lohner_eps[amr_idx], i, j, k);
}

int prj_amr_criteria_need_eosvar(const prj_mesh *mesh)
{
    int amr_idx;

    if (mesh == 0) {
        return 0;
    }

    for (amr_idx = 0; amr_idx < PRJ_AMR_N; ++amr_idx) {
        int estimator;

        if (mesh->amr_criterion_set[amr_idx] == 0) {
            continue;
        }
        estimator = mesh->amr_estimator[amr_idx];
        if (estimator == PRJ_AMR_ESTIMATOR_VELOCITY ||
            estimator == PRJ_AMR_ESTIMATOR_PRESSURE_SCALE_HEIGHT) {
            return 1;
        }
        if (estimator == PRJ_AMR_ESTIMATOR_FRACTIONAL_JUMP &&
            mesh->amr_fractional_jump_var[amr_idx] == PRJ_FRACTIONAL_JUMP_VAR_PRESSURE) {
            return 1;
        }
        if (estimator == PRJ_AMR_ESTIMATOR_LOEHNER &&
            mesh->amr_lohner_var[amr_idx] != PRJ_LOHNER_VAR_DENSITY) {
            return 1;
        }
    }
    return 0;
}

static int prj_has_face_neighbor_coarser_than(const prj_mesh *mesh, const prj_block *b, int min_level)
{
    int n;

    for (n = 0; n < 56; ++n) {
        int id = b->slot[n].id;

        if (id >= 0 && id < mesh->nblocks &&
            prj_is_active_block(&mesh->blocks[id]) &&
            prj_blocks_are_face_neighbors(b, &mesh->blocks[id]) &&
            mesh->blocks[id].level < min_level) {
            return 1;
        }
    }
    return 0;
}

static int prj_block_neighbor_touch_offset(const prj_block *block, const prj_block *neighbor, int axis)
{
    const double tol = 1.0e-12;

    if (block == 0 || neighbor == 0 || axis < 0 || axis >= 3) {
        return 0;
    }
    if (prj_abs_double(block->xmin[axis] - neighbor->xmax[axis]) < tol) {
        return -1;
    }
    if (prj_abs_double(block->xmax[axis] - neighbor->xmin[axis]) < tol) {
        return 1;
    }
    return 0;
}

static void prj_amr_tag_boundary_neighbors(prj_mesh *mesh, const prj_block *block,
    unsigned int boundary_mask, int *local_pos)
{
    int n;

    if (mesh == 0 || block == 0 || local_pos == 0 || boundary_mask == 0U) {
        return;
    }

    for (n = 0; n < 56; ++n) {
        int id = block->slot[n].id;
        prj_block *neighbor;
        int axis;
        int matches = 0;
        int ok = 1;

        if (id < 0 || id >= mesh->nblocks) {
            continue;
        }
        neighbor = &mesh->blocks[id];
        if (!prj_is_active_block(neighbor)) {
            continue;
        }
        if (neighbor->level > block->level) {
            continue;
        }
        for (axis = 0; axis < 3; ++axis) {
            int offset = prj_block_neighbor_touch_offset(block, neighbor, axis);

            if (offset == 0) {
                continue;
            }
            matches = 1;
            if (!prj_amr_boundary_mask_has_side(boundary_mask, axis, offset)) {
                ok = 0;
                break;
            }
        }
        if (ok != 0 && matches != 0) {
            local_pos[id] = 1;
        }
    }
}

static int prj_can_coarsen_parent(const prj_mesh *mesh, int parent_id)
{
    const prj_block *parent;
    int oct;

    if (parent_id < 0 || parent_id >= mesh->nblocks) {
        return 0;
    }
    parent = &mesh->blocks[parent_id];
    if (parent->id < 0 || parent->active == 1) {
        return 0;
    }

    for (oct = 0; oct < 8; ++oct) {
        int child_id = parent->children[oct];
        int n;

        if (child_id < 0 || child_id >= mesh->nblocks ||
            !prj_is_active_block(&mesh->blocks[child_id])) {
            return 0;
        }
        if (mesh->blocks[child_id].refine_flag > 0) {
            return 0;
        }
        if (mesh->blocks[child_id].refine_flag >= 0) {
            return 0;
        }
        for (n = 0; n < 56; ++n) {
            int nid = mesh->blocks[child_id].slot[n].id;

            if (nid >= 0 && nid < mesh->nblocks &&
                prj_is_active_block(&mesh->blocks[nid]) &&
                mesh->blocks[nid].parent != parent_id &&
                mesh->blocks[nid].level > parent->level + 1) {
                return 0;
            }
        }
    }
    return 1;
}

static int prj_amr_has_pending_change(const prj_mesh *mesh)
{
    int i;

    if (mesh == 0) {
        return 0;
    }

    for (i = 0; i < mesh->nblocks; ++i) {
        const prj_block *block = &mesh->blocks[i];

        if (prj_is_active_block(block) &&
            block->refine_flag > 0 && block->can_refine != 0) {
            return 1;
        }
    }
    for (i = 0; i < mesh->nblocks; ++i) {
        const prj_block *parent = &mesh->blocks[i];

        if (parent->id >= 0 && parent->active != 1 &&
            prj_can_coarsen_parent(mesh, i)) {
            return 1;
        }
    }
    return 0;
}

static void prj_sync_primitive_from_conserved(prj_mesh *mesh, prj_eos *eos,
    const prj_mpi *mpi, const int *changed_blocks, double *e_injected)
{
    int bidx;

    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *b = &mesh->blocks[bidx];
        int i;
        int j;
        int k;

        if (!prj_is_local_active_block(mpi, b)) {
            continue;
        }
        /* Only blocks regenerated by this adapt pass (refine children /
         * coarsen parents) need a prim resync and floor; leave others
         * untouched. A NULL mask falls back to syncing every block. */
        if (changed_blocks != 0 && changed_blocks[bidx] == 0) {
            continue;
        }
        for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
            for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    double Uc[PRJ_NVAR_CONS];
                    double Wc[PRJ_NVAR_PRIM];

                    prj_block_load_cons_cell_const(b, i, j, k, Uc);
                    prj_eos_cell_cons2prim(eos, mesh, b, 0, i, j, k, Uc, Wc,
                        PRJ_EOS_CTX_AMR);
                    prj_apply_eint_floor(eos, mesh->E_floor, b->vol, Uc, Wc, e_injected);
                    prj_block_store_cons_cell(b, i, j, k, Uc);
                    prj_block_store_prim_cell(b, 0, i, j, k, Wc);
                }
            }
        }
    }
}

static void prj_sync_conserved_from_primitive(prj_mesh *mesh, prj_eos *eos,
    const prj_mpi *mpi)
{
    int bidx;

    if (mesh == 0 || eos == 0) {
        return;
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        int i;
        int j;
        int k;

        if (!prj_is_local_active_block(mpi, block) || block->W_mhd == 0 ||
            !prj_block_has_cons_storage(block)) {
            continue;
        }
        for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
            for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    double Wc[PRJ_NVAR_PRIM];
                    double Uc[PRJ_NVAR_CONS];

                    prj_block_load_prim_cell_const(block, 0, i, j, k, Wc);
                    prj_eos_cell_prim2cons(eos, mesh, block, 0, i, j, k, Wc, Uc,
                        PRJ_EOS_CTX_AMR);
                    prj_block_store_cons_cell(block, i, j, k, Uc);
                }
            }
        }
    }
}

void prj_amr_enforce_two_to_one(prj_mesh *mesh, const prj_mpi *mpi)
{
    int changed;

    do {
        int i;

        changed = 0;
        for (i = 0; i < mesh->nblocks; ++i) {
            prj_block *b = &mesh->blocks[i];
            int n;

            if (!prj_is_active_block(b) || b->refine_flag <= 0) {
                continue;
            }
            for (n = 0; n < 56; ++n) {
                int id = b->slot[n].id;

                if (id >= 0 && id < mesh->nblocks && prj_is_active_block(&mesh->blocks[id]) &&
                    mesh->blocks[id].level < b->level && mesh->blocks[id].can_refine != 0) {
                    if (mesh->blocks[id].refine_flag != 1) {
                        mesh->blocks[id].refine_flag = 1;
                        changed = 1;
                    }
                }
            }
        }
        prj_amr_sync_refine_flags(mesh, mpi);
    } while (changed != 0);
}

static int prj_amr_refine_marked_blocks_with_dirty(prj_mesh *mesh,
    const prj_mpi *mpi, int *neighbor_dirty, int *changed_blocks)
{
    int i;
    int changed = 0;

    if (mesh == 0) {
        return 0;
    }

    for (i = 0; i < mesh->nblocks; ++i) {
        if (i < mesh->nblocks && prj_is_active_block(&mesh->blocks[i]) &&
            mesh->blocks[i].refine_flag > 0 && mesh->blocks[i].can_refine != 0) {
            int was_active = mesh->blocks[i].active;
            int oct;

            prj_amr_mark_dirty_block_and_neighbors(mesh, i, neighbor_dirty);
            prj_amr_refine_block(mesh, mpi, i);
            if (neighbor_dirty != 0) {
                neighbor_dirty[i] = 1;
            }
            for (oct = 0; oct < 8; ++oct) {
                int child_id = mesh->blocks[i].children[oct];

                if (child_id >= 0 && child_id < mesh->nblocks) {
                    if (neighbor_dirty != 0) {
                        neighbor_dirty[child_id] = 1;
                    }
                    /* Newly created children carry freshly prolongated data
                     * and are the only blocks that need a prim/floor resync. */
                    if (changed_blocks != 0) {
                        changed_blocks[child_id] = 1;
                    }
                }
            }
            if (was_active != mesh->blocks[i].active) {
                changed = 1;
            }
        }
    }
    return changed;
}

int prj_amr_refine_marked_blocks(prj_mesh *mesh, const prj_mpi *mpi)
{
    return prj_amr_refine_marked_blocks_with_dirty(mesh, mpi, 0, 0);
}

static void prj_amr_init_neighbors_with_mask(prj_mesh *mesh,
    const int *rebuild_mask)
{
    int i;

    for (i = 0; i < mesh->nblocks; ++i) {
        if (mesh->blocks[i].id >= 0 &&
            (rebuild_mask == 0 || rebuild_mask[i] != 0)) {
            prj_clear_neighbors(&mesh->blocks[i]);
        }
    }

    if (prj_mesh_rebuild_morton_lookup(mesh) == 0) {
        prj_amr_init_neighbors_morton(mesh, rebuild_mask);
    } else {
        prj_amr_init_neighbors_pairwise(mesh, rebuild_mask);
    }
    prj_mesh_update_cell_derived_mask(mesh);
}

void prj_amr_init_neighbors(prj_mesh *mesh)
{
    prj_amr_init_neighbors_with_mask(mesh, 0);
}

void prj_amr_tag(prj_mesh *mesh, prj_eos *eos, const prj_mpi *mpi)
{
    int i;
    int *local_pos;
    int *local_neg;
    int *global_pos;
    int *global_neg;
    unsigned int *boundary_mask;

    if (mesh == 0 || mesh->nblocks <= 0) {
        return;
    }

    local_pos = (int *)prj_calloc((size_t)mesh->nblocks, sizeof(*local_pos));
    local_neg = (int *)prj_calloc((size_t)mesh->nblocks, sizeof(*local_neg));
    global_pos = (int *)prj_calloc((size_t)mesh->nblocks, sizeof(*global_pos));
    global_neg = (int *)prj_calloc((size_t)mesh->nblocks, sizeof(*global_neg));
    boundary_mask = (unsigned int *)prj_calloc((size_t)mesh->nblocks, sizeof(*boundary_mask));
    if (local_pos == 0 || local_neg == 0 || global_pos == 0 || global_neg == 0 ||
        boundary_mask == 0) {
        free(local_pos);
        free(local_neg);
        free(global_pos);
        free(global_neg);
        free(boundary_mask);
        fprintf(stderr, "prj_amr_tag: allocation failed\n");
        exit(EXIT_FAILURE);
    }

    /* Pass 1: per-block criterion vote, UNCAPPED.
     * Populates local_pos/local_neg/boundary_mask for blocks this rank owns. */
    for (i = 0; i < mesh->nblocks; ++i) {
        prj_block *b = &mesh->blocks[i];
        int refine = 0;
        int derefine = b->parent >= 0 ? 1 : 0;
        int j;
        int k;
        int ii;
        int has_refine_criterion = 0;
        int has_derefine_criterion = 0;
        int init_hook_refine = 0;

        if (!prj_is_local_active_block(mpi, b)) {
            continue;
        }

        if (mesh->amr_init_refine_fn != 0) {
            init_hook_refine = mesh->amr_init_refine_fn(b, mesh->amr_init_refine_userdata);
        }

        for (ii = 0; ii < PRJ_BLOCK_SIZE; ++ii) {
            for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    unsigned int cell_boundary_mask = prj_amr_refine_boundary_mask_for_cell(ii, j, k);
                    double refine_sum = 0.0;
                    int amr_idx;

                    if (refine != 0 && cell_boundary_mask == 0U) {
                        continue;
                    }
                    if (mesh->amr_reach_highest_level_at_density > 0.0 &&
                        b->W_mhd != 0 &&
                        b->W_mhd[WIDX(PRJ_PRIM_RHO, ii, j, k)] >
                        mesh->amr_reach_highest_level_at_density) {
                        has_refine_criterion = 1;
                        refine = 1;
                        derefine = 0;
                        boundary_mask[i] |= cell_boundary_mask;
                    }
                    for (amr_idx = 0; amr_idx < PRJ_AMR_N; ++amr_idx) {
                        if (mesh->amr_criterion_set[amr_idx] == 0) {
                            continue;
                        }
                        if (mesh->amr_refine_thresh[amr_idx] <= 0.0) {
                            continue;
                        }
                        has_refine_criterion = 1;
                        refine_sum += prj_amr_cell_indicator_for_estimator(
                            mesh, b, eos, amr_idx, mesh->amr_estimator[amr_idx], ii, j, k) /
                            mesh->amr_refine_thresh[amr_idx];
                    }

                    if (refine_sum > 1.0) {
                        refine = 1;
                        boundary_mask[i] |= cell_boundary_mask;
                    }
                }
            }
        }

        if (derefine != 0) {
            for (ii = -1; ii <= PRJ_BLOCK_SIZE && derefine != 0; ++ii) {
                for (j = -1; j <= PRJ_BLOCK_SIZE && derefine != 0; ++j) {
                    for (k = -1; k <= PRJ_BLOCK_SIZE; ++k) {
                        double derefine_sum = 0.0;
                        int amr_idx;

                        for (amr_idx = 0; amr_idx < PRJ_AMR_N; ++amr_idx) {
                            if (mesh->amr_criterion_set[amr_idx] == 0) {
                                continue;
                            }
                            if (mesh->amr_derefine_thresh[amr_idx] <= 0.0) {
                                continue;
                            }
                            has_derefine_criterion = 1;
                            derefine_sum += prj_amr_cell_indicator_for_estimator(
                                mesh, b, eos, amr_idx, mesh->amr_estimator[amr_idx], ii, j, k) /
                                mesh->amr_derefine_thresh[amr_idx];
                        }

                        if (derefine_sum >= 1.0) {
                            derefine = 0;
                            break;
                        }
                    }
                }
            }
        }

        if (has_refine_criterion == 0) {
            refine = 0;
        }
        if (has_derefine_criterion == 0) {
            derefine = 0;
        }

        if (init_hook_refine != 0) {
            refine = 1;
            derefine = 0;
        }

        if (refine != 0) {
            local_pos[i] = 1;
        }
        if (derefine != 0) {
            local_neg[i] = -1;
        }
    }

    /* Pass 2: boundary-mask propagation on local_pos.
     * Uses the uncapped vote so a max-level / can_refine=0 / min_dx-capped block
     * can still hand its refine signal to a coarser neighbor. */
    for (i = 0; i < mesh->nblocks; ++i) {
        prj_block *b = &mesh->blocks[i];

        if (boundary_mask[i] != 0U && prj_is_local_active_block(mpi, b)) {
            prj_amr_tag_boundary_neighbors(mesh, b, boundary_mask[i], local_pos);
        }
    }

    /* Pass 3: MPI reduce.
     * refine wins (+1 via MAX), derefine next (-1 via MIN), neutral otherwise.
     * After this every rank has the same global_pos / global_neg view. */
#if defined(PRJ_ENABLE_MPI)
    {
        if (mpi != 0 && mpi->totrank > 1) {
            MPI_Allreduce(local_pos, global_pos, mesh->nblocks, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            MPI_Allreduce(local_neg, global_neg, mesh->nblocks, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        } else {
            for (i = 0; i < mesh->nblocks; ++i) {
                global_pos[i] = local_pos[i];
                global_neg[i] = local_neg[i];
            }
        }
    }
#else
    for (i = 0; i < mesh->nblocks; ++i) {
        global_pos[i] = local_pos[i];
        global_neg[i] = local_neg[i];
    }
#endif

    /* Pass 4: combine + apply caps to produce the final refine_flag.
     * Caps live here so every rank decides the same outcome from the same data. */
    for (i = 0; i < mesh->nblocks; ++i) {
        prj_block *b = &mesh->blocks[i];

        if (!prj_is_active_block(b)) {
            b->refine_flag = 0;
            continue;
        }
        if (global_pos[i] > 0) {
            int allowed = (b->can_refine != 0) &&
                (mesh->max_level < 0 || b->level < mesh->max_level) &&
                (mesh->min_dx <= 0.0 || prj_block_cell_size(b) > mesh->min_dx) &&
                !prj_has_face_neighbor_coarser_than(mesh, b, b->level - 1);

            b->refine_flag = allowed ? 1 : 0;
        } else if (global_neg[i] < 0) {
            b->refine_flag = -1;
        } else {
            b->refine_flag = 0;
        }
    }

    free(local_pos);
    free(local_neg);
    free(global_pos);
    free(global_neg);
    free(boundary_mask);
}

#if PRJ_MHD
static void prj_amr_mhd_fail(const char *message)
{
    fprintf(stderr, "%s\n", message);
    exit(EXIT_FAILURE);
}

static void prj_amr_mhd_check_block(const prj_block *block, const char *caller)
{
    int d;

    if (block == 0 || !prj_block_has_cons_storage(block) || block->W_mhd == 0) {
        prj_amr_mhd_fail(caller);
    }
    for (d = 0; d < 3; ++d) {
        if (block->face_fidelity[d] == 0 || block->Bf[d] == 0) {
            prj_amr_mhd_fail(caller);
        }
    }
}

static int prj_amr_mhd_face_axis_max(int dir, int axis)
{
    return dir == axis ? PRJ_BLOCK_SIZE : PRJ_BLOCK_SIZE - 1;
}

static void prj_amr_mhd_clear_faces(prj_block *block, int use_bf1)
{
    int d;

    prj_amr_mhd_check_block(block, "prj_amr_mhd_clear_faces: missing MHD storage");
    for (d = 0; d < 3; ++d) {
        prj_fill(prj_block_bf_stage(block, d, use_bf1 != 0 ? 1 : 0),
            (size_t)PRJ_BLOCK_NFACES, 0.0);
    }
    for (d = 0; d < 3; ++d) {
        int n;

        for (n = 0; n < PRJ_BLOCK_NFACES; ++n) {
            block->face_fidelity[d][n] = PRJ_MHD_FIDELITY_NONE;
        }
    }
}

static void prj_amr_mhd_mark_active_faces(prj_block *block, int fidelity)
{
    int dir;

    prj_amr_mhd_check_block(block, "prj_amr_mhd_mark_active_faces: missing MHD storage");
    for (dir = 0; dir < 3; ++dir) {
        int i;
        int j;
        int k;

        for (i = 0; i <= prj_amr_mhd_face_axis_max(dir, 0); ++i) {
            for (j = 0; j <= prj_amr_mhd_face_axis_max(dir, 1); ++j) {
                for (k = 0; k <= prj_amr_mhd_face_axis_max(dir, 2); ++k) {
                    block->face_fidelity[dir][FACE_IDX(dir, i, j, k)] = fidelity;
                }
            }
        }
    }
}

static double prj_amr_mhd_cell_sqrt_gamma(const prj_mesh *mesh,
    const prj_block *block, int use_bf1, int i, int j, int k)
{
#if PRJ_DYNAMIC_GR
    if (prj_eos_full_dynamic_gr_enabled(mesh)) {
        prj_z4c_hydro_geom geom;

        if (!prj_z4c_load_hydro_geom(mesh, block,
                prj_stage_slot_from_bf_arg(use_bf1), i, j, k, &geom)) {
            prj_amr_mhd_fail("prj_amr_mhd_set_cons_b_from_bf: failed to load full-GR geometry");
        }
        return geom.sqrt_gamma;
    }
#else
    (void)mesh;
    (void)block;
    (void)use_bf1;
    (void)i;
    (void)j;
    (void)k;
#endif
    return 1.0;
}

static void prj_amr_mhd_set_cons_b_from_bf(const prj_mesh *mesh,
    prj_block *block, int use_bf1)
{
    double *bf[3];
    int d;
    int i;
    int j;
    int k;

    prj_amr_mhd_check_block(block, "prj_amr_mhd_set_cons_b_from_bf: missing MHD storage");
    for (d = 0; d < 3; ++d) {
        bf[d] = prj_block_bf_stage(block, d, use_bf1 != 0 ? 1 : 0);
    }
    for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
        for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
            for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                double b1 = 0.5 * (bf[X1DIR][FACE_IDX(X1DIR, i, j, k)] + bf[X1DIR][FACE_IDX(X1DIR, i + 1, j, k)]);
                double b2 = 0.5 * (bf[X2DIR][FACE_IDX(X2DIR, i, j, k)] + bf[X2DIR][FACE_IDX(X2DIR, i, j + 1, k)]);
                double b3 = 0.5 * (bf[X3DIR][FACE_IDX(X3DIR, i, j, k)] + bf[X3DIR][FACE_IDX(X3DIR, i, j, k + 1)]);
                double sqrt_gamma = prj_amr_mhd_cell_sqrt_gamma(mesh, block,
                    use_bf1, i, j, k);
                double b1_prim = b1 / sqrt_gamma;
                double b2_prim = b2 / sqrt_gamma;
                double b3_prim = b3 / sqrt_gamma;

                if (!isfinite(b1) || !isfinite(b2) || !isfinite(b3) ||
                    !isfinite(b1_prim) || !isfinite(b2_prim) ||
                    !isfinite(b3_prim)) {
                    prj_amr_mhd_fail("prj_amr_mhd_set_cons_b_from_bf: non-finite magnetic field");
                }
                prj_block_set_cons_value(block, PRJ_CONS_B1, i, j, k, b1);
                prj_block_set_cons_value(block, PRJ_CONS_B2, i, j, k, b2);
                prj_block_set_cons_value(block, PRJ_CONS_B3, i, j, k, b3);
                block->W_mhd[WIDX(PRJ_PRIM_B1, i, j, k)] = b1_prim;
                block->W_mhd[WIDX(PRJ_PRIM_B2, i, j, k)] = b2_prim;
                block->W_mhd[WIDX(PRJ_PRIM_B3, i, j, k)] = b3_prim;
                prj_block_set_prim_value(block, 1, PRJ_PRIM_B1, i, j, k, b1_prim);
                prj_block_set_prim_value(block, 1, PRJ_PRIM_B2, i, j, k, b2_prim);
                prj_block_set_prim_value(block, 1, PRJ_PRIM_B3, i, j, k, b3_prim);
            }
        }
    }
}

static void prj_amr_mhd_free_prolong_bf_buffer(double *buf[3])
{
    int dir;

    for (dir = 0; dir < 3; ++dir) {
        free(buf[dir]);
        buf[dir] = 0;
    }
}

static void prj_amr_mhd_pack_prolong_bf_buffer(const prj_block *parent,
    int ci0, int cj0, int ck0, int use_bf1, double *buf[3],
    int buf_lo[3][3], int buf_n[3][3])
{
    int coarse_lo[3];
    int dir;

    coarse_lo[0] = ci0;
    coarse_lo[1] = cj0;
    coarse_lo[2] = ck0;
    for (dir = 0; dir < 3; ++dir) {
        buf[dir] = 0;
    }
    for (dir = 0; dir < 3; ++dir) {
        const double *src = prj_block_bf_stage_const(parent, dir, use_bf1 != 0 ? 1 : 0);
        int count = 1;
        int d;
        int i;
        int j;
        int k;
        int idx;

        for (d = 0; d < 3; ++d) {
            if (d == dir) {
                buf_lo[dir][d] = coarse_lo[d];
                buf_n[dir][d] = PRJ_BLOCK_SIZE / 2 + 1;
            } else {
                buf_lo[dir][d] = coarse_lo[d] - 1;
                buf_n[dir][d] = PRJ_BLOCK_SIZE / 2 + 2;
            }
            count *= buf_n[dir][d];
        }

        buf[dir] = (double *)prj_malloc((size_t)count * sizeof(double));
        if (buf[dir] == 0) {
            prj_amr_mhd_free_prolong_bf_buffer(buf);
            prj_amr_mhd_fail("prj_amr_mhd_pack_prolong_bf_buffer: malloc failed");
        }

        idx = 0;
        for (i = buf_lo[dir][0]; i < buf_lo[dir][0] + buf_n[dir][0]; ++i) {
            for (j = buf_lo[dir][1]; j < buf_lo[dir][1] + buf_n[dir][1]; ++j) {
                for (k = buf_lo[dir][2]; k < buf_lo[dir][2] + buf_n[dir][2]; ++k) {
                    double value = src[FACE_IDX(dir, i, j, k)];

                    if (!isfinite(value)) {
                        prj_amr_mhd_free_prolong_bf_buffer(buf);
                        prj_amr_mhd_fail("prj_amr_mhd_pack_prolong_bf_buffer: non-finite face-centered magnetic field");
                    }
                    buf[dir][idx++] = value;
                }
            }
        }
    }
}

static void prj_amr_mhd_prolongate_bf_one(const prj_mesh *mesh, const prj_mpi *mpi,
    const prj_block *parent, prj_block *child,
    int child_oct, int use_bf1)
{
    int xoct = child_oct & 1;
    int yoct = (child_oct >> 1) & 1;
    int zoct = (child_oct >> 2) & 1;
    int ci0 = xoct * (PRJ_BLOCK_SIZE / 2);
    int cj0 = yoct * (PRJ_BLOCK_SIZE / 2);
    int ck0 = zoct * (PRJ_BLOCK_SIZE / 2);
    int ci;
    int cj;
    int ck;
    double *buf[3] = {0, 0, 0};
    const double *cbuf[3];
    int buf_lo[3][3];
    int buf_n[3][3];
    
    prj_amr_mhd_check_block(parent, "prj_amr_mhd_prolongate_bf_one: missing parent MHD storage");
    prj_amr_mhd_clear_faces(child, use_bf1);


    int n,dir;
    for (n = 0; n < 56; ++n){
        const prj_neighbor *slot = &parent->slot[n];
        if (slot->id<0||slot->type!=PRJ_NEIGHBOR_FACE||slot->rel_level<1) {continue;}
        if (slot->rank == parent->rank) {
            for (dir = 0; dir < 3; ++dir){
                int tan0 = (dir+1)%3;
                int tan1 = (dir+2)%3;
                if (fabs(child->xmin[tan0] - slot->xmin[tan0]) < 1.0e-12*child->dx[tan0]){
                    if (fabs(child->xmin[tan1] - slot->xmin[tan1]) < 1.0e-12*child->dx[tan1]){
                        const prj_block *neighbor = &mesh->blocks[slot->id];
                        const double *src = prj_block_bf_stage_const(neighbor, dir,
                            use_bf1 != 0 ? 1 : 0);
                        int i, j;
                        if (fabs(child->xmax[dir] - slot->xmin[dir]) < 1.0e-12*child->dx[dir]){
                            double buffer[PRJ_BLOCK_SIZE*PRJ_BLOCK_SIZE];
                            int buf_idx = 0;
                            for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
                                for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                                    int it_send[3] = {0,0,0};
                                    it_send[dir] = 0;
                                    it_send[tan0] = i;
                                    it_send[tan1] = j;
                                    buffer[buf_idx++] = src[FACE_IDX(dir, it_send[0], it_send[1], it_send[2])];
                                }
                            }


                            buf_idx = 0;
                            for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
                                for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                                    int it_recv[3] = {0,0,0};
                                    it_recv[dir] = PRJ_BLOCK_SIZE;
                                    it_recv[tan0] = i;
                                    it_recv[tan1] = j;
                                    prj_boundary_write_bf_face(child, use_bf1, dir,
                                                               it_recv[0],
                                                               it_recv[1],
                                                               it_recv[2],
                                                               buffer[buf_idx++], PRJ_MHD_FIDELITY_SAME);
                                }
                            }
                        }
                        if (fabs(child->xmin[dir] - slot->xmax[dir]) < 1.0e-12*child->dx[dir]){
                            double buffer[PRJ_BLOCK_SIZE*PRJ_BLOCK_SIZE];
                            int buf_idx = 0;
                            for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
                                for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                                    int it_send[3] = {0,0,0};
                                    it_send[dir] = PRJ_BLOCK_SIZE;
                                    it_send[tan0] = i;
                                    it_send[tan1] = j;
                                    buffer[buf_idx++] = src[FACE_IDX(dir, it_send[0], it_send[1], it_send[2])];
                                }
                            }

                            buf_idx = 0;
                            for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
                                for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                                    int it_recv[3] = {0,0,0};
                                    it_recv[dir] = 0;
                                    it_recv[tan0] = i;
                                    it_recv[tan1] = j;
                                    prj_boundary_write_bf_face(child, use_bf1, dir,
                                                               it_recv[0],
                                                               it_recv[1],
                                                               it_recv[2],
                                                               buffer[buf_idx++], PRJ_MHD_FIDELITY_SAME);
                                }
                            }
                        }
                    }
                }
            }

        } else {
            for (dir = 0; dir < 3; ++dir) {
                prj_mpi_amr_mhd_prolongate_bf_one(mpi, parent, slot, child,
                    child_oct, use_bf1, dir);
            }
        }
    }


    prj_amr_mhd_pack_prolong_bf_buffer(parent, ci0, cj0, ck0, use_bf1,
        buf, buf_lo, buf_n);
    cbuf[0] = buf[0];
    cbuf[1] = buf[1];
    cbuf[2] = buf[2];

    for (ci = ci0; ci < ci0 + PRJ_BLOCK_SIZE / 2; ++ci) {
        for (cj = cj0; cj < cj0 + PRJ_BLOCK_SIZE / 2; ++cj) {
            for (ck = ck0; ck < ck0 + PRJ_BLOCK_SIZE / 2; ++ck) {
                int fi = 2 * (ci - ci0);
                int fj = 2 * (cj - cj0);
                int fk = 2 * (ck - ck0);

                prj_mhd_prolong_bf_from_buffer(cbuf, buf_lo, buf_n,
                    parent->dx, child, ci, cj, ck, fi, fj, fk, use_bf1,
                    mesh != 0 && mesh->use_BJ != 0);
            }
        }
    }

    prj_amr_mhd_free_prolong_bf_buffer(buf);
}

static void prj_amr_mhd_prolongate_bf(const prj_mesh *mesh, const prj_mpi *mpi,
    const prj_block *parent, prj_block *child, int child_oct)
{
    prj_amr_mhd_prolongate_bf_one(mesh, mpi, parent, child, child_oct, 0);
    prj_amr_mhd_prolongate_bf_one(mesh, mpi, parent, child, child_oct, 1);
    prj_amr_mhd_mark_active_faces(child, PRJ_MHD_FIDELITY_COARSER);
    prj_amr_mhd_set_cons_b_from_bf(mesh, child, 0);
}

static const prj_block *prj_amr_mhd_child_for_face(const prj_block *children[8],
    int bit[3])
{
    int oct = bit[0] + 2 * bit[1] + 4 * bit[2];

    if (oct < 0 || oct >= 8 || children[oct] == 0) {
        prj_amr_mhd_fail("prj_amr_mhd_child_for_face: missing child block");
    }
    return children[oct];
}

static int prj_amr_mhd_face_bit(int global_face_index)
{
    if (global_face_index < PRJ_BLOCK_SIZE) {
        return 0;
    }
    if (global_face_index == PRJ_BLOCK_SIZE) {
        return 0;
    }
    return 1;
}

static int prj_amr_mhd_face_local_index(int global_face_index)
{
    if (global_face_index < PRJ_BLOCK_SIZE) {
        return global_face_index;
    }
    if (global_face_index == PRJ_BLOCK_SIZE) {
        return PRJ_BLOCK_SIZE;
    }
    return global_face_index - PRJ_BLOCK_SIZE;
}

static double prj_amr_mhd_restrict_face_value(const prj_block *children[8],
    int dir, int i, int j, int k, int use_bf1)
{
    int tan0 = (dir + 1) % 3;
    int tan1 = (dir + 2) % 3;
    int coarse_idx[3];
    double sum = 0.0;
    int a;
    int b;

    coarse_idx[0] = i;
    coarse_idx[1] = j;
    coarse_idx[2] = k;
    for (a = 0; a < 2; ++a) {
        for (b = 0; b < 2; ++b) {
            int global[3];
            int local[3];
            int bit[3];
            const prj_block *child;
            const double *src;
            double value;

            global[0] = 2 * coarse_idx[0];
            global[1] = 2 * coarse_idx[1];
            global[2] = 2 * coarse_idx[2];
            global[tan0] += a;
            global[tan1] += b;
            bit[dir] = prj_amr_mhd_face_bit(global[dir]);
            local[dir] = prj_amr_mhd_face_local_index(global[dir]);
            bit[tan0] = global[tan0] / PRJ_BLOCK_SIZE;
            bit[tan1] = global[tan1] / PRJ_BLOCK_SIZE;
            local[tan0] = global[tan0] - bit[tan0] * PRJ_BLOCK_SIZE;
            local[tan1] = global[tan1] - bit[tan1] * PRJ_BLOCK_SIZE;
            child = prj_amr_mhd_child_for_face(children, bit);
            src = prj_block_bf_stage_const(child, dir, use_bf1 != 0 ? 1 : 0);
            if (src == 0) {
                prj_amr_mhd_fail("prj_amr_mhd_restrict_face_value: missing child Bf storage");
            }
            value = src[FACE_IDX(dir, local[0], local[1], local[2])];
            if (!isfinite(value)) {
                prj_amr_mhd_fail("prj_amr_mhd_restrict_face_value: non-finite child Bf");
            }
            sum += value;
        }
    }
    return 0.25 * sum;
}

static void prj_amr_mhd_restrict_bf_one(const prj_block *children[8],
    prj_block *parent, int use_bf1)
{
    int dir;

    prj_amr_mhd_clear_faces(parent, use_bf1);
    for (dir = 0; dir < 3; ++dir) {
        double *dst = prj_block_bf_stage(parent, dir, use_bf1 != 0 ? 1 : 0);
        int i;
        int j;
        int k;

        for (i = 0; i <= prj_amr_mhd_face_axis_max(dir, 0); ++i) {
            for (j = 0; j <= prj_amr_mhd_face_axis_max(dir, 1); ++j) {
                for (k = 0; k <= prj_amr_mhd_face_axis_max(dir, 2); ++k) {
                    dst[FACE_IDX(dir, i, j, k)] = prj_amr_mhd_restrict_face_value(
                        children, dir, i, j, k, use_bf1);
                }
            }
        }
    }
}

static void prj_amr_mhd_restrict_bf(const prj_mesh *mesh,
    const prj_block *children[8], prj_block *parent)
{
    int oct;

    prj_amr_mhd_check_block(parent, "prj_amr_mhd_restrict_bf: missing parent MHD storage");
    for (oct = 0; oct < 8; ++oct) {
        prj_amr_mhd_check_block(children[oct], "prj_amr_mhd_restrict_bf: missing child MHD storage");
    }
    prj_amr_mhd_restrict_bf_one(children, parent, 0);
    prj_amr_mhd_restrict_bf_one(children, parent, 1);
    prj_amr_mhd_mark_active_faces(parent, PRJ_MHD_FIDELITY_FINER);
    prj_amr_mhd_set_cons_b_from_bf(mesh, parent, 0);
}
#endif

void prj_amr_prolongate(const prj_mesh *mesh, const prj_mpi *mpi, const prj_block *parent,
    prj_block *child, int child_oct)
{
    int i;
    int j;
    int k;
    int v;
    int xoct = child_oct & 1;
    int yoct = (child_oct >> 1) & 1;
    int zoct = (child_oct >> 2) & 1;
    int use_BJ = mesh != 0 && mesh->use_BJ != 0;

#if !PRJ_MHD
    (void)mpi;
#endif
    prj_zero_block_arrays(child);

    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
            for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    int gi = xoct * PRJ_BLOCK_SIZE + i;
                    int gj = yoct * PRJ_BLOCK_SIZE + j;
                    int gk = zoct * PRJ_BLOCK_SIZE + k;
                    int ip = gi / 2;
                    int jp = gj / 2;
                    int kp = gk / 2;
                    double stencil[27];
                    double target[3];
                    int di;
                    int dj;
                    int dk;

                    for (di = -1; di <= 1; ++di) {
                        for (dj = -1; dj <= 1; ++dj) {
                            for (dk = -1; dk <= 1; ++dk) {
                                stencil[prj_reconstruct_stencil3_index(di, dj, dk)] =
                                    prj_block_conserved_at(parent, v, ip + di, jp + dj, kp + dk);
                            }
                        }
                    }
                    target[0] = ((gi % 2) == 0) ? -0.25 : 0.25;
                    target[1] = ((gj % 2) == 0) ? -0.25 : 0.25;
                    target[2] = ((gk % 2) == 0) ? -0.25 : 0.25;
                    prj_block_set_cons_value(child, v, i, j, k,
                        prj_reconstruct_cell_for_prolongate(stencil, target, use_BJ));
                }
            }
        }
    }
    /* Guard against a nonpositive prolongated density; the internal-energy
     * floor is applied later in prj_sync_primitive_from_conserved. */
    for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
        for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
            for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                if (prj_block_cons_value_const(child, PRJ_CONS_RHO, i, j, k) <= 0.0) {
                    prj_block_set_cons_value(child, PRJ_CONS_RHO, i, j, k, 1.0e-10);
                }
            }
        }
    }
#if PRJ_MHD
    prj_amr_mhd_prolongate_bf(mesh, mpi, parent, child, child_oct);
#endif
#if PRJ_DYNAMIC_GR
    if (prj_z4c_runtime_enabled(mesh)) {
        prj_z4c_amr_prolongate_child(parent, child, child_oct);
    }
#endif
}

void prj_amr_restrict(const prj_mesh *mesh, const prj_block *children[8], prj_block *parent)
{
    int v;
    int i;
    int j;
    int k;

#if !PRJ_DYNAMIC_GR
    (void)mesh;
#endif
    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
            for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    double sum = 0.0;
                    double parent_cell_vol = parent->dx[0] * parent->dx[1] * parent->dx[2];
                    double parent_xmin = parent->xmin[0] + (double)i * parent->dx[0];
                    double parent_xmax = parent_xmin + parent->dx[0];
                    double parent_ymin = parent->xmin[1] + (double)j * parent->dx[1];
                    double parent_ymax = parent_ymin + parent->dx[1];
                    double parent_zmin = parent->xmin[2] + (double)k * parent->dx[2];
                    double parent_zmax = parent_zmin + parent->dx[2];
                    int oct;

                    for (oct = 0; oct < 8; ++oct) {
                        const prj_block *child = children[oct];
                        int fi;
                        int fj;
                        int fk;

                        for (fi = 0; fi < PRJ_BLOCK_SIZE; ++fi) {
                            double fine_xmin = child->xmin[0] + (double)fi * child->dx[0];
                            double fine_xmax = fine_xmin + child->dx[0];
                            double xover = prj_min_double(parent_xmax, fine_xmax) -
                                prj_max_double(parent_xmin, fine_xmin);

                            if (xover <= 0.0) {
                                continue;
                            }
                            for (fj = 0; fj < PRJ_BLOCK_SIZE; ++fj) {
                                double fine_ymin = child->xmin[1] + (double)fj * child->dx[1];
                                double fine_ymax = fine_ymin + child->dx[1];
                                double yover = prj_min_double(parent_ymax, fine_ymax) -
                                    prj_max_double(parent_ymin, fine_ymin);

                                if (yover <= 0.0) {
                                    continue;
                                }
                                for (fk = 0; fk < PRJ_BLOCK_SIZE; ++fk) {
                                    double fine_zmin = child->xmin[2] + (double)fk * child->dx[2];
                                    double fine_zmax = fine_zmin + child->dx[2];
                                    double zover = prj_min_double(parent_zmax, fine_zmax) -
                                        prj_max_double(parent_zmin, fine_zmin);

                                    if (zover > 0.0) {
                                        double overlap_vol = xover * yover * zover;

                                        sum += prj_block_cons_value_const(child, v, fi, fj, fk) *
                                            overlap_vol;
                                    }
                                }
                            }
                        }
                    }
                    prj_block_set_cons_value(parent, v, i, j, k, sum / parent_cell_vol);
                }
            }
        }
    }
#if PRJ_MHD
    prj_amr_mhd_restrict_bf(mesh, children, parent);
#endif
#if PRJ_DYNAMIC_GR
    if (prj_z4c_runtime_enabled(mesh)) {
        prj_z4c_amr_restrict_parent(children, parent);
    }
#endif
}

void prj_amr_refine_block(prj_mesh *mesh, const prj_mpi *mpi, int block_id)
{
    prj_block *parent;
    int oct;
    double xmid[3];
    int owner_local;

    if (mesh == 0 || block_id < 0 || block_id >= mesh->nblocks) {
        return;
    }
    parent = &mesh->blocks[block_id];
    if (!prj_is_active_block(parent) ||
        (mesh->max_level >= 0 && parent->level >= mesh->max_level)) {
        return;
    }
    if (parent->can_refine == 0) {
        parent->refine_flag = 0;
        return;
    }
    if (mesh->min_dx > 0.0 && prj_block_cell_size(parent) <= mesh->min_dx) {
        parent->refine_flag = 0;
        return;
    }
    owner_local = prj_is_local_block_owner(mpi, parent);

    xmid[0] = 0.5 * (parent->xmin[0] + parent->xmax[0]);
    xmid[1] = 0.5 * (parent->xmin[1] + parent->xmax[1]);
    xmid[2] = 0.5 * (parent->xmin[2] + parent->xmax[2]);

    for (oct = 0; oct < 8; ++oct) {
        int id = prj_find_free_block_slot(mesh);
        prj_block *child;
        int bitx;
        int bity;
        int bitz;

        if (id < 0) {
            return;
        }
        child = &mesh->blocks[id];
        prj_block_free_data(child);
        child->id = id;
        child->rank = parent->rank;
        child->level = parent->level + 1;
        child->active = 1;
        child->refine_flag = 0;
        child->can_refine = 1;
        child->parent = parent->id;
        prj_reset_children(child);
        prj_clear_neighbors(child);
        bitx = oct & 1;
        bity = (oct >> 1) & 1;
        bitz = (oct >> 2) & 1;
        child->xmin[0] = bitx == 0 ? parent->xmin[0] : xmid[0];
        child->xmax[0] = bitx == 0 ? xmid[0] : parent->xmax[0];
        child->xmin[1] = bity == 0 ? parent->xmin[1] : xmid[1];
        child->xmax[1] = bity == 0 ? xmid[1] : parent->xmax[1];
        child->xmin[2] = bitz == 0 ? parent->xmin[2] : xmid[2];
        child->xmax[2] = bitz == 0 ? xmid[2] : parent->xmax[2];
        child->dx[0] = parent->dx[0] * 0.5;
        child->dx[1] = parent->dx[1] * 0.5;
        child->dx[2] = parent->dx[2] * 0.5;
        prj_block_update_can_refine(child, mesh);
        if (owner_local) {
            if (prj_block_alloc_data(child) != 0) {
                child->id = -1;
                child->active = 0;
                return;
            }
            prj_block_setup_geometry(child, &mesh->coord);
            prj_mesh_update_block_r_com(child, mesh);
            prj_amr_prolongate(mesh, mpi, parent, child, oct);
        } else {
            child->W_mhd = 0;
            child->eosvar = 0;
            child->U_mhd = 0;
            child->U_rad = 0;
            child->mhd_rhs = 0;
            child->rad_rhs = 0;
#if TIME_INTEGRATION == PRJ_TIMEINT_IMEX
            child->deriv_ex = 0;
            child->deriv_im = 0;
#endif
            child->flux[0] = 0;
            child->flux[1] = 0;
            child->flux[2] = 0;
            child->v_riemann[0] = 0;
            child->v_riemann[1] = 0;
            child->v_riemann[2] = 0;
            child->kappa_cell = 0;
            child->sigma_cell = 0;
            child->lapse = 0;
            child->grav[0] = 0;
            child->grav[1] = 0;
            child->grav[2] = 0;
            child->r_com = 0;
            child->Ylm = 0;
#if PRJ_USE_RADIATION_FSA && PRJ_USE_RADIAL_FRAME_FSA
            child->rotation_matrix_fsa = 0;
            child->ang_geom_fsa = 0;
#endif
#if PRJ_MHD
            child->deriv_Bf[0] = 0;
            child->deriv_Bf[1] = 0;
            child->deriv_Bf[2] = 0;
#endif
            child->ridx = 0;
            child->fr = 0;
            prj_block_setup_geometry(child, &mesh->coord);
        }
        parent->children[oct] = id;
    }

    parent->active = 0;
    parent->refine_flag = 0;
}

int prj_amr_coarsen_block(prj_mesh *mesh, const prj_mpi *mpi, int parent_id)
{
    prj_block *parent;
    const prj_block *children[8];
    int child_ranks[8];
    int child_counts[8];
    int oct;
    int owner_local;

    if (mesh == 0 || parent_id < 0 || parent_id >= mesh->nblocks) {
        return 0;
    }
    parent = &mesh->blocks[parent_id];
    for (oct = 0; oct < 8; ++oct) {
        int id = parent->children[oct];

        if (id < 0 || id >= mesh->nblocks || !prj_is_active_block(&mesh->blocks[id])) {
            return 0;
        }
        children[oct] = &mesh->blocks[id];
        child_ranks[oct] = mesh->blocks[id].rank;
        child_counts[oct] = 0;
    }

    for (oct = 0; oct < 8; ++oct) {
        int other;

        for (other = 0; other < 8; ++other) {
            if (child_ranks[other] == child_ranks[oct]) {
                child_counts[oct] += 1;
            }
        }
    }
    parent->rank = child_ranks[0];
    for (oct = 1; oct < 8; ++oct) {
        if (child_counts[oct] > child_counts[0]) {
            child_counts[0] = child_counts[oct];
            parent->rank = child_ranks[oct];
        }
    }
    owner_local = prj_is_local_block_owner(mpi, parent);

    prj_amr_move_children_to_parent_rank(mesh, mpi, parent);
    if (owner_local) {
        if (parent->W_mhd == 0) {
            prj_block_alloc_data(parent);
            prj_block_setup_geometry(parent, &mesh->coord);
            prj_mesh_update_block_r_com(parent, mesh);
        }
        prj_amr_restrict(mesh, children, parent);
    }
    parent->active = 1;
    parent->refine_flag = 0;
    for (oct = 0; oct < 8; ++oct) {
        prj_block *child = &mesh->blocks[parent->children[oct]];

        if (owner_local) {
            prj_block_free_data(child);
        }
        child->id = -1;
        child->active = 0;
        child->can_refine = 1;
        child->parent = -1;
        prj_reset_children(child);
        prj_clear_neighbors(child);
        parent->children[oct] = -1;
    }
    return 1;
}

int prj_amr_adapt(prj_mesh *mesh, prj_eos *eos, prj_mpi *mpi)
{
    int i;
    int refined = 0;
    int coarsened = 0;
    int pending_change = 0;
    int *neighbor_dirty = 0;
    int *changed_blocks = 0;

    if (mesh == 0) {
        return 0;
    }
    if (mesh->nblocks_max > 0) {
        neighbor_dirty = (int *)prj_calloc((size_t)mesh->nblocks_max,
            sizeof(*neighbor_dirty));
        changed_blocks = (int *)prj_calloc((size_t)mesh->nblocks_max,
            sizeof(*changed_blocks));
    }

    prj_amr_tag(mesh, eos, mpi);
    prj_amr_enforce_two_to_one(mesh, mpi);
    prj_amr_sync_refine_flags(mesh, mpi);

    pending_change = prj_amr_has_pending_change(mesh);
    if (!pending_change) {
        for (i = 0; i < mesh->nblocks; ++i) {
            if (mesh->blocks[i].id >= 0) {
                mesh->blocks[i].refine_flag = 0;
            }
        }
        free(neighbor_dirty);
        free(changed_blocks);
        return 0;
    }

    prj_sync_conserved_from_primitive(mesh, eos, mpi);

    prj_eos_fill_ghost_cons(mesh, eos, mpi, 1, PRJ_EOS_CTX_AMR);

#if PRJ_MHD
    prj_mpi_exchange_amr_mhd_prolongate_bf(mesh, mpi);
#endif
    refined = prj_amr_refine_marked_blocks_with_dirty(mesh, mpi, neighbor_dirty, changed_blocks);
    if (refined) {
        if (neighbor_dirty != 0) {
            prj_amr_init_neighbors_with_mask(mesh, neighbor_dirty);
            prj_amr_clear_dirty_mask(neighbor_dirty, mesh->nblocks_max);
        } else {
            prj_amr_init_neighbors(mesh);
        }
    }

    for (i = 0; i < mesh->nblocks; ++i) {
        prj_block *parent = &mesh->blocks[i];
        int can_coarsen;

        if (parent->id < 0 || parent->active == 1) {
            continue;
        }
        can_coarsen = prj_can_coarsen_parent(mesh, i);
        if (can_coarsen) {
            prj_amr_mark_coarsen_dirty(mesh, i, neighbor_dirty);
            if (prj_amr_coarsen_block(mesh, mpi, i)) {
                if (neighbor_dirty != 0) {
                    neighbor_dirty[i] = 1;
                }
                /* The newly activated parent holds freshly restricted data
                 * and is the only block from this coarsen that needs a resync. */
                if (changed_blocks != 0) {
                    changed_blocks[i] = 1;
                }
                coarsened = 1;
            }
        }
    }
    if (coarsened) {
        if (neighbor_dirty != 0) {
            prj_amr_init_neighbors_with_mask(mesh, neighbor_dirty);
        } else {
            prj_amr_init_neighbors(mesh);
        }
    }

    if (!refined && !coarsened) {
        for (i = 0; i < mesh->nblocks; ++i) {
            if (mesh->blocks[i].id >= 0) {
                mesh->blocks[i].refine_flag = 0;
            }
        }
        free(neighbor_dirty);
        free(changed_blocks);
        return 0;
    }

    for (i = 0; i < mesh->nblocks; ++i) {
        if (mesh->blocks[i].id >= 0) {
            prj_block_setup_geometry(&mesh->blocks[i], &mesh->coord);
            mesh->blocks[i].refine_flag = 0;
        }
    }
    prj_mesh_update_r_com(mesh);
    prj_mesh_update_max_active_level(mesh);
    {
        double injected_local = 0.0;
        int global_changed;
        double injected_global;

        prj_sync_primitive_from_conserved(mesh, eos, mpi, changed_blocks, &injected_local);

        /* E_injected is only synchronized when the floor actually fired on some
         * rank; the cheap flag reduce decides that, the energy sum runs only then. */
        global_changed = (injected_local != 0.0);
        injected_global = injected_local;
#if defined(PRJ_ENABLE_MPI)
        if (mpi != 0 && mpi->totrank > 1) {
            int local_changed = global_changed;

            MPI_Allreduce(&local_changed, &global_changed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            if (global_changed) {
                MPI_Allreduce(&injected_local, &injected_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            }
        }
#endif
        if (global_changed) {
            eos->E_injected += injected_global;
        }
    }
    free(neighbor_dirty);
    free(changed_blocks);
    return 1;
}
