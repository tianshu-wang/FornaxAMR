#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "prj.h"

#if defined(PRJ_ENABLE_MPI)
#include <mpi.h>
#endif

static prj_mpi *prj_mpi_active = 0;

#define PRJ_MPI_GHOST_NVAR (PRJ_NVAR_PRIM + PRJ_NVAR_EOSVAR)
#define PRJ_MPI_SAMPLE_SAME_LEVEL_OFFSET 4

static int prj_block_is_active(const prj_block *block)
{
    return block != 0 && block->id >= 0 && block->active == 1;
}

static int prj_mpi_block_is_local(const prj_block *block)
{
    return prj_block_is_active(block) &&
        (prj_mpi_active == 0 || block->rank == prj_mpi_active->rank);
}

static double *prj_mpi_stage_array(prj_block *block, int stage)
{
    return stage == 2 ? block->W1 : block->W;
}

static void prj_mpi_buffer_free(prj_mpi_buffer *buffer)
{
    int axis;

    if (buffer == 0) {
        return;
    }
    free(buffer->receiver_blocks);
    free(buffer->cell_data_size_send);
    free(buffer->cell_buffer_send);
    free(buffer->face_data_size_send);
    free(buffer->face_buffer_send);
    free(buffer->cell_data_size_recv);
    free(buffer->cell_buffer_recv);
    free(buffer->face_data_size_recv);
    free(buffer->face_buffer_recv);
    for (axis = 0; axis < 3; ++axis) {
        free(buffer->cell_data_idx_send[axis]);
        free(buffer->face_data_idx_send[axis]);
        free(buffer->cell_data_idx_recv[axis]);
        free(buffer->face_data_idx_recv[axis]);
    }
    memset(buffer, 0, sizeof(*buffer));
}

static void prj_mpi_clear_neighbors(prj_mpi *mpi)
{
    int i;

    if (mpi == 0) {
        return;
    }
    for (i = 0; i < mpi->neighbor_number; ++i) {
        prj_mpi_buffer_free(&mpi->neighbor_buffer[i]);
    }
    free(mpi->neighbor_buffer);
    mpi->neighbor_buffer = 0;
    mpi->neighbor_number = 0;
}

static unsigned long long prj_mpi_spread_bits(unsigned int value)
{
    unsigned long long x = value;

    x = (x | (x << 16)) & 0x0000FFFF0000FFFFULL;
    x = (x | (x << 8)) & 0x00FF00FF00FF00FFULL;
    x = (x | (x << 4)) & 0x0F0F0F0F0F0F0F0FULL;
    x = (x | (x << 2)) & 0x3333333333333333ULL;
    x = (x | (x << 1)) & 0x5555555555555555ULL;
    return x;
}

static unsigned long long prj_mpi_morton_code(const prj_mesh *mesh, const prj_block *block)
{
    double scale_x;
    double scale_y;
    double scale_z;
    unsigned int ix;
    unsigned int iy;
    unsigned int iz;
    unsigned long long level_bits;
    unsigned int level_scale;

    level_scale = 1U << (unsigned int)block->level;
    scale_x = ((double)mesh->root_nx[0] * (double)level_scale) /
        (mesh->coord.x1max - mesh->coord.x1min);
    scale_y = ((double)mesh->root_nx[1] * (double)level_scale) /
        (mesh->coord.x2max - mesh->coord.x2min);
    scale_z = ((double)mesh->root_nx[2] * (double)level_scale) /
        (mesh->coord.x3max - mesh->coord.x3min);
    ix = (unsigned int)(((block->xmin[0] - mesh->coord.x1min) * scale_x) + 0.5);
    iy = (unsigned int)(((block->xmin[1] - mesh->coord.x2min) * scale_y) + 0.5);
    iz = (unsigned int)(((block->xmin[2] - mesh->coord.x3min) * scale_z) + 0.5);
    level_bits = (unsigned long long)(unsigned int)block->level << 48;
    return level_bits |
        (prj_mpi_spread_bits(ix) << 2) |
        (prj_mpi_spread_bits(iy) << 1) |
        prj_mpi_spread_bits(iz);
}

static int prj_mpi_encode_cell_index(int i, int j, int k)
{
    return ((i + PRJ_NGHOST) * PRJ_BS + (j + PRJ_NGHOST)) * PRJ_BS + (k + PRJ_NGHOST);
}

static void prj_mpi_decode_cell_index(int code, int *i, int *j, int *k)
{
    *k = code % PRJ_BS - PRJ_NGHOST;
    code /= PRJ_BS;
    *j = code % PRJ_BS - PRJ_NGHOST;
    code /= PRJ_BS;
    *i = code - PRJ_NGHOST;
}

static int prj_mpi_cell_point_inside(const prj_block *block, double x1, double x2, double x3)
{
    const double tol = 1.0e-12;

    return x1 > block->xmin[0] - tol && x1 < block->xmax[0] + tol &&
        x2 > block->xmin[1] - tol && x2 < block->xmax[1] + tol &&
        x3 > block->xmin[2] - tol && x3 < block->xmax[2] + tol;
}

static int prj_mpi_floor_to_int(double x)
{
    int i = (int)x;

    if ((double)i > x) {
        i -= 1;
    }
    return i;
}

#if PRJ_MHD
static int prj_mpi_min_int(int a, int b)
{
    return a < b ? a : b;
}

static int prj_mpi_max_int(int a, int b)
{
    return a > b ? a : b;
}

static int prj_mpi_valid_restrict_count(int count)
{
    return count == 0 || count == 4 || count == 8;
}
#endif

static int prj_mpi_abs_kind(double x, double target, double tol)
{
    return fabs(x - target) < tol;
}

static int prj_mpi_fraction_case(double frac)
{
    const double tol = 1.0e-2;

    if (prj_mpi_abs_kind(frac, 0.0, tol) || prj_mpi_abs_kind(frac, 1.0, tol)) {
        return 1;
    }
    if (prj_mpi_abs_kind(frac, 0.5, tol)) {
        return 2;
    }
    if (prj_mpi_abs_kind(frac, 0.25, tol) || prj_mpi_abs_kind(frac, 0.75, tol)) {
        return 3;
    }
    return 0;
}

static int prj_mpi_sample_kind(const prj_block *block, double x1, double x2, double x3)
{
    double ox[3];
    double frac[3];
    int cases[3];

    if (block == 0) {
        return -1;
    }
    ox[0] = (x1 - block->xmin[0]) / block->dx[0] - 0.5;
    ox[1] = (x2 - block->xmin[1]) / block->dx[1] - 0.5;
    ox[2] = (x3 - block->xmin[2]) / block->dx[2] - 0.5;
    frac[0] = ox[0] - (double)prj_mpi_floor_to_int(ox[0]);
    frac[1] = ox[1] - (double)prj_mpi_floor_to_int(ox[1]);
    frac[2] = ox[2] - (double)prj_mpi_floor_to_int(ox[2]);
    if (frac[0] < 0.0) {
        frac[0] += 1.0;
    }
    if (frac[1] < 0.0) {
        frac[1] += 1.0;
    }
    if (frac[2] < 0.0) {
        frac[2] += 1.0;
    }

    cases[0] = prj_mpi_fraction_case(frac[0]);
    cases[1] = prj_mpi_fraction_case(frac[1]);
    cases[2] = prj_mpi_fraction_case(frac[2]);
    if (cases[0] == 0 || cases[1] == 0 || cases[2] == 0) {
        return -1;
    }
    if (cases[0] == 1 && cases[1] == 1 && cases[2] == 1) {
        return 0;
    }
    if (cases[0] == 2 && cases[1] == 2 && cases[2] == 2) {
        return 0;
    }
    if (cases[0] == 3 && cases[1] == 3 && cases[2] == 3) {
        return 1;
    }
    fprintf(stderr,
        "prj_mpi_sample_kind: invalid fraction case tuple (%d,%d,%d) "
        "at (x1=%g, x2=%g, x3=%g) -- violates 2:1 AMR constraint\n",
        cases[0], cases[1], cases[2], x1, x2, x3);
    abort();
}

static int prj_mpi_sample_code(int sample_kind, int same_level)
{
    return sample_kind + (same_level != 0 ? PRJ_MPI_SAMPLE_SAME_LEVEL_OFFSET : 0);
}

static int prj_mpi_sample_kind_from_code(int sample_code)
{
    if (sample_code >= PRJ_MPI_SAMPLE_SAME_LEVEL_OFFSET) {
        return sample_code - PRJ_MPI_SAMPLE_SAME_LEVEL_OFFSET;
    }
    return sample_code;
}

static int prj_mpi_sample_is_same_level(int sample_code)
{
    return sample_code >= PRJ_MPI_SAMPLE_SAME_LEVEL_OFFSET;
}

static int prj_mpi_append_double(double **array, int *count, int *capacity, double value);

#if PRJ_MHD
#define PRJ_MPI_BF_PATCH_VALUE_COUNT 54

static int prj_mpi_encode_face_code(int dir, int i, int j, int k)
{
    return dir * PRJ_BLOCK_NCELLS + prj_mpi_encode_cell_index(i, j, k);
}

static void prj_mpi_decode_face_code(int code, int *dir, int *i, int *j, int *k)
{
    if (dir == 0 || i == 0 || j == 0 || k == 0) {
        fprintf(stderr, "prj_mpi_decode_face_code: null output pointer\n");
        abort();
    }
    *dir = code / PRJ_BLOCK_NCELLS;
    prj_mpi_decode_cell_index(code % PRJ_BLOCK_NCELLS, i, j, k);
}

static int prj_mpi_encode_edge_code(int dir, int i, int j, int k)
{
    return dir * PRJ_BLOCK_NCELLS + prj_mpi_encode_cell_index(i, j, k);
}

static void prj_mpi_decode_edge_code(int code, int *dir, int *i, int *j, int *k)
{
    if (dir == 0 || i == 0 || j == 0 || k == 0) {
        fprintf(stderr, "prj_mpi_decode_edge_code: null output pointer\n");
        abort();
    }
    *dir = code / PRJ_BLOCK_NCELLS;
    prj_mpi_decode_cell_index(code % PRJ_BLOCK_NCELLS, i, j, k);
}

static void prj_mpi_coarse_subface_position(const prj_block *block, int dir, int i, int j, int k,
    int aidx, int bidx, double x[3])
{
    if (block == 0 || x == 0) {
        fprintf(stderr, "prj_mpi_coarse_subface_position: null input\n");
        abort();
    }
    if (dir == X1DIR) {
        x[0] = block->xmin[0] + (double)i * block->dx[0];
        x[1] = block->xmin[1] + ((double)j + 0.25 + 0.5 * (double)aidx) * block->dx[1];
        x[2] = block->xmin[2] + ((double)k + 0.25 + 0.5 * (double)bidx) * block->dx[2];
        return;
    }
    if (dir == X2DIR) {
        x[0] = block->xmin[0] + ((double)i + 0.25 + 0.5 * (double)aidx) * block->dx[0];
        x[1] = block->xmin[1] + (double)j * block->dx[1];
        x[2] = block->xmin[2] + ((double)k + 0.25 + 0.5 * (double)bidx) * block->dx[2];
        return;
    }
    if (dir == X3DIR) {
        x[0] = block->xmin[0] + ((double)i + 0.25 + 0.5 * (double)aidx) * block->dx[0];
        x[1] = block->xmin[1] + ((double)j + 0.25 + 0.5 * (double)bidx) * block->dx[1];
        x[2] = block->xmin[2] + (double)k * block->dx[2];
        return;
    }
    fprintf(stderr, "prj_mpi_coarse_subface_position: invalid direction %d\n", dir);
    abort();
}

static void prj_mpi_coarse_subedge_position(const prj_block *block, int dir, int i, int j, int k,
    int subidx, double x[3])
{
    if (block == 0 || x == 0) {
        fprintf(stderr, "prj_mpi_coarse_subedge_position: null input\n");
        abort();
    }
    if (subidx < 0 || subidx > 1) {
        fprintf(stderr, "prj_mpi_coarse_subedge_position: invalid sub-edge index %d\n", subidx);
        abort();
    }
    if (dir == X1DIR) {
        x[0] = block->xmin[0] + ((double)i + 0.25 + 0.5 * (double)subidx) * block->dx[0];
        x[1] = block->xmin[1] + (double)j * block->dx[1];
        x[2] = block->xmin[2] + (double)k * block->dx[2];
        return;
    }
    if (dir == X2DIR) {
        x[0] = block->xmin[0] + (double)i * block->dx[0];
        x[1] = block->xmin[1] + ((double)j + 0.25 + 0.5 * (double)subidx) * block->dx[1];
        x[2] = block->xmin[2] + (double)k * block->dx[2];
        return;
    }
    if (dir == X3DIR) {
        x[0] = block->xmin[0] + (double)i * block->dx[0];
        x[1] = block->xmin[1] + (double)j * block->dx[1];
        x[2] = block->xmin[2] + ((double)k + 0.25 + 0.5 * (double)subidx) * block->dx[2];
        return;
    }
    fprintf(stderr, "prj_mpi_coarse_subedge_position: invalid direction %d\n", dir);
    abort();
}

static int prj_mpi_axis_touches(const prj_block *block, const prj_block *neighbor, int axis)
{
    const double tol = 1.0e-12;

    if (block == 0 || neighbor == 0 || axis < 0 || axis >= 3) {
        fprintf(stderr, "prj_mpi_axis_touches: invalid input\n");
        abort();
    }
    return fabs(block->xmax[axis] - neighbor->xmin[axis]) < tol ||
        fabs(neighbor->xmax[axis] - block->xmin[axis]) < tol;
}

static int prj_mpi_axis_side(const prj_block *block, const prj_block *neighbor, int axis)
{
    const double tol = 1.0e-12;

    if (block == 0 || neighbor == 0 || axis < 0 || axis >= 3) {
        fprintf(stderr, "prj_mpi_axis_side: invalid input\n");
        abort();
    }
    if (fabs(block->xmax[axis] - neighbor->xmin[axis]) < tol) {
        return 1;
    }
    if (fabs(neighbor->xmax[axis] - block->xmin[axis]) < tol) {
        return -1;
    }
    return 0;
}

static void prj_mpi_same_level_face_bounds(const prj_block *block, const prj_block *neighbor,
    const prj_neighbor *slot, int dir, int use_send_start, int start[3], int end[3])
{
    int axis;

    if (block == 0 || neighbor == 0 || slot == 0 || start == 0 || end == 0) {
        fprintf(stderr, "prj_mpi_same_level_face_bounds: invalid input\n");
        abort();
    }
    for (axis = 0; axis < 3; ++axis) {
        int anchor = use_send_start != 0 ? slot->send_loc_start[axis] : slot->recv_loc_start[axis];

        if (prj_mpi_axis_touches(block, neighbor, axis)) {
            start[axis] = anchor == 0 ? -PRJ_NGHOST : PRJ_BLOCK_SIZE;
            end[axis] = anchor == 0 ?
                (dir == axis ? 0 : -1) :
                (PRJ_BLOCK_SIZE + PRJ_NGHOST - 1);
        } else {
            start[axis] = 0;
            end[axis] = dir == axis ? PRJ_BLOCK_SIZE : (PRJ_BLOCK_SIZE - 1);
        }
    }
}

static void prj_mpi_coarse_face_bounds(const prj_block *block, const prj_block *neighbor,
    const prj_neighbor *slot, int dir, int use_send_start, int start[3], int end[3])
{
    int axis;

    if (block == 0 || neighbor == 0 || slot == 0 || start == 0 || end == 0) {
        fprintf(stderr, "prj_mpi_coarse_face_bounds: invalid input\n");
        abort();
    }
    for (axis = 0; axis < 3; ++axis) {
        int side = prj_mpi_axis_side(block, neighbor, axis);
        int anchor = use_send_start != 0 ? slot->send_loc_start[axis] : slot->recv_loc_start[axis];

        if (side != 0) {
            start[axis] = side < 0 ? -PRJ_NGHOST : PRJ_BLOCK_SIZE;
            end[axis] = side < 0 ?
                (dir == axis ? 0 : -1) :
                (PRJ_BLOCK_SIZE + PRJ_NGHOST - 1);
        } else {
            start[axis] = anchor;
            end[axis] = anchor + PRJ_BLOCK_SIZE / 2 - 1 + (dir == axis ? 1 : 0);
        }
    }
}

static void prj_mpi_patch_bounds(const prj_block *block, const prj_block *neighbor,
    const prj_neighbor *slot, int use_send_start, int start[3], int end[3])
{
    int axis;

    if (block == 0 || neighbor == 0 || slot == 0 || start == 0 || end == 0) {
        fprintf(stderr, "prj_mpi_patch_bounds: invalid input\n");
        abort();
    }
    for (axis = 0; axis < 3; ++axis) {
        int anchor = use_send_start != 0 ? slot->send_loc_start[axis] : slot->recv_loc_start[axis];

        if (prj_mpi_axis_touches(block, neighbor, axis)) {
            start[axis] = anchor == 0 ? -1 : (PRJ_BLOCK_SIZE - 1);
            end[axis] = anchor == 0 ? 0 : PRJ_BLOCK_SIZE;
        } else {
            start[axis] = prj_mpi_max_int(-PRJ_NGHOST + 1, anchor - 1);
            end[axis] = prj_mpi_min_int(PRJ_BLOCK_SIZE + PRJ_NGHOST - 2, anchor + PRJ_BLOCK_SIZE / 2);
        }
    }
}

static void prj_mpi_local_restrict_partial(const prj_mesh *mesh, const prj_block *block, int stage, int dir,
    int i, int j, int k, double *sum, int *count)
{
    int n;

    if (mesh == 0 || block == 0 || sum == 0 || count == 0) {
        fprintf(stderr, "prj_mpi_local_restrict_partial: null input\n");
        abort();
    }
    *sum = 0.0;
    *count = 0;
    for (n = 0; n < 56; ++n) {
        int nid = block->slot[n].id;
        const prj_block *neighbor;
        int start[3];
        int end[3];
        int aidx;
        int bidx;

        if (nid < 0 || nid >= mesh->nblocks) {
            continue;
        }
        neighbor = &mesh->blocks[nid];
        if (!prj_mpi_block_is_local(neighbor) || neighbor->level != block->level + 1) {
            continue;
        }
        prj_mpi_coarse_face_bounds(block, neighbor, &block->slot[n], dir, 1, start, end);
        if (i < start[0] || i > end[0] ||
            j < start[1] || j > end[1] ||
            k < start[2] || k > end[2]) {
            continue;
        }
        for (aidx = 0; aidx < 2; ++aidx) {
            for (bidx = 0; bidx < 2; ++bidx) {
                double x[3];
                int ii;
                int jj;
                int kk;

                prj_mpi_coarse_subface_position(block, dir, i, j, k, aidx, bidx, x);
                if (!prj_mhd_block_owns_face_position(neighbor, dir, x, &ii, &jj, &kk)) {
                    continue;
                }
                *sum += prj_mhd_bface_stage_const(neighbor, stage, dir)[IDX(ii, jj, kk)];
                *count += 1;
            }
        }
    }
}

static int prj_mpi_coarse_edge_bounds(const prj_block *block, const prj_block *neighbor,
    const prj_neighbor *slot, int dir, int use_send_start, int start[3], int end[3])
{
    int axis;
    int touch_count = 0;
    int overlap_axis = -1;
    int side[3];

    if (block == 0 || neighbor == 0 || slot == 0 || start == 0 || end == 0) {
        fprintf(stderr, "prj_mpi_coarse_edge_bounds: invalid input\n");
        abort();
    }
    for (axis = 0; axis < 3; ++axis) {
        side[axis] = prj_mpi_axis_side(block, neighbor, axis);
        if (side[axis] != 0) {
            touch_count += 1;
        } else {
            overlap_axis = axis;
        }
    }
    if (touch_count <= 0 || touch_count > 2) {
        return 0;
    }
    if (touch_count == 1) {
        for (axis = 0; axis < 3; ++axis) {
            if (side[axis] != 0 && dir == axis) {
                return 0;
            }
        }
    } else if (dir != overlap_axis) {
        return 0;
    }

    for (axis = 0; axis < 3; ++axis) {
        int anchor = use_send_start != 0 ? slot->send_loc_start[axis] : slot->recv_loc_start[axis];

        if (side[axis] != 0) {
            start[axis] = side[axis] < 0 ? 0 : PRJ_BLOCK_SIZE;
            end[axis] = start[axis];
        } else if (dir == axis) {
            start[axis] = anchor;
            end[axis] = anchor + PRJ_BLOCK_SIZE / 2 - 1;
        } else {
            start[axis] = anchor + (anchor == 0 ? 0 : 1);
            end[axis] = anchor + PRJ_BLOCK_SIZE / 2;
        }
    }
    return 1;
}

static void prj_mpi_local_edge_partial(const prj_mesh *mesh, const prj_block *block, int dir,
    int i, int j, int k, double *sum, int *count)
{
    int n;

    if (mesh == 0 || block == 0 || sum == 0 || count == 0) {
        fprintf(stderr, "prj_mpi_local_edge_partial: null input\n");
        abort();
    }
    *sum = 0.0;
    *count = 0;
    for (n = 0; n < 56; ++n) {
        int nid = block->slot[n].id;
        const prj_block *neighbor;
        int start[3];
        int end[3];
        int subidx;

        if (nid < 0 || nid >= mesh->nblocks) {
            continue;
        }
        neighbor = &mesh->blocks[nid];
        if (block->slot[n].type == PRJ_NEIGHBOR_CORNER ||
            !prj_mpi_block_is_local(neighbor) || neighbor->level != block->level + 1) {
            continue;
        }
        if (!prj_mpi_coarse_edge_bounds(block, neighbor, &block->slot[n], dir, 1, start, end)) {
            continue;
        }
        if (i < start[0] || i > end[0] ||
            j < start[1] || j > end[1] ||
            k < start[2] || k > end[2]) {
            continue;
        }
        for (subidx = 0; subidx < 2; ++subidx) {
            double x[3];
            int ii;
            int jj;
            int kk;

            prj_mpi_coarse_subedge_position(block, dir, i, j, k, subidx, x);
            if (!prj_mhd_block_owns_edge_position(neighbor, dir, x, &ii, &jj, &kk)) {
                continue;
            }
            *sum += neighbor->emf[dir][IDX(ii, jj, kk)];
            *count += 1;
        }
    }
}

static int prj_mpi_patch_overlaps_fine_extended_box(const prj_block *coarse, int i, int j, int k,
    const prj_block *fine)
{
    double patch_min[3];
    double patch_max[3];
    double fine_min[3];
    double fine_max[3];
    int d;

    if (coarse == 0 || fine == 0) {
        fprintf(stderr, "prj_mpi_patch_overlaps_fine_extended_box: null input\n");
        abort();
    }
    for (d = 0; d < 3; ++d) {
        double fine_dx = fine->dx[d];

        patch_min[d] = coarse->xmin[d] + (double)((d == 0) ? i : (d == 1 ? j : k)) * coarse->dx[d];
        patch_max[d] = patch_min[d] + coarse->dx[d];
        fine_min[d] = fine->xmin[d] - (double)PRJ_NGHOST * fine_dx;
        fine_max[d] = fine->xmax[d] + (double)PRJ_NGHOST * fine_dx;
        if (prj_riemann_min_double(patch_max[d], fine_max[d]) -
                prj_riemann_max_double(patch_min[d], fine_min[d]) <= 1.0e-12) {
            return 0;
        }
    }
    return 1;
}

static int prj_mpi_append_bf_patch_values(const double *Bface[3], int i, int j, int k,
    double **buffer, int *count, int *capacity)
{
    int dir;
    int side_slot;
    int aoff;
    int boff;

    for (dir = 0; dir < 3; ++dir) {
        for (side_slot = 0; side_slot < 2; ++side_slot) {
            for (aoff = -1; aoff <= 1; ++aoff) {
                for (boff = -1; boff <= 1; ++boff) {
                    int ii = i;
                    int jj = j;
                    int kk = k;

                    if (dir == X1DIR) {
                        ii = i + side_slot;
                        jj = j + aoff;
                        kk = k + boff;
                    } else if (dir == X2DIR) {
                        ii = i + aoff;
                        jj = j + side_slot;
                        kk = k + boff;
                    } else {
                        ii = i + aoff;
                        jj = j + boff;
                        kk = k + side_slot;
                    }
                    if (prj_mpi_append_double(buffer, count, capacity,
                            Bface[dir][IDX(ii, jj, kk)]) != 0) {
                        return 1;
                    }
                }
            }
        }
    }
    return 0;
}

static void prj_mpi_unpack_bf_patch_values(double *dst[3], int i, int j, int k, const double *buffer)
{
    int dir;
    int side_slot;
    int aoff;
    int boff;
    int pos = 0;

    if (dst == 0 || dst[0] == 0 || dst[1] == 0 || dst[2] == 0 || buffer == 0) {
        fprintf(stderr, "prj_mpi_unpack_bf_patch_values: invalid input\n");
        abort();
    }
    prj_fill(dst[0], PRJ_BLOCK_NCELLS, 0.0);
    prj_fill(dst[1], PRJ_BLOCK_NCELLS, 0.0);
    prj_fill(dst[2], PRJ_BLOCK_NCELLS, 0.0);
    for (dir = 0; dir < 3; ++dir) {
        for (side_slot = 0; side_slot < 2; ++side_slot) {
            for (aoff = -1; aoff <= 1; ++aoff) {
                for (boff = -1; boff <= 1; ++boff) {
                    int ii = i;
                    int jj = j;
                    int kk = k;

                    if (dir == X1DIR) {
                        ii = i + side_slot;
                        jj = j + aoff;
                        kk = k + boff;
                    } else if (dir == X2DIR) {
                        ii = i + aoff;
                        jj = j + side_slot;
                        kk = k + boff;
                    } else {
                        ii = i + aoff;
                        jj = j + boff;
                        kk = k + side_slot;
                    }
                    dst[dir][IDX(ii, jj, kk)] = buffer[pos++];
                }
            }
        }
    }
}
#endif

static int prj_mpi_append_double(double **array, int *count, int *capacity, double value)
{
    double *next;
    int new_capacity;

    if (*count >= *capacity) {
        new_capacity = *capacity == 0 ? 128 : 2 * (*capacity);
        next = (double *)realloc(*array, (size_t)new_capacity * sizeof(**array));
        if (next == 0) {
            return 1;
        }
        *array = next;
        *capacity = new_capacity;
    }
    (*array)[*count] = value;
    *count += 1;
    return 0;
}

static int prj_mpi_append_triplet(int **arrays, int *count, int *capacity, int a, int b, int c)
{
    int new_capacity;
    int *next0;
    int *next1;
    int *next2;

    if (*count >= *capacity) {
        new_capacity = *capacity == 0 ? 32 : 2 * (*capacity);
        next0 = (int *)realloc(arrays[0], (size_t)new_capacity * sizeof(*next0));
        next1 = (int *)realloc(arrays[1], (size_t)new_capacity * sizeof(*next1));
        next2 = (int *)realloc(arrays[2], (size_t)new_capacity * sizeof(*next2));
        if (next0 == 0 || next1 == 0 || next2 == 0) {
            free(next0);
            free(next1);
            free(next2);
            return 1;
        }
        arrays[0] = next0;
        arrays[1] = next1;
        arrays[2] = next2;
        *capacity = new_capacity;
    }

    arrays[0][*count] = a;
    arrays[1][*count] = b;
    arrays[2][*count] = c;
    *count += 1;
    return 0;
}

static void prj_mpi_assign_block_storage(prj_mesh *mesh)
{
    int bidx;

    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (block->id < 0) {
            continue;
        }
        if (prj_mpi_active != 0 && block->rank != prj_mpi_active->rank) {
            if (block->W != 0) {
                prj_block_free_data(block);
            }
            continue;
        }
        if (block->W == 0) {
            prj_block_alloc_data(block);
        }
    }
}

static size_t prj_mpi_block_data_count(void)
{
    size_t prim_count;
    size_t cons_count;

    prim_count = (size_t)PRJ_NVAR_PRIM * (size_t)PRJ_BLOCK_NCELLS;
    cons_count = (size_t)PRJ_NVAR_CONS * (size_t)PRJ_BLOCK_NCELLS;
    return 2U * prim_count + (size_t)PRJ_NVAR_EOSVAR * (size_t)PRJ_BLOCK_NCELLS +
        5U * cons_count + 9U * (size_t)PRJ_BLOCK_NCELLS
#if PRJ_MHD
        + 15U * (size_t)PRJ_BLOCK_NCELLS
#endif
        ;
}

static void prj_mpi_sync_slot_ranks(prj_mesh *mesh)
{
    int i;

    if (mesh == 0) {
        return;
    }
    for (i = 0; i < mesh->nblocks; ++i) {
        int n;

        for (n = 0; n < 56; ++n) {
            int nid = mesh->blocks[i].slot[n].id;

            if (nid >= 0 && nid < mesh->nblocks) {
                mesh->blocks[i].slot[n].rank = mesh->blocks[nid].rank;
            }
        }
    }
}

static void prj_mpi_compute_decomposition(prj_mesh *mesh)
{
    struct prj_order_item {
        unsigned long long key;
        int id;
    };
    struct prj_order_item *items;
    int count;
    int i;
    int cursor;
    int rank;

    if (mesh == 0) {
        return;
    }
    count = 0;
    for (i = 0; i < mesh->nblocks; ++i) {
        if (prj_block_is_active(&mesh->blocks[i])) {
            count += 1;
        }
    }
    items = (struct prj_order_item *)malloc((size_t)count * sizeof(*items));
    if (items == 0) {
        return;
    }
    cursor = 0;
    for (i = 0; i < mesh->nblocks; ++i) {
        if (prj_block_is_active(&mesh->blocks[i])) {
            items[cursor].id = i;
            items[cursor].key = prj_mpi_morton_code(mesh, &mesh->blocks[i]);
            cursor += 1;
        }
    }
    for (i = 0; i < count; ++i) {
        int j;

        for (j = i + 1; j < count; ++j) {
            if (items[j].key < items[i].key ||
                (items[j].key == items[i].key && items[j].id < items[i].id)) {
                struct prj_order_item tmp = items[i];

                items[i] = items[j];
                items[j] = tmp;
            }
        }
    }
    for (i = 0; i < count; ++i) {
        int nrank = prj_mpi_active != 0 ? prj_mpi_active->totrank : 1;

        rank = (int)(((long long)i * (long long)nrank) / (long long)count);
        if (rank >= nrank) {
            rank = nrank - 1;
        }
        mesh->blocks[items[i].id].rank = rank;
    }
    prj_mpi_sync_slot_ranks(mesh);
    free(items);
}

static void prj_mpi_migrate_active_blocks(prj_mesh *mesh, const int *old_ranks)
{
#if defined(PRJ_ENABLE_MPI)
    size_t data_count;
    int bidx;

    if (mesh == 0 || prj_mpi_active == 0 || old_ranks == 0 || prj_mpi_active->totrank <= 1) {
        return;
    }

    data_count = prj_mpi_block_data_count();
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        int old_rank;
        int new_rank;
        int tag;

        if (!prj_block_is_active(block)) {
            continue;
        }
        old_rank = old_ranks[bidx];
        new_rank = block->rank;
        if (old_rank == new_rank) {
            continue;
        }
        tag = 500 + bidx;
        if (prj_mpi_active->rank == new_rank || prj_mpi_active->rank == old_rank) {
            double *sendbuf = 0;
            double *recvbuf = 0;
            int source = MPI_PROC_NULL;
            int dest = MPI_PROC_NULL;
            int recvcount = 0;
            int sendcount = 0;

            if (prj_mpi_active->rank == new_rank) {
                if (block->W == 0 && prj_block_alloc_data(block) != 0) {
                    continue;
                }
                recvbuf = block->W;
                recvcount = (int)data_count;
                source = old_rank;
            }
            if (prj_mpi_active->rank == old_rank && block->W != 0) {
                sendbuf = block->W;
                sendcount = (int)data_count;
                dest = new_rank;
            }

            MPI_Sendrecv(sendbuf, sendcount, MPI_DOUBLE, dest, tag,
                recvbuf, recvcount, MPI_DOUBLE, source, tag,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            if (prj_mpi_active->rank == old_rank && block->W != 0) {
                prj_block_free_data(block);
            }
        }
    }
#else
    (void)mesh;
    (void)old_ranks;
#endif
}

static void prj_mpi_print_balance(const prj_mesh *mesh)
{
#if defined(PRJ_ENABLE_MPI)
    int *counts;
    int local_active;
    int bidx;
    int rank;
    int min_count;
    int max_count;
    int total_active;

    if (mesh == 0 || prj_mpi_active == 0 || prj_mpi_active->totrank <= 1) {
        return;
    }

    counts = (int *)calloc((size_t)prj_mpi_active->totrank, sizeof(*counts));
    if (counts == 0) {
        return;
    }

    local_active = 0;
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        if (prj_block_is_active(&mesh->blocks[bidx]) &&
            mesh->blocks[bidx].rank == prj_mpi_active->rank) {
            local_active += 1;
        }
    }

    MPI_Gather(&local_active, 1, MPI_INT, counts, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (prj_mpi_active->rank == 0) {
        min_count = counts[0];
        max_count = counts[0];
        total_active = counts[0];
        for (rank = 1; rank < prj_mpi_active->totrank; ++rank) {
            int count = counts[rank];

            if (count < min_count) {
                min_count = count;
            }
            if (count > max_count) {
                max_count = count;
            }
            total_active += count;
        }
        fprintf(stderr, "[mpi rebalance] total=%d min=%d max=%d spread=%d\n",
            total_active, min_count, max_count, max_count - min_count);
        fflush(stderr);
    }

    free(counts);
#else
    (void)mesh;
#endif
}

static int prj_mpi_has_rebalanced(const prj_mesh *mesh, const int *old_ranks)
{
    int bidx;

    if (mesh == 0 || old_ranks == 0) {
        return 0;
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        if (!prj_block_is_active(&mesh->blocks[bidx])) {
            continue;
        }
        if (mesh->blocks[bidx].rank != old_ranks[bidx]) {
            return 1;
        }
    }
    return 0;
}

static void prj_mpi_collect_active_counts(const prj_mesh *mesh, int *counts)
{
    int bidx;

    if (mesh == 0 || counts == 0 || prj_mpi_active == 0) {
        return;
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        const prj_block *block = &mesh->blocks[bidx];

        if (!prj_block_is_active(block)) {
            continue;
        }
        if (block->rank >= 0 && block->rank < prj_mpi_active->totrank) {
            counts[block->rank] += 1;
        }
    }
}

static int prj_mpi_counts_changed(const int *before, const int *after, int nrank)
{
    int rank;

    if (before == 0 || after == 0) {
        return 0;
    }
    for (rank = 0; rank < nrank; ++rank) {
        if (before[rank] != after[rank]) {
            return 1;
        }
    }
    return 0;
}

static int prj_mpi_buffer_record_total(const int *sizes, int count)
{
    int total;
    int i;

    total = 0;
    if (sizes == 0) {
        return 0;
    }
    for (i = 0; i < count; ++i) {
        total += sizes[i];
    }
    return total;
}

static int prj_mpi_buffer_recv_count(const prj_mpi_buffer *buffer)
{
    int count;

    if (buffer == 0 || buffer->cell_data_size_recv == 0) {
        return 0;
    }
    count = 0;
    while (buffer->cell_data_size_recv[count] >= 0) {
        count += 1;
    }
    return count;
}

static int prj_mpi_build_ghost_plan_for_neighbor(prj_mesh *mesh, prj_mpi *mpi, prj_mpi_buffer *buffer)
{
#if defined(PRJ_ENABLE_MPI)
    int *send_sizes;
    int send_entries;
    int cell_size_total;
    int axis;
    int bidx;
    int occ;
    int *idx_send[3];
    int idx_capacity;
    int record_count;
    MPI_Request requests[6];

    send_sizes = (int *)calloc((size_t)buffer->number, sizeof(*send_sizes));
    idx_send[0] = 0;
    idx_send[1] = 0;
    idx_send[2] = 0;
    idx_capacity = 0;
    record_count = 0;
    occ = 0;
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        int n;

        if (!prj_block_is_active(block) || block->rank != mpi->rank) {
            continue;
        }
        for (n = 0; n < 56; ++n) {
            int nid = block->slot[n].id;
            const prj_block *neighbor;
            int i;
            int j;
            int k;
            int before;

            if (nid < 0 || nid >= mesh->nblocks || mesh->blocks[nid].rank != buffer->receiver_rank) {
                continue;
            }
            neighbor = &mesh->blocks[nid];
            before = record_count;
            for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
                for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
                    for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                        double x1;
                        double x2;
                        double x3;

                        if (i >= 0 && i < PRJ_BLOCK_SIZE &&
                            j >= 0 && j < PRJ_BLOCK_SIZE &&
                            k >= 0 && k < PRJ_BLOCK_SIZE) {
                            continue;
                        }
                        x1 = neighbor->xmin[0] + ((double)i + 0.5) * neighbor->dx[0];
                        x2 = neighbor->xmin[1] + ((double)j + 0.5) * neighbor->dx[1];
                        x3 = neighbor->xmin[2] + ((double)k + 0.5) * neighbor->dx[2];
                        if (prj_mpi_cell_point_inside(block, x1, x2, x3)) {
                            int sample_kind = prj_mpi_sample_kind(block, x1, x2, x3);
                            int sample_code = prj_mpi_sample_code(sample_kind, block->level == neighbor->level);

                            if (prj_mpi_append_triplet(idx_send, &record_count, &idx_capacity,
                                    nid, prj_mpi_encode_cell_index(i, j, k), sample_code) != 0) {
                                free(send_sizes);
                                for (axis = 0; axis < 3; ++axis) {
                                    free(idx_send[axis]);
                                }
                                return 1;
                            }
                        }
                    }
                }
            }
            send_sizes[occ] = record_count - before;
            occ += 1;
        }
    }

    buffer->cell_data_size_send = send_sizes;
    for (axis = 0; axis < 3; ++axis) {
        buffer->cell_data_idx_send[axis] = idx_send[axis];
    }
    buffer->cell_buffer_send = (double *)calloc((size_t)record_count * PRJ_MPI_GHOST_NVAR, sizeof(*buffer->cell_buffer_send));

    MPI_Sendrecv(&buffer->number, 1, MPI_INT, buffer->receiver_rank, 100,
        &send_entries, 1, MPI_INT, buffer->receiver_rank, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    buffer->cell_data_size_recv = (int *)calloc((size_t)send_entries + 1U, sizeof(*buffer->cell_data_size_recv));
    MPI_Sendrecv(buffer->cell_data_size_send, buffer->number, MPI_INT, buffer->receiver_rank, 101,
        buffer->cell_data_size_recv, send_entries, MPI_INT, buffer->receiver_rank, 101, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if (buffer->cell_data_size_recv != 0) {
        buffer->cell_data_size_recv[send_entries] = -1;
    }
    cell_size_total = prj_mpi_buffer_record_total(buffer->cell_data_size_recv, send_entries);
    for (axis = 0; axis < 3; ++axis) {
        buffer->cell_data_idx_recv[axis] = (int *)calloc((size_t)cell_size_total, sizeof(*buffer->cell_data_idx_recv[axis]));
    }
    buffer->cell_buffer_recv = (double *)calloc((size_t)cell_size_total * PRJ_MPI_GHOST_NVAR, sizeof(*buffer->cell_buffer_recv));

    for (axis = 0; axis < 3; ++axis) {
        MPI_Irecv(buffer->cell_data_idx_recv[axis], cell_size_total, MPI_INT, buffer->receiver_rank, 110 + axis,
            MPI_COMM_WORLD, &requests[2 * axis]);
        MPI_Isend(buffer->cell_data_idx_send[axis], record_count, MPI_INT, buffer->receiver_rank, 110 + axis,
            MPI_COMM_WORLD, &requests[2 * axis + 1]);
    }
    MPI_Waitall(6, requests, MPI_STATUSES_IGNORE);
    return 0;
#else
    (void)mesh;
    (void)mpi;
    (void)buffer;
    return 0;
#endif
}

static void prj_mpi_pack_ghost_values(prj_mesh *mesh, prj_mpi *mpi, prj_mpi_buffer *buffer,
    int stage, int fill_kind)
{
    int bidx;
    int occ;
    int pos;

    if (mesh == 0 || mpi == 0 || buffer == 0 || buffer->cell_buffer_send == 0) {
        return;
    }

    pos = 0;
    occ = 0;
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        int n;

        if (!prj_block_is_active(block) || block->rank != mpi->rank) {
            continue;
        }
        for (n = 0; n < 56; ++n) {
            int nid = block->slot[n].id;
            const prj_block *neighbor;
            int i;
            int j;
            int k;

            if (nid < 0 || nid >= mesh->nblocks || mesh->blocks[nid].rank != buffer->receiver_rank) {
                continue;
            }
            neighbor = &mesh->blocks[nid];
            for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
                for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
                    for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                        double x1;
                        double x2;
                        double x3;

                        if (i >= 0 && i < PRJ_BLOCK_SIZE &&
                            j >= 0 && j < PRJ_BLOCK_SIZE &&
                            k >= 0 && k < PRJ_BLOCK_SIZE) {
                            continue;
                        }
                        x1 = neighbor->xmin[0] + ((double)i + 0.5) * neighbor->dx[0];
                        x2 = neighbor->xmin[1] + ((double)j + 0.5) * neighbor->dx[1];
                        x3 = neighbor->xmin[2] + ((double)k + 0.5) * neighbor->dx[2];
                        if (prj_mpi_cell_point_inside(block, x1, x2, x3)) {
                            double w[PRJ_NVAR_PRIM];
                            double eosv[PRJ_NVAR_EOSVAR];
                            int sample_kind;
                            int v;

                            sample_kind = prj_mpi_sample_kind(block, x1, x2, x3);
                            if (sample_kind < 0) {
                                fprintf(stderr,
                                    "prj_mpi_pack_ghost_values: failed to classify ghost sample "
                                    "src_block=%d dst_block=%d i=%d j=%d k=%d x=(%.17g, %.17g, %.17g)\n",
                                    block->id, nid, i, j, k, x1, x2, x3);
                                exit(EXIT_FAILURE);
                            }
                            if (fill_kind == 2 || sample_kind == fill_kind) {
                                prj_boundary_get_prim(block, stage, x1, x2, x3, w);
                                prj_boundary_get_eosvar(block, x1, x2, x3, eosv);
                            } else {
                                for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
                                    w[v] = 0.0;
                                }
                                for (v = 0; v < PRJ_NVAR_EOSVAR; ++v) {
                                    eosv[v] = 0.0;
                                }
                            }
                            for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
                                buffer->cell_buffer_send[pos++] = w[v];
                            }
                            for (v = 0; v < PRJ_NVAR_EOSVAR; ++v) {
                                buffer->cell_buffer_send[pos++] = eosv[v];
                            }
                        }
                    }
                }
            }
            occ += 1;
        }
    }

    (void)occ;
}

prj_mpi *prj_mpi_current(void)
{
    return prj_mpi_active;
}

void prj_mpi_init(int *argc, char ***argv, prj_mpi *mpi)
{
    if (mpi == 0) {
        return;
    }
    memset(mpi, 0, sizeof(*mpi));
#if defined(PRJ_ENABLE_MPI)
    MPI_Init(argc, argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi->rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi->totrank);
#else
    (void)argc;
    (void)argv;
    mpi->rank = 0;
    mpi->totrank = 1;
#endif
    prj_mpi_active = mpi;
}

void prj_mpi_decompose(prj_mesh *mesh)
{
    prj_mpi_compute_decomposition(mesh);
    prj_mpi_assign_block_storage(mesh);
}

void prj_mpi_prepare(prj_mesh *mesh, prj_mpi *mpi)
{
    int *rank_seen;
    int count;
    int bidx;
    int i;

    if (mesh == 0 || mpi == 0) {
        return;
    }
    prj_mpi_clear_neighbors(mpi);
    rank_seen = (int *)calloc((size_t)mpi->totrank, sizeof(*rank_seen));
    if (rank_seen == 0) {
        return;
    }
    count = 0;
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        const prj_block *block = &mesh->blocks[bidx];
        int n;

        if (!prj_block_is_active(block) || block->rank != mpi->rank) {
            continue;
        }
        for (n = 0; n < 56; ++n) {
            int nid = block->slot[n].id;
            int nrank;

            if (nid < 0 || nid >= mesh->nblocks) {
                continue;
            }
            nrank = mesh->blocks[nid].rank;
            if (nrank != mpi->rank && nrank >= 0 && nrank < mpi->totrank && rank_seen[nrank] == 0) {
                rank_seen[nrank] = 1;
                count += 1;
            }
        }
    }
    mpi->neighbor_number = count;
    mpi->neighbor_buffer = (prj_mpi_buffer *)calloc((size_t)count, sizeof(*mpi->neighbor_buffer));
    if (mpi->neighbor_buffer == 0) {
        mpi->neighbor_number = 0;
        free(rank_seen);
        return;
    }
    count = 0;
    for (i = 0; i < mpi->totrank; ++i) {
        if (rank_seen[i] != 0) {
            mpi->neighbor_buffer[count].receiver_rank = i;
            count += 1;
        }
    }
    for (i = 0; i < mpi->neighbor_number; ++i) {
        prj_mpi_buffer *buffer = &mpi->neighbor_buffer[i];
        int occ = 0;

        for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
            const prj_block *block = &mesh->blocks[bidx];
            int n;

            if (!prj_block_is_active(block) || block->rank != mpi->rank) {
                continue;
            }
            for (n = 0; n < 56; ++n) {
                int nid = block->slot[n].id;

                if (nid >= 0 && nid < mesh->nblocks && mesh->blocks[nid].rank == buffer->receiver_rank) {
                    occ += 1;
                }
            }
        }
        buffer->number = occ;
        buffer->receiver_blocks = (int *)calloc((size_t)occ, sizeof(*buffer->receiver_blocks));
        if (buffer->receiver_blocks != 0) {
            int idx = 0;

            for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
                const prj_block *block = &mesh->blocks[bidx];
                int n;

                if (!prj_block_is_active(block) || block->rank != mpi->rank) {
                    continue;
                }
                for (n = 0; n < 56; ++n) {
                    int nid = block->slot[n].id;

                    if (nid >= 0 && nid < mesh->nblocks && mesh->blocks[nid].rank == buffer->receiver_rank) {
                        buffer->receiver_blocks[idx] = nid;
                        idx += 1;
                    }
                }
            }
        }
        prj_mpi_build_ghost_plan_for_neighbor(mesh, mpi, buffer);
    }
    free(rank_seen);
}

void prj_mpi_exchange_ghosts(prj_mesh *mesh, prj_mpi *mpi, int stage, int fill_kind)
{
#if defined(PRJ_ENABLE_MPI)
    int nb;
    MPI_Request *requests;
    int request_count;
    int request_capacity;

    if (mesh == 0 || mpi == 0 || mpi->totrank <= 1) {
        return;
    }
    requests = 0;
    request_count = 0;
    request_capacity = 0;
    for (nb = 0; nb < mpi->neighbor_number; ++nb) {
        prj_mpi_buffer *buffer = &mpi->neighbor_buffer[nb];
        int record_count;
        int recv_count;
        int cell_size_total;

        record_count = prj_mpi_buffer_record_total(buffer->cell_data_size_send, buffer->number);
        recv_count = prj_mpi_buffer_recv_count(buffer);
        cell_size_total = prj_mpi_buffer_record_total(buffer->cell_data_size_recv, recv_count);
        prj_mpi_pack_ghost_values(mesh, mpi, buffer, stage, fill_kind);

        if (request_count + 8 > request_capacity) {
            int new_capacity = request_capacity == 0 ? 16 : request_capacity * 2;
            MPI_Request *next = (MPI_Request *)realloc(requests, (size_t)new_capacity * sizeof(*next));

            if (next == 0) {
                free(requests);
                return;
            }
            requests = next;
            request_capacity = new_capacity;
        }
        MPI_Irecv(buffer->cell_buffer_recv, cell_size_total * PRJ_MPI_GHOST_NVAR, MPI_DOUBLE, buffer->receiver_rank, 120,
            MPI_COMM_WORLD, &requests[request_count++]);
        MPI_Isend(buffer->cell_buffer_send, record_count * PRJ_MPI_GHOST_NVAR, MPI_DOUBLE, buffer->receiver_rank, 120,
            MPI_COMM_WORLD, &requests[request_count++]);
    }
    if (request_count > 0) {
        MPI_Waitall(request_count, requests, MPI_STATUSES_IGNORE);
    }
    for (nb = 0; nb < mpi->neighbor_number; ++nb) {
        prj_mpi_buffer *buffer = &mpi->neighbor_buffer[nb];
        int total = 0;
        int pos = 0;
        int i;

        for (i = 0; buffer->cell_data_size_recv != 0 && buffer->cell_data_size_recv[i] >= 0; ++i) {
            total += buffer->cell_data_size_recv[i];
        }
        for (i = 0; i < total; ++i) {
            int block_id = buffer->cell_data_idx_recv[0][i];
            int code = buffer->cell_data_idx_recv[1][i];
            int sample_code = buffer->cell_data_idx_recv[2][i];
            int sample_kind = prj_mpi_sample_kind_from_code(sample_code);
            int same_level = prj_mpi_sample_is_same_level(sample_code);
            int ii;
            int jj;
            int kk;
            int v;
            prj_block *block;
            double *dst;

            if (block_id < 0 || block_id >= mesh->nblocks) {
                pos += PRJ_MPI_GHOST_NVAR;
                continue;
            }
            if (sample_kind < 0 || (fill_kind != 2 && sample_kind != fill_kind)) {
                pos += PRJ_MPI_GHOST_NVAR;
                continue;
            }
            block = &mesh->blocks[block_id];
            if (block->rank != mpi->rank || block->W == 0 || block->eosvar == 0) {
                pos += PRJ_MPI_GHOST_NVAR;
                continue;
            }
            dst = prj_mpi_stage_array(block, stage);
            prj_mpi_decode_cell_index(code, &ii, &jj, &kk);
            for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
                dst[VIDX(v, ii, jj, kk)] = buffer->cell_buffer_recv[pos++];
            }
            for (v = 0; v < PRJ_NVAR_EOSVAR; ++v) {
                block->eosvar[EIDX(v, ii, jj, kk)] = buffer->cell_buffer_recv[pos++];
            }
            if (block->eos_done != 0) {
                block->eos_done[IDX(ii, jj, kk)] = same_level ? 1 : 0;
            }
        }
    }
    free(requests);
#else
    (void)mesh;
    (void)mpi;
    (void)stage;
    (void)fill_kind;
#endif
}

void prj_mpi_exchange_bf(prj_mesh *mesh, prj_mpi *mpi, int stage, int fill_kind)
{
#if defined(PRJ_ENABLE_MPI) && PRJ_MHD
    int nb;
    MPI_Request *requests;
    int request_count;
    int request_capacity;
    int value_stride;

    if (mesh == 0 || mpi == 0 || mpi->totrank <= 1) {
        return;
    }
    if (fill_kind != 0 && fill_kind != 1) {
        fprintf(stderr, "prj_mpi_exchange_bf: invalid fill kind %d\n", fill_kind);
        abort();
    }

    value_stride = fill_kind == 0 ? 1 : PRJ_MPI_BF_PATCH_VALUE_COUNT;
    requests = 0;
    request_count = 0;
    request_capacity = 0;
    for (nb = 0; nb < mpi->neighbor_number; ++nb) {
        prj_mpi_buffer *buffer = &mpi->neighbor_buffer[nb];
        int *send_sizes;
        int recv_entries;
        int face_total_recv;
        int axis;
        int bidx;
        int occ;
        int *idx_send[3];
        double *value_send;
        int idx_count;
        int idx_capacity;
        int value_count;
        int value_capacity;

        send_sizes = (int *)calloc((size_t)buffer->number, sizeof(*send_sizes));
        idx_send[0] = 0;
        idx_send[1] = 0;
        idx_send[2] = 0;
        value_send = 0;
        idx_count = 0;
        idx_capacity = 0;
        value_count = 0;
        value_capacity = 0;
        occ = 0;
        for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
            prj_block *block = &mesh->blocks[bidx];
            int n;

            if (!prj_mpi_block_is_local(block)) {
                continue;
            }
            for (n = 0; n < 56; ++n) {
                int nid = block->slot[n].id;
                const prj_block *neighbor;
                int before = idx_count;

                if (nid < 0 || nid >= mesh->nblocks || mesh->blocks[nid].rank != buffer->receiver_rank) {
                    continue;
                }
                neighbor = &mesh->blocks[nid];
                if (fill_kind == 0) {
                    if (block->level == neighbor->level) {
                        int dir;

                        for (dir = 0; dir < 3; ++dir) {
                            const double *srcface = prj_mhd_bface_stage_const(block, stage, dir);
                            int dst_start[3];
                            int dst_end[3];
                            int ii;
                            int jj;
                            int kk;

                            if (srcface == 0) {
                                fprintf(stderr, "prj_mpi_exchange_bf: same-level source face storage is not allocated\n");
                                abort();
                            }
                            prj_mpi_same_level_face_bounds(block, neighbor, &block->slot[n], dir, 0, dst_start, dst_end);
                            for (ii = dst_start[0]; ii <= dst_end[0]; ++ii) {
                                for (jj = dst_start[1]; jj <= dst_end[1]; ++jj) {
                                    for (kk = dst_start[2]; kk <= dst_end[2]; ++kk) {
                                        double x[3];
                                        int i;
                                        int j;
                                        int k;

                                        prj_mhd_face_position(neighbor, dir, ii, jj, kk, x);
                                        if (!prj_mhd_block_owns_face_position(block, dir, x, &i, &j, &k)) {
                                            continue;
                                        }
                                        if (prj_mpi_append_triplet(idx_send, &idx_count, &idx_capacity,
                                                neighbor->id, block->id,
                                                prj_mpi_encode_face_code(dir, ii, jj, kk)) != 0 ||
                                            prj_mpi_append_double(&value_send, &value_count, &value_capacity,
                                                srcface[IDX(i, j, k)]) != 0) {
                                            fprintf(stderr, "prj_mpi_exchange_bf: same-level pack failed\n");
                                            abort();
                                        }
                                    }
                                }
                            }
                        }
                    } else if (block->level == neighbor->level + 1) {
                        int dir;

                        for (dir = 0; dir < 3; ++dir) {
                            const double *srcface = prj_mhd_bface_stage_const(block, stage, dir);
                            int dst_start[3];
                            int dst_end[3];
                            int i;
                            int j;
                            int k;

                            if (srcface == 0) {
                                fprintf(stderr, "prj_mpi_exchange_bf: restriction source face storage is not allocated\n");
                                abort();
                            }
                            prj_mpi_coarse_face_bounds(block, neighbor, &block->slot[n], dir, 0, dst_start, dst_end);
                            for (i = dst_start[0]; i <= dst_end[0]; ++i) {
                                for (j = dst_start[1]; j <= dst_end[1]; ++j) {
                                    for (k = dst_start[2]; k <= dst_end[2]; ++k) {
                                        int aidx;
                                        int bidx2;

                                        if (prj_mhd_face_is_interior(dir, i, j, k)) {
                                            continue;
                                        }
                                        for (aidx = 0; aidx < 2; ++aidx) {
                                            for (bidx2 = 0; bidx2 < 2; ++bidx2) {
                                                double x[3];
                                                int ii;
                                                int jj;
                                                int kk;

                                                prj_mpi_coarse_subface_position(neighbor, dir, i, j, k, aidx, bidx2, x);
                                                if (!prj_mhd_block_owns_face_position(block, dir, x, &ii, &jj, &kk)) {
                                                    continue;
                                                }
                                                if (prj_mpi_append_triplet(idx_send, &idx_count, &idx_capacity,
                                                        neighbor->id, block->id,
                                                        prj_mpi_encode_face_code(dir, i, j, k)) != 0 ||
                                                    prj_mpi_append_double(&value_send, &value_count, &value_capacity,
                                                        srcface[IDX(ii, jj, kk)]) != 0) {
                                                    fprintf(stderr, "prj_mpi_exchange_bf: restriction pack failed\n");
                                                    abort();
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                } else if (block->level == neighbor->level - 1) {
                    const double *srcface[3];
                    int start[3];
                    int end[3];
                    int i;
                    int j;
                    int k;

                    srcface[0] = prj_mhd_bface_stage_const(block, stage, X1DIR);
                    srcface[1] = prj_mhd_bface_stage_const(block, stage, X2DIR);
                    srcface[2] = prj_mhd_bface_stage_const(block, stage, X3DIR);
                    prj_mpi_patch_bounds(block, neighbor, &block->slot[n], 1, start, end);
                    for (i = start[0]; i <= end[0]; ++i) {
                        for (j = start[1]; j <= end[1]; ++j) {
                            for (k = start[2]; k <= end[2]; ++k) {
                                if (!prj_mpi_patch_overlaps_fine_extended_box(block, i, j, k, neighbor)) {
                                    continue;
                                }
                                if (prj_mpi_append_triplet(idx_send, &idx_count, &idx_capacity,
                                        neighbor->id, block->id,
                                        prj_mpi_encode_cell_index(i, j, k)) != 0 ||
                                    prj_mpi_append_bf_patch_values(srcface, i, j, k,
                                        &value_send, &value_count, &value_capacity) != 0) {
                                    fprintf(stderr, "prj_mpi_exchange_bf: recon pack failed\n");
                                    abort();
                                }
                            }
                        }
                    }
                }
                send_sizes[occ++] = idx_count - before;
            }
        }

        free(buffer->face_data_size_send);
        free(buffer->face_buffer_send);
        for (axis = 0; axis < 3; ++axis) {
            free(buffer->face_data_idx_send[axis]);
            free(buffer->face_data_idx_recv[axis]);
            buffer->face_data_idx_send[axis] = idx_send[axis];
            buffer->face_data_idx_recv[axis] = 0;
        }
        free(buffer->face_data_size_recv);
        free(buffer->face_buffer_recv);
        buffer->face_data_size_send = send_sizes;
        buffer->face_buffer_send = value_send;
        buffer->face_data_size_recv = 0;
        buffer->face_buffer_recv = 0;

        MPI_Sendrecv(&buffer->number, 1, MPI_INT, buffer->receiver_rank, 300 + fill_kind,
            &recv_entries, 1, MPI_INT, buffer->receiver_rank, 300 + fill_kind,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        buffer->face_data_size_recv = (int *)calloc((size_t)recv_entries + 1U, sizeof(*buffer->face_data_size_recv));
        MPI_Sendrecv(buffer->face_data_size_send, buffer->number, MPI_INT, buffer->receiver_rank, 310 + fill_kind,
            buffer->face_data_size_recv, recv_entries, MPI_INT, buffer->receiver_rank, 310 + fill_kind,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        buffer->face_data_size_recv[recv_entries] = -1;
        face_total_recv = 0;
        for (axis = 0; axis < recv_entries; ++axis) {
            face_total_recv += buffer->face_data_size_recv[axis];
        }
        for (axis = 0; axis < 3; ++axis) {
            buffer->face_data_idx_recv[axis] = (int *)calloc((size_t)face_total_recv, sizeof(*buffer->face_data_idx_recv[axis]));
        }
        buffer->face_buffer_recv = (double *)calloc((size_t)face_total_recv * (size_t)value_stride,
            sizeof(*buffer->face_buffer_recv));

        if (request_count + 8 > request_capacity) {
            int new_capacity = request_capacity == 0 ? 16 : request_capacity * 2;
            MPI_Request *next = (MPI_Request *)realloc(requests, (size_t)new_capacity * sizeof(*next));

            if (next == 0) {
                free(requests);
                return;
            }
            requests = next;
            request_capacity = new_capacity;
        }
        for (axis = 0; axis < 3; ++axis) {
            MPI_Irecv(buffer->face_data_idx_recv[axis], face_total_recv, MPI_INT, buffer->receiver_rank, 320 + 3 * fill_kind + axis,
                MPI_COMM_WORLD, &requests[request_count++]);
            MPI_Isend(buffer->face_data_idx_send[axis], idx_count, MPI_INT, buffer->receiver_rank, 320 + 3 * fill_kind + axis,
                MPI_COMM_WORLD, &requests[request_count++]);
        }
        MPI_Irecv(buffer->face_buffer_recv, face_total_recv * value_stride, MPI_DOUBLE, buffer->receiver_rank, 330 + fill_kind,
            MPI_COMM_WORLD, &requests[request_count++]);
        MPI_Isend(buffer->face_buffer_send, value_count, MPI_DOUBLE, buffer->receiver_rank, 330 + fill_kind,
            MPI_COMM_WORLD, &requests[request_count++]);
    }
    if (request_count > 0) {
        MPI_Waitall(request_count, requests, MPI_STATUSES_IGNORE);
    }

    if (fill_kind == 0) {
        int *acc_block = 0;
        int *acc_face = 0;
        int *acc_count = 0;
        double *acc_sum = 0;
        int acc_n = 0;
        int acc_cap = 0;

        for (nb = 0; nb < mpi->neighbor_number; ++nb) {
            prj_mpi_buffer *buffer = &mpi->neighbor_buffer[nb];
            int total = 0;
            int pos = 0;
            int i;

            if (buffer->face_data_size_recv == 0) {
                continue;
            }
            for (i = 0; buffer->face_data_size_recv[i] >= 0; ++i) {
                total += buffer->face_data_size_recv[i];
            }
            for (i = 0; i < total; ++i) {
                int block_id = buffer->face_data_idx_recv[0][i];
                int source_id = buffer->face_data_idx_recv[1][i];
                int face_code = buffer->face_data_idx_recv[2][i];
                double value = buffer->face_buffer_recv[pos++];

                if (block_id < 0 || block_id >= mesh->nblocks ||
                    source_id < 0 || source_id >= mesh->nblocks) {
                    continue;
                }
                if (mesh->blocks[block_id].rank != mpi->rank) {
                    continue;
                }
                if (mesh->blocks[source_id].level == mesh->blocks[block_id].level) {
                    int dir;
                    int ii;
                    int jj;
                    int kk;

                    prj_mpi_decode_face_code(face_code, &dir, &ii, &jj, &kk);
                    prj_mhd_bface_stage(&mesh->blocks[block_id], stage, dir)[IDX(ii, jj, kk)] = value;
                } else if (mesh->blocks[source_id].level == mesh->blocks[block_id].level + 1) {
                    int found = -1;
                    int a;

                    for (a = 0; a < acc_n; ++a) {
                        if (acc_block[a] == block_id && acc_face[a] == face_code) {
                            found = a;
                            break;
                        }
                    }
                    if (found < 0) {
                        if (acc_n >= acc_cap) {
                            int new_cap = acc_cap == 0 ? 32 : 2 * acc_cap;
                            int *next_block = (int *)realloc(acc_block, (size_t)new_cap * sizeof(*next_block));
                            int *next_face = (int *)realloc(acc_face, (size_t)new_cap * sizeof(*next_face));
                            int *next_count = (int *)realloc(acc_count, (size_t)new_cap * sizeof(*next_count));
                            double *next_sum = (double *)realloc(acc_sum, (size_t)new_cap * sizeof(*next_sum));

                            if (next_block == 0 || next_face == 0 || next_count == 0 || next_sum == 0) {
                                free(next_sum);
                                free(next_count);
                                free(next_face);
                                free(next_block);
                                free(acc_sum);
                                free(acc_count);
                                free(acc_face);
                                free(acc_block);
                                free(requests);
                                return;
                            }
                            acc_block = next_block;
                            acc_face = next_face;
                            acc_count = next_count;
                            acc_sum = next_sum;
                            acc_cap = new_cap;
                        }
                        found = acc_n++;
                        acc_block[found] = block_id;
                        acc_face[found] = face_code;
                        acc_count[found] = 0;
                        acc_sum[found] = 0.0;
                    }
                    acc_count[found] += 1;
                    acc_sum[found] += value;
                } else {
                    fprintf(stderr, "prj_mpi_exchange_bf: invalid NONRECON level relation\n");
                    abort();
                }
            }
        }
        for (nb = 0; nb < acc_n; ++nb) {
            int dir;
            int ii;
            int jj;
            int kk;
            prj_block *block = &mesh->blocks[acc_block[nb]];
            double local_sum;
            int local_count;

            prj_mpi_decode_face_code(acc_face[nb], &dir, &ii, &jj, &kk);
            prj_mpi_local_restrict_partial(mesh, block, stage, dir, ii, jj, kk, &local_sum, &local_count);
            if (!prj_mpi_valid_restrict_count(acc_count[nb] + local_count)) {
                fprintf(stderr,
                    "prj_mpi_exchange_bf: invalid restriction coverage block=%d dir=%d "
                    "face=(%d,%d,%d) recv_count=%d local_count=%d\n",
                    block->id, dir, ii, jj, kk, acc_count[nb], local_count);
                abort();
            }
            if (acc_count[nb] + local_count >= 4) {
                prj_mhd_bface_stage(block, stage, dir)[IDX(ii, jj, kk)] =
                    (acc_sum[nb] + local_sum) / (double)(acc_count[nb] + local_count);
            }
        }
        free(acc_sum);
        free(acc_count);
        free(acc_face);
        free(acc_block);
    } else {
        double *tmp_face[3];

        tmp_face[0] = (double *)calloc((size_t)PRJ_BLOCK_NCELLS, sizeof(*tmp_face[0]));
        tmp_face[1] = (double *)calloc((size_t)PRJ_BLOCK_NCELLS, sizeof(*tmp_face[1]));
        tmp_face[2] = (double *)calloc((size_t)PRJ_BLOCK_NCELLS, sizeof(*tmp_face[2]));
        if (tmp_face[0] == 0 || tmp_face[1] == 0 || tmp_face[2] == 0) {
            free(tmp_face[2]);
            free(tmp_face[1]);
            free(tmp_face[0]);
            free(requests);
            return;
        }
        for (nb = 0; nb < mpi->neighbor_number; ++nb) {
            prj_mpi_buffer *buffer = &mpi->neighbor_buffer[nb];
            int total = 0;
            int pos = 0;
            int irec;

            if (buffer->face_data_size_recv == 0) {
                continue;
            }
            for (irec = 0; buffer->face_data_size_recv[irec] >= 0; ++irec) {
                total += buffer->face_data_size_recv[irec];
            }
            for (irec = 0; irec < total; ++irec) {
                int block_id = buffer->face_data_idx_recv[0][irec];
                int source_id = buffer->face_data_idx_recv[1][irec];
                int coarse_code = buffer->face_data_idx_recv[2][irec];
                int ii;
                int jj;
                int kk;
                const prj_block *coarse;
                const double *srcface[3];

                if (block_id < 0 || block_id >= mesh->nblocks ||
                    source_id < 0 || source_id >= mesh->nblocks) {
                    pos += PRJ_MPI_BF_PATCH_VALUE_COUNT;
                    continue;
                }
                if (mesh->blocks[block_id].rank != mpi->rank) {
                    pos += PRJ_MPI_BF_PATCH_VALUE_COUNT;
                    continue;
                }
                coarse = &mesh->blocks[source_id];
                if (coarse->level != mesh->blocks[block_id].level - 1) {
                    fprintf(stderr, "prj_mpi_exchange_bf: invalid RECON level relation\n");
                    abort();
                }
                prj_mpi_decode_cell_index(coarse_code, &ii, &jj, &kk);
                srcface[0] = tmp_face[0];
                srcface[1] = tmp_face[1];
                srcface[2] = tmp_face[2];
                prj_mpi_unpack_bf_patch_values(tmp_face, ii, jj, kk, &buffer->face_buffer_recv[pos]);
                pos += PRJ_MPI_BF_PATCH_VALUE_COUNT;
                prj_mhd_apply_bf_patch(mesh, &mesh->blocks[block_id], stage, coarse, srcface, ii, jj, kk);
            }
        }
        free(tmp_face[2]);
        free(tmp_face[1]);
        free(tmp_face[0]);
    }
    free(requests);
#else
    (void)mesh;
    (void)mpi;
    (void)stage;
    (void)fill_kind;
#endif
}

void prj_mpi_exchange_emf(prj_mesh *mesh, prj_mpi *mpi)
{
#if defined(PRJ_ENABLE_MPI) && PRJ_MHD
    int nb;
    MPI_Request *requests;
    int request_count;
    int request_capacity;

    if (mesh == 0 || mpi == 0 || mpi->totrank <= 1) {
        return;
    }

    requests = 0;
    request_count = 0;
    request_capacity = 0;
    for (nb = 0; nb < mpi->neighbor_number; ++nb) {
        prj_mpi_buffer *buffer = &mpi->neighbor_buffer[nb];
        int *send_sizes;
        int recv_entries;
        int recv_total;
        int axis;
        int bidx;
        int occ;
        int *idx_send[3];
        double *value_send;
        int idx_count;
        int idx_capacity;
        int value_count;
        int value_capacity;

        send_sizes = (int *)calloc((size_t)buffer->number, sizeof(*send_sizes));
        idx_send[0] = 0;
        idx_send[1] = 0;
        idx_send[2] = 0;
        value_send = 0;
        idx_count = 0;
        idx_capacity = 0;
        value_count = 0;
        value_capacity = 0;
        occ = 0;
        for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
            prj_block *block = &mesh->blocks[bidx];
            int n;

            if (!prj_mpi_block_is_local(block)) {
                continue;
            }
            for (n = 0; n < 56; ++n) {
                int nid = block->slot[n].id;
                const prj_block *neighbor;
                int before = idx_count;

                if (nid < 0 || nid >= mesh->nblocks || mesh->blocks[nid].rank != buffer->receiver_rank) {
                    continue;
                }
                neighbor = &mesh->blocks[nid];
                if (block->slot[n].type == PRJ_NEIGHBOR_CORNER || block->level != neighbor->level + 1) {
                    send_sizes[occ++] = 0;
                    continue;
                }
                {
                    int dir;

                    for (dir = 0; dir < 3; ++dir) {
                        int start[3];
                        int end[3];
                        int i;
                        int j;
                        int k;

                        if (!prj_mpi_coarse_edge_bounds(block, neighbor, &block->slot[n], dir, 0, start, end)) {
                            continue;
                        }
                        for (i = start[0]; i <= end[0]; ++i) {
                            for (j = start[1]; j <= end[1]; ++j) {
                                for (k = start[2]; k <= end[2]; ++k) {
                                    int subidx;

                                    if (prj_mhd_edge_is_interior(dir, i, j, k)) {
                                        continue;
                                    }
                                    for (subidx = 0; subidx < 2; ++subidx) {
                                        double x[3];
                                        int ii;
                                        int jj;
                                        int kk;

                                        prj_mpi_coarse_subedge_position(neighbor, dir, i, j, k, subidx, x);
                                        if (!prj_mhd_block_owns_edge_position(block, dir, x, &ii, &jj, &kk)) {
                                            continue;
                                        }
                                        if (prj_mpi_append_triplet(idx_send, &idx_count, &idx_capacity,
                                                neighbor->id, block->id,
                                                prj_mpi_encode_edge_code(dir, i, j, k)) != 0 ||
                                            prj_mpi_append_double(&value_send, &value_count, &value_capacity,
                                                block->emf[dir][IDX(ii, jj, kk)]) != 0) {
                                            fprintf(stderr, "prj_mpi_exchange_emf: pack failed\n");
                                            abort();
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                send_sizes[occ++] = idx_count - before;
            }
        }

        free(buffer->face_data_size_send);
        free(buffer->face_buffer_send);
        for (axis = 0; axis < 3; ++axis) {
            free(buffer->face_data_idx_send[axis]);
            free(buffer->face_data_idx_recv[axis]);
            buffer->face_data_idx_send[axis] = idx_send[axis];
            buffer->face_data_idx_recv[axis] = 0;
        }
        free(buffer->face_data_size_recv);
        free(buffer->face_buffer_recv);
        buffer->face_data_size_send = send_sizes;
        buffer->face_buffer_send = value_send;
        buffer->face_data_size_recv = 0;
        buffer->face_buffer_recv = 0;

        MPI_Sendrecv(&buffer->number, 1, MPI_INT, buffer->receiver_rank, 340,
            &recv_entries, 1, MPI_INT, buffer->receiver_rank, 340, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        buffer->face_data_size_recv = (int *)calloc((size_t)recv_entries + 1U, sizeof(*buffer->face_data_size_recv));
        MPI_Sendrecv(buffer->face_data_size_send, buffer->number, MPI_INT, buffer->receiver_rank, 341,
            buffer->face_data_size_recv, recv_entries, MPI_INT, buffer->receiver_rank, 341, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        buffer->face_data_size_recv[recv_entries] = -1;
        recv_total = 0;
        for (axis = 0; axis < recv_entries; ++axis) {
            recv_total += buffer->face_data_size_recv[axis];
        }
        for (axis = 0; axis < 3; ++axis) {
            buffer->face_data_idx_recv[axis] = (int *)calloc((size_t)recv_total, sizeof(*buffer->face_data_idx_recv[axis]));
        }
        buffer->face_buffer_recv = (double *)calloc((size_t)recv_total, sizeof(*buffer->face_buffer_recv));

        if (request_count + 8 > request_capacity) {
            int new_capacity = request_capacity == 0 ? 16 : request_capacity * 2;
            MPI_Request *next = (MPI_Request *)realloc(requests, (size_t)new_capacity * sizeof(*next));

            if (next == 0) {
                free(requests);
                return;
            }
            requests = next;
            request_capacity = new_capacity;
        }
        for (axis = 0; axis < 3; ++axis) {
            MPI_Irecv(buffer->face_data_idx_recv[axis], recv_total, MPI_INT, buffer->receiver_rank, 342 + axis,
                MPI_COMM_WORLD, &requests[request_count++]);
            MPI_Isend(buffer->face_data_idx_send[axis], idx_count, MPI_INT, buffer->receiver_rank, 342 + axis,
                MPI_COMM_WORLD, &requests[request_count++]);
        }
        MPI_Irecv(buffer->face_buffer_recv, recv_total, MPI_DOUBLE, buffer->receiver_rank, 345,
            MPI_COMM_WORLD, &requests[request_count++]);
        MPI_Isend(buffer->face_buffer_send, value_count, MPI_DOUBLE, buffer->receiver_rank, 345,
            MPI_COMM_WORLD, &requests[request_count++]);
    }
    if (request_count > 0) {
        MPI_Waitall(request_count, requests, MPI_STATUSES_IGNORE);
    }
    {
        int *acc_block = 0;
        int *acc_edge = 0;
        int *acc_count = 0;
        double *acc_sum = 0;
        int acc_n = 0;
        int acc_cap = 0;

        for (nb = 0; nb < mpi->neighbor_number; ++nb) {
            prj_mpi_buffer *buffer = &mpi->neighbor_buffer[nb];
            int total = 0;
            int pos = 0;
            int irec;

            if (buffer->face_data_size_recv == 0) {
                continue;
            }
            for (irec = 0; buffer->face_data_size_recv[irec] >= 0; ++irec) {
                total += buffer->face_data_size_recv[irec];
            }
            for (irec = 0; irec < total; ++irec) {
                int block_id = buffer->face_data_idx_recv[0][irec];
                int source_id = buffer->face_data_idx_recv[1][irec];
                int edge_code = buffer->face_data_idx_recv[2][irec];
                double value = buffer->face_buffer_recv[pos++];
                int found = -1;
                int a;

                if (block_id < 0 || block_id >= mesh->nblocks ||
                    source_id < 0 || source_id >= mesh->nblocks) {
                    continue;
                }
                if (mesh->blocks[block_id].rank != mpi->rank) {
                    continue;
                }
                if (mesh->blocks[source_id].level != mesh->blocks[block_id].level + 1) {
                    fprintf(stderr, "prj_mpi_exchange_emf: invalid level relation\n");
                    abort();
                }
                for (a = 0; a < acc_n; ++a) {
                    if (acc_block[a] == block_id && acc_edge[a] == edge_code) {
                        found = a;
                        break;
                    }
                }
                if (found < 0) {
                    if (acc_n >= acc_cap) {
                        int new_cap = acc_cap == 0 ? 32 : 2 * acc_cap;
                        int *next_block = (int *)realloc(acc_block, (size_t)new_cap * sizeof(*next_block));
                        int *next_edge = (int *)realloc(acc_edge, (size_t)new_cap * sizeof(*next_edge));
                        int *next_count = (int *)realloc(acc_count, (size_t)new_cap * sizeof(*next_count));
                        double *next_sum = (double *)realloc(acc_sum, (size_t)new_cap * sizeof(*next_sum));

                        if (next_block == 0 || next_edge == 0 || next_count == 0 || next_sum == 0) {
                            free(next_sum);
                            free(next_count);
                            free(next_edge);
                            free(next_block);
                            free(acc_sum);
                            free(acc_count);
                            free(acc_edge);
                            free(acc_block);
                            free(requests);
                            return;
                        }
                        acc_block = next_block;
                        acc_edge = next_edge;
                        acc_count = next_count;
                        acc_sum = next_sum;
                        acc_cap = new_cap;
                    }
                    found = acc_n++;
                    acc_block[found] = block_id;
                    acc_edge[found] = edge_code;
                    acc_count[found] = 0;
                    acc_sum[found] = 0.0;
                }
                acc_count[found] += 1;
                acc_sum[found] += value;
            }
        }
        for (nb = 0; nb < acc_n; ++nb) {
            int dir;
            int ii;
            int jj;
            int kk;
            prj_block *block = &mesh->blocks[acc_block[nb]];
            double local_sum;
            int local_count;

            prj_mpi_decode_edge_code(acc_edge[nb], &dir, &ii, &jj, &kk);
            prj_mpi_local_edge_partial(mesh, block, dir, ii, jj, kk, &local_sum, &local_count);
            if (acc_count[nb] + local_count == 2) {
                block->emf[dir][IDX(ii, jj, kk)] = 0.5 * (acc_sum[nb] + local_sum);
            }
        }
        free(acc_sum);
        free(acc_count);
        free(acc_edge);
        free(acc_block);
    }
    free(requests);
#else
    (void)mesh;
    (void)mpi;
#endif
}

void prj_mpi_exchange_fluxes(prj_mesh *mesh, prj_mpi *mpi)
{
#if defined(PRJ_ENABLE_MPI)
    int nb;
    MPI_Request *requests;
    int request_count;
    int request_capacity;

    if (mesh == 0 || mpi == 0 || mpi->totrank <= 1) {
        return;
    }
    requests = 0;
    request_count = 0;
    request_capacity = 0;
    for (nb = 0; nb < mpi->neighbor_number; ++nb) {
        prj_mpi_buffer *buffer = &mpi->neighbor_buffer[nb];
        int *send_sizes;
        int recv_entries;
        int face_total_recv;
        int axis;
        int bidx;
        int occ;
        int *idx_send[3];
        double *value_send;
        int idx_count;
        int idx_capacity;
        int value_count;
        int value_capacity;

        send_sizes = (int *)calloc((size_t)buffer->number, sizeof(*send_sizes));
        idx_send[0] = 0;
        idx_send[1] = 0;
        idx_send[2] = 0;
        value_send = 0;
        idx_count = 0;
        idx_capacity = 0;
        value_count = 0;
        value_capacity = 0;
        occ = 0;
        for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
            prj_block *block = &mesh->blocks[bidx];
            int n;

            if (!prj_mpi_block_is_local(block)) {
                continue;
            }
            for (n = 0; n < 56; ++n) {
                int nid = block->slot[n].id;
                const prj_block *neighbor;
                int face_axis;
                int side;
                int d;
                int before = idx_count;

                if (nid < 0 || nid >= mesh->nblocks || mesh->blocks[nid].rank != buffer->receiver_rank) {
                    continue;
                }
                neighbor = &mesh->blocks[nid];
                face_axis = -1;
                side = -1;
                for (d = 0; d < 3; ++d) {
                    if (fabs(block->xmax[d] - neighbor->xmin[d]) < 1.0e-12) {
                        face_axis = d;
                        side = 1;
                        break;
                    }
                    if (fabs(neighbor->xmax[d] - block->xmin[d]) < 1.0e-12) {
                        face_axis = d;
                        side = 0;
                        break;
                    }
                }
                if (face_axis < 0 || neighbor->dx[face_axis] <= block->dx[face_axis]) {
                    send_sizes[occ++] = 0;
                    continue;
                }
                for (d = 0; d < PRJ_BLOCK_SIZE; ++d) {
                    int e;

                    for (e = 0; e < PRJ_BLOCK_SIZE; ++e) {
                        int face[3] = {0, 0, 0};
                        int v;

                        face[face_axis] = side == 1 ? PRJ_BLOCK_SIZE : 0;
                        face[(face_axis + 1) % 3] = d;
                        face[(face_axis + 2) % 3] = e;
                        if (prj_mpi_append_triplet(idx_send, &idx_count, &idx_capacity,
                                nid, prj_mpi_encode_cell_index(face[0], face[1], face[2]),
                                face_axis) != 0) {
                            continue;
                        }
                        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
                            double flux = block->flux[face_axis][VIDX(v, face[0], face[1], face[2])];

                            if (prj_mpi_append_double(&value_send, &value_count, &value_capacity, flux) != 0) {
                                break;
                            }
                        }
                    }
                }
                send_sizes[occ++] = idx_count - before;
            }
        }

        free(buffer->face_data_size_send);
        free(buffer->face_buffer_send);
        for (axis = 0; axis < 3; ++axis) {
            free(buffer->face_data_idx_send[axis]);
            free(buffer->face_data_idx_recv[axis]);
            buffer->face_data_idx_send[axis] = idx_send[axis];
            buffer->face_data_idx_recv[axis] = 0;
        }
        free(buffer->face_data_size_recv);
        free(buffer->face_buffer_recv);
        buffer->face_data_size_send = send_sizes;
        buffer->face_buffer_send = value_send;
        buffer->face_data_size_recv = 0;
        buffer->face_buffer_recv = 0;

        MPI_Sendrecv(&buffer->number, 1, MPI_INT, buffer->receiver_rank, 200,
            &recv_entries, 1, MPI_INT, buffer->receiver_rank, 200, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        buffer->face_data_size_recv = (int *)calloc((size_t)recv_entries, sizeof(*buffer->face_data_size_recv));
        MPI_Sendrecv(buffer->face_data_size_send, buffer->number, MPI_INT, buffer->receiver_rank, 201,
            buffer->face_data_size_recv, recv_entries, MPI_INT, buffer->receiver_rank, 201, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        face_total_recv = 0;
        for (axis = 0; axis < recv_entries; ++axis) {
            face_total_recv += buffer->face_data_size_recv[axis];
        }
        for (axis = 0; axis < 3; ++axis) {
            buffer->face_data_idx_recv[axis] = (int *)calloc((size_t)face_total_recv, sizeof(*buffer->face_data_idx_recv[axis]));
        }
        buffer->face_buffer_recv = (double *)calloc((size_t)face_total_recv * PRJ_NVAR_CONS, sizeof(*buffer->face_buffer_recv));

        if (request_count + 8 > request_capacity) {
            int new_capacity = request_capacity == 0 ? 16 : request_capacity * 2;
            MPI_Request *next = (MPI_Request *)realloc(requests, (size_t)new_capacity * sizeof(*next));

            if (next == 0) {
                free(requests);
                return;
            }
            requests = next;
            request_capacity = new_capacity;
        }
        for (axis = 0; axis < 3; ++axis) {
            MPI_Irecv(buffer->face_data_idx_recv[axis], face_total_recv, MPI_INT, buffer->receiver_rank, 210 + axis,
                MPI_COMM_WORLD, &requests[request_count++]);
            MPI_Isend(buffer->face_data_idx_send[axis], idx_count, MPI_INT, buffer->receiver_rank, 210 + axis,
                MPI_COMM_WORLD, &requests[request_count++]);
        }
        MPI_Irecv(buffer->face_buffer_recv, face_total_recv * PRJ_NVAR_CONS, MPI_DOUBLE, buffer->receiver_rank, 220,
            MPI_COMM_WORLD, &requests[request_count++]);
        MPI_Isend(buffer->face_buffer_send, value_count, MPI_DOUBLE, buffer->receiver_rank, 220,
            MPI_COMM_WORLD, &requests[request_count++]);
    }
    if (request_count > 0) {
        MPI_Waitall(request_count, requests, MPI_STATUSES_IGNORE);
    }
    for (nb = 0; nb < mpi->neighbor_number; ++nb) {
        prj_mpi_buffer *buffer = &mpi->neighbor_buffer[nb];
        int total = 0;
        int pos = 0;
        int i;

        if (buffer->face_data_size_recv == 0) {
            continue;
        }
        for (i = 0; i < buffer->number; ++i) {
            total += buffer->face_data_size_recv[i];
        }
        for (i = 0; i < total; ++i) {
            int block_id = buffer->face_data_idx_recv[0][i];
            int code = buffer->face_data_idx_recv[1][i];
            int face_axis = buffer->face_data_idx_recv[2][i];
            int ii;
            int jj;
            int kk;
            int v;
            prj_block *block;

            if (block_id < 0 || block_id >= mesh->nblocks || face_axis < 0 || face_axis >= 3) {
                pos += PRJ_NVAR_CONS;
                continue;
            }
            block = &mesh->blocks[block_id];
            if (!prj_mpi_block_is_local(block) || block->flux[face_axis] == 0) {
                pos += PRJ_NVAR_CONS;
                continue;
            }
            prj_mpi_decode_cell_index(code, &ii, &jj, &kk);
            for (v = 0; v < PRJ_NVAR_CONS; ++v) {
                block->flux[face_axis][VIDX(v, ii, jj, kk)] = buffer->face_buffer_recv[pos++];
            }
        }
    }
    free(requests);
#else
    (void)mesh;
    (void)mpi;
#endif
}

double prj_mpi_min_dt(double local_dt)
{
#if defined(PRJ_ENABLE_MPI)
    double global_dt = local_dt;

    if (prj_mpi_active != 0 && prj_mpi_active->totrank > 1) {
        MPI_Allreduce(&local_dt, &global_dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    }
    return global_dt;
#else
    return local_dt;
#endif
}

double prj_mpi_global_sum(double local_val)
{
#if defined(PRJ_ENABLE_MPI)
    double global_val = local_val;

    if (prj_mpi_active != 0 && prj_mpi_active->totrank > 1) {
        MPI_Allreduce(&local_val, &global_val, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    return global_val;
#else
    return local_val;
#endif
}

void prj_mpi_rebalance(prj_mesh *mesh)
{
#if defined(PRJ_ENABLE_MPI)
    int *old_ranks;
    int *counts_before;
    int *counts_after;
    int did_rebalance;
    int bidx;

    if (mesh == 0 || prj_mpi_active == 0) {
        return;
    }
    old_ranks = (int *)calloc((size_t)mesh->nblocks, sizeof(*old_ranks));
    counts_before = (int *)calloc((size_t)prj_mpi_active->totrank, sizeof(*counts_before));
    counts_after = (int *)calloc((size_t)prj_mpi_active->totrank, sizeof(*counts_after));
    if (old_ranks == 0 || counts_before == 0 || counts_after == 0) {
        free(counts_after);
        free(counts_before);
        free(old_ranks);
        return;
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        old_ranks[bidx] = mesh->blocks[bidx].rank;
    }
    prj_mpi_collect_active_counts(mesh, counts_before);
    prj_mpi_compute_decomposition(mesh);
    did_rebalance = prj_mpi_has_rebalanced(mesh, old_ranks);
    prj_mpi_collect_active_counts(mesh, counts_after);
    if (did_rebalance) {
        prj_mpi_migrate_active_blocks(mesh, old_ranks);
    }
    prj_mpi_assign_block_storage(mesh);
    prj_mpi_prepare(mesh, prj_mpi_active);
    if (did_rebalance &&
            prj_mpi_counts_changed(counts_before, counts_after, prj_mpi_active->totrank)) {
        prj_mpi_print_balance(mesh);
    }
    free(counts_after);
    free(counts_before);
    free(old_ranks);
#else
    if (mesh != 0) {
        prj_mpi_decompose(mesh);
    }
#endif
}

int prj_mpi_get_neighbor_rank(const prj_mesh *mesh, int neighbor_block_id)
{
    if (mesh == 0 || neighbor_block_id < 0 || neighbor_block_id >= mesh->nblocks) {
        return -1;
    }
    return mesh->blocks[neighbor_block_id].rank;
}

void prj_mpi_finalize(void)
{
    if (prj_mpi_active != 0) {
        prj_mpi_clear_neighbors(prj_mpi_active);
    }
#if defined(PRJ_ENABLE_MPI)
    MPI_Finalize();
#endif
    prj_mpi_active = 0;
}
