#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(PRJ_ENABLE_MPI)
#include <mpi.h>
#endif

#include "prj.h"

#define PRJ_MHD_FOUR_PI 12.56637061435917295385

static void prj_mhd_abort(const char *message)
{
    prj_mpi *mpi = prj_mpi_current();

    fprintf(stderr, "%s\n", message);
#if defined(PRJ_ENABLE_MPI)
    if (mpi != 0 && mpi->totrank > 1) {
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
#endif
    abort();
}

void prj_mhd_face_patch_clear(prj_mhd_face_patch *patch)
{
    if (patch == 0) {
        prj_mhd_abort("prj_mhd_face_patch_clear: patch is null");
    }
    memset(patch, 0, sizeof(*patch));
}

static int prj_mhd_index_near_integer(double value, int *idx)
{
    int rounded;

    if (idx == 0) {
        prj_mhd_abort("prj_mhd_index_near_integer: idx is null");
    }
    rounded = value >= 0.0 ? (int)(value + 0.5) : (int)(value - 0.5);
    if (fabs(value - (double)rounded) > 1.0e-8) {
        return 0;
    }
    *idx = rounded;
    return 1;
}

#if PRJ_MHD
static int prj_mhd_encode_cell_index(int i, int j, int k)
{
    return ((i + PRJ_NGHOST) * PRJ_BS + (j + PRJ_NGHOST)) * PRJ_BS + (k + PRJ_NGHOST);
}

static void prj_mhd_decode_cell_index(int code, int *i, int *j, int *k)
{
    if (i == 0 || j == 0 || k == 0) {
        prj_mhd_abort("prj_mhd_decode_cell_index: null output pointer");
    }
    *k = code % PRJ_BS - PRJ_NGHOST;
    code /= PRJ_BS;
    *j = code % PRJ_BS - PRJ_NGHOST;
    code /= PRJ_BS;
    *i = code - PRJ_NGHOST;
}

static int prj_mhd_encode_edge_code(int dir, int i, int j, int k)
{
    return dir * PRJ_BLOCK_NCELLS + prj_mhd_encode_cell_index(i, j, k);
}

static void prj_mhd_decode_edge_code(int code, int *dir, int *i, int *j, int *k)
{
    if (dir == 0 || i == 0 || j == 0 || k == 0) {
        prj_mhd_abort("prj_mhd_decode_edge_code: null output pointer");
    }
    *dir = code / PRJ_BLOCK_NCELLS;
    prj_mhd_decode_cell_index(code % PRJ_BLOCK_NCELLS, i, j, k);
}

static int prj_mhd_append_double(double **array, int *count, int *capacity, double value)
{
    double *next;
    int new_capacity;

    if (array == 0 || count == 0 || capacity == 0) {
        prj_mhd_abort("prj_mhd_append_double: null input");
    }
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

static int prj_mhd_append_triplet(int **arrays, int *count, int *capacity, int a, int b, int c)
{
    int new_capacity;
    int *next0;
    int *next1;
    int *next2;

    if (arrays == 0 || count == 0 || capacity == 0) {
        prj_mhd_abort("prj_mhd_append_triplet: null input");
    }
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

static int prj_mhd_min_int(int a, int b)
{
    return a < b ? a : b;
}

static int prj_mhd_max_int(int a, int b)
{
    return a > b ? a : b;
}
#endif

double *prj_mhd_bface_stage(prj_block *block, int stage, int dir)
{
    if (block == 0 || dir < 0 || dir >= 3) {
        prj_mhd_abort("prj_mhd_bface_stage: invalid input");
    }
    if (stage == 1) {
        return block->Bf[dir];
    }
    if (stage == 2) {
        return block->Bf1[dir];
    }
    prj_mhd_abort("prj_mhd_bface_stage: invalid RK stage");
    return 0;
}

const double *prj_mhd_bface_stage_const(const prj_block *block, int stage, int dir)
{
    if (block == 0 || dir < 0 || dir >= 3) {
        prj_mhd_abort("prj_mhd_bface_stage_const: invalid input");
    }
    if (stage == 1) {
        return block->Bf[dir];
    }
    if (stage == 2) {
        return block->Bf1[dir];
    }
    prj_mhd_abort("prj_mhd_bface_stage_const: invalid RK stage");
    return 0;
}

int prj_mhd_face_is_interior(int dir, int i, int j, int k)
{
    if (dir == X1DIR) {
        return i > 0 && i < PRJ_BLOCK_SIZE &&
            j >= 0 && j < PRJ_BLOCK_SIZE &&
            k >= 0 && k < PRJ_BLOCK_SIZE;
    }
    if (dir == X2DIR) {
        return j > 0 && j < PRJ_BLOCK_SIZE &&
            i >= 0 && i < PRJ_BLOCK_SIZE &&
            k >= 0 && k < PRJ_BLOCK_SIZE;
    }
    if (dir == X3DIR) {
        return k > 0 && k < PRJ_BLOCK_SIZE &&
            i >= 0 && i < PRJ_BLOCK_SIZE &&
            j >= 0 && j < PRJ_BLOCK_SIZE;
    }
    prj_mhd_abort("prj_mhd_face_is_interior: invalid direction");
    return 0;
}

int prj_mhd_face_is_owned(int dir, int i, int j, int k)
{
    if (dir == X1DIR) {
        return i >= 0 && i <= PRJ_BLOCK_SIZE &&
            j >= 0 && j < PRJ_BLOCK_SIZE &&
            k >= 0 && k < PRJ_BLOCK_SIZE;
    }
    if (dir == X2DIR) {
        return j >= 0 && j <= PRJ_BLOCK_SIZE &&
            i >= 0 && i < PRJ_BLOCK_SIZE &&
            k >= 0 && k < PRJ_BLOCK_SIZE;
    }
    if (dir == X3DIR) {
        return k >= 0 && k <= PRJ_BLOCK_SIZE &&
            i >= 0 && i < PRJ_BLOCK_SIZE &&
            j >= 0 && j < PRJ_BLOCK_SIZE;
    }
    prj_mhd_abort("prj_mhd_face_is_owned: invalid direction");
    return 0;
}

void prj_mhd_face_position(const prj_block *block, int dir, int i, int j, int k, double x[3])
{
    if (block == 0 || x == 0) {
        prj_mhd_abort("prj_mhd_face_position: null input");
    }
    if (dir == X1DIR) {
        x[0] = block->xmin[0] + (double)i * block->dx[0];
        x[1] = block->xmin[1] + ((double)j + 0.5) * block->dx[1];
        x[2] = block->xmin[2] + ((double)k + 0.5) * block->dx[2];
        return;
    }
    if (dir == X2DIR) {
        x[0] = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
        x[1] = block->xmin[1] + (double)j * block->dx[1];
        x[2] = block->xmin[2] + ((double)k + 0.5) * block->dx[2];
        return;
    }
    if (dir == X3DIR) {
        x[0] = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
        x[1] = block->xmin[1] + ((double)j + 0.5) * block->dx[1];
        x[2] = block->xmin[2] + (double)k * block->dx[2];
        return;
    }
    prj_mhd_abort("prj_mhd_face_position: invalid direction");
}

int prj_mhd_face_index_from_position(const prj_block *block, int dir, const double x[3], int *i, int *j, int *k)
{
    double rx;
    double ry;
    double rz;

    if (block == 0 || x == 0 || i == 0 || j == 0 || k == 0) {
        prj_mhd_abort("prj_mhd_face_index_from_position: null input");
    }
    if (dir == X1DIR) {
        rx = (x[0] - block->xmin[0]) / block->dx[0];
        ry = (x[1] - block->xmin[1]) / block->dx[1] - 0.5;
        rz = (x[2] - block->xmin[2]) / block->dx[2] - 0.5;
    } else if (dir == X2DIR) {
        rx = (x[0] - block->xmin[0]) / block->dx[0] - 0.5;
        ry = (x[1] - block->xmin[1]) / block->dx[1];
        rz = (x[2] - block->xmin[2]) / block->dx[2] - 0.5;
    } else if (dir == X3DIR) {
        rx = (x[0] - block->xmin[0]) / block->dx[0] - 0.5;
        ry = (x[1] - block->xmin[1]) / block->dx[1] - 0.5;
        rz = (x[2] - block->xmin[2]) / block->dx[2];
    } else {
        prj_mhd_abort("prj_mhd_face_index_from_position: invalid direction");
        return 0;
    }
    return prj_mhd_index_near_integer(rx, i) &&
        prj_mhd_index_near_integer(ry, j) &&
        prj_mhd_index_near_integer(rz, k);
}

int prj_mhd_block_owns_face_position(const prj_block *block, int dir, const double x[3], int *i, int *j, int *k)
{
    int ii;
    int jj;
    int kk;

    if (block == 0) {
        return 0;
    }
    if (!prj_mhd_face_index_from_position(block, dir, x, &ii, &jj, &kk)) {
        return 0;
    }
    if (!prj_mhd_face_is_owned(dir, ii, jj, kk)) {
        return 0;
    }
    if (i != 0) {
        *i = ii;
    }
    if (j != 0) {
        *j = jj;
    }
    if (k != 0) {
        *k = kk;
    }
    return 1;
}

int prj_mhd_edge_is_interior(int dir, int i, int j, int k)
{
    if (dir == X1DIR) {
        return i >= 0 && i < PRJ_BLOCK_SIZE &&
            j > 0 && j < PRJ_BLOCK_SIZE &&
            k > 0 && k < PRJ_BLOCK_SIZE;
    }
    if (dir == X2DIR) {
        return j >= 0 && j < PRJ_BLOCK_SIZE &&
            i > 0 && i < PRJ_BLOCK_SIZE &&
            k > 0 && k < PRJ_BLOCK_SIZE;
    }
    if (dir == X3DIR) {
        return k >= 0 && k < PRJ_BLOCK_SIZE &&
            i > 0 && i < PRJ_BLOCK_SIZE &&
            j > 0 && j < PRJ_BLOCK_SIZE;
    }
    prj_mhd_abort("prj_mhd_edge_is_interior: invalid direction");
    return 0;
}

int prj_mhd_edge_is_owned(int dir, int i, int j, int k)
{
    if (dir == X1DIR) {
        return i >= 0 && i < PRJ_BLOCK_SIZE &&
            j >= 0 && j <= PRJ_BLOCK_SIZE &&
            k >= 0 && k <= PRJ_BLOCK_SIZE;
    }
    if (dir == X2DIR) {
        return j >= 0 && j < PRJ_BLOCK_SIZE &&
            i >= 0 && i <= PRJ_BLOCK_SIZE &&
            k >= 0 && k <= PRJ_BLOCK_SIZE;
    }
    if (dir == X3DIR) {
        return k >= 0 && k < PRJ_BLOCK_SIZE &&
            i >= 0 && i <= PRJ_BLOCK_SIZE &&
            j >= 0 && j <= PRJ_BLOCK_SIZE;
    }
    prj_mhd_abort("prj_mhd_edge_is_owned: invalid direction");
    return 0;
}

void prj_mhd_edge_position(const prj_block *block, int dir, int i, int j, int k, double x[3])
{
    if (block == 0 || x == 0) {
        prj_mhd_abort("prj_mhd_edge_position: null input");
    }

    if (dir == X1DIR) {
        x[0] = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
        x[1] = block->xmin[1] + (double)j * block->dx[1];
        x[2] = block->xmin[2] + (double)k * block->dx[2];
        return;
    }
    if (dir == X2DIR) {
        x[0] = block->xmin[0] + (double)i * block->dx[0];
        x[1] = block->xmin[1] + ((double)j + 0.5) * block->dx[1];
        x[2] = block->xmin[2] + (double)k * block->dx[2];
        return;
    }
    if (dir == X3DIR) {
        x[0] = block->xmin[0] + (double)i * block->dx[0];
        x[1] = block->xmin[1] + (double)j * block->dx[1];
        x[2] = block->xmin[2] + ((double)k + 0.5) * block->dx[2];
        return;
    }
    prj_mhd_abort("prj_mhd_edge_position: invalid direction");
}

int prj_mhd_edge_index_from_position(const prj_block *block, int dir, const double x[3], int *i, int *j, int *k)
{
    double rx;
    double ry;
    double rz;

    if (block == 0 || x == 0 || i == 0 || j == 0 || k == 0) {
        prj_mhd_abort("prj_mhd_edge_index_from_position: null input");
    }
    if (dir == X1DIR) {
        rx = (x[0] - block->xmin[0]) / block->dx[0] - 0.5;
        ry = (x[1] - block->xmin[1]) / block->dx[1];
        rz = (x[2] - block->xmin[2]) / block->dx[2];
    } else if (dir == X2DIR) {
        rx = (x[0] - block->xmin[0]) / block->dx[0];
        ry = (x[1] - block->xmin[1]) / block->dx[1] - 0.5;
        rz = (x[2] - block->xmin[2]) / block->dx[2];
    } else if (dir == X3DIR) {
        rx = (x[0] - block->xmin[0]) / block->dx[0];
        ry = (x[1] - block->xmin[1]) / block->dx[1];
        rz = (x[2] - block->xmin[2]) / block->dx[2] - 0.5;
    } else {
        prj_mhd_abort("prj_mhd_edge_index_from_position: invalid direction");
        return 0;
    }
    return prj_mhd_index_near_integer(rx, i) &&
        prj_mhd_index_near_integer(ry, j) &&
        prj_mhd_index_near_integer(rz, k);
}

int prj_mhd_block_owns_edge_position(const prj_block *block, int dir, const double x[3], int *i, int *j, int *k)
{
    int ii;
    int jj;
    int kk;

    if (block == 0) {
        return 0;
    }
    if (!prj_mhd_edge_index_from_position(block, dir, x, &ii, &jj, &kk)) {
        return 0;
    }
    if (!prj_mhd_edge_is_owned(dir, ii, jj, kk)) {
        return 0;
    }
    if (i != 0) {
        *i = ii;
    }
    if (j != 0) {
        *j = jj;
    }
    if (k != 0) {
        *k = kk;
    }
    return 1;
}

#if PRJ_MHD
static int prj_mhd_storage_index_valid(int i, int j, int k)
{
    return i >= -PRJ_NGHOST && i < PRJ_BLOCK_SIZE + PRJ_NGHOST &&
        j >= -PRJ_NGHOST && j < PRJ_BLOCK_SIZE + PRJ_NGHOST &&
        k >= -PRJ_NGHOST && k < PRJ_BLOCK_SIZE + PRJ_NGHOST;
}

static int prj_mhd_local_block(const prj_block *block)
{
    prj_mpi *mpi = prj_mpi_current();

    return block != 0 && block->id >= 0 && block->active == 1 &&
        (mpi == 0 || block->rank == mpi->rank);
}

static int prj_mhd_active_block(const prj_block *block)
{
    return block != 0 && block->id >= 0 && block->active == 1;
}

static void prj_mhd_patch_face_position(const prj_block *coarse, int dir, int i, int j, int k,
    int nidx, int aidx, int bidx, double x[3])
{
    if (coarse == 0 || x == 0) {
        prj_mhd_abort("prj_mhd_patch_face_position: null input");
    }
    if (nidx < 0 || nidx > 2 || aidx < 0 || aidx > 1 || bidx < 0 || bidx > 1) {
        prj_mhd_abort("prj_mhd_patch_face_position: invalid patch index");
    }

    if (dir == X1DIR) {
        x[0] = coarse->xmin[0] + ((double)i + 0.5 * (double)nidx) * coarse->dx[0];
        x[1] = coarse->xmin[1] + ((double)j + 0.25 + 0.5 * (double)aidx) * coarse->dx[1];
        x[2] = coarse->xmin[2] + ((double)k + 0.25 + 0.5 * (double)bidx) * coarse->dx[2];
        return;
    }
    if (dir == X2DIR) {
        x[0] = coarse->xmin[0] + ((double)i + 0.25 + 0.5 * (double)aidx) * coarse->dx[0];
        x[1] = coarse->xmin[1] + ((double)j + 0.5 * (double)nidx) * coarse->dx[1];
        x[2] = coarse->xmin[2] + ((double)k + 0.25 + 0.5 * (double)bidx) * coarse->dx[2];
        return;
    }
    if (dir == X3DIR) {
        x[0] = coarse->xmin[0] + ((double)i + 0.25 + 0.5 * (double)aidx) * coarse->dx[0];
        x[1] = coarse->xmin[1] + ((double)j + 0.25 + 0.5 * (double)bidx) * coarse->dx[1];
        x[2] = coarse->xmin[2] + ((double)k + 0.5 * (double)nidx) * coarse->dx[2];
        return;
    }
    prj_mhd_abort("prj_mhd_patch_face_position: invalid direction");
}

static int prj_mhd_face_has_fine_owner(const prj_mesh *mesh, const prj_block *fine, int skip_block_id,
    int dir, const double x[3])
{
    int n;

    if (mesh == 0 || fine == 0 || x == 0) {
        prj_mhd_abort("prj_mhd_face_has_fine_owner: null input");
    }
    if (fine->level < 0) {
        prj_mhd_abort("prj_mhd_face_has_fine_owner: invalid fine level");
    }

    if (fine->id != skip_block_id &&
        prj_mhd_active_block(fine) &&
        prj_mhd_block_owns_face_position(fine, dir, x, 0, 0, 0)) {
        return 1;
    }
    for (n = 0; n < 56; ++n) {
        int nid = fine->slot[n].id;
        const prj_block *neighbor;

        if (nid < 0 || nid >= mesh->nblocks || nid == skip_block_id) {
            continue;
        }
        neighbor = &mesh->blocks[nid];
        if (!prj_mhd_active_block(neighbor) || neighbor->level < fine->level) {
            continue;
        }
        if (prj_mhd_block_owns_face_position(neighbor, dir, x, 0, 0, 0)) {
            return 1;
        }
    }
    return 0;
}

static void prj_mhd_coarse_subedge_position(const prj_block *block, int dir, int i, int j, int k,
    int subidx, double x[3])
{
    if (block == 0 || x == 0) {
        prj_mhd_abort("prj_mhd_coarse_subedge_position: null input");
    }
    if (subidx < 0 || subidx > 1) {
        prj_mhd_abort("prj_mhd_coarse_subedge_position: invalid sub-edge index");
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
    prj_mhd_abort("prj_mhd_coarse_subedge_position: invalid direction");
}

static int prj_mhd_axis_touches(const prj_block *block, const prj_block *neighbor, int axis)
{
    const double tol = 1.0e-12;

    if (block == 0 || neighbor == 0 || axis < 0 || axis >= 3) {
        prj_mhd_abort("prj_mhd_axis_touches: invalid input");
    }
    return fabs(block->xmax[axis] - neighbor->xmin[axis]) < tol ||
        fabs(neighbor->xmax[axis] - block->xmin[axis]) < tol;
}

static int prj_mhd_axis_side(const prj_block *block, const prj_block *neighbor, int axis)
{
    const double tol = 1.0e-12;

    if (block == 0 || neighbor == 0 || axis < 0 || axis >= 3) {
        prj_mhd_abort("prj_mhd_axis_side: invalid input");
    }
    if (fabs(block->xmax[axis] - neighbor->xmin[axis]) < tol) {
        return 1;
    }
    if (fabs(neighbor->xmax[axis] - block->xmin[axis]) < tol) {
        return -1;
    }
    return 0;
}

static void prj_mhd_edge_slot_bounds(const prj_block *block, const prj_block *neighbor,
    const prj_neighbor *slot, int use_send_start, int interior_span, int start[3], int end[3])
{
    int axis;

    if (block == 0 || neighbor == 0 || slot == 0 || start == 0 || end == 0) {
        prj_mhd_abort("prj_mhd_edge_slot_bounds: invalid input");
    }
    for (axis = 0; axis < 3; ++axis) {
        int anchor = use_send_start != 0 ? slot->send_loc_start[axis] : slot->recv_loc_start[axis];

        if (prj_mhd_axis_touches(block, neighbor, axis)) {
            start[axis] = anchor == 0 ? -PRJ_NGHOST :
                prj_mhd_max_int(-PRJ_NGHOST, anchor - 1);
            end[axis] = anchor == 0 ?
                prj_mhd_min_int(PRJ_BLOCK_SIZE + PRJ_NGHOST - 1, PRJ_NGHOST) :
                (PRJ_BLOCK_SIZE + PRJ_NGHOST - 1);
        } else {
            start[axis] = prj_mhd_max_int(-PRJ_NGHOST, anchor - 1);
            end[axis] = prj_mhd_min_int(PRJ_BLOCK_SIZE + PRJ_NGHOST - 1,
                anchor + interior_span + 1);
        }
    }
}

static int prj_mhd_coarse_edge_bounds(const prj_block *block, const prj_block *neighbor,
    const prj_neighbor *slot, int dir, int use_send_start, int start[3], int end[3])
{
    int axis;
    int touch_count = 0;
    int overlap_axis = -1;
    int side[3];

    if (block == 0 || neighbor == 0 || slot == 0 || start == 0 || end == 0) {
        prj_mhd_abort("prj_mhd_coarse_edge_bounds: invalid input");
    }
    for (axis = 0; axis < 3; ++axis) {
        side[axis] = prj_mhd_axis_side(block, neighbor, axis);
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

static int prj_mhd_edge_in_bounds(int i, int j, int k, const int start[3], const int end[3])
{
    return i >= start[0] && i <= end[0] &&
        j >= start[1] && j <= end[1] &&
        k >= start[2] && k <= end[2];
}

static int prj_mhd_edge_is_ct_active(int dir, int i, int j, int k)
{
    if (dir == X1DIR) {
        return i >= 0 && i < PRJ_BLOCK_SIZE &&
            j >= 0 && j <= PRJ_BLOCK_SIZE &&
            k >= 0 && k <= PRJ_BLOCK_SIZE;
    }
    if (dir == X2DIR) {
        return i >= 0 && i <= PRJ_BLOCK_SIZE &&
            j >= 0 && j < PRJ_BLOCK_SIZE &&
            k >= 0 && k <= PRJ_BLOCK_SIZE;
    }
    if (dir == X3DIR) {
        return i >= 0 && i <= PRJ_BLOCK_SIZE &&
            j >= 0 && j <= PRJ_BLOCK_SIZE &&
            k >= 0 && k < PRJ_BLOCK_SIZE;
    }
    return 0;
}

static void prj_mhd_local_edge_partial(const prj_mesh *mesh, const prj_block *block, int dir,
    int i, int j, int k, double *sum, int *count)
{
    int n;

    if (mesh == 0 || block == 0 || sum == 0 || count == 0) {
        prj_mhd_abort("prj_mhd_local_edge_partial: invalid input");
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
            !prj_mhd_local_block(neighbor) || neighbor->level != block->level + 1) {
            continue;
        }
        if (!prj_mhd_coarse_edge_bounds(block, neighbor, &block->slot[n], dir, 1, start, end)) {
            continue;
        }
        if (!prj_mhd_edge_in_bounds(i, j, k, start, end)) {
            continue;
        }
        for (subidx = 0; subidx < 2; ++subidx) {
            double x[3];
            int ii;
            int jj;
            int kk;

            prj_mhd_coarse_subedge_position(block, dir, i, j, k, subidx, x);
            if (!prj_mhd_block_owns_edge_position(neighbor, dir, x, &ii, &jj, &kk)) {
                continue;
            }
            *sum += neighbor->emf[dir][IDX(ii, jj, kk)];
            *count += 1;
        }
    }
}

static void prj_mhd_debug_check_same_level_edges(const prj_mesh *mesh, const prj_block *block, int dir)
{
    int n;

    if (mesh == 0 || block == 0) {
        prj_mhd_abort("prj_mhd_debug_check_same_level_edges: invalid input");
    }
    for (n = 0; n < 56; ++n) {
        int nid = block->slot[n].id;
        const prj_block *neighbor;
        int start[3];
        int end[3];
        int i;
        int j;
        int k;

        if (nid < 0 || nid >= mesh->nblocks) {
            continue;
        }
        neighbor = &mesh->blocks[nid];
        if (!prj_mhd_local_block(neighbor) || neighbor->level != block->level) {
            continue;
        }
        prj_mhd_edge_slot_bounds(block, neighbor, &block->slot[n], 1, PRJ_BLOCK_SIZE, start, end);
        for (i = start[0]; i <= end[0]; ++i) {
            for (j = start[1]; j <= end[1]; ++j) {
                for (k = start[2]; k <= end[2]; ++k) {
                    double x[3];
                    int ii;
                    int jj;
                    int kk;
                    double diff;

                    if (prj_mhd_edge_is_interior(dir, i, j, k)) {
                        continue;
                    }
                    if (!prj_mhd_edge_is_ct_active(dir, i, j, k)) {
                        continue;
                    }
                    prj_mhd_edge_position(block, dir, i, j, k, x);
                    if (!prj_mhd_block_owns_edge_position(neighbor, dir, x, &ii, &jj, &kk)) {
                        continue;
                    }
                    diff = fabs(block->emf[dir][IDX(i, j, k)] - neighbor->emf[dir][IDX(ii, jj, kk)]);
                    if (diff > 1.0e-12 * prj_riemann_max_double(1.0,
                            fabs(block->emf[dir][IDX(i, j, k)]))) {
                        fprintf(stderr,
                            "prj_mhd_debug_check_emf: same-level mismatch block=%d neighbor=%d "
                            "dir=%d idx=(%d,%d,%d) diff=%g\n",
                            block->id, neighbor->id, dir, i, j, k, diff);
                        abort();
                    }
                }
            }
        }
    }
}

static double prj_mhd_field_scale(double B_norm)
{
    return B_norm / sqrt(PRJ_MHD_FOUR_PI);
}

static void prj_mhd_vector_potential(const prj_sim *sim, double x1, double x2, double x3, double Avec[3])
{
    double Bcode;

    if (sim == 0 || Avec == 0) {
        prj_mhd_abort("prj_mhd_vector_potential: null input");
    }

    Bcode = prj_mhd_field_scale(sim->B_norm);
    if (sim->mhd_init_type == PRJ_MHD_INIT_UNIFORM) {
        Avec[0] = -0.5 * x2 * Bcode;
        Avec[1] = 0.5 * x1 * Bcode;
        Avec[2] = 0.0;
        return;
    }
    if (sim->mhd_init_type == PRJ_MHD_INIT_DIPOLE) {
        double r = sqrt(x1 * x1 + x2 * x2 + x3 * x3);
        double factor;

        if (sim->B_scale <= 0.0) {
            prj_mhd_abort("prj_mhd_vector_potential: B_scale must be positive for dipole initialization");
        }
        factor = Bcode / (1.0 + pow(r / sim->B_scale, 3.0));
        Avec[0] = -0.5 * x2 * factor;
        Avec[1] = 0.5 * x1 * factor;
        Avec[2] = 0.0;
        return;
    }

    prj_mhd_abort("prj_mhd_vector_potential: unknown mhd_init_type");
}

static void prj_mhd_edge_coordinates(const prj_block *block, int dir, int i, int j, int k, double x[3])
{
    if (block == 0 || x == 0) {
        prj_mhd_abort("prj_mhd_edge_coordinates: null input");
    }

    if (dir == X1DIR) {
        x[0] = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
        x[1] = block->xmin[1] + (double)j * block->dx[1];
        x[2] = block->xmin[2] + (double)k * block->dx[2];
    } else if (dir == X2DIR) {
        x[0] = block->xmin[0] + (double)i * block->dx[0];
        x[1] = block->xmin[1] + ((double)j + 0.5) * block->dx[1];
        x[2] = block->xmin[2] + (double)k * block->dx[2];
    } else if (dir == X3DIR) {
        x[0] = block->xmin[0] + (double)i * block->dx[0];
        x[1] = block->xmin[1] + (double)j * block->dx[1];
        x[2] = block->xmin[2] + ((double)k + 0.5) * block->dx[2];
    } else {
        prj_mhd_abort("prj_mhd_edge_coordinates: invalid direction");
    }
}

static void prj_mhd_fill_vector_potential(prj_sim *sim, prj_block *block)
{
    int dir;
    int i;
    int j;
    int k;

    if (sim == 0 || block == 0) {
        prj_mhd_abort("prj_mhd_fill_vector_potential: null input");
    }
    if (block->emf[0] == 0 || block->emf[1] == 0 || block->emf[2] == 0) {
        prj_mhd_abort("prj_mhd_fill_vector_potential: emf storage is not allocated");
    }

    for (dir = 0; dir < 3; ++dir) {
        for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
            for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
                for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                    double x[3];
                    double Avec[3];

                    prj_mhd_edge_coordinates(block, dir, i, j, k, x);
                    prj_mhd_vector_potential(sim, x[0], x[1], x[2], Avec);
                    block->emf[dir][IDX(i, j, k)] = Avec[dir];
                }
            }
        }
    }
}

static void prj_mhd_curl_vector_potential(prj_block *block)
{
    int i;
    int j;
    int k;

    if (block == 0) {
        prj_mhd_abort("prj_mhd_curl_vector_potential: null block");
    }
    if (block->Bf[0] == 0 || block->Bf[1] == 0 || block->Bf[2] == 0) {
        prj_mhd_abort("prj_mhd_curl_vector_potential: Bf storage is not allocated");
    }

    prj_fill(block->Bf[0], PRJ_BLOCK_NCELLS, 0.0);
    prj_fill(block->Bf[1], PRJ_BLOCK_NCELLS, 0.0);
    prj_fill(block->Bf[2], PRJ_BLOCK_NCELLS, 0.0);

    for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
        for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST - 1; ++j) {
            for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST - 1; ++k) {
                block->Bf[X1DIR][IDX(i, j, k)] =
                    (block->emf[X3DIR][IDX(i, j + 1, k)] - block->emf[X3DIR][IDX(i, j, k)]) / block->dx[1] -
                    (block->emf[X2DIR][IDX(i, j, k + 1)] - block->emf[X2DIR][IDX(i, j, k)]) / block->dx[2];
            }
        }
    }

    for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST - 1; ++i) {
        for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
            for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST - 1; ++k) {
                block->Bf[X2DIR][IDX(i, j, k)] =
                    (block->emf[X1DIR][IDX(i, j, k + 1)] - block->emf[X1DIR][IDX(i, j, k)]) / block->dx[2] -
                    (block->emf[X3DIR][IDX(i + 1, j, k)] - block->emf[X3DIR][IDX(i, j, k)]) / block->dx[0];
            }
        }
    }

    for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST - 1; ++i) {
        for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST - 1; ++j) {
            for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                block->Bf[X3DIR][IDX(i, j, k)] =
                    (block->emf[X2DIR][IDX(i + 1, j, k)] - block->emf[X2DIR][IDX(i, j, k)]) / block->dx[0] -
                    (block->emf[X1DIR][IDX(i, j + 1, k)] - block->emf[X1DIR][IDX(i, j, k)]) / block->dx[1];
            }
        }
    }
}

static void prj_mhd_clear_emf(prj_block *block)
{
    int dir;

    if (block == 0) {
        prj_mhd_abort("prj_mhd_clear_emf: null block");
    }

    for (dir = 0; dir < 3; ++dir) {
        if (block->emf[dir] == 0) {
            prj_mhd_abort("prj_mhd_clear_emf: emf storage is not allocated");
        }
        prj_fill(block->emf[dir], PRJ_BLOCK_NCELLS, 0.0);
    }
}

static void prj_mhd_copy_face_field(double *dst[3], double *src[3])
{
    int dir;

    for (dir = 0; dir < 3; ++dir) {
        if (dst[dir] == 0 || src[dir] == 0) {
            prj_mhd_abort("prj_mhd_copy_face_field: face-field storage is not allocated");
        }
        prj_fill(dst[dir], PRJ_BLOCK_NCELLS, 0.0);
        for (int idx = 0; idx < PRJ_BLOCK_NCELLS; ++idx) {
            dst[dir][idx] = src[dir][idx];
        }
    }
}

static void prj_mhd_finalize_initialized_blocks(prj_sim *sim)
{
    int bidx;

    if (sim == 0) {
        prj_mhd_abort("prj_mhd_finalize_initialized_blocks: sim is null");
    }

    for (bidx = 0; bidx < sim->mesh.nblocks; ++bidx) {
        prj_block *block = &sim->mesh.blocks[bidx];

        if (!prj_mhd_local_block(block)) {
            continue;
        }
        prj_mhd_curl_vector_potential(block);
        prj_mhd_copy_face_field(block->Bf1, block->Bf);
        prj_mhd_bf2bc(block, 1);
        prj_mhd_bf2bc(block, 2);
        prj_mhd_clear_emf(block);
    }
}

static void prj_mhd_init_impl(prj_sim *sim, int sync_emf_before_curl)
{
    int bidx;

    if (sim == 0) {
        prj_mhd_abort("prj_mhd_init_impl: sim is null");
    }

    for (bidx = 0; bidx < sim->mesh.nblocks; ++bidx) {
        prj_block *block = &sim->mesh.blocks[bidx];

        if (!prj_mhd_local_block(block)) {
            continue;
        }
        prj_mhd_fill_vector_potential(sim, block);
    }

    if (sync_emf_before_curl != 0) {
        prj_mpi *mpi = prj_mpi_current();

        prj_mhd_emf_send(&sim->mesh);
#if defined(PRJ_ENABLE_MPI)
        if (mpi != 0 && mpi->totrank > 1 && mpi->neighbor_number > 0) {
            prj_mpi_exchange_emf(&sim->mesh, mpi);
        }
#else
        (void)mpi;
#endif
#if PRJ_MHD_DEBUG
        prj_mhd_debug_check_emf(&sim->mesh);
#endif
    }

    prj_mhd_finalize_initialized_blocks(sim);
}

static int prj_mhd_sign_to_index(int sign)
{
    return sign > 0 ? 1 : 0;
}

static int prj_mhd_index_to_sign(int idx)
{
    return idx == 0 ? -1 : 1;
}


static double prj_mhd_interpolate_outer_face(const prj_block *coarse, const double *Bface[3],
    int dir, int side_idx, int i, int j, int k, int a_sign, int b_sign)
{
    double base;
    double slope_a;
    double slope_b;
    double stencil_a[3];
    double stencil_b[3];
    double off_a;
    double off_b;

    if (coarse == 0 || Bface == 0 || Bface[0] == 0 || Bface[1] == 0 || Bface[2] == 0) {
        prj_mhd_abort("prj_mhd_interpolate_outer_face: invalid input");
    }
    if (side_idx != 0 && side_idx != 2) {
        prj_mhd_abort("prj_mhd_interpolate_outer_face: invalid side index");
    }

    if (dir == X1DIR) {
        int iface = i + (side_idx == 2 ? 1 : 0);

        base = Bface[X1DIR][IDX(iface, j, k)];
        stencil_a[0] = Bface[X1DIR][IDX(iface, j - 1, k)];
        stencil_a[1] = base;
        stencil_a[2] = Bface[X1DIR][IDX(iface, j + 1, k)];
        stencil_b[0] = Bface[X1DIR][IDX(iface, j, k - 1)];
        stencil_b[1] = base;
        stencil_b[2] = Bface[X1DIR][IDX(iface, j, k + 1)];
        slope_a = prj_reconstruct_slope(stencil_a, coarse->dx[1]);
        slope_b = prj_reconstruct_slope(stencil_b, coarse->dx[2]);
        off_a = 0.25 * (double)a_sign * coarse->dx[1];
        off_b = 0.25 * (double)b_sign * coarse->dx[2];
    } else if (dir == X2DIR) {
        int jface = j + (side_idx == 2 ? 1 : 0);

        base = Bface[X2DIR][IDX(i, jface, k)];
        stencil_a[0] = Bface[X2DIR][IDX(i - 1, jface, k)];
        stencil_a[1] = base;
        stencil_a[2] = Bface[X2DIR][IDX(i + 1, jface, k)];
        stencil_b[0] = Bface[X2DIR][IDX(i, jface, k - 1)];
        stencil_b[1] = base;
        stencil_b[2] = Bface[X2DIR][IDX(i, jface, k + 1)];
        slope_a = prj_reconstruct_slope(stencil_a, coarse->dx[0]);
        slope_b = prj_reconstruct_slope(stencil_b, coarse->dx[2]);
        off_a = 0.25 * (double)a_sign * coarse->dx[0];
        off_b = 0.25 * (double)b_sign * coarse->dx[2];
    } else if (dir == X3DIR) {
        int kface = k + (side_idx == 2 ? 1 : 0);

        base = Bface[X3DIR][IDX(i, j, kface)];
        stencil_a[0] = Bface[X3DIR][IDX(i - 1, j, kface)];
        stencil_a[1] = base;
        stencil_a[2] = Bface[X3DIR][IDX(i + 1, j, kface)];
        stencil_b[0] = Bface[X3DIR][IDX(i, j - 1, kface)];
        stencil_b[1] = base;
        stencil_b[2] = Bface[X3DIR][IDX(i, j + 1, kface)];
        slope_a = prj_reconstruct_slope(stencil_a, coarse->dx[0]);
        slope_b = prj_reconstruct_slope(stencil_b, coarse->dx[1]);
        off_a = 0.25 * (double)a_sign * coarse->dx[0];
        off_b = 0.25 * (double)b_sign * coarse->dx[1];
    } else {
        prj_mhd_abort("prj_mhd_interpolate_outer_face: invalid direction");
        return 0.0;
    }

    return base + slope_a * off_a + slope_b * off_b;
}

static void prj_mhd_fill_outer_patch(const prj_block *coarse, const double *Bface[3],
    int i, int j, int k, prj_mhd_face_patch *patch)
{
    int dir;
    int side_slot;
    int aidx;
    int bidx;

    if (patch == 0) {
        prj_mhd_abort("prj_mhd_fill_outer_patch: patch is null");
    }

    for (dir = 0; dir < 3; ++dir) {
        for (side_slot = 0; side_slot < 2; ++side_slot) {
            int side_idx = side_slot == 0 ? 0 : 2;

            for (aidx = 0; aidx < 2; ++aidx) {
                for (bidx = 0; bidx < 2; ++bidx) {
                    if (patch->fixed[dir][side_idx][aidx][bidx] == 0) {
                        patch->value[dir][side_idx][aidx][bidx] =
                            prj_mhd_interpolate_outer_face(coarse, Bface, dir, side_idx, i, j, k,
                                prj_mhd_index_to_sign(aidx), prj_mhd_index_to_sign(bidx));
                    }
                }
            }
        }
    }
}



/* Directly compute the 12 inner face B values from the outer patch values using the
 * analytic formulas (8)-(12) from Balsara (2001).  The outer faces in patch->value
 * (nidx 0 and 2) must already be filled before calling this function.
 *
 * The notation follows the paper: u,v,w are Bx,By,Bz face-centered fields; indices
 * ±1 label the four transverse sub-faces, and ±2 labels the two outer (coarse-side)
 * faces in the face-normal direction.  The inner faces (at the midplane of the coarse
 * cell, index 0 in the paper) are written to patch->value[dir][1][...].
 */
static void prj_mhd_compute_inner_faces(prj_mhd_face_patch *patch, const double dx[3])
{
    double dx2;
    double dy2;
    double dz2;
    double U_xx;
    double V_yy;
    double W_zz;
    double U_xyz;
    double V_xyz;
    double W_xyz;
    int i;
    int j;
    int k;

    if (patch == 0 || dx == 0) {
        prj_mhd_abort("prj_mhd_compute_inner_faces: null input");
    }

    dx2 = dx[0] * dx[0];
    dy2 = dx[1] * dx[1];
    dz2 = dx[2] * dx[2];
    U_xx  = 0.0;
    V_yy  = 0.0;
    W_zz  = 0.0;
    U_xyz = 0.0;
    V_xyz = 0.0;
    W_xyz = 0.0;

    for (i = -1; i <= 1; i += 2) {
        for (j = -1; j <= 1; j += 2) {
            for (k = -1; k <= 1; k += 2) {
                int ii = prj_mhd_sign_to_index(i);
                int jj = prj_mhd_sign_to_index(j);
                int kk = prj_mhd_sign_to_index(k);
                /* u^{2i,j,k}: outer x-face on the i-side */
                double u_o = patch->value[X1DIR][i > 0 ? 2 : 0][jj][kk];
                /* v^{i,2j,k}: outer y-face on the j-side */
                double v_o = patch->value[X2DIR][j > 0 ? 2 : 0][ii][kk];
                /* w^{i,j,2k}: outer z-face on the k-side */
                double w_o = patch->value[X3DIR][k > 0 ? 2 : 0][ii][jj];

                /* Eq. (11) and cyclic permutations */
                U_xx += (double)(i * j) * v_o + (double)(i * k) * w_o;
                V_yy += (double)(j * k) * w_o + (double)(j * i) * u_o;
                W_zz += (double)(k * i) * u_o + (double)(k * j) * v_o;

                /* Eq. (12) and cyclic permutations */
                U_xyz += (double)(i * j * k) * u_o;
                V_xyz += (double)(i * j * k) * v_o;
                W_xyz += (double)(i * j * k) * w_o;
            }
        }
    }
    U_xx  *= 0.125;
    V_yy  *= 0.125;
    W_zz  *= 0.125;
    U_xyz /= 8.0 * (dy2 + dz2);
    V_xyz /= 8.0 * (dx2 + dz2);
    W_xyz /= 8.0 * (dx2 + dy2);

    /* Eq. (8): u^{0,j,k} for j,k = +-1 */
    for (j = -1; j <= 1; j += 2) {
        for (k = -1; k <= 1; k += 2) {
            int jj = prj_mhd_sign_to_index(j);
            int kk = prj_mhd_sign_to_index(k);

            patch->value[X1DIR][1][jj][kk] =
                0.5 * (patch->value[X1DIR][2][jj][kk] + patch->value[X1DIR][0][jj][kk])
                + U_xx + (double)k * dz2 * V_xyz + (double)j * dy2 * W_xyz;
        }
    }
    /* Eq. (9): v^{i,0,k} for i,k = +-1 */
    for (i = -1; i <= 1; i += 2) {
        for (k = -1; k <= 1; k += 2) {
            int ii = prj_mhd_sign_to_index(i);
            int kk = prj_mhd_sign_to_index(k);

            patch->value[X2DIR][1][ii][kk] =
                0.5 * (patch->value[X2DIR][2][ii][kk] + patch->value[X2DIR][0][ii][kk])
                + V_yy + (double)i * dx2 * W_xyz + (double)k * dz2 * U_xyz;
        }
    }
    /* Eq. (10): w^{i,j,0} for i,j = +-1 */
    for (i = -1; i <= 1; i += 2) {
        for (j = -1; j <= 1; j += 2) {
            int ii = prj_mhd_sign_to_index(i);
            int jj = prj_mhd_sign_to_index(j);

            patch->value[X3DIR][1][ii][jj] =
                0.5 * (patch->value[X3DIR][2][ii][jj] + patch->value[X3DIR][0][ii][jj])
                + W_zz + (double)j * dy2 * U_xyz + (double)i * dx2 * V_xyz;
        }
    }
}

void prj_mhd_bf_prolong_direct(const prj_block *coarse, const double *Bface[3],
    int i, int j, int k, prj_mhd_face_patch *patch)
{
    if (coarse == 0 || Bface == 0 || patch == 0) {
        prj_mhd_abort("prj_mhd_bf_prolong_direct: null input");
    }
    if (Bface[0] == 0 || Bface[1] == 0 || Bface[2] == 0) {
        prj_mhd_abort("prj_mhd_bf_prolong_direct: face-centered field is not allocated");
    }
    if (i < -PRJ_NGHOST + 1 || i > PRJ_BLOCK_SIZE + PRJ_NGHOST - 2 ||
        j < -PRJ_NGHOST + 1 || j > PRJ_BLOCK_SIZE + PRJ_NGHOST - 2 ||
        k < -PRJ_NGHOST + 1 || k > PRJ_BLOCK_SIZE + PRJ_NGHOST - 2) {
        prj_mhd_abort("prj_mhd_bf_prolong_direct: coarse-cell index out of range");
    }
    prj_mhd_fill_outer_patch(coarse, Bface, i, j, k, patch);
    prj_mhd_compute_inner_faces(patch, coarse->dx);
}

void prj_mhd_apply_bf_patch(prj_mesh *mesh, prj_block *fine, int stage,
    const prj_block *coarse, const double *Bface[3], int i, int j, int k)
{
    prj_mhd_face_patch patch;
    double *fine_face[3];
    int dir;
    int nidx;
    int aidx;
    int bidx;
    int has_target_face;

    if (mesh == 0 || fine == 0 || coarse == 0 || Bface == 0) {
        prj_mhd_abort("prj_mhd_apply_bf_patch: null input");
    }
    fine_face[0] = prj_mhd_bface_stage(fine, stage, X1DIR);
    fine_face[1] = prj_mhd_bface_stage(fine, stage, X2DIR);
    fine_face[2] = prj_mhd_bface_stage(fine, stage, X3DIR);
    if (fine_face[0] == 0 || fine_face[1] == 0 || fine_face[2] == 0) {
        prj_mhd_abort("prj_mhd_apply_bf_patch: fine face storage is not allocated");
    }

    prj_mhd_face_patch_clear(&patch);
    has_target_face = 0;
    for (dir = 0; dir < 3; ++dir) {
        for (nidx = 0; nidx < 3; nidx += 2) {
            for (aidx = 0; aidx < 2; ++aidx) {
                for (bidx = 0; bidx < 2; ++bidx) {
                    double x[3];
                    int ii;
                    int jj;
                    int kk;

                    prj_mhd_patch_face_position(coarse, dir, i, j, k, nidx, aidx, bidx, x);
                    if (!prj_mhd_face_index_from_position(fine, dir, x, &ii, &jj, &kk)) {
                        continue;
                    }
                    if (!prj_mhd_storage_index_valid(ii, jj, kk)) {
                        continue;
                    }
                    if (!prj_mhd_face_is_interior(dir, ii, jj, kk)) {
                        has_target_face = 1;
                    }
                    if (prj_mhd_face_has_fine_owner(mesh, fine, coarse->id, dir, x)) {
                        patch.fixed[dir][nidx][aidx][bidx] = 1;
                        patch.value[dir][nidx][aidx][bidx] = fine_face[dir][IDX(ii, jj, kk)];
                    }
                }
            }
        }
    }
    if (has_target_face == 0) {
        return;
    }

    prj_mhd_bf_prolong_direct(coarse, Bface, i, j, k, &patch);
    for (dir = 0; dir < 3; ++dir) {
        for (nidx = 0; nidx < 3; ++nidx) {
            for (aidx = 0; aidx < 2; ++aidx) {
                for (bidx = 0; bidx < 2; ++bidx) {
                    double x[3];
                    int ii;
                    int jj;
                    int kk;

                    prj_mhd_patch_face_position(coarse, dir, i, j, k, nidx, aidx, bidx, x);
                    if (!prj_mhd_face_index_from_position(fine, dir, x, &ii, &jj, &kk)) {
                        continue;
                    }
                    if (!prj_mhd_storage_index_valid(ii, jj, kk) ||
                        prj_mhd_face_is_interior(dir, ii, jj, kk)) {
                        continue;
                    }
                    if ((nidx == 0 || nidx == 2) && patch.fixed[dir][nidx][aidx][bidx] != 0) {
                        continue;
                    }
                    fine_face[dir][IDX(ii, jj, kk)] = patch.value[dir][nidx][aidx][bidx];
#if PRJ_MHD_DEBUG
                    if (dir == X2DIR &&
                        ((fine->id == 1 && ii == 16 && jj == 17) ||
                            (fine->id == 5 && ii == 16 && jj == 1))) {
                        fprintf(stderr,
                            "[mhd-bf-fill-debug] kind=RECON dst_block=%d coarse_block=%d "
                            "coarse_cell=(%d,%d,%d) face=(%d,%d,%d) value=% .17g fixed=%d\n",
                            fine->id, coarse->id, i, j, k, ii, jj, kk,
                            fine_face[dir][IDX(ii, jj, kk)], patch.fixed[dir][nidx][aidx][bidx]);
                    }
                    if (dir == X3DIR &&
                        ((fine->id == 1 && ii == 16 && jj == 16) ||
                            (fine->id == 5 && ii == 16 && jj == 0))) {
                        fprintf(stderr,
                            "[mhd-bf-fill-debug] kind=RECON dst_block=%d coarse_block=%d "
                            "coarse_cell=(%d,%d,%d) face=(%d,%d,%d) value=% .17g fixed=%d\n",
                            fine->id, coarse->id, i, j, k, ii, jj, kk,
                            fine_face[dir][IDX(ii, jj, kk)], patch.fixed[dir][nidx][aidx][bidx]);
                    }
#endif
                }
            }
        }
    }
}
#endif

#if !PRJ_MHD
void prj_mhd_bf_prolong_direct(const prj_block *coarse, const double *Bface[3],
    int i, int j, int k, prj_mhd_face_patch *patch)
{
    (void)coarse;
    (void)Bface;
    (void)i;
    (void)j;
    (void)k;
    (void)patch;
    prj_mhd_abort("prj_mhd_bf_prolong_direct: called with PRJ_MHD disabled");
}
#endif

void prj_mhd_init(prj_sim *sim)
{
#if PRJ_MHD
    prj_mhd_init_impl(sim, 1);
#else
    (void)sim;
#endif
}

void prj_mhd_bf2bc(prj_block *block, int stage)
{
#if PRJ_MHD
    double *Wstage;
    double *Bface[3];
    int i;
    int j;
    int k;

    if (block == 0) {
        prj_mhd_abort("prj_mhd_bf2bc: block is null");
    }
    if (block->U == 0 || block->W == 0 || block->W1 == 0) {
        prj_mhd_abort("prj_mhd_bf2bc: cell data is not allocated");
    }

    if (stage == 1) {
        Wstage = block->W;
        Bface[0] = block->Bf[0];
        Bface[1] = block->Bf[1];
        Bface[2] = block->Bf[2];
    } else if (stage == 2) {
        Wstage = block->W1;
        Bface[0] = block->Bf1[0];
        Bface[1] = block->Bf1[1];
        Bface[2] = block->Bf1[2];
    } else {
        prj_mhd_abort("prj_mhd_bf2bc: invalid RK stage");
        return;
    }

    if (Bface[0] == 0 || Bface[1] == 0 || Bface[2] == 0) {
        prj_mhd_abort("prj_mhd_bf2bc: face-centered magnetic field is not allocated");
    }

    for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
        for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
            for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                double old_b1 = block->U[VIDX(PRJ_CONS_B1, i, j, k)];
                double old_b2 = block->U[VIDX(PRJ_CONS_B2, i, j, k)];
                double old_b3 = block->U[VIDX(PRJ_CONS_B3, i, j, k)];
                double new_b1 = old_b1;
                double new_b2 = old_b2;
                double new_b3 = old_b3;
                double old_mag = 0.5 * (old_b1 * old_b1 + old_b2 * old_b2 + old_b3 * old_b3);
                double new_mag;

                if (i < PRJ_BLOCK_SIZE + PRJ_NGHOST - 1) {
                    new_b1 = 0.5 * (Bface[X1DIR][IDX(i, j, k)] + Bface[X1DIR][IDX(i + 1, j, k)]);
                }
                if (j < PRJ_BLOCK_SIZE + PRJ_NGHOST - 1) {
                    new_b2 = 0.5 * (Bface[X2DIR][IDX(i, j, k)] + Bface[X2DIR][IDX(i, j + 1, k)]);
                }
                if (k < PRJ_BLOCK_SIZE + PRJ_NGHOST - 1) {
                    new_b3 = 0.5 * (Bface[X3DIR][IDX(i, j, k)] + Bface[X3DIR][IDX(i, j, k + 1)]);
                }
                new_mag = 0.5 * (new_b1 * new_b1 + new_b2 * new_b2 + new_b3 * new_b3);

                block->U[VIDX(PRJ_CONS_ETOT, i, j, k)] += new_mag - old_mag;
                block->U[VIDX(PRJ_CONS_B1, i, j, k)] = new_b1;
                block->U[VIDX(PRJ_CONS_B2, i, j, k)] = new_b2;
                block->U[VIDX(PRJ_CONS_B3, i, j, k)] = new_b3;
                Wstage[VIDX(PRJ_PRIM_B1, i, j, k)] = new_b1;
                Wstage[VIDX(PRJ_PRIM_B2, i, j, k)] = new_b2;
                Wstage[VIDX(PRJ_PRIM_B3, i, j, k)] = new_b3;
            }
        }
    }
#else
    (void)block;
    (void)stage;
#endif
}

double prj_mhd_emf_upwind(const double emf_face[4], const double emf_cell[4], const double face_velocity[4])
{
    double emf;

    if (emf_face == 0 || emf_cell == 0 || face_velocity == 0) {
        prj_mhd_abort("prj_mhd_emf_upwind: null input");
    }

    emf = 0.25 * (emf_face[0] + emf_face[1] + emf_face[2] + emf_face[3]);
    if (face_velocity[0] > 0.0) {
        emf += 0.25 * (emf_face[3] - emf_cell[3]);
    } else if (face_velocity[0] < 0.0) {
        emf += 0.25 * (emf_face[1] - emf_cell[0]);
    } else {
        emf += 0.125 * (emf_face[3] - emf_cell[3] + emf_face[1] - emf_cell[0]);
    }
    if (face_velocity[1] > 0.0) {
        emf += 0.25 * (emf_face[2] - emf_cell[1]);
    } else if (face_velocity[1] < 0.0) {
        emf += 0.25 * (emf_face[0] - emf_cell[0]);
    } else {
        emf += 0.125 * (emf_face[2] - emf_cell[1] + emf_face[0] - emf_cell[0]);
    }
    if (face_velocity[2] > 0.0) {
        emf += 0.25 * (emf_face[3] - emf_cell[2]);
    } else if (face_velocity[2] < 0.0) {
        emf += 0.25 * (emf_face[1] - emf_cell[1]);
    } else {
        emf += 0.125 * (emf_face[3] - emf_cell[2] + emf_face[1] - emf_cell[1]);
    }
    if (face_velocity[3] > 0.0) {
        emf += 0.25 * (emf_face[2] - emf_cell[2]);
    } else if (face_velocity[3] < 0.0) {
        emf += 0.25 * (emf_face[0] - emf_cell[3]);
    } else {
        emf += 0.125 * (emf_face[2] - emf_cell[2] + emf_face[0] - emf_cell[3]);
    }
    return emf;
}

void prj_mhd_emf_send(prj_mesh *mesh)
{
#if PRJ_MHD
    int bidx;

    if (mesh == 0) {
        prj_mhd_abort("prj_mhd_emf_send: mesh is null");
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        double *sum_buf[3] = {0, 0, 0};
        int *count_buf[3] = {0, 0, 0};
        int have_fine_neighbor = 0;
        int dir;

        if (!prj_mhd_local_block(block)) {
            continue;
        }
        for (dir = 0; dir < 3; ++dir) {
            if (block->emf[dir] == 0) {
                prj_mhd_abort("prj_mhd_emf_send: emf storage is not allocated");
            }
        }
        {
            int n;

            for (n = 0; n < 56; ++n) {
                int nid = block->slot[n].id;
                const prj_block *neighbor;

                if (nid < 0 || nid >= mesh->nblocks) {
                    continue;
                }
                neighbor = &mesh->blocks[nid];
                if (block->slot[n].type == PRJ_NEIGHBOR_CORNER ||
                    !prj_mhd_local_block(neighbor) || neighbor->level != block->level + 1) {
                    continue;
                }
                if (!have_fine_neighbor) {
                    int axis;

                    for (axis = 0; axis < 3; ++axis) {
                        sum_buf[axis] = (double *)calloc((size_t)PRJ_BLOCK_NCELLS, sizeof(*sum_buf[axis]));
                        count_buf[axis] = (int *)calloc((size_t)PRJ_BLOCK_NCELLS, sizeof(*count_buf[axis]));
                        if (sum_buf[axis] == 0 || count_buf[axis] == 0) {
                            prj_mhd_abort("prj_mhd_emf_send: failed to allocate edge accumulation buffers");
                        }
                    }
                    have_fine_neighbor = 1;
                }
                for (dir = 0; dir < 3; ++dir) {
                    int start[3];
                    int end[3];
                    int i;
                    int j;
                    int k;

                    if (!prj_mhd_coarse_edge_bounds(block, neighbor, &block->slot[n], dir, 1, start, end)) {
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

                                    prj_mhd_coarse_subedge_position(block, dir, i, j, k, subidx, x);
                                    if (!prj_mhd_block_owns_edge_position(neighbor, dir, x, &ii, &jj, &kk)) {
                                        continue;
                                    }
                                    sum_buf[dir][IDX(i, j, k)] += neighbor->emf[dir][IDX(ii, jj, kk)];
                                    count_buf[dir][IDX(i, j, k)] += 1;
                                }
                            }
                        }
                    }
                }
            }
        }
        if (have_fine_neighbor) {
            for (dir = 0; dir < 3; ++dir) {
                int i;
                int j;
                int k;

                for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
                    for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
                        for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                            if (prj_mhd_edge_is_interior(dir, i, j, k)) {
                                continue;
                            }
                            if (count_buf[dir][IDX(i, j, k)] == 2) {
                                block->emf[dir][IDX(i, j, k)] = 0.5 * sum_buf[dir][IDX(i, j, k)];
                            }
                        }
                    }
                }
                free(count_buf[dir]);
                free(sum_buf[dir]);
            }
        }
    }
#else
    (void)mesh;
#endif
}

void prj_mhd_debug_check_emf(const prj_mesh *mesh)
{
#if PRJ_MHD
    int bidx;

    if (mesh == 0) {
        prj_mhd_abort("prj_mhd_debug_check_emf: mesh is null");
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        const prj_block *block = &mesh->blocks[bidx];
        int dir;

        if (!prj_mhd_local_block(block)) {
            continue;
        }
        for (dir = 0; dir < 3; ++dir) {
            int i;
            int j;
            int k;

            prj_mhd_debug_check_same_level_edges(mesh, block, dir);
            for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
                for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
                    for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                        double local_sum = 0.0;
                        int local_count = 0;

                        if (prj_mhd_edge_is_interior(dir, i, j, k)) {
                            continue;
                        }
                        if (!prj_mhd_edge_is_ct_active(dir, i, j, k)) {
                            continue;
                        }
                        prj_mhd_local_edge_partial(mesh, block, dir, i, j, k, &local_sum, &local_count);
                        if (local_count == 2) {
                            double expected = 0.5 * local_sum;
                            double diff = fabs(block->emf[dir][IDX(i, j, k)] - expected);

                            if (diff > 1.0e-12 * prj_riemann_max_double(1.0, fabs(expected))) {
                                fprintf(stderr,
                                    "prj_mhd_debug_check_emf: coarse-fine mismatch block=%d dir=%d "
                                    "idx=(%d,%d,%d) diff=%g\n",
                                    block->id, dir, i, j, k, diff);
                                abort();
                            }
                        }
                    }
                }
            }
        }
    }
#if defined(PRJ_ENABLE_MPI)
    {
        prj_mpi *mpi = prj_mpi_current();

        if (mpi != 0 && mpi->totrank > 1) {
            int nb;
            MPI_Request *requests = 0;
            int request_count = 0;
            int request_capacity = 0;

            for (nb = 0; nb < mpi->neighbor_number; ++nb) {
                prj_mpi_buffer *buffer = &mpi->neighbor_buffer[nb];
                int *send_sizes = (int *)calloc((size_t)buffer->number, sizeof(*send_sizes));
                int *idx_send[3] = {0, 0, 0};
                double *value_send = 0;
                int idx_count = 0;
                int idx_capacity = 0;
                int value_count = 0;
                int value_capacity = 0;
                int recv_entries;
                int recv_total;
                int axis;
                int occ = 0;
                int block_id;

                for (block_id = 0; block_id < mesh->nblocks; ++block_id) {
                    const prj_block *block = &mesh->blocks[block_id];
                    int n;

                    if (!prj_mhd_local_block(block)) {
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
                        if (block->level == neighbor->level) {
                            int dir;

                            for (dir = 0; dir < 3; ++dir) {
                                int start[3];
                                int end[3];
                                int i;
                                int j;
                                int k;

                                prj_mhd_edge_slot_bounds(block, neighbor, &block->slot[n], 0, PRJ_BLOCK_SIZE, start, end);
                                for (i = start[0]; i <= end[0]; ++i) {
                                    for (j = start[1]; j <= end[1]; ++j) {
                                        for (k = start[2]; k <= end[2]; ++k) {
                                            double x[3];
                                            int ii;
                                            int jj;
                                            int kk;

                                            if (prj_mhd_edge_is_interior(dir, i, j, k)) {
                                                continue;
                                            }
                                            if (!prj_mhd_edge_is_ct_active(dir, i, j, k)) {
                                                continue;
                                            }
                                            prj_mhd_edge_position(neighbor, dir, i, j, k, x);
                                            if (!prj_mhd_block_owns_edge_position(block, dir, x, &ii, &jj, &kk)) {
                                                continue;
                                            }
                                            if (prj_mhd_append_triplet(idx_send, &idx_count, &idx_capacity,
                                                    neighbor->id, block->id,
                                                    prj_mhd_encode_edge_code(dir, i, j, k)) != 0 ||
                                                prj_mhd_append_double(&value_send, &value_count, &value_capacity,
                                                    block->emf[dir][IDX(ii, jj, kk)]) != 0) {
                                                prj_mhd_abort("prj_mhd_debug_check_emf: failed to pack same-level MPI records");
                                            }
                                        }
                                    }
                                }
                            }
                        } else if (block->slot[n].type != PRJ_NEIGHBOR_CORNER &&
                                   block->level == neighbor->level + 1) {
                            int dir;

                            for (dir = 0; dir < 3; ++dir) {
                                int start[3];
                                int end[3];
                                int i;
                                int j;
                                int k;

                                if (!prj_mhd_coarse_edge_bounds(block, neighbor, &block->slot[n], dir, 0, start, end)) {
                                    continue;
                                }
                                for (i = start[0]; i <= end[0]; ++i) {
                                    for (j = start[1]; j <= end[1]; ++j) {
                                        for (k = start[2]; k <= end[2]; ++k) {
                                            int subidx;

                                            if (prj_mhd_edge_is_interior(dir, i, j, k)) {
                                                continue;
                                            }
                                            if (!prj_mhd_edge_is_ct_active(dir, i, j, k)) {
                                                continue;
                                            }
                                            for (subidx = 0; subidx < 2; ++subidx) {
                                                double x[3];
                                                int ii;
                                                int jj;
                                                int kk;

                                                prj_mhd_coarse_subedge_position(neighbor, dir, i, j, k, subidx, x);
                                                if (!prj_mhd_block_owns_edge_position(block, dir, x, &ii, &jj, &kk)) {
                                                    continue;
                                                }
                                                if (prj_mhd_append_triplet(idx_send, &idx_count, &idx_capacity,
                                                        neighbor->id, block->id,
                                                        prj_mhd_encode_edge_code(dir, i, j, k)) != 0 ||
                                                    prj_mhd_append_double(&value_send, &value_count, &value_capacity,
                                                        block->emf[dir][IDX(ii, jj, kk)]) != 0) {
                                                    prj_mhd_abort("prj_mhd_debug_check_emf: failed to pack coarse-fine MPI records");
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

                MPI_Sendrecv(&buffer->number, 1, MPI_INT, buffer->receiver_rank, 360,
                    &recv_entries, 1, MPI_INT, buffer->receiver_rank, 360,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                buffer->face_data_size_recv = (int *)calloc((size_t)recv_entries + 1U, sizeof(*buffer->face_data_size_recv));
                MPI_Sendrecv(buffer->face_data_size_send, buffer->number, MPI_INT, buffer->receiver_rank, 361,
                    buffer->face_data_size_recv, recv_entries, MPI_INT, buffer->receiver_rank, 361,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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
                        prj_mhd_abort("prj_mhd_debug_check_emf: failed to grow MPI request array");
                    }
                    requests = next;
                    request_capacity = new_capacity;
                }
                for (axis = 0; axis < 3; ++axis) {
                    MPI_Irecv(buffer->face_data_idx_recv[axis], recv_total, MPI_INT, buffer->receiver_rank, 362 + axis,
                        MPI_COMM_WORLD, &requests[request_count++]);
                    MPI_Isend(buffer->face_data_idx_send[axis], idx_count, MPI_INT, buffer->receiver_rank, 362 + axis,
                        MPI_COMM_WORLD, &requests[request_count++]);
                }
                MPI_Irecv(buffer->face_buffer_recv, recv_total, MPI_DOUBLE, buffer->receiver_rank, 365,
                    MPI_COMM_WORLD, &requests[request_count++]);
                MPI_Isend(buffer->face_buffer_send, value_count, MPI_DOUBLE, buffer->receiver_rank, 365,
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
                    const prj_mpi_buffer *buffer = &mpi->neighbor_buffer[nb];
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
                        int block_id2 = buffer->face_data_idx_recv[0][irec];
                        int source_id = buffer->face_data_idx_recv[1][irec];
                        int edge_code = buffer->face_data_idx_recv[2][irec];
                        double value = buffer->face_buffer_recv[pos++];
                        const prj_block *receiver;
                        const prj_block *source;

                        if (block_id2 < 0 || block_id2 >= mesh->nblocks ||
                            source_id < 0 || source_id >= mesh->nblocks) {
                            continue;
                        }
                        receiver = &mesh->blocks[block_id2];
                        source = &mesh->blocks[source_id];
                        if (receiver->rank != mpi->rank) {
                            continue;
                        }
                        if (source->level == receiver->level) {
                            int dir;
                            int ii;
                            int jj;
                            int kk;
                            double diff;

                            prj_mhd_decode_edge_code(edge_code, &dir, &ii, &jj, &kk);
                            if (!prj_mhd_edge_is_ct_active(dir, ii, jj, kk)) {
                                continue;
                            }
                            diff = fabs(receiver->emf[dir][IDX(ii, jj, kk)] - value);
                            if (diff > 1.0e-12 * prj_riemann_max_double(1.0, fabs(value))) {
                                fprintf(stderr,
                                    "prj_mhd_debug_check_emf: MPI same-level mismatch block=%d src=%d "
                                    "dir=%d idx=(%d,%d,%d) diff=%g\n",
                                    receiver->id, source->id, dir, ii, jj, kk, diff);
                                abort();
                            }
                        } else if (source->level == receiver->level + 1) {
                            int found = -1;
                            int a;

                            for (a = 0; a < acc_n; ++a) {
                                if (acc_block[a] == block_id2 && acc_edge[a] == edge_code) {
                                    found = a;
                                    break;
                                }
                            }
                            if (found < 0) {
                                int *next_block;
                                int *next_edge;
                                int *next_count;
                                double *next_sum;

                                if (acc_n >= acc_cap) {
                                    int new_cap = acc_cap == 0 ? 32 : 2 * acc_cap;

                                    next_block = (int *)realloc(acc_block, (size_t)new_cap * sizeof(*next_block));
                                    next_edge = (int *)realloc(acc_edge, (size_t)new_cap * sizeof(*next_edge));
                                    next_count = (int *)realloc(acc_count, (size_t)new_cap * sizeof(*next_count));
                                    next_sum = (double *)realloc(acc_sum, (size_t)new_cap * sizeof(*next_sum));
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
                                        prj_mhd_abort("prj_mhd_debug_check_emf: failed to grow coarse-fine accumulator");
                                    }
                                    acc_block = next_block;
                                    acc_edge = next_edge;
                                    acc_count = next_count;
                                    acc_sum = next_sum;
                                    acc_cap = new_cap;
                                }
                                found = acc_n++;
                                acc_block[found] = block_id2;
                                acc_edge[found] = edge_code;
                                acc_count[found] = 0;
                                acc_sum[found] = 0.0;
                            }
                            acc_count[found] += 1;
                            acc_sum[found] += value;
                        } else {
                            prj_mhd_abort("prj_mhd_debug_check_emf: invalid MPI edge relation");
                        }
                    }
                }
                for (nb = 0; nb < acc_n; ++nb) {
                    const prj_block *block = &mesh->blocks[acc_block[nb]];
                    int dir;
                    int ii;
                    int jj;
                    int kk;
                    double local_sum = 0.0;
                    int local_count = 0;
                    double diff;
                    double expected;

                    prj_mhd_decode_edge_code(acc_edge[nb], &dir, &ii, &jj, &kk);
                    if (!prj_mhd_edge_is_ct_active(dir, ii, jj, kk)) {
                        continue;
                    }
                    prj_mhd_local_edge_partial(mesh, block, dir, ii, jj, kk, &local_sum, &local_count);
                    if (acc_count[nb] + local_count != 2) {
                        fprintf(stderr,
                            "prj_mhd_debug_check_emf: MPI coarse-fine coverage mismatch block=%d edge=%d "
                            "(remote=%d local=%d)\n",
                            block->id, acc_edge[nb], acc_count[nb], local_count);
                        abort();
                    }
                    expected = 0.5 * (acc_sum[nb] + local_sum);
                    diff = fabs(block->emf[dir][IDX(ii, jj, kk)] - expected);
                    if (diff > 1.0e-12 * prj_riemann_max_double(1.0, fabs(expected))) {
                        fprintf(stderr,
                            "prj_mhd_debug_check_emf: MPI coarse-fine mismatch block=%d dir=%d "
                            "idx=(%d,%d,%d) diff=%g\n",
                            block->id, dir, ii, jj, kk, diff);
                        abort();
                    }
                }
                free(acc_sum);
                free(acc_count);
                free(acc_edge);
                free(acc_block);
            }
            free(requests);
        }
    }
#endif
#else
    (void)mesh;
#endif
}

void prj_mhd_debug_check_divergence(const prj_mesh *mesh, int stage)
{
#if PRJ_MHD
    int bidx;

    if (mesh == 0) {
        prj_mhd_abort("prj_mhd_debug_check_divergence: mesh is null");
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        const prj_block *block = &mesh->blocks[bidx];
        const double *Bface[3];
        int i;
        int j;
        int k;

        if (!prj_mhd_local_block(block)) {
            continue;
        }
        Bface[0] = prj_mhd_bface_stage_const(block, stage, X1DIR);
        Bface[1] = prj_mhd_bface_stage_const(block, stage, X2DIR);
        Bface[2] = prj_mhd_bface_stage_const(block, stage, X3DIR);
        for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
            for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    double divb =
                        (Bface[X1DIR][IDX(i + 1, j, k)] - Bface[X1DIR][IDX(i, j, k)]) / block->dx[0] +
                        (Bface[X2DIR][IDX(i, j + 1, k)] - Bface[X2DIR][IDX(i, j, k)]) / block->dx[1] +
                        (Bface[X3DIR][IDX(i, j, k + 1)] - Bface[X3DIR][IDX(i, j, k)]) / block->dx[2];
                    double scale =
                        fabs(Bface[X1DIR][IDX(i + 1, j, k)]) + fabs(Bface[X1DIR][IDX(i, j, k)]) +
                        fabs(Bface[X2DIR][IDX(i, j + 1, k)]) + fabs(Bface[X2DIR][IDX(i, j, k)]) +
                        fabs(Bface[X3DIR][IDX(i, j, k + 1)]) + fabs(Bface[X3DIR][IDX(i, j, k)]);

                    if (fabs(divb) > 1.0e-12 * prj_riemann_max_double(1.0, scale)) {
                        fprintf(stderr,
                            "prj_mhd_debug_check_divergence: block=%d stage=%d cell=(%d,%d,%d) divB=%g\n",
                            block->id, stage, i, j, k, divb);
                        abort();
                    }
                }
            }
        }
    }
#else
    (void)mesh;
    (void)stage;
#endif
}
