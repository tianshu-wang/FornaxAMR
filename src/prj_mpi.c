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

#if PRJ_MHD
static int prj_mpi_block_is_local(const prj_block *block)
{
    return prj_block_is_active(block) &&
        (prj_mpi_active == 0 || block->rank == prj_mpi_active->rank);
}
#endif

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

#if PRJ_MHD
static int prj_mpi_append_int(int **array, int *count, int *capacity, int value)
{
    int *next;
    int new_capacity;

    if (*count >= *capacity) {
        new_capacity = *capacity == 0 ? 32 : 2 * (*capacity);
        next = (int *)realloc(*array, (size_t)new_capacity * sizeof(**array));
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
#endif /* PRJ_MHD */

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
    size_t eosvar_count;
    size_t cons_count;

    prim_count = (size_t)PRJ_NVAR_PRIM * (size_t)PRJ_BLOCK_NCELLS;
    eosvar_count = (size_t)PRJ_NVAR_EOSVAR * (size_t)PRJ_BLOCK_NCELLS;
    cons_count = (size_t)PRJ_NVAR_CONS * (size_t)PRJ_BLOCK_NCELLS;
    return 2U * prim_count + eosvar_count + 5U * cons_count +
        9U * (size_t)PRJ_BLOCK_NCELLS
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

#if PRJ_MHD
static int prj_mpi_append_unique_int(int **array, int *count, int *capacity, int value)
{
    int i;

    for (i = 0; i < *count; ++i) {
        if ((*array)[i] == value) {
            return 0;
        }
    }
    return prj_mpi_append_int(array, count, capacity, value);
}

static int prj_mpi_collect_bf_source_blocks(const prj_mesh *mesh, const prj_mpi *mpi,
    int receiver_rank, int **ids, int *count)
{
    int capacity = 0;
    int bidx;

    *ids = 0;
    *count = 0;
    if (mesh == 0 || mpi == 0) {
        return 1;
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        const prj_block *block = &mesh->blocks[bidx];
        int n;

        if (!prj_mpi_block_is_local(block)) {
            continue;
        }
        for (n = 0; n < 56; ++n) {
            int nid = block->slot[n].id;

            if (nid >= 0 && nid < mesh->nblocks &&
                mesh->blocks[nid].rank == receiver_rank) {
                if (prj_mpi_append_unique_int(ids, count, &capacity, bidx) != 0) {
                    free(*ids);
                    *ids = 0;
                    *count = 0;
                    return 1;
                }
            }
        }
    }
    return 0;
}

static double *prj_mpi_bf_array(prj_block *block, int dir, int use_bf1)
{
    return use_bf1 != 0 ? block->Bf1[dir] : block->Bf[dir];
}

static const double *prj_mpi_bf_array_const(const prj_block *block, int dir, int use_bf1)
{
    return use_bf1 != 0 ? block->Bf1[dir] : block->Bf[dir];
}

static int prj_mpi_pack_bf_blocks(const prj_mesh *mesh, const int *ids, int count,
    int use_bf1, double **values)
{
    int b;
    int pos = 0;
    size_t nvalue;

    *values = 0;
    if (count <= 0) {
        return 0;
    }
    nvalue = (size_t)count * 3U * (size_t)PRJ_BLOCK_NCELLS;
    *values = (double *)malloc(nvalue * sizeof(**values));
    if (*values == 0) {
        return 1;
    }
    for (b = 0; b < count; ++b) {
        const prj_block *block;
        int dir;

        if (ids[b] < 0 || ids[b] >= mesh->nblocks) {
            free(*values);
            *values = 0;
            return 1;
        }
        block = &mesh->blocks[ids[b]];
        for (dir = 0; dir < 3; ++dir) {
            const double *src = prj_mpi_bf_array_const(block, dir, use_bf1);
            int n;

            if (src == 0) {
                free(*values);
                *values = 0;
                return 1;
            }
            for (n = 0; n < PRJ_BLOCK_NCELLS; ++n) {
                double value = src[n];

                if (!isfinite(value)) {
                    free(*values);
                    *values = 0;
                    return 1;
                }
                (*values)[pos++] = value;
            }
        }
    }
    return 0;
}

static int prj_mpi_install_remote_bf_blocks(prj_mesh *mesh, const int *ids, int count,
    const double *values, int use_bf1, int *old_ranks)
{
    int b;
    int pos = 0;

    if (count <= 0) {
        return 0;
    }
    if (ids == 0 || values == 0 || old_ranks == 0) {
        return 1;
    }
    for (b = 0; b < count; ++b) {
        prj_block *block;
        double *base;
        int dir;
        int n;

        if (ids[b] < 0 || ids[b] >= mesh->nblocks) {
            return 1;
        }
        block = &mesh->blocks[ids[b]];
        if (block->rank == prj_mpi_active->rank || block->Bf[0] != 0 ||
            block->Bf1[0] != 0 || block->face_fidelity[0] != 0 ||
            block->face_fidelity[1] != 0 || block->face_fidelity[2] != 0) {
            return 1;
        }
        base = (double *)calloc(6U * (size_t)PRJ_BLOCK_NCELLS, sizeof(*base));
        for (dir = 0; dir < 3; ++dir) {
            block->face_fidelity[dir] = (int *)calloc((size_t)PRJ_BLOCK_NCELLS,
                sizeof(*block->face_fidelity[dir]));
        }
        if (base == 0 || block->face_fidelity[0] == 0 ||
            block->face_fidelity[1] == 0 || block->face_fidelity[2] == 0) {
            free(base);
            for (dir = 0; dir < 3; ++dir) {
                free(block->face_fidelity[dir]);
                block->face_fidelity[dir] = 0;
            }
            return 1;
        }
        for (dir = 0; dir < 3; ++dir) {
            block->Bf[dir] = base + (size_t)dir * (size_t)PRJ_BLOCK_NCELLS;
            block->Bf1[dir] = base + (size_t)(dir + 3) * (size_t)PRJ_BLOCK_NCELLS;
        }
        for (dir = 0; dir < 3; ++dir) {
            double *dst = prj_mpi_bf_array(block, dir, use_bf1);

            for (n = 0; n < PRJ_BLOCK_NCELLS; ++n) {
                dst[n] = values[pos++];
            }
        }
        for (dir = 0; dir < 3; ++dir) {
            for (n = 0; n < PRJ_BLOCK_NCELLS; ++n) {
                block->face_fidelity[dir][n] = PRJ_MHD_FIDELITY_SAME;
            }
        }
        old_ranks[b] = block->rank;
        block->rank = prj_mpi_active->rank;
    }
    return 0;
}

static void prj_mpi_uninstall_remote_bf_blocks(prj_mesh *mesh, const int *ids, int count,
    const int *old_ranks)
{
    int b;

    if (mesh == 0 || ids == 0 || old_ranks == 0) {
        return;
    }
    for (b = 0; b < count; ++b) {
        prj_block *block;
        int dir;

        if (ids[b] < 0 || ids[b] >= mesh->nblocks) {
            continue;
        }
        block = &mesh->blocks[ids[b]];
        free(block->Bf[0]);
        for (dir = 0; dir < 3; ++dir) {
            free(block->face_fidelity[dir]);
            block->face_fidelity[dir] = 0;
            block->Bf[dir] = 0;
            block->Bf1[dir] = 0;
        }
        block->rank = old_ranks[b];
    }
}

void prj_mpi_exchange_bf(prj_mesh *mesh, prj_mpi *mpi, int use_bf1, int fill_kind)
{
#if defined(PRJ_ENABLE_MPI)
    int nb;

    if (mesh == 0 || mpi == 0 || mpi->totrank <= 1) {
        return;
    }
    for (nb = 0; nb < mpi->neighbor_number; ++nb) {
        prj_mpi_buffer *buffer = &mpi->neighbor_buffer[nb];
        int *send_ids = 0;
        int *recv_ids = 0;
        int send_count = 0;
        int recv_count = 0;
        double *send_values = 0;
        double *recv_values = 0;
        int *old_ranks = 0;
        int b;
        int send_value_count;
        int recv_value_count;

        if (prj_mpi_collect_bf_source_blocks(mesh, mpi, buffer->receiver_rank,
                &send_ids, &send_count) != 0 ||
            prj_mpi_pack_bf_blocks(mesh, send_ids, send_count, use_bf1, &send_values) != 0) {
            fprintf(stderr, "prj_mpi_exchange_bf: failed to pack local Bf blocks\n");
            free(send_ids);
            free(send_values);
            exit(EXIT_FAILURE);
        }
        MPI_Sendrecv(&send_count, 1, MPI_INT, buffer->receiver_rank, 300,
            &recv_count, 1, MPI_INT, buffer->receiver_rank, 300,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        recv_ids = recv_count > 0 ? (int *)calloc((size_t)recv_count, sizeof(*recv_ids)) : 0;
        if (recv_count > 0 && recv_ids == 0) {
            fprintf(stderr, "prj_mpi_exchange_bf: failed to allocate received block ids\n");
            free(send_ids);
            free(send_values);
            exit(EXIT_FAILURE);
        }
        MPI_Sendrecv(send_ids, send_count, MPI_INT, buffer->receiver_rank, 301,
            recv_ids, recv_count, MPI_INT, buffer->receiver_rank, 301,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        send_value_count = send_count * 3 * PRJ_BLOCK_NCELLS;
        recv_value_count = recv_count * 3 * PRJ_BLOCK_NCELLS;
        recv_values = recv_value_count > 0 ?
            (double *)malloc((size_t)recv_value_count * sizeof(*recv_values)) : 0;
        if (recv_value_count > 0 && recv_values == 0) {
            fprintf(stderr, "prj_mpi_exchange_bf: failed to allocate received Bf values\n");
            free(recv_ids);
            free(send_ids);
            free(send_values);
            exit(EXIT_FAILURE);
        }
        MPI_Sendrecv(send_values, send_value_count, MPI_DOUBLE, buffer->receiver_rank, 302,
            recv_values, recv_value_count, MPI_DOUBLE, buffer->receiver_rank, 302,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        old_ranks = recv_count > 0 ? (int *)calloc((size_t)recv_count, sizeof(*old_ranks)) : 0;
        if (recv_count > 0 && old_ranks == 0) {
            fprintf(stderr, "prj_mpi_exchange_bf: failed to allocate rank scratch\n");
            free(recv_values);
            free(recv_ids);
            free(send_ids);
            free(send_values);
            exit(EXIT_FAILURE);
        }
        if (prj_mpi_install_remote_bf_blocks(mesh, recv_ids, recv_count,
                recv_values, use_bf1, old_ranks) != 0) {
            fprintf(stderr, "prj_mpi_exchange_bf: failed to install remote Bf scratch blocks\n");
            prj_mpi_uninstall_remote_bf_blocks(mesh, recv_ids, recv_count, old_ranks);
            free(old_ranks);
            free(recv_values);
            free(recv_ids);
            free(send_ids);
            free(send_values);
            exit(EXIT_FAILURE);
        }
        for (b = 0; b < recv_count; ++b) {
            prj_boundary_send_bf(&mesh->blocks[recv_ids[b]], use_bf1, fill_kind);
        }
        prj_mpi_uninstall_remote_bf_blocks(mesh, recv_ids, recv_count, old_ranks);
        free(old_ranks);
        free(recv_values);
        free(recv_ids);
        free(send_values);
        free(send_ids);
    }
#else
    (void)mesh;
    (void)mpi;
    (void)use_bf1;
    (void)fill_kind;
#endif
}
#endif

#if PRJ_MHD
static int prj_mpi_edge_axis_active_max(int dir, int axis)
{
    return dir == axis ? PRJ_BLOCK_SIZE - 1 : PRJ_BLOCK_SIZE;
}

static double prj_mpi_edge_coord(const prj_block *block, int axis, int dir, int idx)
{
    double offset = axis == dir ? 0.5 : 0.0;

    return block->xmin[axis] + ((double)idx + offset) * block->dx[axis];
}

static int prj_mpi_edge_point_inside(const prj_block *block, const double x[3])
{
    const double tol = 1.0e-12;
    int d;

    if (block == 0) {
        return 0;
    }
    for (d = 0; d < 3; ++d) {
        if (x[d] < block->xmin[d] - tol || x[d] > block->xmax[d] + tol) {
            return 0;
        }
    }
    return 1;
}

static int prj_mpi_storage_index_ok(int i, int j, int k)
{
    return i >= -PRJ_NGHOST && i < PRJ_BLOCK_SIZE + PRJ_NGHOST &&
        j >= -PRJ_NGHOST && j < PRJ_BLOCK_SIZE + PRJ_NGHOST &&
        k >= -PRJ_NGHOST && k < PRJ_BLOCK_SIZE + PRJ_NGHOST;
}

static int prj_mpi_nearest_int(double x)
{
    return x >= 0.0 ? (int)(x + 0.5) : (int)(x - 0.5);
}

static double prj_mpi_restrict_emf_value(const prj_block *fine, int dir,
    const double x[3])
{
    const double tol = 1.0e-8;
    int idx[3] = {0, 0, 0};
    int edge_base = 0;
    double sum = 0.0;
    int n;
    int d;

    if (fine == 0 || dir < 0 || dir >= 3 || fine->emf[dir] == 0) {
        fprintf(stderr, "prj_mpi_restrict_emf_value: invalid fine emf block\n");
        exit(EXIT_FAILURE);
    }
    for (d = 0; d < 3; ++d) {
        double q;

        if (fine->dx[d] <= 0.0) {
            fprintf(stderr, "prj_mpi_restrict_emf_value: invalid fine cell size\n");
            exit(EXIT_FAILURE);
        }
        q = (x[d] - fine->xmin[d]) / fine->dx[d];
        if (d == dir) {
            double r = q - 0.5;

            edge_base = prj_mpi_nearest_int(r - 0.5);
            if (fabs(r - ((double)edge_base + 0.5)) > tol) {
                fprintf(stderr, "prj_mpi_restrict_emf_value: coarse edge is not centered on two fine edges\n");
                exit(EXIT_FAILURE);
            }
            idx[d] = edge_base;
        } else {
            int edge_idx = prj_mpi_nearest_int(q);

            if (fabs(q - (double)edge_idx) > tol) {
                fprintf(stderr, "prj_mpi_restrict_emf_value: coarse edge is not fine-aligned\n");
                exit(EXIT_FAILURE);
            }
            idx[d] = edge_idx;
        }
    }
    for (n = 0; n < 2; ++n) {
        int eidx[3];
        double value;

        eidx[0] = idx[0];
        eidx[1] = idx[1];
        eidx[2] = idx[2];
        eidx[dir] = edge_base + n;
        if (!prj_mpi_storage_index_ok(eidx[0], eidx[1], eidx[2])) {
            fprintf(stderr, "prj_mpi_restrict_emf_value: fine edge index out of storage\n");
            exit(EXIT_FAILURE);
        }
        value = fine->emf[dir][IDX(eidx[0], eidx[1], eidx[2])];
        if (!isfinite(value)) {
            fprintf(stderr, "prj_mpi_restrict_emf_value: non-finite emf\n");
            exit(EXIT_FAILURE);
        }
        sum += value;
    }
    return 0.5 * sum;
}

static int prj_mpi_pack_emf_records(prj_mesh *mesh, prj_mpi *mpi,
    int receiver_rank, int **idx_send, double **value_send, int *count)
{
    int idx_capacity = 0;
    int value_capacity = 0;
    int value_count = 0;
    int bidx;

    idx_send[0] = 0;
    idx_send[1] = 0;
    idx_send[2] = 0;
    *value_send = 0;
    *count = 0;
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        int slot;

        if (!prj_mpi_block_is_local(block)) {
            continue;
        }
        for (slot = 0; slot < 56; ++slot) {
            int nid = block->slot[slot].id;
            const prj_block *coarse;
            int dir;

            if (nid < 0 || nid >= mesh->nblocks || block->slot[slot].rel_level >= 0 ||
                mesh->blocks[nid].rank != receiver_rank) {
                continue;
            }
            coarse = &mesh->blocks[nid];
            for (dir = 0; dir < 3; ++dir) {
                int i;
                int j;
                int k;

                for (i = 0; i <= prj_mpi_edge_axis_active_max(dir, 0); ++i) {
                    for (j = 0; j <= prj_mpi_edge_axis_active_max(dir, 1); ++j) {
                        for (k = 0; k <= prj_mpi_edge_axis_active_max(dir, 2); ++k) {
                            double x[3];
                            double value;

                            x[0] = prj_mpi_edge_coord(coarse, 0, dir, i);
                            x[1] = prj_mpi_edge_coord(coarse, 1, dir, j);
                            x[2] = prj_mpi_edge_coord(coarse, 2, dir, k);
                            if (!prj_mpi_edge_point_inside(block, x)) {
                                continue;
                            }
                            value = prj_mpi_restrict_emf_value(block, dir, x);
                            if (prj_mpi_append_triplet(idx_send, count, &idx_capacity,
                                    nid, prj_mpi_encode_cell_index(i, j, k), dir) != 0 ||
                                prj_mpi_append_double(value_send, &value_count, &value_capacity, value) != 0) {
                                free(idx_send[0]);
                                free(idx_send[1]);
                                free(idx_send[2]);
                                free(*value_send);
                                idx_send[0] = 0;
                                idx_send[1] = 0;
                                idx_send[2] = 0;
                                *value_send = 0;
                                *count = 0;
                                return 1;
                            }
                        }
                    }
                }
            }
        }
    }
    (void)mpi;
    return 0;
}

void prj_mpi_exchange_emf(prj_mesh *mesh, prj_mpi *mpi)
{
#if defined(PRJ_ENABLE_MPI)
    int nb;

    if (mesh == 0 || mpi == 0 || mpi->totrank <= 1) {
        return;
    }
    for (nb = 0; nb < mpi->neighbor_number; ++nb) {
        prj_mpi_buffer *buffer = &mpi->neighbor_buffer[nb];
        int *idx_send[3];
        int *idx_recv[3];
        double *value_send = 0;
        double *value_recv = 0;
        int send_count = 0;
        int recv_count = 0;
        int axis;
        int i;

        if (prj_mpi_pack_emf_records(mesh, mpi, buffer->receiver_rank,
                idx_send, &value_send, &send_count) != 0) {
            fprintf(stderr, "prj_mpi_exchange_emf: failed to pack emf records\n");
            exit(EXIT_FAILURE);
        }
        MPI_Sendrecv(&send_count, 1, MPI_INT, buffer->receiver_rank, 400,
            &recv_count, 1, MPI_INT, buffer->receiver_rank, 400,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (axis = 0; axis < 3; ++axis) {
            idx_recv[axis] = recv_count > 0 ?
                (int *)calloc((size_t)recv_count, sizeof(*idx_recv[axis])) : 0;
            if (recv_count > 0 && idx_recv[axis] == 0) {
                fprintf(stderr, "prj_mpi_exchange_emf: failed to allocate received indices\n");
                exit(EXIT_FAILURE);
            }
        }
        value_recv = recv_count > 0 ?
            (double *)malloc((size_t)recv_count * sizeof(*value_recv)) : 0;
        if (recv_count > 0 && value_recv == 0) {
            fprintf(stderr, "prj_mpi_exchange_emf: failed to allocate received values\n");
            exit(EXIT_FAILURE);
        }
        for (axis = 0; axis < 3; ++axis) {
            MPI_Sendrecv(idx_send[axis], send_count, MPI_INT, buffer->receiver_rank, 401 + axis,
                idx_recv[axis], recv_count, MPI_INT, buffer->receiver_rank, 401 + axis,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        MPI_Sendrecv(value_send, send_count, MPI_DOUBLE, buffer->receiver_rank, 404,
            value_recv, recv_count, MPI_DOUBLE, buffer->receiver_rank, 404,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for (i = 0; i < recv_count; ++i) {
            int block_id = idx_recv[0][i];
            int code = idx_recv[1][i];
            int dir = idx_recv[2][i];
            int ii;
            int jj;
            int kk;
            prj_block *block;
            int flat;

            if (block_id < 0 || block_id >= mesh->nblocks || dir < 0 || dir >= 3) {
                continue;
            }
            block = &mesh->blocks[block_id];
            if (!prj_mpi_block_is_local(block) || block->edge_fidelity[dir] == 0 ||
                block->emf[dir] == 0) {
                continue;
            }
            prj_mpi_decode_cell_index(code, &ii, &jj, &kk);
            if (!prj_mpi_storage_index_ok(ii, jj, kk)) {
                fprintf(stderr, "prj_mpi_exchange_emf: received edge index out of storage\n");
                exit(EXIT_FAILURE);
            }
            flat = IDX(ii, jj, kk);
            if (PRJ_MHD_FIDELITY_FINER < block->edge_fidelity[dir][flat]) {
                continue;
            }
            if (!isfinite(value_recv[i])) {
                fprintf(stderr, "prj_mpi_exchange_emf: received non-finite emf\n");
                exit(EXIT_FAILURE);
            }
            block->emf[dir][flat] = value_recv[i];
            block->edge_fidelity[dir][flat] = PRJ_MHD_FIDELITY_FINER;
        }

        for (axis = 0; axis < 3; ++axis) {
            free(idx_send[axis]);
            free(idx_recv[axis]);
        }
        free(value_send);
        free(value_recv);
    }
#else
    (void)mesh;
    (void)mpi;
#endif
}
#endif

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
