#include <math.h>
#include <limits.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "prj.h"

#if defined(PRJ_ENABLE_MPI)
#include <mpi.h>
#endif

/* Ghost-zone exchange buffer is split into two concatenated streams per
 * neighbor: the front stream covers hydro primitives + EOS vars for every
 * ghost cell, the back stream covers radiation primitives only for cells in
 * the narrower radiation ghost band. */
#define PRJ_MPI_GHOST_NVAR_HE  (PRJ_NHYDRO + PRJ_NVAR_EOSVAR)
#define PRJ_MPI_GHOST_NVAR_RAD (PRJ_NRAD_VAR)
#define PRJ_MPI_GHOST_FILL_KIND_N 6
#define PRJ_MPI_SAMPLE_SAME_LEVEL_OFFSET 4

static int prj_block_is_active(const prj_block *block)
{
    return block != 0 && block->id >= 0 && block->active == 1;
}

#if PRJ_MHD
static int prj_mpi_block_is_local(const prj_mpi *mpi, const prj_block *block)
{
    return prj_block_is_active(block) &&
        (mpi == 0 || block->rank == mpi->rank);
}
#endif

static double *prj_mpi_stage_array(prj_block *block, int stage)
{
    return stage == 2 ? block->W1 : block->W;
}

static void prj_mpi_buffer_free(prj_mpi_buffer *buffer)
{
    int axis;
    int fill_kind;

    if (buffer == 0) {
        return;
    }
    free(buffer->receiver_blocks);
    free(buffer->cell_data_size_send);
    free(buffer->cell_buffer_send);
    for (fill_kind = 0; fill_kind < PRJ_MPI_GHOST_FILL_KIND_N; ++fill_kind) {
        free(buffer->cell_buffer_send_by_kind[fill_kind]);
        free(buffer->cell_buffer_recv_by_kind[fill_kind]);
    }
    free(buffer->face_data_size_send);
    free(buffer->face_buffer_send);
    free(buffer->cell_data_size_recv);
    free(buffer->cell_buffer_recv);
    free(buffer->face_data_size_recv);
    free(buffer->face_buffer_recv);
    free(buffer->cell_data_size_send_rad);
    free(buffer->cell_data_size_recv_rad);
    free(buffer->flux_idx_send);
    free(buffer->flux_value_send);
    free(buffer->flux_idx_recv);
    free(buffer->flux_value_recv);
    for (axis = 0; axis < 3; ++axis) {
        free(buffer->cell_data_idx_send[axis]);
        free(buffer->face_data_idx_send[axis]);
        free(buffer->cell_data_idx_recv[axis]);
        free(buffer->face_data_idx_recv[axis]);
        free(buffer->cell_data_idx_send_rad[axis]);
        free(buffer->cell_data_idx_recv_rad[axis]);
    }
#if PRJ_MHD
    free(buffer->bf_headers_send);
    free(buffer->bf_values_send);
    free(buffer->bf_headers_recv);
    free(buffer->bf_values_recv);
    free(buffer->emf_value_send);
    free(buffer->emf_value_recv);
    free(buffer->emf_src_block);
    free(buffer->emf_src_dir);
    free(buffer->amr_bf_values_send);
    free(buffer->amr_bf_headers_recv);
    free(buffer->amr_bf_values_recv);
    for (axis = 0; axis < 3; ++axis) {
        free(buffer->emf_idx_send[axis]);
        free(buffer->emf_idx_recv[axis]);
        free(buffer->emf_src_idx[axis]);
    }
#endif
    memset(buffer, 0, sizeof(*buffer));
}

static void prj_mpi_clear_shared_buffers(prj_mpi *mpi)
{
    if (mpi == 0) {
        return;
    }
    free(mpi->request_buffer);
    mpi->request_buffer = 0;
    mpi->request_capacity = 0;
#if PRJ_MHD
    free(mpi->amr_bf_headers);
    free(mpi->amr_bf_values);
    mpi->amr_bf_headers = 0;
    mpi->amr_bf_values = 0;
    mpi->amr_bf_record_count = 0;
    mpi->amr_bf_value_count = 0;
    mpi->amr_bf_record_capacity = 0;
    mpi->amr_bf_value_capacity = 0;
#endif
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
    prj_mpi_clear_shared_buffers(mpi);
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

typedef struct {
    unsigned long long key;
    int id;
} prj_mpi_order_item;

static int prj_mpi_order_item_less(const prj_mpi_order_item *a,
    const prj_mpi_order_item *b)
{
    if (a->key != b->key) {
        return a->key < b->key;
    }
    return a->id < b->id;
}

static void prj_mpi_order_merge(prj_mpi_order_item *items,
    prj_mpi_order_item *scratch, int left, int mid, int right)
{
    int i = left;
    int j = mid;
    int k = left;

    while (i < mid && j < right) {
        if (!prj_mpi_order_item_less(&items[j], &items[i])) {
            scratch[k++] = items[i++];
        } else {
            scratch[k++] = items[j++];
        }
    }
    while (i < mid) {
        scratch[k++] = items[i++];
    }
    while (j < right) {
        scratch[k++] = items[j++];
    }
    for (i = left; i < right; ++i) {
        items[i] = scratch[i];
    }
}

static void prj_mpi_order_sort_range(prj_mpi_order_item *items,
    prj_mpi_order_item *scratch, int left, int right)
{
    int mid;

    if (right - left <= 1) {
        return;
    }
    mid = left + (right - left) / 2;
    prj_mpi_order_sort_range(items, scratch, left, mid);
    prj_mpi_order_sort_range(items, scratch, mid, right);
    if (!prj_mpi_order_item_less(&items[mid], &items[mid - 1])) {
        return;
    }
    prj_mpi_order_merge(items, scratch, left, mid, right);
}

static int prj_mpi_order_sort(prj_mpi_order_item *items, int count)
{
    prj_mpi_order_item *scratch;

    if (count <= 1) {
        return 0;
    }
    scratch = (prj_mpi_order_item *)malloc((size_t)count * sizeof(*scratch));
    if (scratch == 0) {
        return 1;
    }
    prj_mpi_order_sort_range(items, scratch, 0, count);
    free(scratch);
    return 0;
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

static int prj_mpi_sample_kind_from_code(int sample_code)
{
    return sample_code;
}

static void prj_mpi_assign_block_storage(prj_mesh *mesh, const prj_mpi *mpi)
{
    int bidx;

    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (block->id < 0) {
            continue;
        }
        if (mpi != 0 && block->rank != mpi->rank) {
            if (block->W != 0) {
                prj_block_free_data(block);
            }
            continue;
        }
        if (block->W == 0) {
            prj_block_alloc_data(block);
        }
    }
    prj_mesh_update_cell_derived_mask(mesh);
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
        + 5U * (size_t)PRJ_BLOCK_NCELLS
        + (size_t)(LMAX*LMAX) * (size_t)PRJ_BLOCK_NCELLS
#if PRJ_MHD
        + 6U * (size_t)PRJ_BLOCK_NFACES + 6U * (size_t)PRJ_BLOCK_NCELLS + 3U * (size_t)PRJ_BLOCK_NEDGES
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

static void prj_mpi_compute_decomposition(prj_mesh *mesh, const prj_mpi *mpi)
{
    prj_mpi_order_item *items;
    int count;
    int i;
    int cursor;
    int rank;
    int nrank;

    if (mesh == 0) {
        return;
    }
    count = 0;
    for (i = 0; i < mesh->nblocks; ++i) {
        if (prj_block_is_active(&mesh->blocks[i])) {
            count += 1;
        }
    }
    items = (prj_mpi_order_item *)malloc((size_t)count * sizeof(*items));
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
    if (prj_mpi_order_sort(items, count) != 0) {
        free(items);
        return;
    }
    nrank = mpi != 0 ? mpi->totrank : 1;
    for (i = 0; i < count; ++i) {
        rank = (int)(((long long)i * (long long)nrank) / (long long)count);
        if (rank >= nrank) {
            rank = nrank - 1;
        }
        mesh->blocks[items[i].id].rank = rank;
    }
    prj_mpi_sync_slot_ranks(mesh);
    free(items);
}

static void prj_mpi_migrate_active_blocks(prj_mesh *mesh, const prj_mpi *mpi, const int *old_ranks)
{
#if defined(PRJ_ENABLE_MPI)
    size_t data_count;
    int bidx;

    if (mesh == 0 || mpi == 0 || old_ranks == 0 || mpi->totrank <= 1) {
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
        if (mpi->rank == new_rank || mpi->rank == old_rank) {
            double *sendbuf = 0;
            double *recvbuf = 0;
            int source = MPI_PROC_NULL;
            int dest = MPI_PROC_NULL;
            int recvcount = 0;
            int sendcount = 0;

            if (mpi->rank == new_rank) {
                if (block->W == 0 && prj_block_alloc_data(block) != 0) {
                    continue;
                }
                recvbuf = block->W;
                recvcount = (int)data_count;
                source = old_rank;
            }
            if (mpi->rank == old_rank && block->W != 0) {
                sendbuf = block->W;
                sendcount = (int)data_count;
                dest = new_rank;
            }

            MPI_Sendrecv(sendbuf, sendcount, MPI_DOUBLE, dest, tag,
                recvbuf, recvcount, MPI_DOUBLE, source, tag,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            if (mpi->rank == old_rank && block->W != 0) {
                prj_block_free_data(block);
            }
        }
    }
#else
    (void)mesh;
    (void)mpi;
    (void)old_ranks;
#endif
}

static void prj_mpi_print_balance(const prj_mesh *mesh, const prj_mpi *mpi)
{
#if defined(PRJ_ENABLE_MPI)
    int *counts;
    int local_active;
    int bidx;
    int rank;
    int min_count;
    int max_count;
    int total_active;

    if (mesh == 0 || mpi == 0 || mpi->totrank <= 1) {
        return;
    }

    counts = (int *)calloc((size_t)mpi->totrank, sizeof(*counts));
    if (counts == 0) {
        return;
    }

    local_active = 0;
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        if (prj_block_is_active(&mesh->blocks[bidx]) &&
            mesh->blocks[bidx].rank == mpi->rank) {
            local_active += 1;
        }
    }

    MPI_Gather(&local_active, 1, MPI_INT, counts, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (mpi->rank == 0) {
        min_count = counts[0];
        max_count = counts[0];
        total_active = counts[0];
        for (rank = 1; rank < mpi->totrank; ++rank) {
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
    (void)mpi;
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

static void prj_mpi_collect_active_counts(const prj_mesh *mesh, const prj_mpi *mpi, int *counts)
{
    int bidx;

    if (mesh == 0 || counts == 0 || mpi == 0) {
        return;
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        const prj_block *block = &mesh->blocks[bidx];

        if (!prj_block_is_active(block)) {
            continue;
        }
        if (block->rank >= 0 && block->rank < mpi->totrank) {
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

static int prj_mpi_ghost_fill_kind_ok(int fill_kind)
{
    return fill_kind >= 0 && fill_kind < PRJ_MPI_GHOST_FILL_KIND_N;
}

static int prj_mpi_ghost_fill_kind_from_rel_level(int rel_level)
{
    if (rel_level == 0) {
        return PRJ_BOUNDARY_FILL_SAME_LEVEL;
    }
    if (rel_level < 0) {
        return PRJ_BOUNDARY_FILL_RESTRICTION;
    }
    return PRJ_BOUNDARY_FILL_PROLONGATION;
}

static int prj_mpi_alloc_ghost_value_buffer(double **buffer,
    int cell_count, int cell_count_rad)
{
    size_t value_count;

    if (buffer == 0 || cell_count < 0 || cell_count_rad < 0) {
        return 1;
    }
    value_count = (size_t)cell_count * (size_t)PRJ_MPI_GHOST_NVAR_HE +
        (size_t)cell_count_rad * (size_t)PRJ_MPI_GHOST_NVAR_RAD;
    if (value_count == 0) {
        *buffer = 0;
        return 0;
    }
    *buffer = (double *)calloc(value_count, sizeof(**buffer));
    return *buffer == 0 ? 1 : 0;
}

static void prj_mpi_count_recv_ghosts_by_kind(prj_mpi_buffer *buffer,
    int cell_count, int cell_count_rad)
{
    int i;

    if (buffer == 0) {
        return;
    }
    memset(buffer->cell_recv_count_by_kind, 0,
        sizeof(buffer->cell_recv_count_by_kind));
    memset(buffer->cell_recv_count_rad_by_kind, 0,
        sizeof(buffer->cell_recv_count_rad_by_kind));
    for (i = 0; i < cell_count; ++i) {
        int fill_kind = prj_mpi_sample_kind_from_code(buffer->cell_data_idx_recv[2][i]);

        if (prj_mpi_ghost_fill_kind_ok(fill_kind)) {
            buffer->cell_recv_count_by_kind[fill_kind] += 1;
        }
    }
    for (i = 0; i < cell_count_rad; ++i) {
        int fill_kind = prj_mpi_sample_kind_from_code(buffer->cell_data_idx_recv_rad[2][i]);

        if (prj_mpi_ghost_fill_kind_ok(fill_kind)) {
            buffer->cell_recv_count_rad_by_kind[fill_kind] += 1;
        }
    }
}

static int prj_mpi_build_ghost_plan_for_neighbor(prj_mesh *mesh, prj_mpi *mpi, prj_mpi_buffer *buffer)
{
#if defined(PRJ_ENABLE_MPI)
    int *send_sizes;
    int *send_sizes_rad;
    int send_entries;
    int cell_size_total;
    int cell_size_total_rad;
    int axis;
    int fill_kind;
    int bidx;
    int occ;
    int *idx_send[3];
    int *idx_send_rad[3];
    int record_count;
    int record_count_rad;
    int record_count_by_kind[PRJ_MPI_GHOST_FILL_KIND_N];
    int record_count_rad_by_kind[PRJ_MPI_GHOST_FILL_KIND_N];
    int pos;
    int pos_rad;
    MPI_Request requests[12];
    int req_n = 0;

    send_sizes = buffer->number > 0 ?
        (int *)calloc((size_t)buffer->number, sizeof(*send_sizes)) : 0;
    send_sizes_rad = buffer->number > 0 ?
        (int *)calloc((size_t)buffer->number, sizeof(*send_sizes_rad)) : 0;
    if (buffer->number > 0 && (send_sizes == 0 || send_sizes_rad == 0)) {
        free(send_sizes);
        free(send_sizes_rad);
        return 1;
    }
    idx_send[0] = 0;
    idx_send[1] = 0;
    idx_send[2] = 0;
    idx_send_rad[0] = 0;
    idx_send_rad[1] = 0;
    idx_send_rad[2] = 0;
    record_count = 0;
    record_count_rad = 0;
    memset(record_count_by_kind, 0, sizeof(record_count_by_kind));
    memset(record_count_rad_by_kind, 0, sizeof(record_count_rad_by_kind));
    occ = 0;
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        int n;

        if (!prj_block_is_active(block) || block->rank != mpi->rank) {
            continue;
        }
        for (n = 0; n < 56; ++n) {
            const prj_neighbor *slot = &block->slot[n];
            int nid = slot->id;
            int i;
            int j;
            int k;
            int before;
            int before_rad;
            int sample_kind;

            if (nid < 0 || nid >= mesh->nblocks || mesh->blocks[nid].rank != buffer->receiver_rank) {
                continue;
            }
            sample_kind = prj_mpi_ghost_fill_kind_from_rel_level(slot->rel_level);
            before = record_count;
            before_rad = record_count_rad;

            for (i = slot->recv_loc_start[0]; i < slot->recv_loc_end[0]; ++i) {
                for (j = slot->recv_loc_start[1]; j < slot->recv_loc_end[1]; ++j) {
                    for (k = slot->recv_loc_start[2]; k < slot->recv_loc_end[2]; ++k) {
                        record_count += 1;
                        record_count_by_kind[sample_kind] += 1;
                    }
                }
            }
            for (i = slot->recv_loc_start_rad[0]; i < slot->recv_loc_end_rad[0]; ++i) {
                for (j = slot->recv_loc_start_rad[1]; j < slot->recv_loc_end_rad[1]; ++j) {
                    for (k = slot->recv_loc_start_rad[2]; k < slot->recv_loc_end_rad[2]; ++k) {
                        record_count_rad += 1;
                        record_count_rad_by_kind[sample_kind] += 1;
                    }
                }
            }
            if (occ >= buffer->number) {
                free(send_sizes);
                free(send_sizes_rad);
                return 1;
            }
            send_sizes[occ] = record_count - before;
            send_sizes_rad[occ] = record_count_rad - before_rad;
            occ += 1;
        }
    }
    if (occ != buffer->number) {
        free(send_sizes);
        free(send_sizes_rad);
        return 1;
    }

    if (record_count > 0) {
        for (axis = 0; axis < 3; ++axis) {
            idx_send[axis] = (int *)calloc((size_t)record_count,
                sizeof(*idx_send[axis]));
        }
        if (idx_send[0] == 0 || idx_send[1] == 0 || idx_send[2] == 0) {
            free(send_sizes);
            free(send_sizes_rad);
            for (axis = 0; axis < 3; ++axis) {
                free(idx_send[axis]);
            }
            return 1;
        }
    }
    if (record_count_rad > 0) {
        for (axis = 0; axis < 3; ++axis) {
            idx_send_rad[axis] = (int *)calloc((size_t)record_count_rad,
                sizeof(*idx_send_rad[axis]));
        }
        if (idx_send_rad[0] == 0 || idx_send_rad[1] == 0 || idx_send_rad[2] == 0) {
            free(send_sizes);
            free(send_sizes_rad);
            for (axis = 0; axis < 3; ++axis) {
                free(idx_send[axis]);
                free(idx_send_rad[axis]);
            }
            return 1;
        }
    }
    for (fill_kind = 0; fill_kind < PRJ_MPI_GHOST_FILL_KIND_N; ++fill_kind) {
        buffer->cell_send_count_by_kind[fill_kind] = record_count_by_kind[fill_kind];
        buffer->cell_send_count_rad_by_kind[fill_kind] = record_count_rad_by_kind[fill_kind];
        if (prj_mpi_alloc_ghost_value_buffer(
                &buffer->cell_buffer_send_by_kind[fill_kind],
                record_count_by_kind[fill_kind],
                record_count_rad_by_kind[fill_kind]) != 0) {
            free(send_sizes);
            free(send_sizes_rad);
            for (axis = 0; axis < 3; ++axis) {
                free(idx_send[axis]);
                free(idx_send_rad[axis]);
            }
            return 1;
        }
    }

    pos = 0;
    pos_rad = 0;
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        int n;

        if (!prj_block_is_active(block) || block->rank != mpi->rank) {
            continue;
        }
        for (n = 0; n < 56; ++n) {
            const prj_neighbor *slot = &block->slot[n];
            int nid = slot->id;
            int i;
            int j;
            int k;
            int sample_kind;
            int sample_code;

            if (nid < 0 || nid >= mesh->nblocks || mesh->blocks[nid].rank != buffer->receiver_rank) {
                continue;
            }
            if (slot->rel_level == 0) {
                sample_kind = PRJ_BOUNDARY_FILL_SAME_LEVEL;
            } else if (slot->rel_level < 0) {
                sample_kind = PRJ_BOUNDARY_FILL_RESTRICTION;
            } else {
                sample_kind = PRJ_BOUNDARY_FILL_PROLONGATION;
            }
            sample_code = sample_kind;

            for (i = slot->recv_loc_start[0]; i < slot->recv_loc_end[0]; ++i) {
                for (j = slot->recv_loc_start[1]; j < slot->recv_loc_end[1]; ++j) {
                    for (k = slot->recv_loc_start[2]; k < slot->recv_loc_end[2]; ++k) {
                        if (pos >= record_count) {
                            free(send_sizes);
                            free(send_sizes_rad);
                            free(buffer->cell_buffer_send);
                            buffer->cell_buffer_send = 0;
                            for (axis = 0; axis < 3; ++axis) {
                                free(idx_send[axis]);
                                free(idx_send_rad[axis]);
                            }
                            return 1;
                        }
                        idx_send[0][pos] = nid;
                        idx_send[1][pos] = prj_mpi_encode_cell_index(i, j, k);
                        idx_send[2][pos] = sample_code;
                        pos += 1;
                    }
                }
            }
            for (i = slot->recv_loc_start_rad[0]; i < slot->recv_loc_end_rad[0]; ++i) {
                for (j = slot->recv_loc_start_rad[1]; j < slot->recv_loc_end_rad[1]; ++j) {
                    for (k = slot->recv_loc_start_rad[2]; k < slot->recv_loc_end_rad[2]; ++k) {
                        if (pos_rad >= record_count_rad) {
                            free(send_sizes);
                            free(send_sizes_rad);
                            free(buffer->cell_buffer_send);
                            buffer->cell_buffer_send = 0;
                            for (axis = 0; axis < 3; ++axis) {
                                free(idx_send[axis]);
                                free(idx_send_rad[axis]);
                            }
                            return 1;
                        }
                        idx_send_rad[0][pos_rad] = nid;
                        idx_send_rad[1][pos_rad] = prj_mpi_encode_cell_index(i, j, k);
                        idx_send_rad[2][pos_rad] = sample_code;
                        pos_rad += 1;
                    }
                }
            }
        }
    }
    if (pos != record_count || pos_rad != record_count_rad) {
        free(send_sizes);
        free(send_sizes_rad);
        free(buffer->cell_buffer_send);
        buffer->cell_buffer_send = 0;
        for (axis = 0; axis < 3; ++axis) {
            free(idx_send[axis]);
            free(idx_send_rad[axis]);
        }
        return 1;
    }

    buffer->cell_data_size_send = send_sizes;
    buffer->cell_data_size_send_rad = send_sizes_rad;
    for (axis = 0; axis < 3; ++axis) {
        buffer->cell_data_idx_send[axis] = idx_send[axis];
        buffer->cell_data_idx_send_rad[axis] = idx_send_rad[axis];
    }

    MPI_Sendrecv(&buffer->number, 1, MPI_INT, buffer->receiver_rank, 100,
        &send_entries, 1, MPI_INT, buffer->receiver_rank, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    buffer->cell_data_size_recv = (int *)calloc((size_t)send_entries + 1U,
        sizeof(*buffer->cell_data_size_recv));
    buffer->cell_data_size_recv_rad = (int *)calloc((size_t)send_entries + 1U,
        sizeof(*buffer->cell_data_size_recv_rad));
    if (buffer->cell_data_size_recv == 0 || buffer->cell_data_size_recv_rad == 0) {
        return 1;
    }
    MPI_Sendrecv(buffer->cell_data_size_send, buffer->number, MPI_INT, buffer->receiver_rank, 101,
        buffer->cell_data_size_recv, send_entries, MPI_INT, buffer->receiver_rank, 101, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    buffer->cell_data_size_recv[send_entries] = -1;
    MPI_Sendrecv(buffer->cell_data_size_send_rad, buffer->number, MPI_INT, buffer->receiver_rank, 102,
        buffer->cell_data_size_recv_rad, send_entries, MPI_INT, buffer->receiver_rank, 102, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    buffer->cell_data_size_recv_rad[send_entries] = -1;
    cell_size_total = prj_mpi_buffer_record_total(buffer->cell_data_size_recv, send_entries);
    cell_size_total_rad = prj_mpi_buffer_record_total(buffer->cell_data_size_recv_rad, send_entries);
    if (cell_size_total > 0) {
        for (axis = 0; axis < 3; ++axis) {
            buffer->cell_data_idx_recv[axis] = (int *)calloc(
                (size_t)cell_size_total,
                sizeof(*buffer->cell_data_idx_recv[axis]));
            if (buffer->cell_data_idx_recv[axis] == 0) {
                return 1;
            }
        }
    }
    if (cell_size_total_rad > 0) {
        for (axis = 0; axis < 3; ++axis) {
            buffer->cell_data_idx_recv_rad[axis] = (int *)calloc(
                (size_t)cell_size_total_rad,
                sizeof(*buffer->cell_data_idx_recv_rad[axis]));
            if (buffer->cell_data_idx_recv_rad[axis] == 0) {
                return 1;
            }
        }
    }
    for (axis = 0; axis < 3; ++axis) {
        MPI_Irecv(buffer->cell_data_idx_recv[axis], cell_size_total, MPI_INT, buffer->receiver_rank, 110 + axis,
            MPI_COMM_WORLD, &requests[req_n++]);
        MPI_Isend(buffer->cell_data_idx_send[axis], record_count, MPI_INT, buffer->receiver_rank, 110 + axis,
            MPI_COMM_WORLD, &requests[req_n++]);
        MPI_Irecv(buffer->cell_data_idx_recv_rad[axis], cell_size_total_rad, MPI_INT, buffer->receiver_rank, 113 + axis,
            MPI_COMM_WORLD, &requests[req_n++]);
        MPI_Isend(buffer->cell_data_idx_send_rad[axis], record_count_rad, MPI_INT, buffer->receiver_rank, 113 + axis,
            MPI_COMM_WORLD, &requests[req_n++]);
    }
    MPI_Waitall(req_n, requests, MPI_STATUSES_IGNORE);
    prj_mpi_count_recv_ghosts_by_kind(buffer, cell_size_total,
        cell_size_total_rad);
    for (fill_kind = 0; fill_kind < PRJ_MPI_GHOST_FILL_KIND_N; ++fill_kind) {
        if (prj_mpi_alloc_ghost_value_buffer(
                &buffer->cell_buffer_recv_by_kind[fill_kind],
                buffer->cell_recv_count_by_kind[fill_kind],
                buffer->cell_recv_count_rad_by_kind[fill_kind]) != 0) {
            return 1;
        }
    }
    return 0;
#else
    (void)mesh;
    (void)mpi;
    (void)buffer;
    return 0;
#endif
}

static double prj_mpi_read_cell_value(const double *src, int var,
    int i, int j, int k, int is_eosvar)
{
    return is_eosvar != 0 ? src[EIDX(var, i, j, k)] : src[VIDX(var, i, j, k)];
}

static double prj_mpi_prolongate_cell_value(const double *src, int var,
    int i, int j, int k, int is_eosvar, const double target[3], int use_BJ)
{
    double stencil[27];
    int di;
    int dj;
    int dk;

    for (di = -1; di <= 1; ++di) {
        for (dj = -1; dj <= 1; ++dj) {
            for (dk = -1; dk <= 1; ++dk) {
                stencil[prj_reconstruct_stencil3_index(di, dj, dk)] =
                    prj_mpi_read_cell_value(src, var, i + di, j + dj, k + dk, is_eosvar);
            }
        }
    }
    return prj_reconstruct_cell_for_prolongate(stencil, target, use_BJ);
}

static void prj_mpi_pack_ghost_values(prj_mesh *mesh, prj_mpi *mpi, prj_mpi_buffer *buffer,
    int stage, int fill_kind)
{
    double *send_buffer;
    int bidx;
    int pos;
    int pos_rad;
    int rad_offset;
    int record_count_total;

    if (mesh == 0 || mpi == 0 || buffer == 0 ||
        !prj_mpi_ghost_fill_kind_ok(fill_kind)) {
        return;
    }
    send_buffer = buffer->cell_buffer_send_by_kind[fill_kind];
    if (send_buffer == 0 &&
        (buffer->cell_send_count_by_kind[fill_kind] > 0 ||
         buffer->cell_send_count_rad_by_kind[fill_kind] > 0)) {
        return;
    }

    /* Radiation stream lives directly after the hydro+EOS stream in the same buffer. */
    record_count_total = buffer->cell_send_count_by_kind[fill_kind];
    rad_offset = record_count_total * PRJ_MPI_GHOST_NVAR_HE;

    pos = 0;
    pos_rad = rad_offset;
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        int n;

        if (!prj_block_is_active(block) || block->rank != mpi->rank) {
            continue;
        }
        double *W_send = stage == 2 ? block->W1 : block->W;
        double *eos_send = block->eosvar;
        for (n = 0; n < 56; ++n) {
            const prj_neighbor *slot = &block->slot[n];
            int nid = slot->id;
            int i;
            int j;
            int k;
            int v;

            int sample_kind = prj_mpi_ghost_fill_kind_from_rel_level(slot->rel_level);

            if (nid < 0 || nid >= mesh->nblocks || mesh->blocks[nid].rank != buffer->receiver_rank) {
                continue;
            }
            if (sample_kind != fill_kind) {
                continue;
            }

            /* Hydro primitives + EOS vars: iterate the full recv box. */
            for (i = 0; i < slot->recv_loc_end[0]-slot->recv_loc_start[0]; ++i) {
                for (j = 0; j < slot->recv_loc_end[1]-slot->recv_loc_start[1]; ++j) {
                    for (k = 0; k < slot->recv_loc_end[2]-slot->recv_loc_start[2]; ++k) {
                        if (slot->rel_level==0){
                            // Same level
                            for (v = 0; v < PRJ_NHYDRO; ++v) {
                                send_buffer[pos++] =
                                    W_send[VIDX(v, i+slot->send_loc_start[0], j+slot->send_loc_start[1], k+slot->send_loc_start[2])];
                            }
                            for (v = 0; v < PRJ_NVAR_EOSVAR; ++v) {
                                send_buffer[pos++] =
                                    eos_send[EIDX(v, i+slot->send_loc_start[0], j+slot->send_loc_start[1], k+slot->send_loc_start[2])];
                            }
                        } else if (slot->rel_level==-1) {
                            // Neighbor is coarser, restriction
                            for (v = 0; v < PRJ_NHYDRO; ++v) {
                                send_buffer[pos++] =
                                    0.125*
                                    (W_send[VIDX(v, 2*i+slot->send_loc_start[0],
                                                    2*j+slot->send_loc_start[1],
                                                    2*k+slot->send_loc_start[2])]
                                    +W_send[VIDX(v, 2*i+slot->send_loc_start[0],
                                                    2*j+slot->send_loc_start[1],
                                                    2*k+slot->send_loc_start[2]+1)]
                                    +W_send[VIDX(v, 2*i+slot->send_loc_start[0],
                                                    2*j+slot->send_loc_start[1]+1,
                                                    2*k+slot->send_loc_start[2])]
                                    +W_send[VIDX(v, 2*i+slot->send_loc_start[0],
                                                    2*j+slot->send_loc_start[1]+1,
                                                    2*k+slot->send_loc_start[2]+1)]
                                    +W_send[VIDX(v, 2*i+slot->send_loc_start[0]+1,
                                                    2*j+slot->send_loc_start[1],
                                                    2*k+slot->send_loc_start[2])]
                                    +W_send[VIDX(v, 2*i+slot->send_loc_start[0]+1,
                                                    2*j+slot->send_loc_start[1],
                                                    2*k+slot->send_loc_start[2]+1)]
                                    +W_send[VIDX(v, 2*i+slot->send_loc_start[0]+1,
                                                    2*j+slot->send_loc_start[1]+1,
                                                    2*k+slot->send_loc_start[2])]
                                    +W_send[VIDX(v, 2*i+slot->send_loc_start[0]+1,
                                                    2*j+slot->send_loc_start[1]+1,
                                                    2*k+slot->send_loc_start[2]+1)]);
                            }
                            for (v = 0; v < PRJ_NVAR_EOSVAR; ++v) {
                                send_buffer[pos++] =
                                    0.125*
                                    (eos_send[EIDX(v, 2*i+slot->send_loc_start[0],
                                                    2*j+slot->send_loc_start[1],
                                                    2*k+slot->send_loc_start[2])]
                                    +eos_send[EIDX(v, 2*i+slot->send_loc_start[0],
                                                    2*j+slot->send_loc_start[1],
                                                    2*k+slot->send_loc_start[2]+1)]
                                    +eos_send[EIDX(v, 2*i+slot->send_loc_start[0],
                                                    2*j+slot->send_loc_start[1]+1,
                                                    2*k+slot->send_loc_start[2])]
                                    +eos_send[EIDX(v, 2*i+slot->send_loc_start[0],
                                                    2*j+slot->send_loc_start[1]+1,
                                                    2*k+slot->send_loc_start[2]+1)]
                                    +eos_send[EIDX(v, 2*i+slot->send_loc_start[0]+1,
                                                    2*j+slot->send_loc_start[1],
                                                    2*k+slot->send_loc_start[2])]
                                    +eos_send[EIDX(v, 2*i+slot->send_loc_start[0]+1,
                                                    2*j+slot->send_loc_start[1],
                                                    2*k+slot->send_loc_start[2]+1)]
                                    +eos_send[EIDX(v, 2*i+slot->send_loc_start[0]+1,
                                                    2*j+slot->send_loc_start[1]+1,
                                                    2*k+slot->send_loc_start[2])]
                                    +eos_send[EIDX(v, 2*i+slot->send_loc_start[0]+1,
                                                    2*j+slot->send_loc_start[1]+1,
                                                    2*k+slot->send_loc_start[2]+1)]);
                            }
                        } else if (slot->rel_level==1) {
                            // Neighbor is finer, prolongation
                            double target[3];

                            target[0] = (i % 2 == 0) ? 0.25 : -0.25;
                            target[1] = (j % 2 == 0) ? 0.25 : -0.25;
                            target[2] = (k % 2 == 0) ? 0.25 : -0.25;
                            for (v = 0; v < PRJ_NHYDRO; ++v) {
                                send_buffer[pos++] =
                                    prj_mpi_prolongate_cell_value(W_send, v,
                                        i/2+slot->send_loc_start[0],
                                        j/2+slot->send_loc_start[1],
                                        k/2+slot->send_loc_start[2],
                                        0, target, mesh->use_BJ);
                            }
                            for (v = 0; v < PRJ_NVAR_EOSVAR; ++v) {
                                send_buffer[pos++] =
                                    prj_mpi_prolongate_cell_value(eos_send, v,
                                        i/2+slot->send_loc_start[0],
                                        j/2+slot->send_loc_start[1],
                                        k/2+slot->send_loc_start[2],
                                        1, target, mesh->use_BJ);
                            }
                        } else {
                          fprintf(stderr,"slot->rel_level unrecognized: %d\n", slot->rel_level);
                          exit(1);
                        }
                    }
                }
            }

#if PRJ_NRAD > 0
            /* Radiation primitives: narrower rad-clipped recv box. */
            for (i = 0; i < slot->recv_loc_end_rad[0]-slot->recv_loc_start_rad[0]; ++i) {
                for (j = 0; j < slot->recv_loc_end_rad[1]-slot->recv_loc_start_rad[1]; ++j) {
                    for (k = 0; k < slot->recv_loc_end_rad[2]-slot->recv_loc_start_rad[2]; ++k) {
                        if (sample_kind != fill_kind) {
                            for (v = PRJ_NHYDRO; v < PRJ_NVAR_PRIM; ++v) {
                                send_buffer[pos_rad++] = 0.0;
                            }
                            continue;
                        }
                        if (slot->rel_level==0){
                            for (v = PRJ_NHYDRO; v < PRJ_NVAR_PRIM; ++v) {
                                send_buffer[pos_rad++] =
                                    W_send[VIDX(v, i+slot->send_loc_start_rad[0], j+slot->send_loc_start_rad[1], k+slot->send_loc_start_rad[2])];
                            }
                        } else if (slot->rel_level==-1) {
                            for (v = PRJ_NHYDRO; v < PRJ_NVAR_PRIM; ++v) {
                                send_buffer[pos_rad++] =
                                    0.125*
                                    (W_send[VIDX(v, 2*i+slot->send_loc_start_rad[0],
                                                    2*j+slot->send_loc_start_rad[1],
                                                    2*k+slot->send_loc_start_rad[2])]
                                    +W_send[VIDX(v, 2*i+slot->send_loc_start_rad[0],
                                                    2*j+slot->send_loc_start_rad[1],
                                                    2*k+slot->send_loc_start_rad[2]+1)]
                                    +W_send[VIDX(v, 2*i+slot->send_loc_start_rad[0],
                                                    2*j+slot->send_loc_start_rad[1]+1,
                                                    2*k+slot->send_loc_start_rad[2])]
                                    +W_send[VIDX(v, 2*i+slot->send_loc_start_rad[0],
                                                    2*j+slot->send_loc_start_rad[1]+1,
                                                    2*k+slot->send_loc_start_rad[2]+1)]
                                    +W_send[VIDX(v, 2*i+slot->send_loc_start_rad[0]+1,
                                                    2*j+slot->send_loc_start_rad[1],
                                                    2*k+slot->send_loc_start_rad[2])]
                                    +W_send[VIDX(v, 2*i+slot->send_loc_start_rad[0]+1,
                                                    2*j+slot->send_loc_start_rad[1],
                                                    2*k+slot->send_loc_start_rad[2]+1)]
                                    +W_send[VIDX(v, 2*i+slot->send_loc_start_rad[0]+1,
                                                    2*j+slot->send_loc_start_rad[1]+1,
                                                    2*k+slot->send_loc_start_rad[2])]
                                    +W_send[VIDX(v, 2*i+slot->send_loc_start_rad[0]+1,
                                                    2*j+slot->send_loc_start_rad[1]+1,
                                                    2*k+slot->send_loc_start_rad[2]+1)]);
                            }
                        } else if (slot->rel_level==1) {
                            double target[3];
                            int ai = i + slot->recv_loc_start_rad[0];
                            int aj = j + slot->recv_loc_start_rad[1];
                            int ak = k + slot->recv_loc_start_rad[2];

                            target[0] = (ai % 2 == 0) ? 0.25 : -0.25;
                            target[1] = (aj % 2 == 0) ? 0.25 : -0.25;
                            target[2] = (ak % 2 == 0) ? 0.25 : -0.25;
                            for (v = PRJ_NHYDRO; v < PRJ_NVAR_PRIM; ++v) {
                                send_buffer[pos_rad++] =
                                    prj_mpi_prolongate_cell_value(W_send, v,
                                        i/2+slot->send_loc_start_rad[0],
                                        j/2+slot->send_loc_start_rad[1],
                                        k/2+slot->send_loc_start_rad[2],
                                        0, target, mesh->use_BJ);
                            }
                        }
                    }
                }
            }
#endif
        }
    }

    (void)pos;
    (void)pos_rad;
}

static void prj_mpi_unpack_ghost_values(prj_mesh *mesh, prj_mpi *mpi,
    prj_mpi_buffer *buffer, int stage, int fill_kind)
{
    double *recv_buffer;
    int total;
    int total_rad;
    int pos;
    int pos_rad;
    int i;

    if (mesh == 0 || mpi == 0 || buffer == 0 ||
        !prj_mpi_ghost_fill_kind_ok(fill_kind)) {
        return;
    }
    recv_buffer = buffer->cell_buffer_recv_by_kind[fill_kind];
    if (recv_buffer == 0 &&
        (buffer->cell_recv_count_by_kind[fill_kind] > 0 ||
         buffer->cell_recv_count_rad_by_kind[fill_kind] > 0)) {
        return;
    }
    total = prj_mpi_buffer_recv_count(buffer);
    total = prj_mpi_buffer_record_total(buffer->cell_data_size_recv, total);
    total_rad = prj_mpi_buffer_recv_count(buffer);
    total_rad = prj_mpi_buffer_record_total(buffer->cell_data_size_recv_rad,
        total_rad);
    pos = 0;
    pos_rad = buffer->cell_recv_count_by_kind[fill_kind] *
        PRJ_MPI_GHOST_NVAR_HE;

    for (i = 0; i < total; ++i) {
        int block_id = buffer->cell_data_idx_recv[0][i];
        int code = buffer->cell_data_idx_recv[1][i];
        int sample_kind = prj_mpi_sample_kind_from_code(buffer->cell_data_idx_recv[2][i]);
        int ii;
        int jj;
        int kk;
        int v;
        prj_block *block;
        double *dst;

        if (sample_kind != fill_kind) {
            continue;
        }
        if (block_id < 0 || block_id >= mesh->nblocks) {
            pos += PRJ_MPI_GHOST_NVAR_HE;
            continue;
        }
        block = &mesh->blocks[block_id];
        if (block->rank != mpi->rank || block->W == 0 || block->eosvar == 0) {
            pos += PRJ_MPI_GHOST_NVAR_HE;
            continue;
        }
        dst = prj_mpi_stage_array(block, stage);
        prj_mpi_decode_cell_index(code, &ii, &jj, &kk);
        for (v = 0; v < PRJ_NHYDRO; ++v) {
            dst[VIDX(v, ii, jj, kk)] = recv_buffer[pos++];
        }
        for (v = 0; v < PRJ_NVAR_EOSVAR; ++v) {
            block->eosvar[EIDX(v, ii, jj, kk)] = recv_buffer[pos++];
        }
    }
#if PRJ_NRAD > 0
    for (i = 0; i < total_rad; ++i) {
        int block_id = buffer->cell_data_idx_recv_rad[0][i];
        int code = buffer->cell_data_idx_recv_rad[1][i];
        int sample_kind = prj_mpi_sample_kind_from_code(buffer->cell_data_idx_recv_rad[2][i]);
        int ii;
        int jj;
        int kk;
        int v;
        prj_block *block;
        double *dst;

        if (sample_kind != fill_kind) {
            continue;
        }
        if (block_id < 0 || block_id >= mesh->nblocks) {
            pos_rad += PRJ_MPI_GHOST_NVAR_RAD;
            continue;
        }
        block = &mesh->blocks[block_id];
        if (block->rank != mpi->rank || block->W == 0) {
            pos_rad += PRJ_MPI_GHOST_NVAR_RAD;
            continue;
        }
        dst = prj_mpi_stage_array(block, stage);
        prj_mpi_decode_cell_index(code, &ii, &jj, &kk);
        for (v = PRJ_NHYDRO; v < PRJ_NVAR_PRIM; ++v) {
            dst[VIDX(v, ii, jj, kk)] = recv_buffer[pos_rad++];
        }
    }
#else
    (void)pos_rad;
    (void)total_rad;
#endif
}

#if PRJ_MHD
enum {
    PRJ_MPI_BF_HDR_DST_ID = 0,
    PRJ_MPI_BF_HDR_REL_LEVEL = 1,
    PRJ_MPI_BF_HDR_VALUE_COUNT = 2,
    PRJ_MPI_BF_HDR_RECV_START0 = 3,
    PRJ_MPI_BF_HDR_RECV_START1 = 4,
    PRJ_MPI_BF_HDR_RECV_START2 = 5,
    PRJ_MPI_BF_HDR_RECV_END0 = 6,
    PRJ_MPI_BF_HDR_RECV_END1 = 7,
    PRJ_MPI_BF_HDR_RECV_END2 = 8,
    PRJ_MPI_BF_HDR_SEND_START0 = 9,
    PRJ_MPI_BF_HDR_SEND_START1 = 10,
    PRJ_MPI_BF_HDR_SEND_START2 = 11,
    PRJ_MPI_BF_HDR_SEND_END0 = 12,
    PRJ_MPI_BF_HDR_SEND_END1 = 13,
    PRJ_MPI_BF_HDR_SEND_END2 = 14,
    PRJ_MPI_BF_HEADER_NINT = 15
};

#define PRJ_MPI_BF_FILL_N 6

static double *prj_mpi_bf_array(prj_block *block, int dir, int use_bf1)
{
    return use_bf1 != 0 ? block->Bf1[dir] : block->Bf[dir];
}

static const double *prj_mpi_bf_array_const(const prj_block *block, int dir, int use_bf1)
{
    return use_bf1 != 0 ? block->Bf1[dir] : block->Bf[dir];
}

static void prj_mpi_check_bf_storage(const prj_block *block, const char *label)
{
    int d;

    if (block == 0) {
        fprintf(stderr, "%s: missing block\n", label);
        exit(EXIT_FAILURE);
    }
    for (d = 0; d < 3; ++d) {
        if (block->face_fidelity[d] == 0) {
            fprintf(stderr, "%s: missing face fidelity storage\n", label);
            exit(EXIT_FAILURE);
        }
        if (block->Bf[d] == 0 || block->Bf1[d] == 0) {
            fprintf(stderr, "%s: missing face-centered magnetic field storage\n", label);
            exit(EXIT_FAILURE);
        }
    }
}

static int prj_mpi_bf_axis_active_max(int dir, int axis)
{
    return dir == axis ? PRJ_BLOCK_SIZE : PRJ_BLOCK_SIZE - 1;
}

static int prj_mpi_bf_face_active(int dir, int i, int j, int k)
{
    int idx[3];
    int d;

    idx[0] = i;
    idx[1] = j;
    idx[2] = k;
    for (d = 0; d < 3; ++d) {
        if (idx[d] < 0 || idx[d] > prj_mpi_bf_axis_active_max(dir, d)) {
            return 0;
        }
    }
    return 1;
}

static int prj_mpi_bf_slot_matches(int fill_kind, int rel_level)
{
    switch (fill_kind) {
    case PRJ_BOUNDARY_FILL_SAME_LEVEL:
        return rel_level == 0;
    case PRJ_BOUNDARY_FILL_RESTRICTION:
        return rel_level < 0;
    case PRJ_BOUNDARY_FILL_PROLONGATION:
        return rel_level > 0;
    case PRJ_BOUNDARY_FILL_NONRECON:
        return rel_level <= 0;
    case PRJ_BOUNDARY_FILL_RECON:
        return rel_level > 0;
    case PRJ_BOUNDARY_FILL_ALL:
        return 1;
    default:
        return 0;
    }
}

static int prj_mpi_bf_append_value(double *values, int *count, int capacity,
    double value)
{
    if (!isfinite(value)) {
        return 1;
    }
    if (*count >= capacity) {
        return 1;
    }
    values[*count] = value;
    *count += 1;
    return 0;
}

static int prj_mpi_bf_append_record(int *headers, int *record_count,
    int record_capacity, const prj_neighbor *slot, int value_count)
{
    int offset;
    int d;

    if (*record_count >= record_capacity) {
        return 1;
    }
    offset = (*record_count) * PRJ_MPI_BF_HEADER_NINT;
    headers[offset + PRJ_MPI_BF_HDR_DST_ID] = slot->id;
    headers[offset + PRJ_MPI_BF_HDR_REL_LEVEL] = slot->rel_level;
    headers[offset + PRJ_MPI_BF_HDR_VALUE_COUNT] = value_count;
    for (d = 0; d < 3; ++d) {
        headers[offset + PRJ_MPI_BF_HDR_RECV_START0 + d] = slot->recv_loc_start[d];
        headers[offset + PRJ_MPI_BF_HDR_RECV_END0 + d] = slot->recv_loc_end[d];
        headers[offset + PRJ_MPI_BF_HDR_SEND_START0 + d] = slot->send_loc_start[d];
        headers[offset + PRJ_MPI_BF_HDR_SEND_END0 + d] = slot->send_loc_end[d];
    }
    *record_count += 1;
    return 0;
}

static void prj_mpi_bf_header_to_slot(const int *header, prj_neighbor *slot)
{
    int d;

    memset(slot, 0, sizeof(*slot));
    slot->id = header[PRJ_MPI_BF_HDR_DST_ID];
    slot->rel_level = header[PRJ_MPI_BF_HDR_REL_LEVEL];
    for (d = 0; d < 3; ++d) {
        slot->recv_loc_start[d] = header[PRJ_MPI_BF_HDR_RECV_START0 + d];
        slot->recv_loc_end[d] = header[PRJ_MPI_BF_HDR_RECV_END0 + d];
        slot->send_loc_start[d] = header[PRJ_MPI_BF_HDR_SEND_START0 + d];
        slot->send_loc_end[d] = header[PRJ_MPI_BF_HDR_SEND_END0 + d];
    }
}

static int prj_mpi_bf_same_level_value_count(const prj_neighbor *slot)
{
    int dir;
    int count = 0;

    for (dir = 0; dir < 3; ++dir) {
        int i;
        int j;
        int k;

        for (i = 0; i < slot->recv_loc_end[0] - slot->recv_loc_start[0] + (dir == 0 ? 1 : 0); ++i) {
            for (j = 0; j < slot->recv_loc_end[1] - slot->recv_loc_start[1] + (dir == 1 ? 1 : 0); ++j) {
                for (k = 0; k < slot->recv_loc_end[2] - slot->recv_loc_start[2] + (dir == 2 ? 1 : 0); ++k) {
                    if (prj_mpi_bf_face_active(dir,
                            i + slot->recv_loc_start[0],
                            j + slot->recv_loc_start[1],
                            k + slot->recv_loc_start[2])) {
                        continue;
                    }
                    count += 1;
                }
            }
        }
    }
    return count;
}

static int prj_mpi_bf_restriction_value_count(const prj_neighbor *slot)
{
    int dir;
    int count = 0;

    for (dir = 0; dir < 3; ++dir) {
        int di = slot->send_loc_end[0] - slot->send_loc_start[0] + (dir == 0 ? 1 : 0) + 1;
        int dj = slot->send_loc_end[1] - slot->send_loc_start[1] + (dir == 1 ? 1 : 0) + 1;
        int dk = slot->send_loc_end[2] - slot->send_loc_start[2] + (dir == 2 ? 1 : 0) + 1;
        count += di * dj * dk;
    }
    return count;
}

static int prj_mpi_bf_prolongation_value_count(const prj_neighbor *slot)
{
    int dir;
    int count = 0;

    for (dir = 0; dir < 3; ++dir) {
        int d;
        int n[3];

        for (d = 0; d < 3; ++d) {
            if (d == dir) {
                n[d] = slot->send_loc_end[d] - slot->send_loc_start[d] + 2;
            } else {
                n[d] = slot->send_loc_end[d] - slot->send_loc_start[d] + 3;
            }
        }
        count += n[0] * n[1] * n[2];
    }
    return count;
}

static int prj_mpi_bf_record_value_count(int fill_kind, const prj_neighbor *slot)
{
    if (!prj_mpi_bf_slot_matches(fill_kind, slot->rel_level)) {
        return 0;
    }
    if (slot->rel_level == 0) {
        return prj_mpi_bf_same_level_value_count(slot);
    }
    if (slot->rel_level < 0) {
        return prj_mpi_bf_restriction_value_count(slot);
    }
    return prj_mpi_bf_prolongation_value_count(slot);
}

static int prj_mpi_pack_bf_same_level_values(const prj_block *block,
    int use_bf1, const prj_neighbor *slot, double *values, int *value_count,
    int value_capacity)
{
    int dir;

    prj_mpi_check_bf_storage(block, "prj_mpi_pack_bf_same_level_values");
    for (dir = 0; dir < 3; ++dir) {
        const double *src = prj_mpi_bf_array_const(block, dir, use_bf1);
        int i;
        int j;
        int k;

        for (i = 0; i < slot->recv_loc_end[0] - slot->recv_loc_start[0] + (dir == 0 ? 1 : 0); ++i) {
            for (j = 0; j < slot->recv_loc_end[1] - slot->recv_loc_start[1] + (dir == 1 ? 1 : 0); ++j) {
                for (k = 0; k < slot->recv_loc_end[2] - slot->recv_loc_start[2] + (dir == 2 ? 1 : 0); ++k) {
                    if (prj_mpi_bf_face_active(dir,
                            i + slot->recv_loc_start[0],
                            j + slot->recv_loc_start[1],
                            k + slot->recv_loc_start[2])) {
                        continue;
                    }
                    if (prj_mpi_bf_append_value(values, value_count, value_capacity,
                            src[FACE_IDX(dir,
                                i + slot->send_loc_start[0],
                                j + slot->send_loc_start[1],
                                k + slot->send_loc_start[2])]) != 0) {
                        return 1;
                    }
                }
            }
        }
    }
    return 0;
}

static int prj_mpi_pack_bf_restriction_values(const prj_block *block,
    int use_bf1, const prj_neighbor *slot, double *values, int *value_count,
    int value_capacity)
{
    int dir;

    prj_mpi_check_bf_storage(block, "prj_mpi_pack_bf_restriction_values");
    for (dir = 0; dir < 3; ++dir) {
        const double *src = prj_mpi_bf_array_const(block, dir, use_bf1);
        int si;
        int sj;
        int sk;

        for (si = slot->send_loc_start[0]; si <= slot->send_loc_end[0] + (dir == 0 ? 1 : 0); ++si) {
            for (sj = slot->send_loc_start[1]; sj <= slot->send_loc_end[1] + (dir == 1 ? 1 : 0); ++sj) {
                for (sk = slot->send_loc_start[2]; sk <= slot->send_loc_end[2] + (dir == 2 ? 1 : 0); ++sk) {
                    if (prj_mpi_bf_append_value(values, value_count, value_capacity,
                            src[FACE_IDX(dir, si, sj, sk)]) != 0) {
                        return 1;
                    }
                }
            }
        }
    }
    return 0;
}

static int prj_mpi_pack_bf_prolongation_values(const prj_block *block,
    int use_bf1, const prj_neighbor *slot, double *values, int *value_count,
    int value_capacity)
{
    int dir;

    prj_mpi_check_bf_storage(block, "prj_mpi_pack_bf_prolongation_values");
    for (dir = 0; dir < 3; ++dir) {
        const double *src = prj_mpi_bf_array_const(block, dir, use_bf1);
        int lo[3];
        int hi[3];
        int d;
        int si;
        int sj;
        int sk;

        for (d = 0; d < 3; ++d) {
            if (d == dir) {
                lo[d] = slot->send_loc_start[d];
                hi[d] = slot->send_loc_end[d] + 1;
            } else {
                lo[d] = slot->send_loc_start[d] - 1;
                hi[d] = slot->send_loc_end[d] + 1;
            }
        }
        for (si = lo[0]; si <= hi[0]; ++si) {
            for (sj = lo[1]; sj <= hi[1]; ++sj) {
                for (sk = lo[2]; sk <= hi[2]; ++sk) {
                    if (prj_mpi_bf_append_value(values, value_count, value_capacity,
                            src[FACE_IDX(dir, si, sj, sk)]) != 0) {
                        return 1;
                    }
                }
            }
        }
    }
    return 0;
}

static int prj_mpi_pack_bf_record_values(const prj_block *block, int use_bf1,
    const prj_neighbor *slot, double *values, int *value_count,
    int value_capacity)
{
    if (slot->rel_level == 0) {
        return prj_mpi_pack_bf_same_level_values(block, use_bf1, slot,
            values, value_count, value_capacity);
    }
    if (slot->rel_level < 0) {
        return prj_mpi_pack_bf_restriction_values(block, use_bf1, slot,
            values, value_count, value_capacity);
    }
    return prj_mpi_pack_bf_prolongation_values(block, use_bf1, slot,
        values, value_count, value_capacity);
}

static int prj_mpi_pack_bf_records(const prj_mesh *mesh, const prj_mpi *mpi,
    int receiver_rank, int use_bf1, int fill_kind, int *headers,
    int record_capacity, int *record_count, double *values,
    int value_capacity, int *value_count)
{
    int bidx;

    *record_count = 0;
    *value_count = 0;
    if (mesh == 0 || mpi == 0) {
        return 1;
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        const prj_block *block = &mesh->blocks[bidx];
        int n;

        if (!prj_mpi_block_is_local(mpi, block)) {
            continue;
        }
        for (n = 0; n < 56; ++n) {
            const prj_neighbor *slot = &block->slot[n];
            int nid = slot->id;
            int before;
            int expected;

            if (nid < 0 || nid >= mesh->nblocks || mesh->blocks[nid].rank != receiver_rank ||
                !prj_mpi_bf_slot_matches(fill_kind, slot->rel_level)) {
                continue;
            }
            expected = prj_mpi_bf_record_value_count(fill_kind, slot);
            if (expected <= 0) {
                continue;
            }
            if (prj_mpi_bf_append_record(headers, record_count, record_capacity,
                    slot, expected) != 0) {
                return 1;
            }
            before = *value_count;
            if (prj_mpi_pack_bf_record_values(block, use_bf1, slot,
                    values, value_count, value_capacity) != 0 ||
                *value_count - before != expected) {
                return 1;
            }
        }
    }
    return 0;
}

static int prj_mpi_pack_bf_values_only(const prj_mesh *mesh, const prj_mpi *mpi,
    int receiver_rank, int use_bf1, int fill_kind, double *values,
    int value_capacity, int *value_count)
{
    int bidx;

    *value_count = 0;
    if (mesh == 0 || mpi == 0) {
        return 1;
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        const prj_block *block = &mesh->blocks[bidx];
        int n;

        if (!prj_mpi_block_is_local(mpi, block)) {
            continue;
        }
        for (n = 0; n < 56; ++n) {
            const prj_neighbor *slot = &block->slot[n];
            int nid = slot->id;
            int expected;
            int before;

            if (nid < 0 || nid >= mesh->nblocks ||
                mesh->blocks[nid].rank != receiver_rank ||
                !prj_mpi_bf_slot_matches(fill_kind, slot->rel_level)) {
                continue;
            }
            expected = prj_mpi_bf_record_value_count(fill_kind, slot);
            if (expected <= 0) {
                continue;
            }
            before = *value_count;
            if (prj_mpi_pack_bf_record_values(block, use_bf1, slot,
                    values, value_count, value_capacity) != 0 ||
                *value_count - before != expected) {
                return 1;
            }
        }
    }
    return 0;
}

static inline double prj_mpi_bf_buf_read(const double *buf,
    const int buf_lo[3], const int buf_n[3], int i, int j, int k)
{
    return buf[((i - buf_lo[0]) * buf_n[1] + (j - buf_lo[1])) * buf_n[2] +
               (k - buf_lo[2])];
}

static inline double prj_mpi_bf_buf_sign_half(int bit)
{
    return bit == 0 ? -1.0 : 1.0;
}

static inline double prj_mpi_bf_interp_x1_buf(const double *buf,
    const int buf_lo[3], const int buf_n[3], const double dx[3],
    int i, int j, int k, int fine_j, int fine_k, double area, int use_BJ)
{
    double stencil[9];
    double target[2];
    double value;
    int dj;
    int dk;

    (void)dx;
    for (dj = -1; dj <= 1; ++dj) {
        for (dk = -1; dk <= 1; ++dk) {
            stencil[prj_reconstruct_stencil2_index(dj, dk)] =
                prj_mpi_bf_buf_read(buf, buf_lo, buf_n, i, j + dj, k + dk);
        }
    }
    target[0] = 0.25 * prj_mpi_bf_buf_sign_half(fine_j);
    target[1] = 0.25 * prj_mpi_bf_buf_sign_half(fine_k);
    value = prj_reconstruct_face_for_prolongate(stencil, target, use_BJ);

    return value * area;
}

static inline double prj_mpi_bf_interp_x2_buf(const double *buf,
    const int buf_lo[3], const int buf_n[3], const double dx[3],
    int i, int j, int k, int fine_i, int fine_k, double area, int use_BJ)
{
    double stencil[9];
    double target[2];
    double value;
    int di;
    int dk;

    (void)dx;
    for (di = -1; di <= 1; ++di) {
        for (dk = -1; dk <= 1; ++dk) {
            stencil[prj_reconstruct_stencil2_index(di, dk)] =
                prj_mpi_bf_buf_read(buf, buf_lo, buf_n, i + di, j, k + dk);
        }
    }
    target[0] = 0.25 * prj_mpi_bf_buf_sign_half(fine_i);
    target[1] = 0.25 * prj_mpi_bf_buf_sign_half(fine_k);
    value = prj_reconstruct_face_for_prolongate(stencil, target, use_BJ);

    return value * area;
}

static inline double prj_mpi_bf_interp_x3_buf(const double *buf,
    const int buf_lo[3], const int buf_n[3], const double dx[3],
    int i, int j, int k, int fine_i, int fine_j, double area, int use_BJ)
{
    double stencil[9];
    double target[2];
    double value;
    int di;
    int dj;

    (void)dx;
    for (di = -1; di <= 1; ++di) {
        for (dj = -1; dj <= 1; ++dj) {
            stencil[prj_reconstruct_stencil2_index(di, dj)] =
                prj_mpi_bf_buf_read(buf, buf_lo, buf_n, i + di, j + dj, k);
        }
    }
    target[0] = 0.25 * prj_mpi_bf_buf_sign_half(fine_i);
    target[1] = 0.25 * prj_mpi_bf_buf_sign_half(fine_j);
    value = prj_reconstruct_face_for_prolongate(stencil, target, use_BJ);

    return value * area;
}

static inline void prj_mpi_bf_compute_inner_fluxes(double u[3][2][2],
    double v[2][3][2], double w[2][2][3], double dx1, double dx2,
    double dx3)
{
    double Uxx = 0.0;
    double Vyy = 0.0;
    double Wzz = 0.0;
    double Uxyz = 0.0;
    double Vxyz = 0.0;
    double Wxyz = 0.0;
    double dx1s = dx1 * dx1;
    double dx2s = dx2 * dx2;
    double dx3s = dx3 * dx3;
    int i;
    int j;
    int k;

    for (i = 0; i < 2; ++i) {
        double si = prj_mpi_bf_buf_sign_half(i);
        for (j = 0; j < 2; ++j) {
            double sj = prj_mpi_bf_buf_sign_half(j);
            for (k = 0; k < 2; ++k) {
                double sk = prj_mpi_bf_buf_sign_half(k);
                Uxx += si * sj * v[i][2 * j][k] + si * sk * w[i][j][2 * k];
                Vyy += sj * sk * w[i][j][2 * k] + sj * si * u[2 * i][j][k];
                Wzz += sk * si * u[2 * i][j][k] + sk * sj * v[i][2 * j][k];
                Uxyz += si * sj * sk * u[2 * i][j][k] / (dx2s + dx3s);
                Vxyz += si * sj * sk * v[i][2 * j][k] / (dx1s + dx3s);
                Wxyz += si * sj * sk * w[i][j][2 * k] / (dx1s + dx2s);
            }
        }
    }
    Uxx *= 0.125;
    Vyy *= 0.125;
    Wzz *= 0.125;
    Uxyz *= 0.125;
    Vxyz *= 0.125;
    Wxyz *= 0.125;

    for (j = 0; j < 2; ++j) {
        double sj = prj_mpi_bf_buf_sign_half(j);
        for (k = 0; k < 2; ++k) {
            double sk = prj_mpi_bf_buf_sign_half(k);
            u[1][j][k] = 0.5 * (u[2][j][k] + u[0][j][k]) +
                Uxx + sk * dx3s * Vxyz + sj * dx2s * Wxyz;
        }
    }
    for (i = 0; i < 2; ++i) {
        double si = prj_mpi_bf_buf_sign_half(i);
        for (k = 0; k < 2; ++k) {
            double sk = prj_mpi_bf_buf_sign_half(k);
            v[i][1][k] = 0.5 * (v[i][2][k] + v[i][0][k]) +
                Vyy + si * dx1s * Wxyz + sk * dx3s * Uxyz;
        }
    }
    for (i = 0; i < 2; ++i) {
        double si = prj_mpi_bf_buf_sign_half(i);
        for (j = 0; j < 2; ++j) {
            double sj = prj_mpi_bf_buf_sign_half(j);
            w[i][j][1] = 0.5 * (w[i][j][2] + w[i][j][0]) +
                Wzz + sj * dx2s * Uxyz + si * dx1s * Vxyz;
        }
    }
}

static void prj_mpi_bf_prolong_from_buffer(const double *buf[3],
    int buf_lo[3][3], int buf_n[3][3], const double coarse_dx[3],
    prj_block *fine, int ci, int cj, int ck, int fi, int fj, int fk,
    int use_bf1, int use_BJ)
{
    double u[3][2][2];
    double v[2][3][2];
    double w[2][2][3];
    double *dst[3];
    double area_u;
    double area_v;
    double area_w;
    int i;
    int j;
    int k;
    int d;

    prj_mpi_check_bf_storage(fine, "prj_mpi_bf_prolong_from_buffer");
    for (d = 0; d < 3; ++d) {
        dst[d] = prj_mpi_bf_array(fine, d, use_bf1);
    }
    area_u = fine->dx[1] * fine->dx[2];
    area_v = fine->dx[0] * fine->dx[2];
    area_w = fine->dx[0] * fine->dx[1];

    for (j = 0; j < 2; ++j) {
        for (k = 0; k < 2; ++k) {
            double flux;
            int idx;

            flux = prj_mpi_bf_interp_x1_buf(buf[X1DIR], buf_lo[X1DIR],
                buf_n[X1DIR], coarse_dx, ci, cj, ck, j, k, area_u, use_BJ);
            idx = FACE_IDX(X1DIR, fi, fj + j, fk + k);
            if (fine->face_fidelity[X1DIR][idx] > PRJ_MHD_FIDELITY_COARSER) {
                u[0][j][k] = dst[X1DIR][idx] * area_u;
            } else {
                u[0][j][k] = flux;
                prj_boundary_write_bf_face(fine, use_bf1, X1DIR,
                    fi, fj + j, fk + k, flux / area_u, PRJ_MHD_FIDELITY_COARSER);
            }

            flux = prj_mpi_bf_interp_x1_buf(buf[X1DIR], buf_lo[X1DIR],
                buf_n[X1DIR], coarse_dx, ci + 1, cj, ck, j, k, area_u, use_BJ);
            idx = FACE_IDX(X1DIR, fi + 2, fj + j, fk + k);
            if (fine->face_fidelity[X1DIR][idx] > PRJ_MHD_FIDELITY_COARSER) {
                u[2][j][k] = dst[X1DIR][idx] * area_u;
            } else {
                u[2][j][k] = flux;
                prj_boundary_write_bf_face(fine, use_bf1, X1DIR,
                    fi + 2, fj + j, fk + k, flux / area_u, PRJ_MHD_FIDELITY_COARSER);
            }
        }
    }

    for (i = 0; i < 2; ++i) {
        for (k = 0; k < 2; ++k) {
            double flux;
            int idx;

            flux = prj_mpi_bf_interp_x2_buf(buf[X2DIR], buf_lo[X2DIR],
                buf_n[X2DIR], coarse_dx, ci, cj, ck, i, k, area_v, use_BJ);
            idx = FACE_IDX(X2DIR, fi + i, fj, fk + k);
            if (fine->face_fidelity[X2DIR][idx] > PRJ_MHD_FIDELITY_COARSER) {
                v[i][0][k] = dst[X2DIR][idx] * area_v;
            } else {
                v[i][0][k] = flux;
                prj_boundary_write_bf_face(fine, use_bf1, X2DIR,
                    fi + i, fj, fk + k, flux / area_v, PRJ_MHD_FIDELITY_COARSER);
            }

            flux = prj_mpi_bf_interp_x2_buf(buf[X2DIR], buf_lo[X2DIR],
                buf_n[X2DIR], coarse_dx, ci, cj + 1, ck, i, k, area_v, use_BJ);
            idx = FACE_IDX(X2DIR, fi + i, fj + 2, fk + k);
            if (fine->face_fidelity[X2DIR][idx] > PRJ_MHD_FIDELITY_COARSER) {
                v[i][2][k] = dst[X2DIR][idx] * area_v;
            } else {
                v[i][2][k] = flux;
                prj_boundary_write_bf_face(fine, use_bf1, X2DIR,
                    fi + i, fj + 2, fk + k, flux / area_v, PRJ_MHD_FIDELITY_COARSER);
            }
        }
    }

    for (i = 0; i < 2; ++i) {
        for (j = 0; j < 2; ++j) {
            double flux;
            int idx;

            flux = prj_mpi_bf_interp_x3_buf(buf[X3DIR], buf_lo[X3DIR],
                buf_n[X3DIR], coarse_dx, ci, cj, ck, i, j, area_w, use_BJ);
            idx = FACE_IDX(X3DIR, fi + i, fj + j, fk);
            if (fine->face_fidelity[X3DIR][idx] > PRJ_MHD_FIDELITY_COARSER) {
                w[i][j][0] = dst[X3DIR][idx] * area_w;
            } else {
                w[i][j][0] = flux;
                prj_boundary_write_bf_face(fine, use_bf1, X3DIR,
                    fi + i, fj + j, fk, flux / area_w, PRJ_MHD_FIDELITY_COARSER);
            }

            flux = prj_mpi_bf_interp_x3_buf(buf[X3DIR], buf_lo[X3DIR],
                buf_n[X3DIR], coarse_dx, ci, cj, ck + 1, i, j, area_w, use_BJ);
            idx = FACE_IDX(X3DIR, fi + i, fj + j, fk + 2);
            if (fine->face_fidelity[X3DIR][idx] > PRJ_MHD_FIDELITY_COARSER) {
                w[i][j][2] = dst[X3DIR][idx] * area_w;
            } else {
                w[i][j][2] = flux;
                prj_boundary_write_bf_face(fine, use_bf1, X3DIR,
                    fi + i, fj + j, fk + 2, flux / area_w, PRJ_MHD_FIDELITY_COARSER);
            }
        }
    }

    prj_mpi_bf_compute_inner_fluxes(u, v, w,
        fine->dx[0], fine->dx[1], fine->dx[2]);

    for (j = 0; j < 2; ++j) {
        for (k = 0; k < 2; ++k) {
            prj_boundary_write_bf_face(fine, use_bf1, X1DIR,
                fi + 1, fj + j, fk + k,
                u[1][j][k] / area_u, PRJ_MHD_FIDELITY_COARSER);
        }
    }
    for (i = 0; i < 2; ++i) {
        for (k = 0; k < 2; ++k) {
            prj_boundary_write_bf_face(fine, use_bf1, X2DIR,
                fi + i, fj + 1, fk + k,
                v[i][1][k] / area_v, PRJ_MHD_FIDELITY_COARSER);
        }
    }
    for (i = 0; i < 2; ++i) {
        for (j = 0; j < 2; ++j) {
            prj_boundary_write_bf_face(fine, use_bf1, X3DIR,
                fi + i, fj + j, fk + 1,
                w[i][j][1] / area_w, PRJ_MHD_FIDELITY_COARSER);
        }
    }
}

static int prj_mpi_apply_bf_same_level_record(prj_block *dst, int use_bf1,
    const prj_neighbor *slot, const double *values, int value_count)
{
    int dir;
    int pos = 0;

    for (dir = 0; dir < 3; ++dir) {
        int i;
        int j;
        int k;

        for (i = 0; i < slot->recv_loc_end[0] - slot->recv_loc_start[0] + (dir == 0 ? 1 : 0); ++i) {
            for (j = 0; j < slot->recv_loc_end[1] - slot->recv_loc_start[1] + (dir == 1 ? 1 : 0); ++j) {
                for (k = 0; k < slot->recv_loc_end[2] - slot->recv_loc_start[2] + (dir == 2 ? 1 : 0); ++k) {
                    if (prj_mpi_bf_face_active(dir,
                            i + slot->recv_loc_start[0],
                            j + slot->recv_loc_start[1],
                            k + slot->recv_loc_start[2])) {
                        continue;
                    }
                    if (pos >= value_count) {
                        return 1;
                    }
                    prj_boundary_write_bf_face(dst, use_bf1, dir,
                        i + slot->recv_loc_start[0],
                        j + slot->recv_loc_start[1],
                        k + slot->recv_loc_start[2],
                        values[pos++], PRJ_MHD_FIDELITY_SAME);
                }
            }
        }
    }
    return pos == value_count ? 0 : 1;
}

static int prj_mpi_apply_bf_restriction_record(prj_block *coarse, int use_bf1,
    const prj_neighbor *slot, const double *values, int value_count)
{
    int dir;
    int pos = 0;

    for (dir = 0; dir < 3; ++dir) {
        int dj = slot->send_loc_end[1] - slot->send_loc_start[1] + (dir == 1 ? 1 : 0) + 1;
        int dk = slot->send_loc_end[2] - slot->send_loc_start[2] + (dir == 2 ? 1 : 0) + 1;
        int dir_count =
            (slot->send_loc_end[0] - slot->send_loc_start[0] + (dir == 0 ? 1 : 0) + 1) *
            dj * dk;
        const double *buf;
        int tan0 = (dir + 1) % 3;
        int tan1 = (dir + 2) % 3;
        int i;
        int j;
        int k;

        if (pos + dir_count > value_count) {
            return 1;
        }
        buf = values + pos;
        for (i = 0; i < slot->recv_loc_end[0] - slot->recv_loc_start[0] + (dir == 0 ? 1 : 0); ++i) {
            for (j = 0; j < slot->recv_loc_end[1] - slot->recv_loc_start[1] + (dir == 1 ? 1 : 0); ++j) {
                for (k = 0; k < slot->recv_loc_end[2] - slot->recv_loc_start[2] + (dir == 2 ? 1 : 0); ++k) {
                    int it_send0[3] = {0, 0, 0};
                    int it_send1[3] = {0, 0, 0};
                    int it_send2[3] = {0, 0, 0};
                    int it_send3[3] = {0, 0, 0};
                    int b0, b1, b2, b3;
                    double value;

                    it_send0[0] = 2 * i + slot->send_loc_start[0];
                    it_send0[1] = 2 * j + slot->send_loc_start[1];
                    it_send0[2] = 2 * k + slot->send_loc_start[2];

                    it_send1[dir] = it_send0[dir];
                    it_send1[tan0] = it_send0[tan0];
                    it_send1[tan1] = it_send0[tan1] + 1;

                    it_send2[dir] = it_send0[dir];
                    it_send2[tan0] = it_send0[tan0] + 1;
                    it_send2[tan1] = it_send0[tan1];

                    it_send3[dir] = it_send0[dir];
                    it_send3[tan0] = it_send0[tan0] + 1;
                    it_send3[tan1] = it_send0[tan1] + 1;

                    b0 = ((it_send0[0] - slot->send_loc_start[0]) * dj +
                          (it_send0[1] - slot->send_loc_start[1])) * dk +
                         (it_send0[2] - slot->send_loc_start[2]);
                    b1 = ((it_send1[0] - slot->send_loc_start[0]) * dj +
                          (it_send1[1] - slot->send_loc_start[1])) * dk +
                         (it_send1[2] - slot->send_loc_start[2]);
                    b2 = ((it_send2[0] - slot->send_loc_start[0]) * dj +
                          (it_send2[1] - slot->send_loc_start[1])) * dk +
                         (it_send2[2] - slot->send_loc_start[2]);
                    b3 = ((it_send3[0] - slot->send_loc_start[0]) * dj +
                          (it_send3[1] - slot->send_loc_start[1])) * dk +
                         (it_send3[2] - slot->send_loc_start[2]);

                    if (b0 < 0 || b1 < 0 || b2 < 0 || b3 < 0 ||
                        b0 >= dir_count || b1 >= dir_count ||
                        b2 >= dir_count || b3 >= dir_count) {
                        return 1;
                    }
                    value = 0.25 * (buf[b0] + buf[b1] + buf[b2] + buf[b3]);
                    prj_boundary_write_bf_face(coarse, use_bf1, dir,
                        i + slot->recv_loc_start[0],
                        j + slot->recv_loc_start[1],
                        k + slot->recv_loc_start[2],
                        value, PRJ_MHD_FIDELITY_FINER);
                }
            }
        }
        pos += dir_count;
    }
    return pos == value_count ? 0 : 1;
}

static int prj_mpi_apply_bf_prolongation_record(prj_block *fine, int use_bf1,
    const prj_neighbor *slot, const double *values, int value_count, int use_BJ)
{
    const double *buf[3];
    int buf_lo[3][3];
    int buf_n[3][3];
    double coarse_dx[3];
    int pos = 0;
    int dir;
    int i;
    int j;
    int k;

    for (dir = 0; dir < 3; ++dir) {
        int d;
        int dir_count;

        for (d = 0; d < 3; ++d) {
            if (d == dir) {
                buf_lo[dir][d] = slot->send_loc_start[d];
                buf_n[dir][d] = slot->send_loc_end[d] - slot->send_loc_start[d] + 2;
            } else {
                buf_lo[dir][d] = slot->send_loc_start[d] - 1;
                buf_n[dir][d] = slot->send_loc_end[d] - slot->send_loc_start[d] + 3;
            }
        }
        dir_count = buf_n[dir][0] * buf_n[dir][1] * buf_n[dir][2];
        if (pos + dir_count > value_count) {
            return 1;
        }
        buf[dir] = values + pos;
        pos += dir_count;
    }
    if (pos != value_count) {
        return 1;
    }
    coarse_dx[0] = 2.0 * fine->dx[0];
    coarse_dx[1] = 2.0 * fine->dx[1];
    coarse_dx[2] = 2.0 * fine->dx[2];

    for (i = 0; i < slot->recv_loc_end[0] - slot->recv_loc_start[0]; i += 2) {
        for (j = 0; j < slot->recv_loc_end[1] - slot->recv_loc_start[1]; j += 2) {
            for (k = 0; k < slot->recv_loc_end[2] - slot->recv_loc_start[2]; k += 2) {
                int fi = i + slot->recv_loc_start[0];
                int fj = j + slot->recv_loc_start[1];
                int fk = k + slot->recv_loc_start[2];
                int ci = i / 2 + slot->send_loc_start[0];
                int cj = j / 2 + slot->send_loc_start[1];
                int ck = k / 2 + slot->send_loc_start[2];

                prj_mpi_bf_prolong_from_buffer(buf, buf_lo, buf_n,
                    coarse_dx, fine, ci, cj, ck, fi, fj, fk, use_bf1, use_BJ);
            }
        }
    }
    return 0;
}

static int prj_mpi_apply_bf_records(prj_mesh *mesh, const prj_mpi *mpi, int use_bf1,
    int fill_kind, const int *headers, int record_count,
    const double *values, int value_count)
{
    int r;
    int pos = 0;

    for (r = 0; r < record_count; ++r) {
        const int *header = headers + r * PRJ_MPI_BF_HEADER_NINT;
        int block_id = header[PRJ_MPI_BF_HDR_DST_ID];
        int rel_level = header[PRJ_MPI_BF_HDR_REL_LEVEL];
        int nvalue = header[PRJ_MPI_BF_HDR_VALUE_COUNT];
        prj_neighbor slot;
        prj_block *dst;
        int err = 0;

        if (nvalue < 0 || pos + nvalue > value_count) {
            return 1;
        }
        if (!prj_mpi_bf_slot_matches(fill_kind, rel_level)) {
            pos += nvalue;
            continue;
        }
        if (block_id < 0 || block_id >= mesh->nblocks) {
            pos += nvalue;
            continue;
        }
        dst = &mesh->blocks[block_id];
        if (!prj_mpi_block_is_local(mpi, dst)) {
            pos += nvalue;
            continue;
        }
        prj_mpi_check_bf_storage(dst, "prj_mpi_apply_bf_records");
        prj_mpi_bf_header_to_slot(header, &slot);
        if (rel_level == 0) {
            err = prj_mpi_apply_bf_same_level_record(dst, use_bf1,
                &slot, values + pos, nvalue);
        } else if (rel_level < 0) {
            err = prj_mpi_apply_bf_restriction_record(dst, use_bf1,
                &slot, values + pos, nvalue);
        } else {
            err = prj_mpi_apply_bf_prolongation_record(dst, use_bf1,
                &slot, values + pos, nvalue, mesh->use_BJ);
        }
        if (err != 0) {
            return 1;
        }
        pos += nvalue;
    }
    return pos == value_count ? 0 : 1;
}

static int prj_mpi_apply_bf_values_from_sender(prj_mesh *mesh, const prj_mpi *mpi,
    int sender_rank, int use_bf1, int fill_kind, const double *values,
    int value_count)
{
    int sbidx;
    int pos = 0;

    if (mesh == 0 || mpi == 0) {
        return 1;
    }
    for (sbidx = 0; sbidx < mesh->nblocks; ++sbidx) {
        const prj_block *sblock = &mesh->blocks[sbidx];
        int n;

        if (!prj_block_is_active(sblock) || sblock->rank != sender_rank) {
            continue;
        }
        for (n = 0; n < 56; ++n) {
            const prj_neighbor *slot = &sblock->slot[n];
            int nid = slot->id;
            int expected;
            prj_block *dst;
            int err = 0;

            if (nid < 0 || nid >= mesh->nblocks ||
                mesh->blocks[nid].rank != mpi->rank ||
                !prj_mpi_bf_slot_matches(fill_kind, slot->rel_level)) {
                continue;
            }
            expected = prj_mpi_bf_record_value_count(fill_kind, slot);
            if (expected <= 0) {
                continue;
            }
            if (pos + expected > value_count) {
                return 1;
            }
            dst = &mesh->blocks[nid];
            if (!prj_mpi_block_is_local(mpi, dst)) {
                pos += expected;
                continue;
            }
            prj_mpi_check_bf_storage(dst, "prj_mpi_apply_bf_values_from_sender");
            if (slot->rel_level == 0) {
                err = prj_mpi_apply_bf_same_level_record(dst, use_bf1,
                    slot, values + pos, expected);
            } else if (slot->rel_level < 0) {
                err = prj_mpi_apply_bf_restriction_record(dst, use_bf1,
                    slot, values + pos, expected);
            } else {
                err = prj_mpi_apply_bf_prolongation_record(dst, use_bf1,
                    slot, values + pos, expected, mesh->use_BJ);
            }
            if (err != 0) {
                return 1;
            }
            pos += expected;
        }
    }
    return pos == value_count ? 0 : 1;
}

static void prj_mpi_count_bf_records(const prj_mesh *mesh, const prj_mpi *mpi,
    int receiver_rank, int fill_kind, int *record_count, int *value_count)
{
    int bidx;

    if (record_count != 0) {
        *record_count = 0;
    }
    if (value_count != 0) {
        *value_count = 0;
    }
    if (mesh == 0 || mpi == 0 || record_count == 0 || value_count == 0) {
        return;
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        const prj_block *block = &mesh->blocks[bidx];
        int n;

        if (!prj_mpi_block_is_local(mpi, block)) {
            continue;
        }
        for (n = 0; n < 56; ++n) {
            const prj_neighbor *slot = &block->slot[n];
            int nid = slot->id;
            int expected;

            if (nid < 0 || nid >= mesh->nblocks ||
                mesh->blocks[nid].rank != receiver_rank ||
                !prj_mpi_bf_slot_matches(fill_kind, slot->rel_level)) {
                continue;
            }
            expected = prj_mpi_bf_record_value_count(fill_kind, slot);
            if (expected <= 0) {
                continue;
            }
            *record_count += 1;
            *value_count += expected;
        }
    }
}

static int prj_mpi_build_bf_plan_for_neighbor(const prj_mesh *mesh,
    const prj_mpi *mpi, prj_mpi_buffer *buffer)
{
#if defined(PRJ_ENABLE_MPI)
    int fill_kind;
    int max_send_records = 0;
    int max_send_values = 0;
    int max_recv_records = 0;
    int max_recv_values = 0;

    if (mesh == 0 || mpi == 0 || buffer == 0) {
        return 1;
    }
    for (fill_kind = 0; fill_kind < PRJ_MPI_BF_FILL_N; ++fill_kind) {
        int send_counts[2];
        int recv_counts[2];

        prj_mpi_count_bf_records(mesh, mpi, buffer->receiver_rank, fill_kind,
            &buffer->bf_send_record_count[fill_kind],
            &buffer->bf_send_value_count[fill_kind]);
        send_counts[0] = buffer->bf_send_record_count[fill_kind];
        send_counts[1] = buffer->bf_send_value_count[fill_kind];
        MPI_Sendrecv(send_counts, 2, MPI_INT, buffer->receiver_rank,
            310 + fill_kind, recv_counts, 2, MPI_INT,
            buffer->receiver_rank, 310 + fill_kind, MPI_COMM_WORLD,
            MPI_STATUS_IGNORE);
        if (recv_counts[0] < 0 || recv_counts[1] < 0) {
            return 1;
        }
        buffer->bf_recv_record_count[fill_kind] = recv_counts[0];
        buffer->bf_recv_value_count[fill_kind] = recv_counts[1];
        if (send_counts[0] > max_send_records) {
            max_send_records = send_counts[0];
        }
        if (send_counts[1] > max_send_values) {
            max_send_values = send_counts[1];
        }
        if (recv_counts[0] > max_recv_records) {
            max_recv_records = recv_counts[0];
        }
        if (recv_counts[1] > max_recv_values) {
            max_recv_values = recv_counts[1];
        }
    }
    buffer->bf_send_record_capacity = max_send_records;
    buffer->bf_send_value_capacity = max_send_values;
    buffer->bf_recv_record_capacity = max_recv_records;
    buffer->bf_recv_value_capacity = max_recv_values;
    if (max_send_records > 0) {
        buffer->bf_headers_send = (int *)calloc(
            (size_t)max_send_records * (size_t)PRJ_MPI_BF_HEADER_NINT,
            sizeof(*buffer->bf_headers_send));
    }
    if (max_send_values > 0) {
        buffer->bf_values_send = (double *)calloc((size_t)max_send_values,
            sizeof(*buffer->bf_values_send));
    }
    if (max_recv_records > 0) {
        buffer->bf_headers_recv = (int *)calloc(
            (size_t)max_recv_records * (size_t)PRJ_MPI_BF_HEADER_NINT,
            sizeof(*buffer->bf_headers_recv));
    }
    if (max_recv_values > 0) {
        buffer->bf_values_recv = (double *)calloc((size_t)max_recv_values,
            sizeof(*buffer->bf_values_recv));
    }
    if ((max_send_records > 0 && buffer->bf_headers_send == 0) ||
        (max_send_values > 0 && buffer->bf_values_send == 0) ||
        (max_recv_records > 0 && buffer->bf_headers_recv == 0) ||
        (max_recv_values > 0 && buffer->bf_values_recv == 0)) {
        return 1;
    }
    return 0;
#else
    (void)mesh;
    (void)mpi;
    (void)buffer;
    return 0;
#endif
}

void prj_mpi_exchange_bf(prj_mesh *mesh, prj_mpi *mpi, int use_bf1, int fill_kind)
{
#if defined(PRJ_ENABLE_MPI)
    int nb;

    if (mesh == 0 || mpi == 0 || mpi->totrank <= 1) {
        return;
    }
    if (fill_kind < 0 || fill_kind >= PRJ_MPI_BF_FILL_N) {
        fprintf(stderr, "prj_mpi_exchange_bf: invalid fill_kind=%d\n",
            fill_kind);
        exit(EXIT_FAILURE);
    }
    for (nb = 0; nb < mpi->neighbor_number; ++nb) {
        prj_mpi_buffer *buffer = &mpi->neighbor_buffer[nb];
        int send_record_count = 0;
        int send_value_count = 0;
        int recv_record_count = buffer->bf_recv_record_count[fill_kind];
        int recv_value_count = buffer->bf_recv_value_count[fill_kind];

        if (prj_mpi_pack_bf_records(mesh, mpi, buffer->receiver_rank, use_bf1,
                fill_kind, buffer->bf_headers_send,
                buffer->bf_send_record_capacity, &send_record_count,
                buffer->bf_values_send, buffer->bf_send_value_capacity,
                &send_value_count) != 0 ||
            send_record_count != buffer->bf_send_record_count[fill_kind] ||
            send_value_count != buffer->bf_send_value_count[fill_kind]) {
            fprintf(stderr, "prj_mpi_exchange_bf: failed to pack Bf records\n");
            exit(EXIT_FAILURE);
        }
        MPI_Sendrecv(buffer->bf_headers_send,
            send_record_count * PRJ_MPI_BF_HEADER_NINT, MPI_INT,
            buffer->receiver_rank, 301,
            buffer->bf_headers_recv,
            recv_record_count * PRJ_MPI_BF_HEADER_NINT, MPI_INT,
            buffer->receiver_rank, 301,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Sendrecv(buffer->bf_values_send, send_value_count, MPI_DOUBLE,
            buffer->receiver_rank, 302, buffer->bf_values_recv,
            recv_value_count, MPI_DOUBLE,
            buffer->receiver_rank, 302,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        if (prj_mpi_apply_bf_records(mesh, mpi, use_bf1, fill_kind,
                buffer->bf_headers_recv, recv_record_count,
                buffer->bf_values_recv,
                recv_value_count) != 0) {
            fprintf(stderr, "prj_mpi_exchange_bf: failed to unpack Bf records\n");
            exit(EXIT_FAILURE);
        }
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
enum {
    PRJ_MPI_AMR_BF_HDR_PARENT_ID = 0,
    PRJ_MPI_AMR_BF_HDR_CHILD_OCT = 1,
    PRJ_MPI_AMR_BF_HDR_NEIGHBOR_ID = 2,
    PRJ_MPI_AMR_BF_HDR_USE_BF1 = 3,
    PRJ_MPI_AMR_BF_HDR_DIR = 4,
    PRJ_MPI_AMR_BF_HDR_RECV_FACE = 5,
    PRJ_MPI_AMR_BF_HEADER_NINT = 6
};

#define PRJ_MPI_AMR_BF_FACE_NVALUE (PRJ_BLOCK_SIZE * PRJ_BLOCK_SIZE)

static void prj_mpi_amr_bf_clear_cache(prj_mpi *mpi)
{
    if (mpi == 0) {
        return;
    }
    mpi->amr_bf_record_count = 0;
    mpi->amr_bf_value_count = 0;
}

static int prj_mpi_amr_bf_cache_append(prj_mpi *mpi, const int *headers,
    int record_count, const double *values, int value_count)
{
    if (mpi == 0 || record_count < 0 || value_count < 0 ||
        value_count != record_count * PRJ_MPI_AMR_BF_FACE_NVALUE) {
        return 1;
    }
    if (record_count == 0) {
        return 0;
    }
    if (headers == 0 || values == 0 ||
        mpi->amr_bf_record_count + record_count > mpi->amr_bf_record_capacity ||
        mpi->amr_bf_value_count + value_count > mpi->amr_bf_value_capacity) {
        return 1;
    }
    memcpy(mpi->amr_bf_headers +
            (size_t)mpi->amr_bf_record_count * (size_t)PRJ_MPI_AMR_BF_HEADER_NINT,
        headers, (size_t)record_count * (size_t)PRJ_MPI_AMR_BF_HEADER_NINT * sizeof(*headers));
    memcpy(mpi->amr_bf_values + mpi->amr_bf_value_count,
        values, (size_t)value_count * sizeof(*values));
    mpi->amr_bf_record_count += record_count;
    mpi->amr_bf_value_count += value_count;
    return 0;
}

static void prj_mpi_amr_bf_child_bounds(const prj_block *parent, int child_oct,
    double child_xmin[3], double child_xmax[3], double child_dx[3])
{
    int d;

    for (d = 0; d < 3; ++d) {
        double xmid = 0.5 * (parent->xmin[d] + parent->xmax[d]);
        int bit = (child_oct >> d) & 1;

        child_xmin[d] = bit == 0 ? parent->xmin[d] : xmid;
        child_xmax[d] = bit == 0 ? xmid : parent->xmax[d];
        child_dx[d] = 0.5 * parent->dx[d];
    }
}

static int prj_mpi_amr_bf_slot_face_for_child(const prj_block *parent,
    const prj_neighbor *slot, int child_oct, int dir, int *recv_face,
    int *send_face)
{
    double child_xmin[3];
    double child_xmax[3];
    double child_dx[3];
    int tan0;
    int tan1;

    if (parent == 0 || slot == 0 || child_oct < 0 || child_oct >= 8 ||
        dir < 0 || dir >= 3) {
        return 0;
    }
    prj_mpi_amr_bf_child_bounds(parent, child_oct, child_xmin, child_xmax, child_dx);
    tan0 = (dir + 1) % 3;
    tan1 = (dir + 2) % 3;
    if (fabs(child_xmin[tan0] - slot->xmin[tan0]) >= 1.0e-12 * child_dx[tan0] ||
        fabs(child_xmin[tan1] - slot->xmin[tan1]) >= 1.0e-12 * child_dx[tan1]) {
        return 0;
    }
    if (fabs(child_xmax[dir] - slot->xmin[dir]) < 1.0e-12 * child_dx[dir]) {
        if (recv_face != 0) {
            *recv_face = PRJ_BLOCK_SIZE;
        }
        if (send_face != 0) {
            *send_face = 0;
        }
        return 1;
    }
    if (fabs(child_xmin[dir] - slot->xmax[dir]) < 1.0e-12 * child_dx[dir]) {
        if (recv_face != 0) {
            *recv_face = 0;
        }
        if (send_face != 0) {
            *send_face = PRJ_BLOCK_SIZE;
        }
        return 1;
    }
    return 0;
}

static int prj_mpi_amr_bf_pack_face_values(const prj_block *neighbor,
    int use_bf1, int dir, int send_face, double *values, int value_capacity,
    int *value_count)
{
    const double *src;
    int tan0 = (dir + 1) % 3;
    int tan1 = (dir + 2) % 3;
    int i;
    int j;

    prj_mpi_check_bf_storage(neighbor, "prj_mpi_amr_bf_pack_face_values");
    src = prj_mpi_bf_array_const(neighbor, dir, use_bf1);
    for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
        for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
            int it_send[3] = {0, 0, 0};

            it_send[dir] = send_face;
            it_send[tan0] = i;
            it_send[tan1] = j;
            if (prj_mpi_bf_append_value(values, value_count, value_capacity,
                    src[FACE_IDX(dir, it_send[0], it_send[1], it_send[2])]) != 0) {
                return 1;
            }
        }
    }
    return 0;
}

static int prj_mpi_amr_bf_pack_records(const prj_mesh *mesh, const prj_mpi *mpi,
    int receiver_rank,
    int require_refine_flag, double *values, int value_capacity, int *value_count)
{
    int parent_idx;

    *value_count = 0;
    if (mesh == 0) {
        return 1;
    }
    for (parent_idx = 0; parent_idx < mesh->nblocks; ++parent_idx) {
        const prj_block *parent = &mesh->blocks[parent_idx];
        int n;

        if (!prj_block_is_active(parent) || parent->rank != receiver_rank ||
            (require_refine_flag != 0 && parent->refine_flag <= 0)) {
            continue;
        }
        for (n = 0; n < 56; ++n) {
            const prj_neighbor *slot = &parent->slot[n];
            const prj_block *neighbor;
            int child_oct;

            if (slot->id < 0 || slot->id >= mesh->nblocks ||
                slot->type != PRJ_NEIGHBOR_FACE || slot->rel_level < 1) {
                continue;
            }
            neighbor = &mesh->blocks[slot->id];
            if (!prj_mpi_block_is_local(mpi, neighbor)) {
                continue;
            }
            for (child_oct = 0; child_oct < 8; ++child_oct) {
                int dir;

                for (dir = 0; dir < 3; ++dir) {
                    int recv_face;
                    int send_face;
                    int use_bf1;

                    if (!prj_mpi_amr_bf_slot_face_for_child(parent, slot,
                            child_oct, dir, &recv_face, &send_face)) {
                        continue;
                    }
                    (void)recv_face;
                    for (use_bf1 = 0; use_bf1 < 2; ++use_bf1) {
                        if (prj_mpi_amr_bf_pack_face_values(neighbor, use_bf1,
                                dir, send_face, values, value_capacity,
                                value_count) != 0) {
                            *value_count = 0;
                            return 1;
                        }
                    }
                }
            }
        }
    }
    return 0;
}

static int prj_mpi_amr_bf_count_records_pair(const prj_mesh *mesh,
    int sender_rank, int receiver_rank, int require_refine_flag)
{
    int count = 0;
    int parent_idx;

    if (mesh == 0) {
        return 0;
    }
    for (parent_idx = 0; parent_idx < mesh->nblocks; ++parent_idx) {
        const prj_block *parent = &mesh->blocks[parent_idx];
        int n;

        if (!prj_block_is_active(parent) || parent->rank != receiver_rank ||
            (require_refine_flag != 0 && parent->refine_flag <= 0)) {
            continue;
        }
        for (n = 0; n < 56; ++n) {
            const prj_neighbor *slot = &parent->slot[n];
            const prj_block *neighbor;
            int child_oct;

            if (slot->id < 0 || slot->id >= mesh->nblocks ||
                slot->type != PRJ_NEIGHBOR_FACE || slot->rel_level < 1) {
                continue;
            }
            neighbor = &mesh->blocks[slot->id];
            if (!prj_block_is_active(neighbor) || neighbor->rank != sender_rank) {
                continue;
            }
            for (child_oct = 0; child_oct < 8; ++child_oct) {
                int dir;

                for (dir = 0; dir < 3; ++dir) {
                    if (prj_mpi_amr_bf_slot_face_for_child(parent, slot,
                            child_oct, dir, 0, 0)) {
                        count += 2;
                    }
                }
            }
        }
    }
    return count;
}

static int prj_mpi_amr_bf_count_records(const prj_mesh *mesh, const prj_mpi *mpi, int receiver_rank)
{
    int sender_rank = mpi != 0 ? mpi->rank : 0;

    return prj_mpi_amr_bf_count_records_pair(mesh, sender_rank, receiver_rank, 0);
}

/*
 * Reproduce the iteration order used by prj_mpi_amr_bf_pack_records so the
 * receiver can derive the exact header sequence the sender will pack — the
 * walk is over receiver-owned parents (with refine_flag set) whose finer
 * face neighbors belong to the sender, and the header fields are pure
 * mesh-topology values both ranks already share.
 */
static int prj_mpi_amr_bf_derive_headers_pair(const prj_mesh *mesh,
    int sender_rank, int receiver_rank, int require_refine_flag,
    int *headers, int header_capacity, int *record_count)
{
    int parent_idx;

    *record_count = 0;
    if (mesh == 0) {
        return 1;
    }
    for (parent_idx = 0; parent_idx < mesh->nblocks; ++parent_idx) {
        const prj_block *parent = &mesh->blocks[parent_idx];
        int n;

        if (!prj_block_is_active(parent) || parent->rank != receiver_rank ||
            (require_refine_flag != 0 && parent->refine_flag <= 0)) {
            continue;
        }
        for (n = 0; n < 56; ++n) {
            const prj_neighbor *slot = &parent->slot[n];
            const prj_block *neighbor;
            int child_oct;

            if (slot->id < 0 || slot->id >= mesh->nblocks ||
                slot->type != PRJ_NEIGHBOR_FACE || slot->rel_level < 1) {
                continue;
            }
            neighbor = &mesh->blocks[slot->id];
            if (!prj_block_is_active(neighbor) || neighbor->rank != sender_rank) {
                continue;
            }
            for (child_oct = 0; child_oct < 8; ++child_oct) {
                int dir;

                for (dir = 0; dir < 3; ++dir) {
                    int recv_face;
                    int send_face;
                    int use_bf1;

                    if (!prj_mpi_amr_bf_slot_face_for_child(parent, slot,
                            child_oct, dir, &recv_face, &send_face)) {
                        continue;
                    }
                    (void)send_face;
                    for (use_bf1 = 0; use_bf1 < 2; ++use_bf1) {
                        int *header;

                        if (*record_count >= header_capacity) {
                            *record_count = 0;
                            return 1;
                        }
                        header = headers + (size_t)(*record_count) *
                            (size_t)PRJ_MPI_AMR_BF_HEADER_NINT;
                        header[PRJ_MPI_AMR_BF_HDR_PARENT_ID] = parent->id;
                        header[PRJ_MPI_AMR_BF_HDR_CHILD_OCT] = child_oct;
                        header[PRJ_MPI_AMR_BF_HDR_NEIGHBOR_ID] = slot->id;
                        header[PRJ_MPI_AMR_BF_HDR_USE_BF1] = use_bf1;
                        header[PRJ_MPI_AMR_BF_HDR_DIR] = dir;
                        header[PRJ_MPI_AMR_BF_HDR_RECV_FACE] = recv_face;
                        *record_count += 1;
                    }
                }
            }
        }
    }
    return 0;
}

static int prj_mpi_build_amr_bf_plan_for_neighbor(const prj_mesh *mesh, const prj_mpi *mpi,
    prj_mpi_buffer *buffer)
{
#if defined(PRJ_ENABLE_MPI)
    int send_counts[2];
    int recv_counts[2];

    if (mesh == 0 || buffer == 0) {
        return 1;
    }
    send_counts[0] = prj_mpi_amr_bf_count_records(mesh, mpi, buffer->receiver_rank);
    send_counts[1] = send_counts[0] * PRJ_MPI_AMR_BF_FACE_NVALUE;
    MPI_Sendrecv(send_counts, 2, MPI_INT, buffer->receiver_rank, 503,
        recv_counts, 2, MPI_INT, buffer->receiver_rank, 503,
        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if (recv_counts[0] < 0 || recv_counts[1] < 0 ||
        recv_counts[1] != recv_counts[0] * PRJ_MPI_AMR_BF_FACE_NVALUE) {
        return 1;
    }
    buffer->amr_bf_send_value_capacity = send_counts[1];
    buffer->amr_bf_recv_record_capacity = recv_counts[0];
    buffer->amr_bf_recv_value_capacity = recv_counts[1];
    if (send_counts[1] > 0) {
        buffer->amr_bf_values_send = (double *)calloc(
            (size_t)send_counts[1], sizeof(*buffer->amr_bf_values_send));
    }
    if (recv_counts[0] > 0) {
        buffer->amr_bf_headers_recv = (int *)calloc(
            (size_t)recv_counts[0] * (size_t)PRJ_MPI_AMR_BF_HEADER_NINT,
            sizeof(*buffer->amr_bf_headers_recv));
    }
    if (recv_counts[1] > 0) {
        buffer->amr_bf_values_recv = (double *)calloc(
            (size_t)recv_counts[1], sizeof(*buffer->amr_bf_values_recv));
    }
    if ((send_counts[1] > 0 && buffer->amr_bf_values_send == 0) ||
        (recv_counts[0] > 0 && buffer->amr_bf_headers_recv == 0) ||
        (recv_counts[1] > 0 && buffer->amr_bf_values_recv == 0)) {
        return 1;
    }
    return 0;
#else
    (void)mesh;
    (void)buffer;
    return 0;
#endif
}

static int prj_mpi_build_amr_bf_cache(prj_mpi *mpi)
{
    int nb;
    int record_capacity = 0;
    int value_capacity = 0;

    if (mpi == 0) {
        return 1;
    }
    for (nb = 0; nb < mpi->neighbor_number; ++nb) {
        record_capacity += mpi->neighbor_buffer[nb].amr_bf_recv_record_capacity;
        value_capacity += mpi->neighbor_buffer[nb].amr_bf_recv_value_capacity;
    }
    mpi->amr_bf_record_capacity = record_capacity;
    mpi->amr_bf_value_capacity = value_capacity;
    mpi->amr_bf_record_count = 0;
    mpi->amr_bf_value_count = 0;
    if (record_capacity > 0) {
        mpi->amr_bf_headers = (int *)calloc(
            (size_t)record_capacity * (size_t)PRJ_MPI_AMR_BF_HEADER_NINT,
            sizeof(*mpi->amr_bf_headers));
    }
    if (value_capacity > 0) {
        mpi->amr_bf_values = (double *)calloc((size_t)value_capacity,
            sizeof(*mpi->amr_bf_values));
    }
    if ((record_capacity > 0 && mpi->amr_bf_headers == 0) ||
        (value_capacity > 0 && mpi->amr_bf_values == 0)) {
        return 1;
    }
    return 0;
}

void prj_mpi_exchange_amr_mhd_prolongate_bf(const prj_mesh *mesh, prj_mpi *mpi)
{
    prj_mpi_amr_bf_clear_cache(mpi);
#if defined(PRJ_ENABLE_MPI)
    if (mesh == 0 || mpi == 0 || mpi->totrank <= 1) {
        return;
    }
    {
        MPI_Request *requests = (MPI_Request *)mpi->request_buffer;
        int nb;
        int request_count = 0;

        if (mpi->neighbor_number > 0 &&
            (requests == 0 || mpi->request_capacity < 2 * mpi->neighbor_number)) {
            fprintf(stderr,
                "prj_mpi_exchange_amr_mhd_prolongate_bf: insufficient MPI request buffer\n");
            exit(EXIT_FAILURE);
        }

        /*
         * Loop 1: derive counts, pack send-side values, post non-blocking
         * sends. Both ranks share mesh topology + synchronized refine_flag,
         * so each rank derives the peer's record counts locally.
         */
        for (nb = 0; nb < mpi->neighbor_number; ++nb) {
            prj_mpi_buffer *buffer = &mpi->neighbor_buffer[nb];
            int send_records = prj_mpi_amr_bf_count_records_pair(mesh,
                mpi->rank, buffer->receiver_rank, 1);
            int recv_records = prj_mpi_amr_bf_count_records_pair(mesh,
                buffer->receiver_rank, mpi->rank, 1);
            int send_value_count = send_records * PRJ_MPI_AMR_BF_FACE_NVALUE;
            int recv_value_count = recv_records * PRJ_MPI_AMR_BF_FACE_NVALUE;

            if (recv_records > buffer->amr_bf_recv_record_capacity ||
                recv_value_count > buffer->amr_bf_recv_value_capacity ||
                send_value_count > buffer->amr_bf_send_value_capacity) {
                fprintf(stderr,
                    "prj_mpi_exchange_amr_mhd_prolongate_bf: derived Bf counts exceed capacity\n");
                exit(EXIT_FAILURE);
            }
            buffer->amr_bf_send_value_count = send_value_count;
            buffer->amr_bf_recv_record_count = recv_records;
            buffer->amr_bf_recv_value_count = recv_value_count;

            if (send_records > 0) {
                int packed_value_count = 0;

                if (prj_mpi_amr_bf_pack_records(mesh, mpi, buffer->receiver_rank,
                        1, buffer->amr_bf_values_send,
                        buffer->amr_bf_send_value_capacity,
                        &packed_value_count) != 0 ||
                    packed_value_count != send_value_count) {
                    fprintf(stderr,
                        "prj_mpi_exchange_amr_mhd_prolongate_bf: packed Bf value count %d does not match derived %d\n",
                        packed_value_count, send_value_count);
                    exit(EXIT_FAILURE);
                }
                MPI_Isend(buffer->amr_bf_values_send, send_value_count,
                    MPI_DOUBLE, buffer->receiver_rank, 502,
                    MPI_COMM_WORLD, &requests[request_count++]);
            }
        }

        /* Loop 2: post non-blocking receives for every active peer. */
        for (nb = 0; nb < mpi->neighbor_number; ++nb) {
            prj_mpi_buffer *buffer = &mpi->neighbor_buffer[nb];

            if (buffer->amr_bf_recv_value_count > 0) {
                MPI_Irecv(buffer->amr_bf_values_recv,
                    buffer->amr_bf_recv_value_count, MPI_DOUBLE,
                    buffer->receiver_rank, 502, MPI_COMM_WORLD,
                    &requests[request_count++]);
            }
        }

        if (request_count > 0) {
            MPI_Waitall(request_count, requests, MPI_STATUSES_IGNORE);
        }

        /* Loop 3: derive recv headers (locally) and append into the cache. */
        for (nb = 0; nb < mpi->neighbor_number; ++nb) {
            prj_mpi_buffer *buffer = &mpi->neighbor_buffer[nb];
            int recv_records = buffer->amr_bf_recv_record_count;
            int derived_recv_records = 0;

            if (recv_records == 0) {
                continue;
            }
            if (prj_mpi_amr_bf_derive_headers_pair(mesh,
                    buffer->receiver_rank, mpi->rank, 1,
                    buffer->amr_bf_headers_recv,
                    buffer->amr_bf_recv_record_capacity,
                    &derived_recv_records) != 0 ||
                derived_recv_records != recv_records) {
                fprintf(stderr,
                    "prj_mpi_exchange_amr_mhd_prolongate_bf: failed to derive Bf recv headers (%d vs %d)\n",
                    derived_recv_records, recv_records);
                exit(EXIT_FAILURE);
            }
            if (prj_mpi_amr_bf_cache_append(mpi,
                    buffer->amr_bf_headers_recv, recv_records,
                    buffer->amr_bf_values_recv,
                    buffer->amr_bf_recv_value_count) != 0) {
                fprintf(stderr,
                    "prj_mpi_exchange_amr_mhd_prolongate_bf: failed to cache Bf records\n");
                exit(EXIT_FAILURE);
            }
        }
    }
#else
    (void)mesh;
    (void)mpi;
#endif
}

int prj_mpi_amr_mhd_prolongate_bf_one(const prj_mpi *mpi,
    const prj_block *parent,
    const prj_neighbor *slot, prj_block *child, int child_oct, int use_bf1,
    int dir)
{
#if defined(PRJ_ENABLE_MPI)
    int recv_face;
    int send_face;
    int r;

    if (mpi == 0 || mpi->totrank <= 1) {
        return 0;
    }
    if (parent == 0 || slot == 0 || child == 0 || dir < 0 || dir >= 3) {
        return 0;
    }
    if (!prj_mpi_amr_bf_slot_face_for_child(parent, slot, child_oct, dir,
            &recv_face, &send_face)) {
        return 0;
    }
    (void)send_face;
    if (mpi->amr_bf_value_count !=
        mpi->amr_bf_record_count * PRJ_MPI_AMR_BF_FACE_NVALUE) {
        fprintf(stderr,
            "prj_mpi_amr_mhd_prolongate_bf_one: invalid cached MPI Bf counts\n");
        exit(EXIT_FAILURE);
    }
    if (mpi->amr_bf_record_count > 0 &&
        (mpi->amr_bf_headers == 0 || mpi->amr_bf_values == 0)) {
        fprintf(stderr,
            "prj_mpi_amr_mhd_prolongate_bf_one: missing cached MPI Bf buffers\n");
        exit(EXIT_FAILURE);
    }
    for (r = 0; r < mpi->amr_bf_record_count; ++r) {
        const int *header = mpi->amr_bf_headers +
            (size_t)r * (size_t)PRJ_MPI_AMR_BF_HEADER_NINT;

        if (header[PRJ_MPI_AMR_BF_HDR_PARENT_ID] == parent->id &&
            header[PRJ_MPI_AMR_BF_HDR_CHILD_OCT] == child_oct &&
            header[PRJ_MPI_AMR_BF_HDR_NEIGHBOR_ID] == slot->id &&
            header[PRJ_MPI_AMR_BF_HDR_USE_BF1] == (use_bf1 != 0 ? 1 : 0) &&
            header[PRJ_MPI_AMR_BF_HDR_DIR] == dir &&
            header[PRJ_MPI_AMR_BF_HDR_RECV_FACE] == recv_face) {
            const double *values = mpi->amr_bf_values +
                (size_t)r * (size_t)PRJ_MPI_AMR_BF_FACE_NVALUE;
            int tan0 = (dir + 1) % 3;
            int tan1 = (dir + 2) % 3;
            int i;
            int j;
            int pos = 0;

            for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
                for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                    int it_recv[3] = {0, 0, 0};

                    it_recv[dir] = recv_face;
                    it_recv[tan0] = i;
                    it_recv[tan1] = j;
                    prj_boundary_write_bf_face(child, use_bf1, dir,
                        it_recv[0], it_recv[1], it_recv[2],
                        values[pos++], PRJ_MHD_FIDELITY_SAME);
                }
            }
            return 1;
        }
    }
    fprintf(stderr,
        "prj_mpi_amr_mhd_prolongate_bf_one: missing MPI Bf face "
        "parent=%d child_oct=%d neighbor=%d use_bf1=%d dir=%d recv_face=%d\n",
        parent->id, child_oct, slot->id, use_bf1 != 0 ? 1 : 0, dir, recv_face);
    exit(EXIT_FAILURE);
#else
    (void)mpi;
    (void)parent;
    (void)slot;
    (void)child;
    (void)child_oct;
    (void)use_bf1;
    (void)dir;
    return 0;
#endif
}

static int prj_mpi_edge_storage_index_ok(int dir, int i, int j, int k)
{
    int idx[3] = {i, j, k};
    int d;

    for (d = 0; d < 3; ++d) {
        int lo = -PRJ_NGHOST;
        int hi = dir == d ? PRJ_BLOCK_SIZE + PRJ_NGHOST - 1 : PRJ_BLOCK_SIZE + PRJ_NGHOST;
        if (idx[d] < lo || idx[d] > hi) {
            return 0;
        }
    }
    return 1;
}

static int prj_mpi_store_emf_edge_plan(const prj_block *fine, int coarse_id,
    int dir, const int send_idx[3], const int recv_idx[3],
    prj_mpi_buffer *buffer, int *count, int store)
{
    if (fine == 0 || dir < 0 || dir >= 3 || fine->emf[dir] == 0 ||
        fine->edge_fidelity[dir] == 0) {
        fprintf(stderr, "prj_mpi_store_emf_edge_plan: invalid fine emf block\n");
        exit(EXIT_FAILURE);
    }
    if (!prj_mpi_edge_storage_index_ok(dir, recv_idx[0], recv_idx[1], recv_idx[2])) {
        fprintf(stderr, "prj_mpi_store_emf_edge_plan: coarse edge index out of storage\n");
        exit(EXIT_FAILURE);
    }
    if (!prj_mpi_edge_storage_index_ok(dir, send_idx[0], send_idx[1], send_idx[2]) ||
        !prj_mpi_edge_storage_index_ok(dir,
            send_idx[0] + (dir == 0 ? 1 : 0),
            send_idx[1] + (dir == 1 ? 1 : 0),
            send_idx[2] + (dir == 2 ? 1 : 0))) {
        fprintf(stderr, "prj_mpi_store_emf_edge_plan: fine edge index out of storage\n");
        exit(EXIT_FAILURE);
    }
    if (store != 0) {
        if (*count >= buffer->emf_send_count) {
            return 1;
        }
        buffer->emf_idx_send[0][*count] = coarse_id;
        buffer->emf_idx_send[1][*count] =
            prj_mpi_encode_cell_index(recv_idx[0], recv_idx[1], recv_idx[2]);
        buffer->emf_idx_send[2][*count] = dir;
        buffer->emf_src_block[*count] = fine->id;
        buffer->emf_src_dir[*count] = dir;
        buffer->emf_src_idx[0][*count] = send_idx[0];
        buffer->emf_src_idx[1][*count] = send_idx[1];
        buffer->emf_src_idx[2][*count] = send_idx[2];
    }
    *count += 1;
    return 0;
}

static int prj_mpi_build_emf_restriction_plan(const prj_block *fine,
    const prj_block *coarse, prj_mpi_buffer *buffer, int *count, int store)
{
    int axis[3] = {0, 0, 0};
    int it_recv[3] = {-100, -100, -100};
    int it_send[3] = {-100, -100, -100};
    int touch = 0;
    int dir;
    int d;

    for (d = 0; d < 3; ++d) {
        if (fabs(fine->xmax[d] - coarse->xmin[d]) < 1.0e-12 * fine->dx[d]) {
            axis[d] = 1;
            touch += 1;
        } else if (fabs(coarse->xmax[d] - fine->xmin[d]) < 1.0e-12 * fine->dx[d]) {
            axis[d] = -1;
            touch += 1;
        }
    }
    if (touch == 0 || touch > 2) {
        return 0;
    }

    for (d = 0; d < 3; ++d) {
        if (axis[d] == 1) {
            it_send[d] = PRJ_BLOCK_SIZE;
            it_recv[d] = 0;
        } else if (axis[d] == -1) {
            it_send[d] = 0;
            it_recv[d] = PRJ_BLOCK_SIZE;
        }
    }

    for (dir = 0; dir < 3; ++dir) {
        int tan0;
        int tan1;
        int i_offset;
        int i;

        if (axis[dir] != 0) {
            continue;
        }
        i_offset = ((fine->xmin[dir] + fine->xmax[dir]) <
            (coarse->xmin[dir] + coarse->xmax[dir])) ? 0 : PRJ_BLOCK_SIZE / 2;
        tan0 = (dir + 1) % 3;
        tan1 = (dir + 2) % 3;

        if (axis[tan0] != 0 && axis[tan1] != 0) {
            for (i = 0; i < PRJ_BLOCK_SIZE; i += 2) {
                it_send[dir] = i;
                it_recv[dir] = i / 2 + i_offset;
                if (prj_mpi_store_emf_edge_plan(fine, coarse->id, dir,
                        it_send, it_recv, buffer, count, store) != 0) {
                    return 1;
                }
            }
        } else if (axis[tan0] != 0) {
            int j;
            int j_offset = ((fine->xmin[tan1] + fine->xmax[tan1]) <
                (coarse->xmin[tan1] + coarse->xmax[tan1])) ? 0 : PRJ_BLOCK_SIZE / 2;

            for (i = 0; i < PRJ_BLOCK_SIZE; i += 2) {
                for (j = 0; j <= PRJ_BLOCK_SIZE; j += 2) {
                    it_send[tan1] = j;
                    it_recv[tan1] = j / 2 + j_offset;
                    it_send[dir] = i;
                    it_recv[dir] = i / 2 + i_offset;
                    if (prj_mpi_store_emf_edge_plan(fine, coarse->id, dir,
                            it_send, it_recv, buffer, count, store) != 0) {
                        return 1;
                    }
                }
            }
        } else if (axis[tan1] != 0) {
            int j;
            int j_offset = ((fine->xmin[tan0] + fine->xmax[tan0]) <
                (coarse->xmin[tan0] + coarse->xmax[tan0])) ? 0 : PRJ_BLOCK_SIZE / 2;

            for (i = 0; i < PRJ_BLOCK_SIZE; i += 2) {
                for (j = 0; j <= PRJ_BLOCK_SIZE; j += 2) {
                    it_send[tan0] = j;
                    it_recv[tan0] = j / 2 + j_offset;
                    it_send[dir] = i;
                    it_recv[dir] = i / 2 + i_offset;
                    if (prj_mpi_store_emf_edge_plan(fine, coarse->id, dir,
                            it_send, it_recv, buffer, count, store) != 0) {
                        return 1;
                    }
                }
            }
        } else {
            fprintf(stderr, "prj_mpi_pack_emf_restriction: unknown emf edge type\n");
            exit(EXIT_FAILURE);
        }
    }
    return 0;
}

static int prj_mpi_build_emf_plan_pass(prj_mesh *mesh, prj_mpi *mpi,
    prj_mpi_buffer *buffer, int store, int *count)
{
    int bidx;

    *count = 0;
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        int slot;

        if (!prj_mpi_block_is_local(mpi, block)) {
            continue;
        }
        for (slot = 0; slot < 56; ++slot) {
            int nid = block->slot[slot].id;
            const prj_block *coarse;

            if (nid < 0 || nid >= mesh->nblocks ||
                block->slot[slot].rel_level >= 0 ||
                mesh->blocks[nid].rank != buffer->receiver_rank) {
                continue;
            }
            coarse = &mesh->blocks[nid];
            if (prj_mpi_build_emf_restriction_plan(block, coarse, buffer,
                    count, store) != 0) {
                return 1;
            }
        }
    }
    (void)mpi;
    return 0;
}

static int prj_mpi_build_emf_plan_for_neighbor(prj_mesh *mesh, prj_mpi *mpi,
    prj_mpi_buffer *buffer)
{
#if defined(PRJ_ENABLE_MPI)
    int count = 0;
    int axis;

    if (mesh == 0 || mpi == 0 || buffer == 0) {
        return 1;
    }
    if (prj_mpi_build_emf_plan_pass(mesh, mpi, buffer, 0, &count) != 0) {
        return 1;
    }
    buffer->emf_send_count = count;
    MPI_Sendrecv(&buffer->emf_send_count, 1, MPI_INT,
        buffer->receiver_rank, 400, &buffer->emf_recv_count, 1, MPI_INT,
        buffer->receiver_rank, 400, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if (buffer->emf_send_count < 0 || buffer->emf_recv_count < 0) {
        return 1;
    }
    for (axis = 0; axis < 3; ++axis) {
        if (buffer->emf_send_count > 0) {
            buffer->emf_idx_send[axis] = (int *)calloc(
                (size_t)buffer->emf_send_count,
                sizeof(*buffer->emf_idx_send[axis]));
            buffer->emf_src_idx[axis] = (int *)calloc(
                (size_t)buffer->emf_send_count,
                sizeof(*buffer->emf_src_idx[axis]));
        }
        if (buffer->emf_recv_count > 0) {
            buffer->emf_idx_recv[axis] = (int *)calloc(
                (size_t)buffer->emf_recv_count,
                sizeof(*buffer->emf_idx_recv[axis]));
        }
    }
    if (buffer->emf_send_count > 0) {
        buffer->emf_value_send = (double *)calloc(
            (size_t)buffer->emf_send_count, sizeof(*buffer->emf_value_send));
        buffer->emf_src_block = (int *)calloc(
            (size_t)buffer->emf_send_count, sizeof(*buffer->emf_src_block));
        buffer->emf_src_dir = (int *)calloc(
            (size_t)buffer->emf_send_count, sizeof(*buffer->emf_src_dir));
    }
    if (buffer->emf_recv_count > 0) {
        buffer->emf_value_recv = (double *)calloc(
            (size_t)buffer->emf_recv_count, sizeof(*buffer->emf_value_recv));
    }
    for (axis = 0; axis < 3; ++axis) {
        if ((buffer->emf_send_count > 0 &&
                (buffer->emf_idx_send[axis] == 0 ||
                 buffer->emf_src_idx[axis] == 0)) ||
            (buffer->emf_recv_count > 0 &&
                buffer->emf_idx_recv[axis] == 0)) {
            return 1;
        }
    }
    if ((buffer->emf_send_count > 0 &&
            (buffer->emf_value_send == 0 || buffer->emf_src_block == 0 ||
             buffer->emf_src_dir == 0)) ||
        (buffer->emf_recv_count > 0 && buffer->emf_value_recv == 0)) {
        return 1;
    }
    count = 0;
    if (prj_mpi_build_emf_plan_pass(mesh, mpi, buffer, 1, &count) != 0 ||
        count != buffer->emf_send_count) {
        return 1;
    }
    return 0;
#else
    (void)mesh;
    (void)mpi;
    (void)buffer;
    return 0;
#endif
}

static int prj_mpi_pack_emf_values_for_neighbor(prj_mesh *mesh,
    prj_mpi_buffer *buffer)
{
    int i;

    if (mesh == 0 || buffer == 0) {
        return 1;
    }
    for (i = 0; i < buffer->emf_send_count; ++i) {
        int block_id = buffer->emf_src_block[i];
        int dir = buffer->emf_src_dir[i];
        int idx0[3];
        int idx1[3];
        const prj_block *fine;
        double v0;
        double v1;

        if (block_id < 0 || block_id >= mesh->nblocks ||
            dir < 0 || dir >= 3) {
            return 1;
        }
        fine = &mesh->blocks[block_id];
        if (fine->emf[dir] == 0) {
            return 1;
        }
        idx0[0] = buffer->emf_src_idx[0][i];
        idx0[1] = buffer->emf_src_idx[1][i];
        idx0[2] = buffer->emf_src_idx[2][i];
        idx1[0] = idx0[0];
        idx1[1] = idx0[1];
        idx1[2] = idx0[2];
        idx1[dir] += 1;
        if (!prj_mpi_edge_storage_index_ok(dir, idx0[0], idx0[1], idx0[2]) ||
            !prj_mpi_edge_storage_index_ok(dir, idx1[0], idx1[1], idx1[2])) {
            return 1;
        }
        v0 = fine->emf[dir][EDGE_IDX(dir, idx0[0], idx0[1], idx0[2])];
        v1 = fine->emf[dir][EDGE_IDX(dir, idx1[0], idx1[1], idx1[2])];
        if (!isfinite(v0) || !isfinite(v1)) {
            return 1;
        }
        buffer->emf_value_send[i] = 0.5 * (v0 + v1);
    }
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
        int axis;
        int i;

        if (prj_mpi_pack_emf_values_for_neighbor(mesh, buffer) != 0) {
            fprintf(stderr, "prj_mpi_exchange_emf: failed to pack emf records\n");
            exit(EXIT_FAILURE);
        }
        for (axis = 0; axis < 3; ++axis) {
            MPI_Sendrecv(buffer->emf_idx_send[axis], buffer->emf_send_count,
                MPI_INT, buffer->receiver_rank, 401 + axis,
                buffer->emf_idx_recv[axis], buffer->emf_recv_count, MPI_INT,
                buffer->receiver_rank, 401 + axis,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        MPI_Sendrecv(buffer->emf_value_send, buffer->emf_send_count,
            MPI_DOUBLE, buffer->receiver_rank, 404, buffer->emf_value_recv,
            buffer->emf_recv_count, MPI_DOUBLE, buffer->receiver_rank, 404,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for (i = 0; i < buffer->emf_recv_count; ++i) {
            int block_id = buffer->emf_idx_recv[0][i];
            int code = buffer->emf_idx_recv[1][i];
            int dir = buffer->emf_idx_recv[2][i];
            int ii;
            int jj;
            int kk;
            prj_block *block;
            int flat;

            if (block_id < 0 || block_id >= mesh->nblocks || dir < 0 || dir >= 3) {
                continue;
            }
            block = &mesh->blocks[block_id];
            if (!prj_mpi_block_is_local(mpi, block) || block->edge_fidelity[dir] == 0 ||
                block->emf[dir] == 0) {
                continue;
            }
            prj_mpi_decode_cell_index(code, &ii, &jj, &kk);
            if (!prj_mpi_edge_storage_index_ok(dir, ii, jj, kk)) {
                fprintf(stderr, "prj_mpi_exchange_emf: received edge index out of storage\n");
                exit(EXIT_FAILURE);
            }
            flat = EDGE_IDX(dir, ii, jj, kk);
            if (PRJ_MHD_FIDELITY_FINER < block->edge_fidelity[dir][flat]) {
                continue;
            }
            if (!isfinite(buffer->emf_value_recv[i])) {
                fprintf(stderr, "prj_mpi_exchange_emf: received non-finite emf\n");
                exit(EXIT_FAILURE);
            }
            block->emf[dir][flat] = buffer->emf_value_recv[i];
            block->edge_fidelity[dir][flat] = PRJ_MHD_FIDELITY_FINER;
        }
    }
#else
    (void)mesh;
    (void)mpi;
#endif
}
#endif

static int prj_mpi_flux_record_count_for_neighbor(prj_mesh *mesh,
    const prj_mpi *mpi, int receiver_rank)
{
    int count = 0;
    int bidx;

    if (mesh == 0 || mpi == 0) {
        return 0;
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        int n;

        if (block->id < 0 || block->active != 1 || block->rank != mpi->rank) {
            continue;
        }
        for (n = 0; n < 56; ++n) {
            const prj_neighbor *slot = &block->slot[n];
            int nid = slot->id;
            prj_block *neighbor;
            int axis;
            int side;
            int tan0;
            int tan1;

            if (slot->type != PRJ_NEIGHBOR_FACE || slot->rel_level >= 0 ||
                slot->rank != receiver_rank || nid < 0 || nid >= mesh->nblocks) {
                continue;
            }
            neighbor = &mesh->blocks[nid];
            axis = prj_riemann_face_axis(block, neighbor->xmin, neighbor->xmax,
                &side);
            if (axis < 0) {
                continue;
            }
            tan0 = (axis + 1) % 3;
            tan1 = (axis + 2) % 3;
            count += (slot->recv_loc_end[tan0] - slot->recv_loc_start[tan0]) *
                (slot->recv_loc_end[tan1] - slot->recv_loc_start[tan1]);
        }
    }
    return count;
}

static int prj_mpi_build_flux_plan_for_neighbor(prj_mesh *mesh,
    const prj_mpi *mpi, prj_mpi_buffer *buffer)
{
#if defined(PRJ_ENABLE_MPI)
    if (mesh == 0 || mpi == 0 || buffer == 0) {
        return 1;
    }
    buffer->flux_send_count = prj_mpi_flux_record_count_for_neighbor(mesh,
        mpi, buffer->receiver_rank);
    MPI_Sendrecv(&buffer->flux_send_count, 1, MPI_INT,
        buffer->receiver_rank, 200, &buffer->flux_recv_count, 1, MPI_INT,
        buffer->receiver_rank, 200, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if (buffer->flux_send_count < 0 || buffer->flux_recv_count < 0) {
        return 1;
    }
    if (buffer->flux_send_count > 0) {
        buffer->flux_idx_send = (int *)calloc(
            (size_t)5 * (size_t)buffer->flux_send_count,
            sizeof(*buffer->flux_idx_send));
        buffer->flux_value_send = (double *)calloc(
            (size_t)buffer->flux_send_count * (size_t)PRJ_NVAR_CONS,
            sizeof(*buffer->flux_value_send));
    }
    if (buffer->flux_recv_count > 0) {
        buffer->flux_idx_recv = (int *)calloc(
            (size_t)5 * (size_t)buffer->flux_recv_count,
            sizeof(*buffer->flux_idx_recv));
        buffer->flux_value_recv = (double *)calloc(
            (size_t)buffer->flux_recv_count * (size_t)PRJ_NVAR_CONS,
            sizeof(*buffer->flux_value_recv));
    }
    if ((buffer->flux_send_count > 0 &&
            (buffer->flux_idx_send == 0 || buffer->flux_value_send == 0)) ||
        (buffer->flux_recv_count > 0 &&
            (buffer->flux_idx_recv == 0 || buffer->flux_value_recv == 0))) {
        return 1;
    }
    return 0;
#else
    (void)mesh;
    (void)mpi;
    (void)buffer;
    return 0;
#endif
}

static int prj_mpi_pack_flux_records_for_neighbor(prj_mesh *mesh,
    const prj_mpi *mpi, prj_mpi_buffer *buffer)
{
    int bidx;
    int count = 0;

    if (mesh == 0 || mpi == 0 || buffer == 0) {
        return 1;
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        int n;

        if (block->id < 0 || block->active != 1 || block->rank != mpi->rank) {
            continue;
        }
        for (n = 0; n < 56; ++n) {
            const prj_neighbor *slot = &block->slot[n];
            int nid;
            prj_block *neighbor;
            int axis;
            int side;
            int tan0;
            int tan1;
            int it0;
            int it1;

            if (slot->type != PRJ_NEIGHBOR_FACE || slot->rel_level >= 0 ||
                slot->rank != buffer->receiver_rank) {
                continue;
            }
            nid = slot->id;
            if (nid < 0 || nid >= mesh->nblocks) {
                continue;
            }
            neighbor = &mesh->blocks[nid];
            axis = prj_riemann_face_axis(block, neighbor->xmin, neighbor->xmax,
                &side);
            if (axis < 0) {
                continue;
            }
            tan0 = (axis + 1) % 3;
            tan1 = (axis + 2) % 3;
            for (it0 = 0; it0 < slot->recv_loc_end[tan0] - slot->recv_loc_start[tan0]; ++it0) {
                for (it1 = 0; it1 < slot->recv_loc_end[tan1] - slot->recv_loc_start[tan1]; ++it1) {
                    int it_recv[3] = {0, 0, 0};
                    int it_send0[3] = {0, 0, 0};
                    int it_send1[3] = {0, 0, 0};
                    int it_send2[3] = {0, 0, 0};
                    int it_send3[3] = {0, 0, 0};
                    int v;

                    if (count >= buffer->flux_send_count) {
                        return 1;
                    }
                    it_recv[axis] = side == 1 ? 0 : PRJ_BLOCK_SIZE;
                    it_recv[tan0] = it0 + slot->recv_loc_start[tan0];
                    it_recv[tan1] = it1 + slot->recv_loc_start[tan1];

                    it_send0[axis] = side == 1 ? PRJ_BLOCK_SIZE : 0;
                    it_send0[tan0] = 2 * it0 + slot->send_loc_start[tan0];
                    it_send0[tan1] = 2 * it1 + slot->send_loc_start[tan1];
                    it_send1[axis] = it_send0[axis];
                    it_send1[tan0] = it_send0[tan0] + 1;
                    it_send1[tan1] = it_send0[tan1];
                    it_send2[axis] = it_send0[axis];
                    it_send2[tan0] = it_send0[tan0];
                    it_send2[tan1] = it_send0[tan1] + 1;
                    it_send3[axis] = it_send0[axis];
                    it_send3[tan0] = it_send0[tan0] + 1;
                    it_send3[tan1] = it_send0[tan1] + 1;

                    buffer->flux_idx_send[5 * count + 0] = nid;
                    buffer->flux_idx_send[5 * count + 1] = it_recv[0];
                    buffer->flux_idx_send[5 * count + 2] = it_recv[1];
                    buffer->flux_idx_send[5 * count + 3] = it_recv[2];
                    buffer->flux_idx_send[5 * count + 4] = axis;
                    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
                        buffer->flux_value_send[(size_t)count * (size_t)PRJ_NVAR_CONS + (size_t)v] =
                            0.25 * (block->flux[axis][VIDX(v, it_send0[0], it_send0[1], it_send0[2])] +
                                     block->flux[axis][VIDX(v, it_send1[0], it_send1[1], it_send1[2])] +
                                     block->flux[axis][VIDX(v, it_send2[0], it_send2[1], it_send2[2])] +
                                     block->flux[axis][VIDX(v, it_send3[0], it_send3[1], it_send3[2])]);
                    }
                    count += 1;
                }
            }
        }
    }
    return count == buffer->flux_send_count ? 0 : 1;
}

static int prj_mpi_prepare_request_buffer(prj_mpi *mpi)
{
#if defined(PRJ_ENABLE_MPI)
    int capacity;

    if (mpi == 0) {
        return 1;
    }
    if (mpi->neighbor_number <= 0) {
        return 0;
    }
    if (mpi->neighbor_number > INT_MAX / 12) {
        return 1;
    }
    capacity = 12 * mpi->neighbor_number;
    mpi->request_buffer = calloc((size_t)capacity, sizeof(MPI_Request));
    if (mpi->request_buffer == 0) {
        return 1;
    }
    mpi->request_capacity = capacity;
    return 0;
#else
    (void)mpi;
    return 0;
#endif
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
}

void prj_mpi_decompose(prj_mesh *mesh, const prj_mpi *mpi)
{
    prj_mpi_compute_decomposition(mesh, mpi);
    prj_mpi_assign_block_storage(mesh, mpi);
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
    mpi->neighbor_buffer = count > 0 ?
        (prj_mpi_buffer *)calloc((size_t)count, sizeof(*mpi->neighbor_buffer)) : 0;
    if (count > 0 && mpi->neighbor_buffer == 0) {
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
        if (prj_mpi_build_ghost_plan_for_neighbor(mesh, mpi, buffer) != 0 ||
            prj_mpi_build_flux_plan_for_neighbor(mesh, mpi, buffer) != 0
#if PRJ_MHD
            || prj_mpi_build_bf_plan_for_neighbor(mesh, mpi, buffer) != 0
            || prj_mpi_build_emf_plan_for_neighbor(mesh, mpi, buffer) != 0
            || prj_mpi_build_amr_bf_plan_for_neighbor(mesh, mpi, buffer) != 0
#endif
            ) {
            free(rank_seen);
            fprintf(stderr, "prj_mpi_prepare: failed to build communication buffers\n");
            exit(EXIT_FAILURE);
        }
    }
    if (prj_mpi_prepare_request_buffer(mpi) != 0) {
        free(rank_seen);
        fprintf(stderr, "prj_mpi_prepare: failed to allocate MPI request buffer\n");
        exit(EXIT_FAILURE);
    }
#if PRJ_MHD
    if (prj_mpi_build_amr_bf_cache(mpi) != 0) {
        free(rank_seen);
        fprintf(stderr, "prj_mpi_prepare: failed to allocate MPI AMR Bf cache\n");
        exit(EXIT_FAILURE);
    }
#endif
    free(rank_seen);
}

void prj_mpi_exchange_ghosts(prj_mesh *mesh, prj_mpi *mpi, int stage, int fill_kind)
{
#if defined(PRJ_ENABLE_MPI)
    int nb;
    MPI_Request *requests;
    int request_count;

    if (mesh == 0 || mpi == 0 || !prj_mpi_ghost_fill_kind_ok(fill_kind) ||
        mpi->totrank <= 1 ||
        mpi->neighbor_number == 0) {
        return;
    }
    requests = (MPI_Request *)mpi->request_buffer;
    request_count = 0;
    if (requests == 0 || mpi->request_capacity < 2 * mpi->neighbor_number) {
        fprintf(stderr, "prj_mpi_exchange_ghosts: missing MPI request buffer\n");
        exit(EXIT_FAILURE);
    }
    for (nb = 0; nb < mpi->neighbor_number; ++nb) {
        prj_mpi_buffer *buffer = &mpi->neighbor_buffer[nb];
        int record_count = buffer->cell_send_count_by_kind[fill_kind];
        int record_count_rad = buffer->cell_send_count_rad_by_kind[fill_kind];
        int cell_size_total = buffer->cell_recv_count_by_kind[fill_kind];
        int cell_size_total_rad = buffer->cell_recv_count_rad_by_kind[fill_kind];
        int send_total;
        int recv_total;

        send_total = record_count * PRJ_MPI_GHOST_NVAR_HE + record_count_rad * PRJ_MPI_GHOST_NVAR_RAD;
        recv_total = cell_size_total * PRJ_MPI_GHOST_NVAR_HE + cell_size_total_rad * PRJ_MPI_GHOST_NVAR_RAD;
        prj_mpi_pack_ghost_values(mesh, mpi, buffer, stage, fill_kind);

        MPI_Irecv(buffer->cell_buffer_recv_by_kind[fill_kind], recv_total, MPI_DOUBLE, buffer->receiver_rank, 120,
            MPI_COMM_WORLD, &requests[request_count++]);
        MPI_Isend(buffer->cell_buffer_send_by_kind[fill_kind], send_total, MPI_DOUBLE, buffer->receiver_rank, 120,
            MPI_COMM_WORLD, &requests[request_count++]);
    }
    if (request_count > 0) {
        MPI_Waitall(request_count, requests, MPI_STATUSES_IGNORE);
    }
    for (nb = 0; nb < mpi->neighbor_number; ++nb) {
        prj_mpi_buffer *buffer = &mpi->neighbor_buffer[nb];
        prj_mpi_unpack_ghost_values(mesh, mpi, buffer, stage, fill_kind);
    }
#else
    (void)mesh;
    (void)mpi;
    (void)stage;
    (void)fill_kind;
#endif
}

void prj_mpi_exchange_ghosts_and_bf(prj_mesh *mesh, prj_mpi *mpi, int stage, int fill_kind, int use_bf1)
{
#if defined(PRJ_ENABLE_MPI)
    int nb;
    MPI_Request *requests;
    int request_count;

    if (mesh == 0 || mpi == 0 || !prj_mpi_ghost_fill_kind_ok(fill_kind) ||
        mpi->totrank <= 1 || mpi->neighbor_number == 0) {
        return;
    }
    requests = (MPI_Request *)mpi->request_buffer;
    request_count = 0;
    if (requests == 0 || mpi->request_capacity < 6 * mpi->neighbor_number) {
        fprintf(stderr, "prj_mpi_exchange_ghosts_and_bf: missing MPI request buffer\n");
        exit(EXIT_FAILURE);
    }
    for (nb = 0; nb < mpi->neighbor_number; ++nb) {
        prj_mpi_buffer *buffer = &mpi->neighbor_buffer[nb];
        int record_count = buffer->cell_send_count_by_kind[fill_kind];
        int record_count_rad = buffer->cell_send_count_rad_by_kind[fill_kind];
        int cell_size_total = buffer->cell_recv_count_by_kind[fill_kind];
        int cell_size_total_rad = buffer->cell_recv_count_rad_by_kind[fill_kind];
        int send_total = record_count * PRJ_MPI_GHOST_NVAR_HE + record_count_rad * PRJ_MPI_GHOST_NVAR_RAD;
        int recv_total = cell_size_total * PRJ_MPI_GHOST_NVAR_HE + cell_size_total_rad * PRJ_MPI_GHOST_NVAR_RAD;

        if (send_total > 0) {
            prj_mpi_pack_ghost_values(mesh, mpi, buffer, stage, fill_kind);
            MPI_Isend(buffer->cell_buffer_send_by_kind[fill_kind], send_total, MPI_DOUBLE,
                buffer->receiver_rank, 120, MPI_COMM_WORLD, &requests[request_count++]);
        }
        if (recv_total > 0) {
            MPI_Irecv(buffer->cell_buffer_recv_by_kind[fill_kind], recv_total, MPI_DOUBLE,
                buffer->receiver_rank, 120, MPI_COMM_WORLD, &requests[request_count++]);
        }
#if PRJ_MHD
        {
            int expected_send = buffer->bf_send_value_count[fill_kind];
            int recv_value_count = buffer->bf_recv_value_count[fill_kind];

            if (expected_send > 0) {
                int send_value_count = 0;

                if (prj_mpi_pack_bf_values_only(mesh, mpi, buffer->receiver_rank,
                        use_bf1, fill_kind, buffer->bf_values_send,
                        buffer->bf_send_value_capacity, &send_value_count) != 0 ||
                    send_value_count != expected_send) {
                    fprintf(stderr, "prj_mpi_exchange_ghosts_and_bf: failed to pack Bf values\n");
                    exit(EXIT_FAILURE);
                }
                MPI_Isend(buffer->bf_values_send, send_value_count, MPI_DOUBLE,
                    buffer->receiver_rank, 302, MPI_COMM_WORLD, &requests[request_count++]);
            }
            if (recv_value_count > 0) {
                MPI_Irecv(buffer->bf_values_recv, recv_value_count, MPI_DOUBLE,
                    buffer->receiver_rank, 302, MPI_COMM_WORLD, &requests[request_count++]);
            }
        }
#else
        (void)use_bf1;
#endif
    }
    if (request_count > 0) {
        MPI_Waitall(request_count, requests, MPI_STATUSES_IGNORE);
    }
    for (nb = 0; nb < mpi->neighbor_number; ++nb) {
        prj_mpi_buffer *buffer = &mpi->neighbor_buffer[nb];
        int cell_size_total = buffer->cell_recv_count_by_kind[fill_kind];
        int cell_size_total_rad = buffer->cell_recv_count_rad_by_kind[fill_kind];
        int recv_total = cell_size_total * PRJ_MPI_GHOST_NVAR_HE + cell_size_total_rad * PRJ_MPI_GHOST_NVAR_RAD;

        if (recv_total > 0) {
            prj_mpi_unpack_ghost_values(mesh, mpi, buffer, stage, fill_kind);
        }
#if PRJ_MHD
        {
            int recv_value_count = buffer->bf_recv_value_count[fill_kind];

            if (recv_value_count > 0 &&
                prj_mpi_apply_bf_values_from_sender(mesh, mpi,
                    buffer->receiver_rank, use_bf1, fill_kind,
                    buffer->bf_values_recv, recv_value_count) != 0) {
                fprintf(stderr, "prj_mpi_exchange_ghosts_and_bf: failed to apply Bf values\n");
                exit(EXIT_FAILURE);
            }
        }
#endif
    }
#else
    (void)mesh;
    (void)mpi;
    (void)stage;
    (void)fill_kind;
    (void)use_bf1;
#endif
}

void prj_mpi_exchange_fluxes(prj_mesh *mesh, prj_mpi *mpi)
{
#if defined(PRJ_ENABLE_MPI)
    MPI_Request *reqs;
    int req_count;
    int nb;

    if (mesh == 0 || mpi == 0 || mpi->totrank <= 1 || mpi->neighbor_number == 0) {
        return;
    }

    reqs = (MPI_Request *)mpi->request_buffer;
    req_count = 0;
    if (reqs == 0 || mpi->request_capacity < 4 * mpi->neighbor_number) {
        fprintf(stderr, "prj_mpi_exchange_fluxes: missing MPI request buffer\n");
        exit(EXIT_FAILURE);
    }

    for (nb = 0; nb < mpi->neighbor_number; ++nb) {
        prj_mpi_buffer *buffer = &mpi->neighbor_buffer[nb];
        int rk = buffer->receiver_rank;

        if (prj_mpi_pack_flux_records_for_neighbor(mesh, mpi, buffer) != 0) {
            fprintf(stderr, "prj_mpi_exchange_fluxes: failed to pack flux records\n");
            exit(EXIT_FAILURE);
        }
        if (buffer->flux_recv_count > 0) {
            MPI_Irecv(buffer->flux_idx_recv, 5 * buffer->flux_recv_count,
                MPI_INT, rk, 201, MPI_COMM_WORLD, &reqs[req_count++]);
            MPI_Irecv(buffer->flux_value_recv,
                buffer->flux_recv_count * PRJ_NVAR_CONS, MPI_DOUBLE, rk,
                202, MPI_COMM_WORLD, &reqs[req_count++]);
        }
        if (buffer->flux_send_count > 0) {
            MPI_Isend(buffer->flux_idx_send, 5 * buffer->flux_send_count,
                MPI_INT, rk, 201, MPI_COMM_WORLD, &reqs[req_count++]);
            MPI_Isend(buffer->flux_value_send,
                buffer->flux_send_count * PRJ_NVAR_CONS, MPI_DOUBLE, rk,
                202, MPI_COMM_WORLD, &reqs[req_count++]);
        }
    }
    if (req_count > 0) {
        MPI_Waitall(req_count, reqs, MPI_STATUSES_IGNORE);
    }

    for (nb = 0; nb < mpi->neighbor_number; ++nb) {
        prj_mpi_buffer *buffer = &mpi->neighbor_buffer[nb];
        int e;

        for (e = 0; e < buffer->flux_recv_count; ++e) {
            int blk_id = buffer->flux_idx_recv[5*e+0];
            int ci     = buffer->flux_idx_recv[5*e+1];
            int cj     = buffer->flux_idx_recv[5*e+2];
            int ck     = buffer->flux_idx_recv[5*e+3];
            int faxis  = buffer->flux_idx_recv[5*e+4];
            prj_block *coarse;
            int v;

            if (blk_id < 0 || blk_id >= mesh->nblocks) continue;
            if (faxis < 0 || faxis >= 3) continue;
            coarse = &mesh->blocks[blk_id];
            if (coarse->id < 0 || coarse->active != 1 || coarse->rank != mpi->rank) continue;
            if (coarse->flux[faxis] == 0) continue;
            for (v = 0; v < PRJ_NVAR_CONS; ++v)
                coarse->flux[faxis][VIDX(v, ci, cj, ck)] =
                    buffer->flux_value_recv[(size_t)e * PRJ_NVAR_CONS + (size_t)v];
        }
    }
#else
    (void)mesh;
    (void)mpi;
#endif
}

void prj_mpi_exchange_fluxes_and_emf(prj_mesh *mesh, prj_mpi *mpi)
{
#if defined(PRJ_ENABLE_MPI)
    MPI_Request *reqs;
    int req_count;
    int nb;

    if (mesh == 0 || mpi == 0 || mpi->totrank <= 1 || mpi->neighbor_number == 0) {
        return;
    }
    reqs = (MPI_Request *)mpi->request_buffer;
    req_count = 0;
    if (reqs == 0 || mpi->request_capacity < 12 * mpi->neighbor_number) {
        fprintf(stderr, "prj_mpi_exchange_fluxes_and_emf: missing MPI request buffer\n");
        exit(EXIT_FAILURE);
    }


    for (nb = 0; nb < mpi->neighbor_number; ++nb) {
        prj_mpi_buffer *buffer = &mpi->neighbor_buffer[nb];
        int rk = buffer->receiver_rank;

        if (prj_mpi_pack_flux_records_for_neighbor(mesh, mpi, buffer) != 0) {
            fprintf(stderr, "prj_mpi_exchange_fluxes_and_emf: failed to pack flux records\n");
            exit(EXIT_FAILURE);
        }

        if (buffer->flux_recv_count > 0) {
            MPI_Irecv(buffer->flux_idx_recv, 5 * buffer->flux_recv_count,
                MPI_INT, rk, 201, MPI_COMM_WORLD, &reqs[req_count++]);
            MPI_Irecv(buffer->flux_value_recv,
                buffer->flux_recv_count * PRJ_NVAR_CONS, MPI_DOUBLE, rk,
                202, MPI_COMM_WORLD, &reqs[req_count++]);
        }
        if (buffer->flux_send_count > 0) {
            MPI_Isend(buffer->flux_idx_send, 5 * buffer->flux_send_count,
                MPI_INT, rk, 201, MPI_COMM_WORLD, &reqs[req_count++]);
            MPI_Isend(buffer->flux_value_send,
                buffer->flux_send_count * PRJ_NVAR_CONS, MPI_DOUBLE, rk,
                202, MPI_COMM_WORLD, &reqs[req_count++]);
        }

#if PRJ_MHD
        {
            int axis;

            if (prj_mpi_pack_emf_values_for_neighbor(mesh, buffer) != 0) {
                fprintf(stderr, "prj_mpi_exchange_fluxes_and_emf: failed to pack emf records\n");
                exit(EXIT_FAILURE);
            }

            for (axis = 0; axis < 3; ++axis) {
                if (buffer->emf_recv_count > 0) {
                    MPI_Irecv(buffer->emf_idx_recv[axis], buffer->emf_recv_count,
                        MPI_INT, rk, 401 + axis, MPI_COMM_WORLD, &reqs[req_count++]);
                }
                if (buffer->emf_send_count > 0) {
                    MPI_Isend(buffer->emf_idx_send[axis], buffer->emf_send_count,
                        MPI_INT, rk, 401 + axis, MPI_COMM_WORLD, &reqs[req_count++]);
                }
            }
            if (buffer->emf_recv_count > 0) {
                MPI_Irecv(buffer->emf_value_recv, buffer->emf_recv_count,
                    MPI_DOUBLE, rk, 404, MPI_COMM_WORLD, &reqs[req_count++]);
            }
            if (buffer->emf_send_count > 0) {
                MPI_Isend(buffer->emf_value_send, buffer->emf_send_count,
                    MPI_DOUBLE, rk, 404, MPI_COMM_WORLD, &reqs[req_count++]);
            }
        }
#endif
    }

    if (req_count > 0) {
        MPI_Waitall(req_count, reqs, MPI_STATUSES_IGNORE);
    }

    for (nb = 0; nb < mpi->neighbor_number; ++nb) {
        prj_mpi_buffer *buffer = &mpi->neighbor_buffer[nb];
        int e;

        for (e = 0; e < buffer->flux_recv_count; ++e) {
            int blk_id = buffer->flux_idx_recv[5*e+0];
            int ci     = buffer->flux_idx_recv[5*e+1];
            int cj     = buffer->flux_idx_recv[5*e+2];
            int ck     = buffer->flux_idx_recv[5*e+3];
            int faxis  = buffer->flux_idx_recv[5*e+4];
            prj_block *coarse;
            int v;

            if (blk_id < 0 || blk_id >= mesh->nblocks) continue;
            if (faxis < 0 || faxis >= 3) continue;
            coarse = &mesh->blocks[blk_id];
            if (coarse->id < 0 || coarse->active != 1 || coarse->rank != mpi->rank) continue;
            if (coarse->flux[faxis] == 0) continue;
            for (v = 0; v < PRJ_NVAR_CONS; ++v)
                coarse->flux[faxis][VIDX(v, ci, cj, ck)] =
                    buffer->flux_value_recv[(size_t)e * PRJ_NVAR_CONS + (size_t)v];
        }
    }

#if PRJ_MHD
    for (nb = 0; nb < mpi->neighbor_number; ++nb) {
        prj_mpi_buffer *buffer = &mpi->neighbor_buffer[nb];
        int i;

        for (i = 0; i < buffer->emf_recv_count; ++i) {
            int block_id = buffer->emf_idx_recv[0][i];
            int code = buffer->emf_idx_recv[1][i];
            int dir = buffer->emf_idx_recv[2][i];
            int ii;
            int jj;
            int kk;
            prj_block *block;
            int flat;

            if (block_id < 0 || block_id >= mesh->nblocks || dir < 0 || dir >= 3) {
                continue;
            }
            block = &mesh->blocks[block_id];
            if (!prj_mpi_block_is_local(mpi, block) || block->edge_fidelity[dir] == 0 ||
                block->emf[dir] == 0) {
                continue;
            }
            prj_mpi_decode_cell_index(code, &ii, &jj, &kk);
            if (!prj_mpi_edge_storage_index_ok(dir, ii, jj, kk)) {
                fprintf(stderr, "prj_mpi_exchange_fluxes_and_emf: received edge index out of storage\n");
                exit(EXIT_FAILURE);
            }
            flat = EDGE_IDX(dir, ii, jj, kk);
            if (PRJ_MHD_FIDELITY_FINER < block->edge_fidelity[dir][flat]) {
                continue;
            }
            if (!isfinite(buffer->emf_value_recv[i])) {
                fprintf(stderr, "prj_mpi_exchange_fluxes_and_emf: received non-finite emf\n");
                exit(EXIT_FAILURE);
            }
            block->emf[dir][flat] = buffer->emf_value_recv[i];
            block->edge_fidelity[dir][flat] = PRJ_MHD_FIDELITY_FINER;
        }
    }
#if PRJ_MHD_DEBUG
    prj_mhd_debug_check_emf(mesh, mpi);
#endif
#endif

#else
    (void)mesh;
    (void)mpi;
#endif
}


void prj_mpi_barrier(const prj_mpi *mpi)
{
#if defined(PRJ_ENABLE_MPI)
    if (mpi != 0 && mpi->totrank > 1) {
        MPI_Barrier(MPI_COMM_WORLD);
    }
#else
    (void)mpi;
#endif
}

double prj_mpi_min_dt(const prj_mpi *mpi, double local_dt)
{
#if defined(PRJ_ENABLE_MPI)
    double global_dt = local_dt;

    if (mpi != 0 && mpi->totrank > 1) {
        MPI_Allreduce(&local_dt, &global_dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    }
    return global_dt;
#else
    (void)mpi;
    return local_dt;
#endif
}

double prj_mpi_global_sum(const prj_mpi *mpi, double local_val)
{
#if defined(PRJ_ENABLE_MPI)
    double global_val = local_val;

    if (mpi != 0 && mpi->totrank > 1) {
        MPI_Allreduce(&local_val, &global_val, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    return global_val;
#else
    (void)mpi;
    return local_val;
#endif
}

void prj_mpi_rebalance(prj_mesh *mesh, prj_mpi *mpi)
{
#if defined(PRJ_ENABLE_MPI)
    int *old_ranks;
    int *counts_before;
    int *counts_after;
    int did_rebalance;
    int bidx;

    if (mesh == 0 || mpi == 0) {
        return;
    }
    old_ranks = (int *)calloc((size_t)mesh->nblocks, sizeof(*old_ranks));
    counts_before = (int *)calloc((size_t)mpi->totrank, sizeof(*counts_before));
    counts_after = (int *)calloc((size_t)mpi->totrank, sizeof(*counts_after));
    if (old_ranks == 0 || counts_before == 0 || counts_after == 0) {
        free(counts_after);
        free(counts_before);
        free(old_ranks);
        return;
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        old_ranks[bidx] = mesh->blocks[bidx].rank;
    }
    prj_mpi_collect_active_counts(mesh, mpi, counts_before);
    prj_mpi_compute_decomposition(mesh, mpi);
    did_rebalance = prj_mpi_has_rebalanced(mesh, old_ranks);
    prj_mpi_collect_active_counts(mesh, mpi, counts_after);
    if (did_rebalance) {
        prj_mpi_migrate_active_blocks(mesh, mpi, old_ranks);
    }
    prj_mpi_assign_block_storage(mesh, mpi);
    prj_mpi_prepare(mesh, mpi);
    if (did_rebalance &&
            prj_mpi_counts_changed(counts_before, counts_after, mpi->totrank)) {
        prj_mpi_print_balance(mesh, mpi);
    }
    free(counts_after);
    free(counts_before);
    free(old_ranks);
#else
    if (mesh != 0) {
        prj_mpi_decompose(mesh, mpi);
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

void prj_mpi_finalize(prj_mpi *mpi)
{
    if (mpi != 0) {
        prj_mpi_clear_neighbors(mpi);
    }
#if defined(PRJ_ENABLE_MPI)
    MPI_Finalize();
#endif
}
