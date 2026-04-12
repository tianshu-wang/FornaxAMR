#include <math.h>
#include <stddef.h>
#include <stdlib.h>

#include "prj.h"

#if defined(PRJ_ENABLE_MPI)
#include <mpi.h>
#endif

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

static int prj_is_local_active_block(const prj_block *b)
{
    prj_mpi *mpi = prj_mpi_current();

    return prj_is_active_block(b) && (mpi == 0 || b->rank == mpi->rank);
}

static int prj_is_local_block_owner(const prj_block *b)
{
    prj_mpi *mpi = prj_mpi_current();

    return b != 0 && b->id >= 0 && (mpi == 0 || b->rank == mpi->rank);
}

static void prj_amr_sync_refine_flags(prj_mesh *mesh)
{
#if defined(PRJ_ENABLE_MPI)
    prj_mpi *mpi = prj_mpi_current();
    int *local_pos;
    int *global_pos;
    int *local_neg;
    int *global_neg;
    int i;

    if (mesh == 0 || mpi == 0 || mpi->totrank <= 1) {
        return;
    }
    local_pos = (int *)calloc((size_t)mesh->nblocks, sizeof(*local_pos));
    global_pos = (int *)calloc((size_t)mesh->nblocks, sizeof(*global_pos));
    local_neg = (int *)calloc((size_t)mesh->nblocks, sizeof(*local_neg));
    global_neg = (int *)calloc((size_t)mesh->nblocks, sizeof(*global_neg));
    if (local_pos == 0 || global_pos == 0 || local_neg == 0 || global_neg == 0) {
        free(local_pos);
        free(global_pos);
        free(local_neg);
        free(global_neg);
        return;
    }
    for (i = 0; i < mesh->nblocks; ++i) {
        if (!prj_is_local_active_block(&mesh->blocks[i])) {
            continue;
        }
        if (mesh->blocks[i].refine_flag > 0) {
            local_pos[i] = 1;
        } else if (mesh->blocks[i].refine_flag < 0) {
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
        if (global_pos[i] > 0) {
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
#endif
}

static void prj_zero_block_arrays(prj_block *b)
{
    size_t n;
    size_t total;

    if (b == 0 || b->W == 0) {
        return;
    }

    total = (size_t)2U * (size_t)PRJ_NVAR_PRIM * (size_t)PRJ_BLOCK_NCELLS +
        (size_t)PRJ_NVAR_EOSVAR * (size_t)PRJ_BLOCK_NCELLS +
        (size_t)5U * (size_t)PRJ_NVAR_CONS * (size_t)PRJ_BLOCK_NCELLS;
    for (n = 0; n < total; ++n) {
        b->W[n] = 0.0;
    }
}

static size_t prj_block_data_count(void)
{
    size_t prim_count;
    size_t cons_count;

    prim_count = (size_t)PRJ_NVAR_PRIM * (size_t)PRJ_BLOCK_NCELLS;
    cons_count = (size_t)PRJ_NVAR_CONS * (size_t)PRJ_BLOCK_NCELLS;
    return 2U * prim_count + (size_t)PRJ_NVAR_EOSVAR * (size_t)PRJ_BLOCK_NCELLS + 5U * cons_count;
}

static void prj_amr_move_children_to_parent_rank(prj_mesh *mesh, prj_block *parent)
{
#if defined(PRJ_ENABLE_MPI)
    prj_mpi *mpi;
    size_t data_count;
    int parent_rank;
    int oct;

    if (mesh == 0 || parent == 0) {
        return;
    }
    mpi = prj_mpi_current();
    if (mpi == 0 || mpi->totrank <= 1) {
        return;
    }

    data_count = prj_block_data_count();
    parent_rank = parent->rank;
    for (oct = 0; oct < 8; ++oct) {
        int child_id;
        prj_block *child;
        int source_rank;
        int tag;

        child_id = parent->children[oct];
        if (child_id < 0 || child_id >= mesh->nblocks) {
            continue;
        }
        child = &mesh->blocks[child_id];
        source_rank = child->rank;
        if (source_rank == parent_rank) {
            continue;
        }

        tag = 300 + child_id;
        if (mpi->rank == parent_rank) {
            if (child->W == 0 && prj_block_alloc_data(child) != 0) {
                continue;
            }
            MPI_Recv(child->W, (int)data_count, MPI_DOUBLE, source_rank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else if (mpi->rank == source_rank && child->W != 0) {
            MPI_Send(child->W, (int)data_count, MPI_DOUBLE, parent_rank, tag, MPI_COMM_WORLD);
            prj_block_free_data(child);
        }
        child->rank = parent_rank;
    }
#else
    (void)mesh;
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
        for (d = 0; d < 3; ++d) {
            b->slot[n].xmin[d] = 0.0;
            b->slot[n].xmax[d] = 0.0;
            b->slot[n].dx[d] = 0.0;
        }
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
    return -1;
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
            return 0;
        }
    }
    return 1;
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

static double prj_block_pressure_at(const prj_block *b, prj_eos *eos, int i, int j, int k)
{
    (void)eos;
    i = prj_amr_clamp_storage_index(i);
    j = prj_amr_clamp_storage_index(j);
    k = prj_amr_clamp_storage_index(k);
    return b->eosvar[EIDX(PRJ_EOSVAR_PRESSURE, i, j, k)];
}

static double prj_block_primitive_at(const prj_block *b, int v, int i, int j, int k)
{
    i = prj_amr_clamp_storage_index(i);
    j = prj_amr_clamp_storage_index(j);
    k = prj_amr_clamp_storage_index(k);
    return b->W[VIDX(v, i, j, k)];
}

static double prj_block_sound_speed_at(const prj_block *b, prj_eos *eos, int i, int j, int k)
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
    pressure = b->eosvar[EIDX(PRJ_EOSVAR_PRESSURE, i, j, k)];
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
    return b->U[VIDX(v, i, j, k)];
}

static void prj_apply_eint_floor(double E_floor, double *U, double *W)
{
    double rho;
    double v1;
    double v2;
    double v3;
    double kinetic;

    if (E_floor <= 0.0 || U == 0 || W == 0) {
        return;
    }

    rho = W[PRJ_PRIM_RHO];
    if (rho <= 0.0) {
        return;
    }

    if (W[PRJ_PRIM_EINT] >= E_floor) {
        return;
    }

    v1 = W[PRJ_PRIM_V1];
    v2 = W[PRJ_PRIM_V2];
    v3 = W[PRJ_PRIM_V3];
    kinetic = 0.5 * (v1 * v1 + v2 * v2 + v3 * v3);
    W[PRJ_PRIM_EINT] = E_floor;
    U[PRJ_CONS_ETOT] = rho * (E_floor + kinetic);
}

static double prj_loehner_cell_indicator(const prj_mesh *mesh, const prj_block *b, prj_eos *eos, int i, int j, int k)
{
    const double small = 1.0e-14;
    double eps = mesh != 0 ? mesh->amr_eps : 0.1;
    double p0 = prj_block_pressure_at(b, eos, i, j, k);
    double pxm = prj_block_pressure_at(b, eos, i - 1, j, k);
    double pxp = prj_block_pressure_at(b, eos, i + 1, j, k);
    double pym = prj_block_pressure_at(b, eos, i, j - 1, k);
    double pyp = prj_block_pressure_at(b, eos, i, j + 1, k);
    double pzm = prj_block_pressure_at(b, eos, i, j, k - 1);
    double pzp = prj_block_pressure_at(b, eos, i, j, k + 1);
    double pxpyp = prj_block_pressure_at(b, eos, i + 1, j + 1, k);
    double pxmyp = prj_block_pressure_at(b, eos, i - 1, j + 1, k);
    double pxpym = prj_block_pressure_at(b, eos, i + 1, j - 1, k);
    double pxmym = prj_block_pressure_at(b, eos, i - 1, j - 1, k);
    double pxpzp = prj_block_pressure_at(b, eos, i + 1, j, k + 1);
    double pxmzp = prj_block_pressure_at(b, eos, i - 1, j, k + 1);
    double pxpzm = prj_block_pressure_at(b, eos, i + 1, j, k - 1);
    double pxmzm = prj_block_pressure_at(b, eos, i - 1, j, k - 1);
    double pypzp = prj_block_pressure_at(b, eos, i, j + 1, k + 1);
    double pymzp = prj_block_pressure_at(b, eos, i, j - 1, k + 1);
    double pypzm = prj_block_pressure_at(b, eos, i, j + 1, k - 1);
    double pymzm = prj_block_pressure_at(b, eos, i, j - 1, k - 1);
    double d2xx = pxp - 2.0 * p0 + pxm;
    double d2yy = pyp - 2.0 * p0 + pym;
    double d2zz = pzp - 2.0 * p0 + pzm;
    double d2xy = 0.25 * (pxpyp - pxmyp - pxpym + pxmym);
    double d2xz = 0.25 * (pxpzp - pxmzp - pxpzm + pxmzm);
    double d2yz = 0.25 * (pypzp - pymzp - pypzm + pymzm);
    double grad_x = prj_abs_double(pxp - p0) + prj_abs_double(p0 - pxm);
    double grad_y = prj_abs_double(pyp - p0) + prj_abs_double(p0 - pym);
    double grad_z = prj_abs_double(pzp - p0) + prj_abs_double(p0 - pzm);
    double numerator = d2xx * d2xx + d2yy * d2yy + d2zz * d2zz +
        2.0 * (d2xy * d2xy + d2xz * d2xz + d2yz * d2yz);
    double denominator = grad_x + grad_y + grad_z +
        eps * (prj_abs_double(d2xx) + prj_abs_double(d2yy) + prj_abs_double(d2zz) +
            2.0 * (prj_abs_double(d2xy) + prj_abs_double(d2xz) + prj_abs_double(d2yz))) +
        small;

    return prj_sqrt_double(numerator) / denominator;
}

static double prj_velocity_cell_indicator(const prj_block *b, prj_eos *eos, int i, int j, int k)
{
    static const int offset[6][3] = {
        {-1, 0, 0},
        {1, 0, 0},
        {0, -1, 0},
        {0, 1, 0},
        {0, 0, -1},
        {0, 0, 1}
    };
    const double small = 1.0e-14;
    double v1;
    double v2;
    double v3;
    double cs;
    double max_indicator = 0.0;
    int n;

    v1 = prj_block_primitive_at(b, PRJ_PRIM_V1, i, j, k);
    v2 = prj_block_primitive_at(b, PRJ_PRIM_V2, i, j, k);
    v3 = prj_block_primitive_at(b, PRJ_PRIM_V3, i, j, k);
    cs = prj_block_sound_speed_at(b, eos, i, j, k);
    for (n = 0; n < 6; ++n) {
        double dv1 = v1 - prj_block_primitive_at(b, PRJ_PRIM_V1, i + offset[n][0], j + offset[n][1], k + offset[n][2]);
        double dv2 = v2 - prj_block_primitive_at(b, PRJ_PRIM_V2, i + offset[n][0], j + offset[n][1], k + offset[n][2]);
        double dv3 = v3 - prj_block_primitive_at(b, PRJ_PRIM_V3, i + offset[n][0], j + offset[n][1], k + offset[n][2]);
        double dv = prj_sqrt_double(dv1 * dv1 + dv2 * dv2 + dv3 * dv3) / (cs + small);

        max_indicator = prj_max_double(max_indicator, dv);
    }
    return max_indicator;
}

static double prj_amr_cell_indicator(const prj_mesh *mesh, const prj_block *b, prj_eos *eos, int i, int j, int k)
{
    if (mesh != 0 && mesh->amr_estimator == PRJ_AMR_ESTIMATOR_VELOCITY) {
        return prj_velocity_cell_indicator(b, eos, i, j, k);
    }
    return prj_loehner_cell_indicator(mesh, b, eos, i, j, k);
}

static int prj_has_active_finer_neighbor_than_level(const prj_mesh *mesh, const prj_block *b, int level)
{
    int n;

    for (n = 0; n < 56; ++n) {
        int id = b->slot[n].id;

        if (id >= 0 && id < mesh->nblocks && prj_is_active_block(&mesh->blocks[id]) &&
            mesh->blocks[id].level > level) {
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

static void prj_sync_primitive_from_conserved(prj_mesh *mesh, prj_eos *eos)
{
    int bidx;

    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *b = &mesh->blocks[bidx];
        int i;
        int j;
        int k;

        if (!prj_is_local_active_block(b)) {
            continue;
        }
        for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
            for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    double Uc[PRJ_NVAR_CONS];
                    double Wc[PRJ_NVAR_PRIM];
                    int v;

                    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
                        Uc[v] = b->U[VIDX(v, i, j, k)];
                    }
                    prj_eos_cons2prim(eos, Uc, Wc);
                    prj_apply_eint_floor(mesh->E_floor, Uc, Wc);
                    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
                        b->U[VIDX(v, i, j, k)] = Uc[v];
                    }
                    for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
                        b->W[VIDX(v, i, j, k)] = Wc[v];
                    }
                }
            }
        }
    }
}

static void prj_enforce_two_to_one(prj_mesh *mesh)
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
                    mesh->blocks[id].level < b->level) {
                    if (mesh->blocks[id].refine_flag != 1) {
                        mesh->blocks[id].refine_flag = 1;
                        changed = 1;
                    }
                }
            }
        }
        prj_amr_sync_refine_flags(mesh);
    } while (changed != 0);
}

void prj_amr_init_neighbors(prj_mesh *mesh)
{
    int i;
    int j;

    for (i = 0; i < mesh->nblocks; ++i) {
        if (mesh->blocks[i].id >= 0) {
            prj_clear_neighbors(&mesh->blocks[i]);
        }
    }

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
            if (prj_boxes_touch_face_edge_corner(a, b)) {
                prj_add_neighbor(a, b);
                prj_add_neighbor(b, a);
            }
        }
    }
}

void prj_amr_tag(prj_mesh *mesh, prj_eos *eos)
{
    int i;

    for (i = 0; i < mesh->nblocks; ++i) {
        prj_block *b = &mesh->blocks[i];
        int refine = 0;
        int derefine = b->parent >= 0 ? 1 : 0;
        int j;
        int k;
        int ii;

        if (!prj_is_local_active_block(b)) {
            continue;
        }
        if (b->refine_flag != 0) {
            continue;
        }

        for (ii = 0; ii < PRJ_BLOCK_SIZE && refine == 0; ++ii) {
            for (j = 0; j < PRJ_BLOCK_SIZE && refine == 0; ++j) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    if (prj_amr_cell_indicator(mesh, b, eos, ii, j, k) > mesh->amr_refine_thresh) {
                        refine = 1;
                        break;
                    }
                }
            }
        }

        if (derefine != 0) {
            for (ii = -1; ii <= PRJ_BLOCK_SIZE && derefine != 0; ++ii) {
                for (j = -1; j <= PRJ_BLOCK_SIZE && derefine != 0; ++j) {
                    for (k = -1; k <= PRJ_BLOCK_SIZE; ++k) {
                        if (prj_amr_cell_indicator(mesh, b, eos, ii, j, k) >= mesh->amr_derefine_thresh) {
                            derefine = 0;
                            break;
                        }
                    }
                }
            }
        }

        if (refine != 0 && b->level < mesh->max_level) {
            if (prj_has_face_neighbor_coarser_than(mesh, b, b->level - 1)) {
                b->refine_flag = 0;
            } else {
                b->refine_flag = 1;
            }
        } else if (derefine != 0) {
            b->refine_flag = -1;
        } else {
            b->refine_flag = 0;
        }
    }
    prj_amr_sync_refine_flags(mesh);
}

void prj_amr_prolongate(const prj_block *parent, prj_block *child, int child_oct, double E_floor)
{
    int i;
    int j;
    int k;
    int v;
    int xoct = child_oct & 1;
    int yoct = (child_oct >> 1) & 1;
    int zoct = (child_oct >> 2) & 1;

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
                    double sx;
                    double sy;
                    double sz;
                    double stx[3];
                    double sty[3];
                    double stz[3];
                    double offx;
                    double offy;
                    double offz;
                    double base;

                    stx[0] = prj_block_conserved_at(parent, v, ip - 1, jp, kp);
                    stx[1] = prj_block_conserved_at(parent, v, ip, jp, kp);
                    stx[2] = prj_block_conserved_at(parent, v, ip + 1, jp, kp);
                    sty[0] = prj_block_conserved_at(parent, v, ip, jp - 1, kp);
                    sty[1] = prj_block_conserved_at(parent, v, ip, jp, kp);
                    sty[2] = prj_block_conserved_at(parent, v, ip, jp + 1, kp);
                    stz[0] = prj_block_conserved_at(parent, v, ip, jp, kp - 1);
                    stz[1] = prj_block_conserved_at(parent, v, ip, jp, kp);
                    stz[2] = prj_block_conserved_at(parent, v, ip, jp, kp + 1);
                    sx = prj_reconstruct_slope(stx, parent->dx[0]);
                    sy = prj_reconstruct_slope(sty, parent->dx[1]);
                    sz = prj_reconstruct_slope(stz, parent->dx[2]);
                    offx = ((gi % 2) == 0 ? -0.25 : 0.25) * parent->dx[0];
                    offy = ((gj % 2) == 0 ? -0.25 : 0.25) * parent->dx[1];
                    offz = ((gk % 2) == 0 ? -0.25 : 0.25) * parent->dx[2];
                    base = parent->U[VIDX(v, ip, jp, kp)];
                    child->U[VIDX(v, i, j, k)] = base + sx * offx + sy * offy + sz * offz;
                }
            }
        }
    }
    for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
        for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
            for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                double rho = child->U[VIDX(PRJ_CONS_RHO, i, j, k)];
                double mom1 = child->U[VIDX(PRJ_CONS_MOM1, i, j, k)];
                double mom2 = child->U[VIDX(PRJ_CONS_MOM2, i, j, k)];
                double mom3 = child->U[VIDX(PRJ_CONS_MOM3, i, j, k)];
                double kinetic;
                double min_etot;

                if (rho <= 0.0) {
                    rho = 1.0e-10;
                    child->U[VIDX(PRJ_CONS_RHO, i, j, k)] = rho;
                }
                kinetic = 0.5 * (mom1 * mom1 + mom2 * mom2 + mom3 * mom3) / rho;
                min_etot = child->U[VIDX(PRJ_CONS_ETOT, i, j, k)];
                if (E_floor > 0.0) {
                    min_etot = kinetic + rho * E_floor;
                }
                if (child->U[VIDX(PRJ_CONS_ETOT, i, j, k)] < min_etot) {
                    child->U[VIDX(PRJ_CONS_ETOT, i, j, k)] = min_etot;
                }
            }
        }
    }
}

void prj_amr_restrict(const prj_block *children[8], prj_block *parent)
{
    int v;
    int i;
    int j;
    int k;

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

                                        sum += child->U[VIDX(v, fi, fj, fk)] * overlap_vol;
                                    }
                                }
                            }
                        }
                    }
                    parent->U[VIDX(v, i, j, k)] = sum / parent_cell_vol;
                }
            }
        }
    }
}

void prj_amr_refine_block(prj_mesh *mesh, int block_id)
{
    prj_block *parent;
    int oct;
    double xmid[3];
    int owner_local;

    if (mesh == 0 || block_id < 0 || block_id >= mesh->nblocks) {
        return;
    }
    parent = &mesh->blocks[block_id];
    if (!prj_is_active_block(parent) || parent->level >= mesh->max_level) {
        return;
    }
    owner_local = prj_is_local_block_owner(parent);

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
        if (owner_local) {
            if (prj_block_alloc_data(child) != 0) {
                child->id = -1;
                child->active = 0;
                return;
            }
            prj_block_setup_geometry(child, &mesh->coord);
            prj_amr_prolongate(parent, child, oct, mesh->E_floor);
        } else {
            child->W = 0;
            child->W1 = 0;
            child->eosvar = 0;
            child->U = 0;
            child->dUdt = 0;
            child->flux[0] = 0;
            child->flux[1] = 0;
            child->flux[2] = 0;
            prj_block_setup_geometry(child, &mesh->coord);
        }
        parent->children[oct] = id;
    }

    parent->active = 0;
    parent->refine_flag = 0;
    prj_amr_init_neighbors(mesh);
}

void prj_amr_coarsen_block(prj_mesh *mesh, int parent_id)
{
    prj_block *parent;
    const prj_block *children[8];
    int oct;
    int owner_local;

    if (mesh == 0 || parent_id < 0 || parent_id >= mesh->nblocks) {
        return;
    }
    parent = &mesh->blocks[parent_id];
    owner_local = prj_is_local_block_owner(parent);
    for (oct = 0; oct < 8; ++oct) {
        int id = parent->children[oct];

        if (id < 0 || id >= mesh->nblocks || !prj_is_active_block(&mesh->blocks[id])) {
            return;
        }
        children[oct] = &mesh->blocks[id];
    }

    prj_amr_move_children_to_parent_rank(mesh, parent);
    if (owner_local) {
        prj_amr_restrict(children, parent);
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
        child->parent = -1;
        prj_reset_children(child);
        prj_clear_neighbors(child);
        parent->children[oct] = -1;
    }
    prj_amr_init_neighbors(mesh);
}

void prj_amr_adapt(prj_mesh *mesh, prj_eos *eos)
{
    int i;

    if (mesh == 0) {
        return;
    }

    prj_amr_init_neighbors(mesh);
    prj_eos_fill_mesh(mesh, eos, 1);
    prj_amr_tag(mesh, eos);
    prj_amr_sync_refine_flags(mesh);
    prj_enforce_two_to_one(mesh);
    prj_amr_sync_refine_flags(mesh);

    for (i = 0; i < mesh->nblocks; ++i) {
        if (i < mesh->nblocks && prj_is_active_block(&mesh->blocks[i]) && mesh->blocks[i].refine_flag > 0) {
            prj_amr_refine_block(mesh, i);
        }
    }

    for (i = 0; i < mesh->nblocks; ++i) {
        prj_block *parent = &mesh->blocks[i];
        int can_coarsen;

        if (parent->id < 0 || parent->active == 1) {
            continue;
        }
        can_coarsen = prj_can_coarsen_parent(mesh, i);
        if (can_coarsen && !prj_has_active_finer_neighbor_than_level(mesh, parent, parent->level + 1)) {
            prj_amr_coarsen_block(mesh, i);
        }
    }

    for (i = 0; i < mesh->nblocks; ++i) {
        if (mesh->blocks[i].id >= 0) {
            prj_block_setup_geometry(&mesh->blocks[i], &mesh->coord);
            mesh->blocks[i].refine_flag = 0;
        }
    }
    prj_amr_init_neighbors(mesh);
    prj_sync_primitive_from_conserved(mesh, eos);
    prj_eos_fill_mesh(mesh, eos, 1);
}
