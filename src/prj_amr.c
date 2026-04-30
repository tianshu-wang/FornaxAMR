#include <math.h>
#include <stddef.h>
#include <stdio.h>
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
    const int boundary_buffer = 2;
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
        if (mesh->blocks[i].refine_flag > 0) {
            if (prj_is_active_block(&mesh->blocks[i])) {
                local_pos[i] = 1;
            }
        } else if (mesh->blocks[i].refine_flag < 0 &&
            prj_is_local_active_block(&mesh->blocks[i])) {
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
        (size_t)5U * (size_t)PRJ_NVAR_CONS * (size_t)PRJ_BLOCK_NCELLS +
        9U * (size_t)PRJ_BLOCK_NCELLS
#if PRJ_MHD
        + 15U * (size_t)PRJ_BLOCK_NCELLS
#endif
        ;
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
    return 2U * prim_count + (size_t)PRJ_NVAR_EOSVAR * (size_t)PRJ_BLOCK_NCELLS +
        5U * cons_count + 9U * (size_t)PRJ_BLOCK_NCELLS
#if PRJ_MHD
        + 15U * (size_t)PRJ_BLOCK_NCELLS
#endif
        ;
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
        b->slot[n].rel_level = 0;
        b->slot[n].type = PRJ_NEIGHBOR_NONE;
        for (d = 0; d < 3; ++d) {
            b->slot[n].xmin[d] = 0.0;
            b->slot[n].xmax[d] = 0.0;
            b->slot[n].dx[d] = 0.0;
            b->slot[n].send_loc_start[d] = 0;
            b->slot[n].recv_loc_start[d] = 0;
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
            prj_neighbor_compute_geometry(a, b, &a->slot[n]);
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

    if (W[PRJ_PRIM_EINT] >= E_floor) {
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
    W[PRJ_PRIM_EINT] = E_floor;
    U[PRJ_CONS_ETOT] = rho * (E_floor + kinetic)
#if PRJ_MHD
        + magnetic
#endif
        ;
}

static double prj_loehner_cell_value(const prj_mesh *mesh, const prj_block *b, int lohner_var, int i, int j, int k)
{
    (void)mesh;

    if (b == 0) {
        return 0.0;
    }
    if (lohner_var == PRJ_LOHNER_VAR_LOG_DENSITY) {
        double rho = prj_block_primitive_at(b, PRJ_PRIM_RHO, i, j, k);

        if (rho <= 0.0) {
            return 0.0;
        }
        return log(rho);
    }
    if (lohner_var == PRJ_LOHNER_VAR_DENSITY) {
        return prj_block_primitive_at(b, PRJ_PRIM_RHO, i, j, k);
    }
    if (lohner_var == PRJ_LOHNER_VAR_TEMPERATURE) {
        return b->eosvar[EIDX(PRJ_EOSVAR_TEMPERATURE, i, j, k)];
    }
    return b->eosvar[EIDX(PRJ_EOSVAR_PRESSURE, i, j, k)];
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
    den_x = prj_abs_double(uxp - u0) + prj_abs_double(u0 - uxm) +
        alpha * (prj_abs_double(uxp) + 2.0 * prj_abs_double(u0) + prj_abs_double(uxm));
    den_y = prj_abs_double(uyp - u0) + prj_abs_double(u0 - uym) +
        alpha * (prj_abs_double(uyp) + 2.0 * prj_abs_double(u0) + prj_abs_double(uym));
    den_z = prj_abs_double(uzp - u0) + prj_abs_double(u0 - uzm) +
        alpha * (prj_abs_double(uzp) + 2.0 * prj_abs_double(u0) + prj_abs_double(uzm));
    numerator = num_x * num_x + num_y * num_y + num_z * num_z;
    denominator = den_x * den_x + den_y * den_y + den_z * den_z + small;

    return prj_sqrt_double(numerator / denominator);
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

static double prj_pressure_scale_height_cell_indicator(const prj_block *b, int i, int j, int k)
{
    double rho;
    double pressure;
    double accel;
    double cell_size;
    double Hp;

    if (b == 0) {
        return 0.0;
    }

    rho = prj_block_primitive_at(b, PRJ_PRIM_RHO, i, j, k);
    pressure = b->eosvar[EIDX(PRJ_EOSVAR_PRESSURE, i, j, k)];
    accel = prj_abs_double(prj_gravity_block_accel_at(b, i, j, k));
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

static double prj_density_jump_cell_indicator(const prj_block *b, int i, int j, int k)
{
    const double small = 1.0e-14;
    double r0;
    double max_indicator = 0.0;
    int di;
    int dj;
    int dk;

    if (b == 0 || b->W == 0) {
        return 0.0;
    }

    r0 = b->W[VIDX(PRJ_PRIM_RHO, i, j, k)];
    if (r0 <= 0.0) {
        return 0.0;
    }

    for (di = -1; di <= 1; ++di) {
        for (dj = -1; dj <= 1; ++dj) {
            for (dk = -1; dk <= 1; ++dk) {
                double rnei;
                double denom;
                double jump;

                if (di == 0 && dj == 0 && dk == 0) {
                    continue;
                }
                rnei = b->W[VIDX(PRJ_PRIM_RHO, i + di, j + dj, k + dk)];
                if (rnei <= 0.0) {
                    continue;
                }
                denom = r0 < rnei ? r0 : rnei;
                jump = prj_abs_double(r0 - rnei) / (denom + small);
                max_indicator = prj_max_double(max_indicator, jump);
            }
        }
    }

    return max_indicator;
}

static double prj_amr_cell_indicator_for_estimator(
    const prj_mesh *mesh, const prj_block *b, prj_eos *eos, int amr_idx, int estimator, int i, int j, int k)
{
    if (mesh != 0 && estimator == PRJ_AMR_ESTIMATOR_DENSITY_JUMP) {
        return prj_density_jump_cell_indicator(b, i, j, k);
    }
    if (mesh != 0 && estimator == PRJ_AMR_ESTIMATOR_PRESSURE_SCALE_HEIGHT) {
        return prj_pressure_scale_height_cell_indicator(b, i, j, k);
    }
    if (mesh != 0 && estimator == PRJ_AMR_ESTIMATOR_VELOCITY) {
        return prj_velocity_cell_indicator(b, eos, i, j, k);
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
        if (estimator == PRJ_AMR_ESTIMATOR_LOEHNER &&
            mesh->amr_lohner_var[amr_idx] != PRJ_LOHNER_VAR_DENSITY &&
            mesh->amr_lohner_var[amr_idx] != PRJ_LOHNER_VAR_LOG_DENSITY) {
            return 1;
        }
    }
    return 0;
}

int prj_amr_criteria_need_gravity(const prj_mesh *mesh)
{
    int amr_idx;

    if (mesh == 0) {
        return 0;
    }

    for (amr_idx = 0; amr_idx < PRJ_AMR_N; ++amr_idx) {
        if (mesh->amr_criterion_set[amr_idx] == 0) {
            continue;
        }
        if (mesh->amr_estimator[amr_idx] == PRJ_AMR_ESTIMATOR_PRESSURE_SCALE_HEIGHT) {
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

static void prj_amr_tag_boundary_neighbors(prj_mesh *mesh, const prj_block *block, unsigned int boundary_mask)
{
    int n;

    if (mesh == 0 || block == 0 || boundary_mask == 0U) {
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
        if (mesh->min_dx > 0.0 && prj_block_cell_size(neighbor) < mesh->min_dx) {
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
            neighbor->refine_flag = 1;
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
        if (mesh->blocks[child_id].base_block != 0) {
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

void prj_amr_enforce_two_to_one(prj_mesh *mesh)
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

int prj_amr_refine_marked_blocks(prj_mesh *mesh)
{
    int i;
    int changed = 0;

    if (mesh == 0) {
        return 0;
    }

    for (i = 0; i < mesh->nblocks; ++i) {
        if (i < mesh->nblocks && prj_is_active_block(&mesh->blocks[i]) &&
            mesh->blocks[i].refine_flag > 0) {
            int was_active = mesh->blocks[i].active;

            prj_amr_refine_block(mesh, i);
            if (was_active != mesh->blocks[i].active) {
                changed = 1;
            }
        }
    }
    return changed;
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
        unsigned int boundary_refine_mask = 0U;
        int derefine = b->parent >= 0 ? 1 : 0;
        int j;
        int k;
        int ii;
        int has_refine_criterion = 0;
        int has_derefine_criterion = 0;

        if (!prj_is_local_active_block(b)) {
            continue;
        }
        if (b->refine_flag != 0) {
            continue;
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
                        boundary_refine_mask |= cell_boundary_mask;
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

        if (refine != 0 && (mesh->max_level < 0 || b->level < mesh->max_level) &&
            (mesh->min_dx <= 0.0 || prj_block_cell_size(b) >= mesh->min_dx)) {
            if (prj_has_face_neighbor_coarser_than(mesh, b, b->level - 1)) {
                b->refine_flag = 0;
            } else {
                b->refine_flag = 1;
            }
            prj_amr_tag_boundary_neighbors(mesh, b, boundary_refine_mask);
        } else if (derefine != 0) {
            b->refine_flag = -1;
        } else {
            b->refine_flag = 0;
        }
    }
    prj_amr_sync_refine_flags(mesh);
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

    if (block == 0 || block->face_fidelity == 0 || block->U == 0 || block->W == 0 || block->W1 == 0) {
        prj_amr_mhd_fail(caller);
    }
    for (d = 0; d < 3; ++d) {
        if (block->Bf[d] == 0 || block->Bf1[d] == 0) {
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
    double **bf;
    int d;

    prj_amr_mhd_check_block(block, "prj_amr_mhd_clear_faces: missing MHD storage");
    bf = use_bf1 != 0 ? block->Bf1 : block->Bf;
    for (d = 0; d < 3; ++d) {
        prj_fill(bf[d], (size_t)PRJ_BLOCK_NCELLS, 0.0);
    }
    for (d = 0; d < PRJ_BLOCK_NCELLS; ++d) {
        block->face_fidelity[d] = PRJ_MHD_FIDELITY_NONE;
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
                    block->face_fidelity[IDX(i, j, k)] = fidelity;
                }
            }
        }
    }
}

static void prj_amr_mhd_set_cons_b_from_bf(prj_block *block, int use_bf1)
{
    double **bf;
    int i;
    int j;
    int k;

    prj_amr_mhd_check_block(block, "prj_amr_mhd_set_cons_b_from_bf: missing MHD storage");
    bf = use_bf1 != 0 ? block->Bf1 : block->Bf;
    for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
        for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
            for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                double b1 = 0.5 * (bf[X1DIR][IDX(i, j, k)] + bf[X1DIR][IDX(i + 1, j, k)]);
                double b2 = 0.5 * (bf[X2DIR][IDX(i, j, k)] + bf[X2DIR][IDX(i, j + 1, k)]);
                double b3 = 0.5 * (bf[X3DIR][IDX(i, j, k)] + bf[X3DIR][IDX(i, j, k + 1)]);

                if (!isfinite(b1) || !isfinite(b2) || !isfinite(b3)) {
                    prj_amr_mhd_fail("prj_amr_mhd_set_cons_b_from_bf: non-finite magnetic field");
                }
                block->U[VIDX(PRJ_CONS_B1, i, j, k)] = b1;
                block->U[VIDX(PRJ_CONS_B2, i, j, k)] = b2;
                block->U[VIDX(PRJ_CONS_B3, i, j, k)] = b3;
                block->W[VIDX(PRJ_PRIM_B1, i, j, k)] = b1;
                block->W[VIDX(PRJ_PRIM_B2, i, j, k)] = b2;
                block->W[VIDX(PRJ_PRIM_B3, i, j, k)] = b3;
                block->W1[VIDX(PRJ_PRIM_B1, i, j, k)] = b1;
                block->W1[VIDX(PRJ_PRIM_B2, i, j, k)] = b2;
                block->W1[VIDX(PRJ_PRIM_B3, i, j, k)] = b3;
            }
        }
    }
}

static void prj_amr_mhd_prolongate_bf_one(const prj_block *parent, prj_block *child,
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

    prj_amr_mhd_check_block(parent, "prj_amr_mhd_prolongate_bf_one: missing parent MHD storage");
    prj_amr_mhd_clear_faces(child, use_bf1);
    for (ci = ci0; ci < ci0 + PRJ_BLOCK_SIZE / 2; ++ci) {
        for (cj = cj0; cj < cj0 + PRJ_BLOCK_SIZE / 2; ++cj) {
            for (ck = ck0; ck < ck0 + PRJ_BLOCK_SIZE / 2; ++ck) {
                int fi = 2 * (ci - ci0);
                int fj = 2 * (cj - cj0);
                int fk = 2 * (ck - ck0);

                prj_mhd_bf_prolongate(parent, child, ci, cj, ck, fi, fj, fk, use_bf1);
            }
        }
    }
}

static void prj_amr_mhd_prolongate_bf(const prj_block *parent, prj_block *child, int child_oct)
{
    prj_amr_mhd_prolongate_bf_one(parent, child, child_oct, 0);
    prj_amr_mhd_prolongate_bf_one(parent, child, child_oct, 1);
    prj_amr_mhd_mark_active_faces(child, PRJ_MHD_FIDELITY_COARSER);
    prj_amr_mhd_set_cons_b_from_bf(child, 0);
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
            double *src;
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
            src = use_bf1 != 0 ? child->Bf1[dir] : child->Bf[dir];
            if (src == 0) {
                prj_amr_mhd_fail("prj_amr_mhd_restrict_face_value: missing child Bf storage");
            }
            value = src[IDX(local[0], local[1], local[2])];
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
    double **dst;
    int dir;

    prj_amr_mhd_clear_faces(parent, use_bf1);
    dst = use_bf1 != 0 ? parent->Bf1 : parent->Bf;
    for (dir = 0; dir < 3; ++dir) {
        int i;
        int j;
        int k;

        for (i = 0; i <= prj_amr_mhd_face_axis_max(dir, 0); ++i) {
            for (j = 0; j <= prj_amr_mhd_face_axis_max(dir, 1); ++j) {
                for (k = 0; k <= prj_amr_mhd_face_axis_max(dir, 2); ++k) {
                    dst[dir][IDX(i, j, k)] = prj_amr_mhd_restrict_face_value(
                        children, dir, i, j, k, use_bf1);
                }
            }
        }
    }
}

static void prj_amr_mhd_restrict_bf(const prj_block *children[8], prj_block *parent)
{
    int oct;

    prj_amr_mhd_check_block(parent, "prj_amr_mhd_restrict_bf: missing parent MHD storage");
    for (oct = 0; oct < 8; ++oct) {
        prj_amr_mhd_check_block(children[oct], "prj_amr_mhd_restrict_bf: missing child MHD storage");
    }
    prj_amr_mhd_restrict_bf_one(children, parent, 0);
    prj_amr_mhd_restrict_bf_one(children, parent, 1);
    prj_amr_mhd_mark_active_faces(parent, PRJ_MHD_FIDELITY_FINER);
    prj_amr_mhd_set_cons_b_from_bf(parent, 0);
}
#endif

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
                double magnetic = 0.0;
                double min_etot;

                if (rho <= 0.0) {
                    rho = 1.0e-10;
                    child->U[VIDX(PRJ_CONS_RHO, i, j, k)] = rho;
                }
                kinetic = 0.5 * (mom1 * mom1 + mom2 * mom2 + mom3 * mom3) / rho;
#if PRJ_MHD
                magnetic = 0.5 * (child->U[VIDX(PRJ_CONS_B1, i, j, k)] *
                    child->U[VIDX(PRJ_CONS_B1, i, j, k)] +
                    child->U[VIDX(PRJ_CONS_B2, i, j, k)] *
                    child->U[VIDX(PRJ_CONS_B2, i, j, k)] +
                    child->U[VIDX(PRJ_CONS_B3, i, j, k)] *
                    child->U[VIDX(PRJ_CONS_B3, i, j, k)]);
#endif
                min_etot = child->U[VIDX(PRJ_CONS_ETOT, i, j, k)];
                if (E_floor > 0.0) {
                    min_etot = kinetic + rho * E_floor + magnetic;
                }
                if (child->U[VIDX(PRJ_CONS_ETOT, i, j, k)] < min_etot) {
                    child->U[VIDX(PRJ_CONS_ETOT, i, j, k)] = min_etot;
                }
            }
        }
    }
#if PRJ_MHD
    prj_amr_mhd_prolongate_bf(parent, child, child_oct);
#endif
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
#if PRJ_MHD
    prj_amr_mhd_restrict_bf(children, parent);
#endif
}

static int prj_amr_all_cells_meet_angle_resolution_limit(const prj_mesh *mesh, const prj_block *block)
{
    int i;
    int j;
    int k;

    if (mesh == 0 || block == 0 || mesh->use_amr_angle_resolution == 0 ||
        mesh->amr_angle_resolution_limit <= 0.0) {
        return 0;
    }

    double cell_size = prj_block_cell_size(block);

    for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
        for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
            for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                double x = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
                double y = block->xmin[1] + ((double)j + 0.5) * block->dx[1];
                double z = block->xmin[2] + ((double)k + 0.5) * block->dx[2];
                double radius = sqrt(x * x + y * y + z * z);

                if (radius <= 0.0) {
                    return 0;
                }
                if (cell_size / radius >= mesh->amr_angle_resolution_limit) {
                    return 0;
                }
            }
        }
    }

    return 1;
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
    if (!prj_is_active_block(parent) ||
        (mesh->max_level >= 0 && parent->level >= mesh->max_level)) {
        return;
    }
    if (mesh->min_dx > 0.0 && prj_block_cell_size(parent) < mesh->min_dx) {
        parent->refine_flag = 0;
        return;
    }
    if (prj_amr_all_cells_meet_angle_resolution_limit(mesh, parent)) {
        parent->refine_flag = 0;
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
        child->base_block = 0;
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
            child->v_riemann[0] = 0;
            child->v_riemann[1] = 0;
            child->v_riemann[2] = 0;
            child->ridx = 0;
            child->fr = 0;
            prj_block_setup_geometry(child, &mesh->coord);
        }
        parent->children[oct] = id;
    }

    parent->active = 0;
    parent->refine_flag = 0;
}

int prj_amr_coarsen_block(prj_mesh *mesh, int parent_id)
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
    owner_local = prj_is_local_block_owner(parent);

    prj_amr_move_children_to_parent_rank(mesh, parent);
    if (owner_local) {
        if (parent->W == 0) {
            prj_block_alloc_data(parent);
            prj_block_setup_geometry(parent, &mesh->coord);
        }
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
        child->base_block = 0;
        child->parent = -1;
        prj_reset_children(child);
        prj_clear_neighbors(child);
        parent->children[oct] = -1;
    }
    return 1;
}

void prj_amr_adapt(prj_mesh *mesh, prj_eos *eos)
{
    int i;
    int refined = 0;
    int coarsened = 0;

    if (mesh == 0) {
        return;
    }

    prj_amr_tag(mesh, eos);
    prj_amr_sync_refine_flags(mesh);
    prj_amr_enforce_two_to_one(mesh);
    prj_amr_sync_refine_flags(mesh);

    refined = prj_amr_refine_marked_blocks(mesh);
    if (refined) {
        prj_amr_init_neighbors(mesh);
    }

    for (i = 0; i < mesh->nblocks; ++i) {
        prj_block *parent = &mesh->blocks[i];
        int can_coarsen;

        if (parent->id < 0 || parent->active == 1) {
            continue;
        }
        can_coarsen = prj_can_coarsen_parent(mesh, i);
        if (can_coarsen) {
            if (prj_amr_coarsen_block(mesh, i)) {
                coarsened = 1;
            }
        }
    }
    if (coarsened) {
        prj_amr_init_neighbors(mesh);
    }

    if (!refined && !coarsened) {
        for (i = 0; i < mesh->nblocks; ++i) {
            if (mesh->blocks[i].id >= 0) {
                mesh->blocks[i].refine_flag = 0;
            }
        }
        return;
    }

    for (i = 0; i < mesh->nblocks; ++i) {
        if (mesh->blocks[i].id >= 0) {
            prj_block_setup_geometry(&mesh->blocks[i], &mesh->coord);
            mesh->blocks[i].refine_flag = 0;
        }
    }
    prj_mesh_update_max_active_level(mesh);
    prj_sync_primitive_from_conserved(mesh, eos);
}
