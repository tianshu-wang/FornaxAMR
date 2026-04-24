#include <stdio.h>
#include <stdlib.h>

#include "prj.h"

static prj_mesh *prj_boundary_active_mesh = 0;

static int prj_boundary_idx_outside(int idx);

static double prj_abs_double(double x)
{
    return x < 0.0 ? -x : x;
}

#if PRJ_MHD
static int prj_min_int(int a, int b)
{
    return a < b ? a : b;
}

static int prj_max_int(int a, int b)
{
    return a > b ? a : b;
}

static int prj_boundary_valid_restrict_count(int count)
{
    return count == 0 || count == 4 || count == 8;
}
#endif

static int prj_floor_to_int(double x)
{
    int i = (int)x;

    if ((double)i > x) {
        i -= 1;
    }
    return i;
}

static double *prj_boundary_stage_array(prj_block *block, int stage)
{
    return stage == 2 ? block->W1 : block->W;
}

static const double *prj_boundary_stage_array_const(const prj_block *block, int stage)
{
    return stage == 2 ? block->W1 : block->W;
}

static int prj_boundary_active_block(const prj_block *block)
{
    prj_mpi *mpi = prj_mpi_current();

    return block != 0 && block->id >= 0 && block->active == 1 &&
        (mpi == 0 || block->rank == mpi->rank);
}

static int prj_boundary_point_inside(const prj_block *block, double x1, double x2, double x3)
{
    const double tol = 1.0e-12;

    return x1 > block->xmin[0] - tol && x1 < block->xmax[0] + tol &&
        x2 > block->xmin[1] - tol && x2 < block->xmax[1] + tol &&
        x3 > block->xmin[2] - tol && x3 < block->xmax[2] + tol;
}

static int prj_boundary_fraction_case(double frac)
{
    const double tol = 1.0e-2;

    if (prj_abs_double(frac - 0.0) < tol || prj_abs_double(frac - 1.0) < tol) {
        return 1;
    }
    if (prj_abs_double(frac - 0.5) < tol) {
        return 2;
    }
    if (prj_abs_double(frac - 0.25) < tol || prj_abs_double(frac - 0.75) < tol) {
        return 3;
    }
    return 0;
}

enum {
    PRJ_BOUNDARY_FILL_NONRECON = 0,
    PRJ_BOUNDARY_FILL_RECON = 1,
    PRJ_BOUNDARY_FILL_ALL = 2
};

enum {
    PRJ_BOUNDARY_PHYS_FACE_ONLY = 0,
    PRJ_BOUNDARY_PHYS_ALL = 1
};

static int prj_boundary_sample_kind(const prj_block *block, double x1, double x2, double x3)
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
    frac[0] = ox[0] - (double)prj_floor_to_int(ox[0]);
    frac[1] = ox[1] - (double)prj_floor_to_int(ox[1]);
    frac[2] = ox[2] - (double)prj_floor_to_int(ox[2]);
    if (frac[0] < 0.0) {
        frac[0] += 1.0;
    }
    if (frac[1] < 0.0) {
        frac[1] += 1.0;
    }
    if (frac[2] < 0.0) {
        frac[2] += 1.0;
    }

    cases[0] = prj_boundary_fraction_case(frac[0]);
    cases[1] = prj_boundary_fraction_case(frac[1]);
    cases[2] = prj_boundary_fraction_case(frac[2]);
    if (cases[0] == 0 || cases[1] == 0 || cases[2] == 0) {
        return -1;
    }
    if (cases[0] == 1 && cases[1] == 1 && cases[2] == 1) {
        return PRJ_BOUNDARY_FILL_NONRECON;
    }
    if (cases[0] == 2 && cases[1] == 2 && cases[2] == 2) {
        return PRJ_BOUNDARY_FILL_NONRECON;
    }
    if (cases[0] == 3 && cases[1] == 3 && cases[2] == 3) {
        return PRJ_BOUNDARY_FILL_RECON;
    }
    fprintf(stderr,
        "prj_boundary_sample_kind: invalid fraction case tuple (%d,%d,%d) "
        "at (x1=%g, x2=%g, x3=%g) -- violates 2:1 AMR constraint\n",
        cases[0], cases[1], cases[2], x1, x2, x3);
    abort();
}

static double prj_boundary_read_value(const double *src, int var, int i, int j, int k, int is_eosvar)
{
    if (i < -PRJ_NGHOST || i >= PRJ_BLOCK_SIZE + PRJ_NGHOST ||
        j < -PRJ_NGHOST || j >= PRJ_BLOCK_SIZE + PRJ_NGHOST ||
        k < -PRJ_NGHOST || k >= PRJ_BLOCK_SIZE + PRJ_NGHOST) {
        fprintf(stderr,
            "prj_boundary_read_value: out-of-range access var=%d i=%d j=%d k=%d "
            "(valid [%d, %d])\n",
            var, i, j, k, -PRJ_NGHOST, PRJ_BLOCK_SIZE + PRJ_NGHOST - 1);
        exit(EXIT_FAILURE);
    }
    return src[is_eosvar != 0 ? EIDX(var, i, j, k) : VIDX(var, i, j, k)];
}

static void prj_boundary_get_values(const prj_block *block, const double *src, int nvar, int is_eosvar,
    const char *label, double x1, double x2, double x3, double *dst)
{
    double ox[3];
    double frac[3];
    int cases[3];
    int v;

    ox[0] = (x1 - block->xmin[0]) / block->dx[0] - 0.5;
    ox[1] = (x2 - block->xmin[1]) / block->dx[1] - 0.5;
    ox[2] = (x3 - block->xmin[2]) / block->dx[2] - 0.5;
    frac[0] = ox[0] - (double)prj_floor_to_int(ox[0]);
    frac[1] = ox[1] - (double)prj_floor_to_int(ox[1]);
    frac[2] = ox[2] - (double)prj_floor_to_int(ox[2]);
    if (frac[0] < 0.0) {
        frac[0] += 1.0;
    }
    if (frac[1] < 0.0) {
        frac[1] += 1.0;
    }
    if (frac[2] < 0.0) {
        frac[2] += 1.0;
    }

    cases[0] = prj_boundary_fraction_case(frac[0]);
    cases[1] = prj_boundary_fraction_case(frac[1]);
    cases[2] = prj_boundary_fraction_case(frac[2]);
    if (cases[0] == 0 || cases[1] == 0 || cases[2] == 0) {
        fprintf(stderr, "%s: unsupported sample location (%g, %g, %g)\n", label, x1, x2, x3);
        exit(EXIT_FAILURE);
    }

    for (v = 0; v < nvar; ++v) {
        if (cases[0] == 1 && cases[1] == 1 && cases[2] == 1) {
            int i = (int)(ox[0] >= 0.0 ? ox[0] + 0.5 : ox[0] - 0.5);
            int j = (int)(ox[1] >= 0.0 ? ox[1] + 0.5 : ox[1] - 0.5);
            int k = (int)(ox[2] >= 0.0 ? ox[2] + 0.5 : ox[2] - 0.5);

            dst[v] = prj_boundary_read_value(src, v, i, j, k, is_eosvar);
        } else if (cases[0] == 2 && cases[1] == 2 && cases[2] == 2) {
            int i = prj_floor_to_int(ox[0]);
            int j = prj_floor_to_int(ox[1]);
            int k = prj_floor_to_int(ox[2]);
            int di;
            int dj;
            int dk;
            double sum = 0.0;

            for (di = 0; di < 2; ++di) {
                for (dj = 0; dj < 2; ++dj) {
                    for (dk = 0; dk < 2; ++dk) {
                        sum += prj_boundary_read_value(src, v, i + di, j + dj, k + dk, is_eosvar);
                    }
                }
            }
            dst[v] = 0.125 * sum;
        } else if ((cases[0] == 1 || cases[0] == 3) &&
                   (cases[1] == 1 || cases[1] == 3) &&
                   (cases[2] == 1 || cases[2] == 3)) {
            int i = prj_floor_to_int(ox[0] + 0.5);
            int j = prj_floor_to_int(ox[1] + 0.5);
            int k = prj_floor_to_int(ox[2] + 0.5);
            double stx[3];
            double sty[3];
            double stz[3];
            double sx;
            double sy;
            double sz;
            double xcenter = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
            double ycenter = block->xmin[1] + ((double)j + 0.5) * block->dx[1];
            double zcenter = block->xmin[2] + ((double)k + 0.5) * block->dx[2];
            double base = prj_boundary_read_value(src, v, i, j, k, is_eosvar);

            stx[0] = prj_boundary_read_value(src, v, i - 1, j, k, is_eosvar);
            stx[1] = base;
            stx[2] = prj_boundary_read_value(src, v, i + 1, j, k, is_eosvar);
            sty[0] = prj_boundary_read_value(src, v, i, j - 1, k, is_eosvar);
            sty[1] = base;
            sty[2] = prj_boundary_read_value(src, v, i, j + 1, k, is_eosvar);
            stz[0] = prj_boundary_read_value(src, v, i, j, k - 1, is_eosvar);
            stz[1] = base;
            stz[2] = prj_boundary_read_value(src, v, i, j, k + 1, is_eosvar);
            sx = prj_reconstruct_slope(stx, block->dx[0]);
            sy = prj_reconstruct_slope(sty, block->dx[1]);
            sz = prj_reconstruct_slope(stz, block->dx[2]);
            dst[v] = base +
                sx * (x1 - xcenter) +
                sy * (x2 - ycenter) +
                sz * (x3 - zcenter);
        } else {
            fprintf(stderr, "%s: mixed unsupported sample location (%g, %g, %g)\n", label, x1, x2, x3);
            exit(EXIT_FAILURE);
        }
    }
}

void prj_boundary_get_prim(const prj_block *block, int stage, double x1, double x2, double x3, double *w)
{
    prj_boundary_get_values(block, prj_boundary_stage_array_const(block, stage),
        PRJ_NVAR_PRIM, 0, "prj_boundary_get_prim", x1, x2, x3, w);
}

void prj_boundary_get_eosvar(const prj_block *block, double x1, double x2, double x3, double *eosv)
{
    prj_boundary_get_values(block, block != 0 ? block->eosvar : 0,
        PRJ_NVAR_EOSVAR, 1, "prj_boundary_get_eosvar", x1, x2, x3, eosv);
}

void prj_boundary_send(prj_block *block, int stage, int fill_kind)
{
    int n;

    for (n = 0; n < 56; ++n) {
        int id = block->slot[n].id;

        if (id >= 0 && prj_boundary_active_mesh != 0 && id < prj_boundary_active_mesh->nblocks) {
            prj_block *neighbor = &prj_boundary_active_mesh->blocks[id];
            int i;
            int j;
            int k;

            if (!prj_boundary_active_block(neighbor)) {
                continue;
            }
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

                        if (prj_boundary_point_inside(block, x1, x2, x3)) {
                            double w[PRJ_NVAR_PRIM];
                            int sample_kind;
                            int v;

                            sample_kind = prj_boundary_sample_kind(block, x1, x2, x3);
                            if (sample_kind < 0) {
                                fprintf(stderr,
                                    "prj_boundary_send: failed to classify ghost sample "
                                    "src_block=%d dst_block=%d i=%d j=%d k=%d x=(%.17g, %.17g, %.17g)\n",
                                    block->id, neighbor->id, i, j, k, x1, x2, x3);
                                exit(EXIT_FAILURE);
                            }
                            if (fill_kind != PRJ_BOUNDARY_FILL_ALL && sample_kind != fill_kind) {
                                continue;
                            }
                            if (neighbor->rank == block->rank) {
                                double eosv[PRJ_NVAR_EOSVAR];
                                double *dst = prj_boundary_stage_array(neighbor, stage);
                                int same_level = block->level == neighbor->level;

                                prj_boundary_get_prim(block, stage, x1, x2, x3, w);
                                prj_boundary_get_eosvar(block, x1, x2, x3, eosv);
                                for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
                                    dst[VIDX(v, i, j, k)] = w[v];
                                }
                                if (neighbor->eosvar != 0) {
                                    for (v = 0; v < PRJ_NVAR_EOSVAR; ++v) {
                                        neighbor->eosvar[EIDX(v, i, j, k)] = eosv[v];
                                    }
                                }
                                if (neighbor->eos_done != 0) {
                                    neighbor->eos_done[IDX(i, j, k)] = same_level ? 1 : 0;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

#if PRJ_MHD
static void prj_boundary_coarse_subface_position(const prj_block *block, int dir, int i, int j, int k,
    int aidx, int bidx, double x[3])
{
    if (block == 0 || x == 0) {
        fprintf(stderr, "prj_boundary_coarse_subface_position: null input\n");
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
    fprintf(stderr, "prj_boundary_coarse_subface_position: invalid direction %d\n", dir);
    abort();
}

static int prj_boundary_axis_touches(const prj_block *block, const prj_block *neighbor, int axis)
{
    const double tol = 1.0e-12;

    if (block == 0 || neighbor == 0 || axis < 0 || axis >= 3) {
        fprintf(stderr, "prj_boundary_axis_touches: invalid input\n");
        abort();
    }
    return prj_abs_double(block->xmax[axis] - neighbor->xmin[axis]) < tol ||
        prj_abs_double(neighbor->xmax[axis] - block->xmin[axis]) < tol;
}

static int prj_boundary_axis_side(const prj_block *block, const prj_block *neighbor, int axis)
{
    const double tol = 1.0e-12;

    if (block == 0 || neighbor == 0 || axis < 0 || axis >= 3) {
        fprintf(stderr, "prj_boundary_axis_side: invalid input\n");
        abort();
    }
    if (prj_abs_double(block->xmax[axis] - neighbor->xmin[axis]) < tol) {
        return 1;
    }
    if (prj_abs_double(neighbor->xmax[axis] - block->xmin[axis]) < tol) {
        return -1;
    }
    return 0;
}

static void prj_boundary_same_level_face_bounds(const prj_block *block, const prj_block *neighbor,
    const prj_neighbor *slot, int dir, int use_send_start, int start[3], int end[3])
{
    int axis;

    if (block == 0 || neighbor == 0 || slot == 0 || start == 0 || end == 0) {
        fprintf(stderr, "prj_boundary_same_level_face_bounds: invalid input\n");
        abort();
    }
    for (axis = 0; axis < 3; ++axis) {
        int anchor = use_send_start != 0 ? slot->send_loc_start[axis] : slot->recv_loc_start[axis];

        if (prj_boundary_axis_touches(block, neighbor, axis)) {
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

static void prj_boundary_coarse_face_bounds(const prj_block *block, const prj_block *neighbor,
    const prj_neighbor *slot, int dir, int use_send_start, int start[3], int end[3])
{
    int axis;

    if (block == 0 || neighbor == 0 || slot == 0 || start == 0 || end == 0) {
        fprintf(stderr, "prj_boundary_coarse_face_bounds: invalid input\n");
        abort();
    }
    for (axis = 0; axis < 3; ++axis) {
        int side = prj_boundary_axis_side(block, neighbor, axis);
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

static void prj_boundary_patch_bounds(const prj_block *block, const prj_block *neighbor,
    const prj_neighbor *slot, int use_send_start, int start[3], int end[3])
{
    int axis;

    if (block == 0 || neighbor == 0 || slot == 0 || start == 0 || end == 0) {
        fprintf(stderr, "prj_boundary_patch_bounds: invalid input\n");
        abort();
    }
    for (axis = 0; axis < 3; ++axis) {
        int anchor = use_send_start != 0 ? slot->send_loc_start[axis] : slot->recv_loc_start[axis];

        if (prj_boundary_axis_touches(block, neighbor, axis)) {
            start[axis] = anchor == 0 ? -1 : (PRJ_BLOCK_SIZE - 1);
            end[axis] = anchor == 0 ? 0 : PRJ_BLOCK_SIZE;
        } else {
            start[axis] = prj_max_int(-PRJ_NGHOST + 1, anchor - 1);
            end[axis] = prj_min_int(PRJ_BLOCK_SIZE + PRJ_NGHOST - 2, anchor + PRJ_BLOCK_SIZE / 2);
        }
    }
}

static void prj_boundary_apply_axis_bf(prj_block *block, int stage, int axis, int side, int face_only)
{
    double *Bface[3];
    int dir;

    if (block == 0) {
        fprintf(stderr, "prj_boundary_apply_axis_bf: null block\n");
        abort();
    }
    Bface[0] = prj_mhd_bface_stage(block, stage, X1DIR);
    Bface[1] = prj_mhd_bface_stage(block, stage, X2DIR);
    Bface[2] = prj_mhd_bface_stage(block, stage, X3DIR);
    if (Bface[0] == 0 || Bface[1] == 0 || Bface[2] == 0) {
        fprintf(stderr,
            "prj_boundary_apply_axis_bf: missing face storage block=%d rank=%d stage=%d "
            "Bf=(%p,%p,%p)\n",
            block->id, block->rank, stage, (void *)Bface[0], (void *)Bface[1], (void *)Bface[2]);
        abort();
    }
    for (dir = 0; dir < 3; ++dir) {
        int i;
        int j;
        int k;

        for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
            for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
                for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                    int idx[3];
                    int src_idx[3];

                    idx[0] = i;
                    idx[1] = j;
                    idx[2] = k;
                    if (prj_mhd_face_is_interior(dir, idx[0], idx[1], idx[2])) {
                        continue;
                    }
                    if (face_only != 0) {
                        int axis1 = (axis + 1) % 3;
                        int axis2 = (axis + 2) % 3;

                        if (prj_boundary_idx_outside(idx[axis1]) ||
                            prj_boundary_idx_outside(idx[axis2])) {
                            continue;
                        }
                    }
                    if (side == 0) {
                        if (idx[axis] >= 0) {
                            continue;
                        }
                        src_idx[axis] = axis == dir ? 0 : 0;
                    } else {
                        int outer = axis == dir ? PRJ_BLOCK_SIZE : PRJ_BLOCK_SIZE - 1;

                        if (axis == dir) {
                            if (idx[axis] <= PRJ_BLOCK_SIZE) {
                                continue;
                            }
                        } else if (idx[axis] < PRJ_BLOCK_SIZE) {
                            continue;
                        }
                        src_idx[axis] = outer;
                    }
                    src_idx[(axis + 1) % 3] = idx[(axis + 1) % 3];
                    src_idx[(axis + 2) % 3] = idx[(axis + 2) % 3];
                    Bface[dir][IDX(idx[0], idx[1], idx[2])] = Bface[dir][IDX(src_idx[0], src_idx[1], src_idx[2])];
                }
            }
        }
    }
}

static void prj_boundary_physical_bf(const prj_bc *bc, prj_block *block, int stage, int mode)
{
    const double tol = 1.0e-12;
    int pass;
    int npass = mode == PRJ_BOUNDARY_PHYS_FACE_ONLY ? 1 : 3;
    int face_only = mode == PRJ_BOUNDARY_PHYS_FACE_ONLY ? 1 : 0;

    if (bc == 0 || block == 0 || prj_boundary_active_mesh == 0) {
        fprintf(stderr, "prj_boundary_physical_bf: invalid input\n");
        abort();
    }

    for (pass = 0; pass < npass; ++pass) {
        if (prj_abs_double(block->xmin[0] - prj_boundary_active_mesh->coord.x1min) < tol) {
            if (bc->bc_x1_inner != PRJ_BC_OUTFLOW) {
                fprintf(stderr, "prj_boundary_physical_bf: only outflow BC is implemented for MHD\n");
                abort();
            }
            prj_boundary_apply_axis_bf(block, stage, 0, 0, face_only);
        }
        if (prj_abs_double(block->xmax[0] - prj_boundary_active_mesh->coord.x1max) < tol) {
            if (bc->bc_x1_outer != PRJ_BC_OUTFLOW) {
                fprintf(stderr, "prj_boundary_physical_bf: only outflow BC is implemented for MHD\n");
                abort();
            }
            prj_boundary_apply_axis_bf(block, stage, 0, 1, face_only);
        }
        if (prj_abs_double(block->xmin[1] - prj_boundary_active_mesh->coord.x2min) < tol) {
            if (bc->bc_x2_inner != PRJ_BC_OUTFLOW) {
                fprintf(stderr, "prj_boundary_physical_bf: only outflow BC is implemented for MHD\n");
                abort();
            }
            prj_boundary_apply_axis_bf(block, stage, 1, 0, face_only);
        }
        if (prj_abs_double(block->xmax[1] - prj_boundary_active_mesh->coord.x2max) < tol) {
            if (bc->bc_x2_outer != PRJ_BC_OUTFLOW) {
                fprintf(stderr, "prj_boundary_physical_bf: only outflow BC is implemented for MHD\n");
                abort();
            }
            prj_boundary_apply_axis_bf(block, stage, 1, 1, face_only);
        }
        if (prj_abs_double(block->xmin[2] - prj_boundary_active_mesh->coord.x3min) < tol) {
            if (bc->bc_x3_inner != PRJ_BC_OUTFLOW) {
                fprintf(stderr, "prj_boundary_physical_bf: only outflow BC is implemented for MHD\n");
                abort();
            }
            prj_boundary_apply_axis_bf(block, stage, 2, 0, face_only);
        }
        if (prj_abs_double(block->xmax[2] - prj_boundary_active_mesh->coord.x3max) < tol) {
            if (bc->bc_x3_outer != PRJ_BC_OUTFLOW) {
                fprintf(stderr, "prj_boundary_physical_bf: only outflow BC is implemented for MHD\n");
                abort();
            }
            prj_boundary_apply_axis_bf(block, stage, 2, 1, face_only);
        }
    }
}
#endif

static int prj_boundary_idx_outside(int idx)
{
    return idx < 0 || idx >= PRJ_BLOCK_SIZE;
}

static void prj_boundary_apply_axis(double *dst, int axis, int side, int bc_type, int face_only)
{
    int i;
    int j;
    int k;

    for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
        for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
            for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                int idx[3];
                int src_idx[3];
                int v;
                int normal_v;

                idx[0] = i;
                idx[1] = j;
                idx[2] = k;
                if (face_only != 0) {
                    int axis1 = (axis + 1) % 3;
                    int axis2 = (axis + 2) % 3;

                    if (prj_boundary_idx_outside(idx[axis1]) ||
                        prj_boundary_idx_outside(idx[axis2])) {
                        continue;
                    }
                }
                if (side == 0) {
                    if (idx[axis] >= 0) {
                        continue;
                    }
                    src_idx[axis] = -idx[axis] - 1;
                } else {
                    if (idx[axis] < PRJ_BLOCK_SIZE) {
                        continue;
                    }
                    src_idx[axis] = 2 * PRJ_BLOCK_SIZE - 1 - idx[axis];
                }
                src_idx[(axis + 1) % 3] = idx[(axis + 1) % 3];
                src_idx[(axis + 2) % 3] = idx[(axis + 2) % 3];
                for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
                    dst[VIDX(v, idx[0], idx[1], idx[2])] = dst[VIDX(v, src_idx[0], src_idx[1], src_idx[2])];
                }
                if (bc_type == PRJ_BC_REFLECT) {
                    normal_v = axis == 0 ? PRJ_PRIM_V1 : (axis == 1 ? PRJ_PRIM_V2 : PRJ_PRIM_V3);
                    dst[VIDX(normal_v, idx[0], idx[1], idx[2])] = -dst[VIDX(normal_v, src_idx[0], src_idx[1], src_idx[2])];
                } else if (bc_type == PRJ_BC_USER) {
                    if (axis == 0 && side == 0) {
                        dst[VIDX(PRJ_PRIM_V1, idx[0], idx[1], idx[2])] = 5.0;
                    }
                }
            }
        }
    }
}

void prj_boundary_physical(const prj_bc *bc, prj_block *block, int stage, int mode)
{
    double *dst = prj_boundary_stage_array(block, stage);
    const double tol = 1.0e-12;
    int pass;
    int npass = mode == PRJ_BOUNDARY_PHYS_FACE_ONLY ? 1 : 3;
    int face_only = mode == PRJ_BOUNDARY_PHYS_FACE_ONLY ? 1 : 0;

    for (pass = 0; pass < npass; ++pass) {
        if (prj_abs_double(block->xmin[0] - prj_boundary_active_mesh->coord.x1min) < tol) {
            prj_boundary_apply_axis(dst, 0, 0, bc->bc_x1_inner, face_only);
        }
        if (prj_abs_double(block->xmax[0] - prj_boundary_active_mesh->coord.x1max) < tol) {
            prj_boundary_apply_axis(dst, 0, 1, bc->bc_x1_outer, face_only);
        }
        if (prj_abs_double(block->xmin[1] - prj_boundary_active_mesh->coord.x2min) < tol) {
            prj_boundary_apply_axis(dst, 1, 0, bc->bc_x2_inner, face_only);
        }
        if (prj_abs_double(block->xmax[1] - prj_boundary_active_mesh->coord.x2max) < tol) {
            prj_boundary_apply_axis(dst, 1, 1, bc->bc_x2_outer, face_only);
        }
        if (prj_abs_double(block->xmin[2] - prj_boundary_active_mesh->coord.x3min) < tol) {
            prj_boundary_apply_axis(dst, 2, 0, bc->bc_x3_inner, face_only);
        }
        if (prj_abs_double(block->xmax[2] - prj_boundary_active_mesh->coord.x3max) < tol) {
            prj_boundary_apply_axis(dst, 2, 1, bc->bc_x3_outer, face_only);
        }
    }
}

void prj_boundary_send_bf(prj_block *block, int stage, int fill_kind)
{
#if PRJ_MHD
    if (block == 0) {
        fprintf(stderr, "prj_boundary_send_bf: block is null\n");
        abort();
    }
    if (fill_kind == PRJ_BOUNDARY_FILL_NONRECON) {
        double *restrict_sum[3] = {0, 0, 0};
        int *restrict_count[3] = {0, 0, 0};
        int have_restrict = 0;
        int n;

        for (n = 0; n < 56; ++n) {
            int nid = block->slot[n].id;
            const prj_block *neighbor;

            if (nid < 0 || prj_boundary_active_mesh == 0 || nid >= prj_boundary_active_mesh->nblocks) {
                continue;
            }
            neighbor = &prj_boundary_active_mesh->blocks[nid];
            if (!prj_boundary_active_block(neighbor)) {
                continue;
            }
            if (neighbor->level == block->level) {
                int dir;

                for (dir = 0; dir < 3; ++dir) {
                    double *dst = prj_mhd_bface_stage(block, stage, dir);
                    const double *src = prj_mhd_bface_stage_const(neighbor, stage, dir);
                    int dst_start[3];
                    int dst_end[3];
                    int i;
                    int j;
                    int k;

                    if (dst == 0 || src == 0) {
                        fprintf(stderr, "prj_boundary_send_bf: face storage is not allocated\n");
                        abort();
                    }
                    prj_boundary_same_level_face_bounds(block, neighbor, &block->slot[n], dir, 1, dst_start, dst_end);
                    for (i = dst_start[0]; i <= dst_end[0]; ++i) {
                        for (j = dst_start[1]; j <= dst_end[1]; ++j) {
                            for (k = dst_start[2]; k <= dst_end[2]; ++k) {
                                double x[3];
                                int ii;
                                int jj;
                                int kk;

                                prj_mhd_face_position(block, dir, i, j, k, x);
                                if (!prj_mhd_block_owns_face_position(neighbor, dir, x, &ii, &jj, &kk)) {
                                    continue;
                                }
                                dst[IDX(i, j, k)] = src[IDX(ii, jj, kk)];
#if PRJ_MHD_DEBUG
                                if (dir == X2DIR &&
                                    ((block->id == 1 && i == 16 && j == 17) ||
                                        (block->id == 5 && i == 16 && j == 1))) {
                                    fprintf(stderr,
                                        "[mhd-bf-fill-debug] kind=NONRECON_SAME dst_block=%d src_block=%d "
                                        "face=(%d,%d,%d) src_face=(%d,%d,%d) value=% .17g\n",
                                        block->id, neighbor->id, i, j, k, ii, jj, kk, dst[IDX(i, j, k)]);
                                }
                                if (dir == X3DIR &&
                                    ((block->id == 1 && i == 16 && j == 16) ||
                                        (block->id == 5 && i == 16 && j == 0))) {
                                    fprintf(stderr,
                                        "[mhd-bf-fill-debug] kind=NONRECON_SAME dst_block=%d src_block=%d "
                                        "face=(%d,%d,%d) src_face=(%d,%d,%d) value=% .17g\n",
                                        block->id, neighbor->id, i, j, k, ii, jj, kk, dst[IDX(i, j, k)]);
                                }
#endif
                            }
                        }
                    }
                }
            } else if (neighbor->level == block->level + 1) {
                int dir;

                if (!have_restrict) {
                    int axis;

                    for (axis = 0; axis < 3; ++axis) {
                        restrict_sum[axis] = (double *)calloc((size_t)PRJ_BLOCK_NCELLS, sizeof(*restrict_sum[axis]));
                        restrict_count[axis] = (int *)calloc((size_t)PRJ_BLOCK_NCELLS, sizeof(*restrict_count[axis]));
                        if (restrict_sum[axis] == 0 || restrict_count[axis] == 0) {
                            fprintf(stderr, "prj_boundary_send_bf: failed to allocate restriction buffers\n");
                            abort();
                        }
                    }
                    have_restrict = 1;
                }
                for (dir = 0; dir < 3; ++dir) {
                    const double *src = prj_mhd_bface_stage_const(neighbor, stage, dir);
                    int dst_start[3];
                    int dst_end[3];
                    int i;
                    int j;
                    int k;

                    if (src == 0) {
                        fprintf(stderr, "prj_boundary_send_bf: fine face storage is not allocated\n");
                        abort();
                    }
                    prj_boundary_coarse_face_bounds(block, neighbor, &block->slot[n], dir, 1, dst_start, dst_end);
                    for (i = dst_start[0]; i <= dst_end[0]; ++i) {
                        for (j = dst_start[1]; j <= dst_end[1]; ++j) {
                            for (k = dst_start[2]; k <= dst_end[2]; ++k) {
                                int aidx;
                                int bidx;

                                if (prj_mhd_face_is_interior(dir, i, j, k)) {
                                    continue;
                                }
                                for (aidx = 0; aidx < 2; ++aidx) {
                                    for (bidx = 0; bidx < 2; ++bidx) {
                                        double x[3];
                                        int ii;
                                        int jj;
                                        int kk;

                                        prj_boundary_coarse_subface_position(block, dir, i, j, k, aidx, bidx, x);
                                        if (!prj_mhd_block_owns_face_position(neighbor, dir, x, &ii, &jj, &kk)) {
                                            continue;
                                        }
                                        restrict_sum[dir][IDX(i, j, k)] += src[IDX(ii, jj, kk)];
                                        restrict_count[dir][IDX(i, j, k)] += 1;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        if (have_restrict) {
            int dir;

            for (dir = 0; dir < 3; ++dir) {
                double *dst = prj_mhd_bface_stage(block, stage, dir);
                int i;
                int j;
                int k;

                for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
                    for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
                        for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                            if (prj_mhd_face_is_interior(dir, i, j, k)) {
                                continue;
                            }
                            if (!prj_boundary_valid_restrict_count(restrict_count[dir][IDX(i, j, k)])) {
                                fprintf(stderr,
                                    "prj_boundary_send_bf: invalid restriction coverage block=%d dir=%d "
                                    "face=(%d,%d,%d) count=%d\n",
                                    block->id, dir, i, j, k, restrict_count[dir][IDX(i, j, k)]);
                                abort();
                            }
                            if (restrict_count[dir][IDX(i, j, k)] >= 4) {
                                dst[IDX(i, j, k)] =
                                    restrict_sum[dir][IDX(i, j, k)] / (double)restrict_count[dir][IDX(i, j, k)];
#if PRJ_MHD_DEBUG
                                if (dir == X2DIR &&
                                    ((block->id == 1 && i == 16 && j == 17) ||
                                        (block->id == 5 && i == 16 && j == 1))) {
                                    fprintf(stderr,
                                        "[mhd-bf-fill-debug] kind=NONRECON_RESTRICT dst_block=%d "
                                        "face=(%d,%d,%d) count=%d value=% .17g\n",
                                        block->id, i, j, k, restrict_count[dir][IDX(i, j, k)],
                                        dst[IDX(i, j, k)]);
                                }
                                if (dir == X3DIR &&
                                    ((block->id == 1 && i == 16 && j == 16) ||
                                        (block->id == 5 && i == 16 && j == 0))) {
                                    fprintf(stderr,
                                        "[mhd-bf-fill-debug] kind=NONRECON_RESTRICT dst_block=%d "
                                        "face=(%d,%d,%d) count=%d value=% .17g\n",
                                        block->id, i, j, k, restrict_count[dir][IDX(i, j, k)],
                                        dst[IDX(i, j, k)]);
                                }
#endif
                            }
                        }
                    }
                }
                free(restrict_count[dir]);
                free(restrict_sum[dir]);
            }
        }
        return;
    }
    if (fill_kind == PRJ_BOUNDARY_FILL_RECON) {
        int n;

        for (n = 0; n < 56; ++n) {
            int nid = block->slot[n].id;
            const prj_block *neighbor;
            const double *srcface[3];
            int i;
            int j;
            int k;

            if (nid < 0 || prj_boundary_active_mesh == 0 || nid >= prj_boundary_active_mesh->nblocks) {
                continue;
            }
            neighbor = &prj_boundary_active_mesh->blocks[nid];
            if (!prj_boundary_active_block(neighbor) || neighbor->level != block->level - 1) {
                continue;
            }
            srcface[0] = prj_mhd_bface_stage_const(neighbor, stage, X1DIR);
            srcface[1] = prj_mhd_bface_stage_const(neighbor, stage, X2DIR);
            srcface[2] = prj_mhd_bface_stage_const(neighbor, stage, X3DIR);
            if (srcface[0] == 0 || srcface[1] == 0 || srcface[2] == 0) {
                fprintf(stderr, "prj_boundary_send_bf: coarse face storage is not allocated\n");
                abort();
            }
            {
                int start[3];
                int end[3];

                prj_boundary_patch_bounds(block, neighbor, &block->slot[n], 0, start, end);
                for (i = start[0]; i <= end[0]; ++i) {
                    for (j = start[1]; j <= end[1]; ++j) {
                        for (k = start[2]; k <= end[2]; ++k) {
                            prj_mhd_apply_bf_patch(prj_boundary_active_mesh, block, stage, neighbor, srcface, i, j, k);
                        }
                    }
                }
            }
        }
        return;
    }
    fprintf(stderr, "prj_boundary_send_bf: invalid fill kind %d\n", fill_kind);
    abort();
#else
    (void)block;
    (void)stage;
    (void)fill_kind;
#endif
}

void prj_boundary_mpi_recv(prj_mesh *mesh, int stage, int fill_kind)
{
    prj_mpi *mpi = prj_mpi_current();

    if (mpi == 0) {
        return;
    }
    prj_mpi_exchange_ghosts(mesh, mpi, stage, fill_kind);
}

void prj_boundary_fill_ghosts(prj_mesh *mesh, const prj_bc *bc, int stage)
{
    int i;

    prj_boundary_active_mesh = mesh;
    for (i = 0; i < mesh->nblocks; ++i) {
        if (prj_boundary_active_block(&mesh->blocks[i])) {
            prj_boundary_physical(bc, &mesh->blocks[i], stage, PRJ_BOUNDARY_PHYS_FACE_ONLY);
        }
    }
    for (i = 0; i < mesh->nblocks; ++i) {
        if (prj_boundary_active_block(&mesh->blocks[i])) {
            prj_boundary_send(&mesh->blocks[i], stage, PRJ_BOUNDARY_FILL_NONRECON);
        }
    }
    prj_boundary_mpi_recv(mesh, stage, PRJ_BOUNDARY_FILL_NONRECON);
    for (i = 0; i < mesh->nblocks; ++i) {
        if (prj_boundary_active_block(&mesh->blocks[i])) {
            prj_boundary_send(&mesh->blocks[i], stage, PRJ_BOUNDARY_FILL_RECON);
        }
    }
    prj_boundary_mpi_recv(mesh, stage, PRJ_BOUNDARY_FILL_RECON);
    for (i = 0; i < mesh->nblocks; ++i) {
        if (prj_boundary_active_block(&mesh->blocks[i])) {
            prj_boundary_physical(bc, &mesh->blocks[i], stage, PRJ_BOUNDARY_PHYS_ALL);
        }
    }
}

void prj_boundary_fill_ghosts_bf(prj_mesh *mesh, const prj_bc *bc, int stage)
{
#if PRJ_MHD
    int i;
    prj_mpi *mpi = prj_mpi_current();

    if (mesh == 0 || bc == 0) {
        fprintf(stderr, "prj_boundary_fill_ghosts_bf: invalid input\n");
        abort();
    }

    prj_boundary_active_mesh = mesh;
    for (i = 0; i < mesh->nblocks; ++i) {
        if (prj_boundary_active_block(&mesh->blocks[i])) {
            prj_boundary_physical_bf(bc, &mesh->blocks[i], stage, PRJ_BOUNDARY_PHYS_FACE_ONLY);
        }
    }
    for (i = 0; i < mesh->nblocks; ++i) {
        if (prj_boundary_active_block(&mesh->blocks[i])) {
            prj_boundary_send_bf(&mesh->blocks[i], stage, PRJ_BOUNDARY_FILL_NONRECON);
        }
    }
    if (mpi != 0) {
        prj_mpi_exchange_bf(mesh, mpi, stage, PRJ_BOUNDARY_FILL_NONRECON);
    }
    /* A first NONRECON pass applies coarse-fine restriction. Run it once more so
     * same-level ghost duplicates see any owned-face updates before RECON. */
    for (i = 0; i < mesh->nblocks; ++i) {
        if (prj_boundary_active_block(&mesh->blocks[i])) {
            prj_boundary_send_bf(&mesh->blocks[i], stage, PRJ_BOUNDARY_FILL_NONRECON);
        }
    }
    if (mpi != 0) {
        prj_mpi_exchange_bf(mesh, mpi, stage, PRJ_BOUNDARY_FILL_NONRECON);
    }
    for (i = 0; i < mesh->nblocks; ++i) {
        if (prj_boundary_active_block(&mesh->blocks[i])) {
            prj_boundary_send_bf(&mesh->blocks[i], stage, PRJ_BOUNDARY_FILL_RECON);
        }
    }
    if (mpi != 0) {
        prj_mpi_exchange_bf(mesh, mpi, stage, PRJ_BOUNDARY_FILL_RECON);
    }
    /* RECON can change owned faces seen by same-level duplicates, so do one
     * final NONRECON sync before converting back to cell-centred fields. */
    for (i = 0; i < mesh->nblocks; ++i) {
        if (prj_boundary_active_block(&mesh->blocks[i])) {
            prj_boundary_send_bf(&mesh->blocks[i], stage, PRJ_BOUNDARY_FILL_NONRECON);
        }
    }
    if (mpi != 0) {
        prj_mpi_exchange_bf(mesh, mpi, stage, PRJ_BOUNDARY_FILL_NONRECON);
    }
    for (i = 0; i < mesh->nblocks; ++i) {
        if (prj_boundary_active_block(&mesh->blocks[i])) {
            prj_boundary_physical_bf(bc, &mesh->blocks[i], stage, PRJ_BOUNDARY_PHYS_ALL);
            prj_mhd_bf2bc(&mesh->blocks[i], stage);
        }
    }
#else
    (void)mesh;
    (void)bc;
    (void)stage;
#endif
}
