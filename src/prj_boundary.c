#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "prj.h"

static prj_mesh *prj_boundary_active_mesh = 0;

static double prj_abs_double(double x)
{
    return x < 0.0 ? -x : x;
}

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

#if PRJ_MHD
static double *prj_boundary_bf_array(prj_block *block, int dir, int use_bf1)
{
    return use_bf1 != 0 ? block->Bf1[dir] : block->Bf[dir];
}

static const double *prj_boundary_bf_array_const(const prj_block *block, int dir, int use_bf1)
{
    return use_bf1 != 0 ? block->Bf1[dir] : block->Bf[dir];
}

static void prj_boundary_check_bf_storage(const prj_block *block, const char *label)
{
    int d;

    if (block == 0) {
        fprintf(stderr, "%s: missing face fidelity storage\n", label);
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

static int prj_boundary_storage_index_ok(int i, int j, int k)
{
    return i >= -PRJ_NGHOST && i < PRJ_BLOCK_SIZE + PRJ_NGHOST &&
        j >= -PRJ_NGHOST && j < PRJ_BLOCK_SIZE + PRJ_NGHOST &&
        k >= -PRJ_NGHOST && k < PRJ_BLOCK_SIZE + PRJ_NGHOST;
}

static int prj_boundary_nearest_int(double x)
{
    return x >= 0.0 ? (int)(x + 0.5) : (int)(x - 0.5);
}

static int prj_boundary_bf_axis_active_max(int dir, int axis)
{
    return dir == axis ? PRJ_BLOCK_SIZE : PRJ_BLOCK_SIZE - 1;
}

static int prj_boundary_bf_face_active(int dir, int i, int j, int k)
{
    int idx[3];
    int d;

    idx[0] = i;
    idx[1] = j;
    idx[2] = k;
    for (d = 0; d < 3; ++d) {
        if (idx[d] < 0 || idx[d] > prj_boundary_bf_axis_active_max(dir, d)) {
            return 0;
        }
    }
    return 1;
}

static double prj_boundary_bf_face_coord(const prj_block *block, int dir, int axis, int idx)
{
    double offset = dir == axis ? 0.0 : 0.5;

    return block->xmin[axis] + ((double)idx + offset) * block->dx[axis];
}

static int prj_boundary_bf_point_inside(const prj_block *block, const double x[3])
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

static int prj_boundary_bf_index_from_coord(const prj_block *block, int dir,
    const double x[3], int idx[3])
{
    const double tol = 1.0e-8;
    int d;

    if (block == 0 || idx == 0) {
        return 0;
    }
    for (d = 0; d < 3; ++d) {
        double offset = dir == d ? 0.0 : 0.5;
        double q;
        int n;

        if (block->dx[d] <= 0.0) {
            return 0;
        }
        q = (x[d] - block->xmin[d]) / block->dx[d] - offset;
        n = prj_boundary_nearest_int(q);
        if (fabs(q - (double)n) > tol) {
            return 0;
        }
        idx[d] = n;
    }
    return prj_boundary_storage_index_ok(idx[0], idx[1], idx[2]);
}

static int prj_boundary_aligned_cell_index(const prj_block *block, int axis, double xlo)
{
    const double tol = 1.0e-8;
    double q;
    int idx;

    if (block == 0 || block->dx[axis] <= 0.0) {
        return -1000000;
    }
    q = (xlo - block->xmin[axis]) / block->dx[axis];
    idx = prj_boundary_nearest_int(q);
    if (fabs(q - (double)idx) > tol) {
        return -1000000;
    }
    return idx;
}

static void prj_boundary_write_bf_face(prj_block *block, int use_bf1, int dir,
    int i, int j, int k, double value, int fidelity)
{
    int idx;
    double *dst;

    prj_boundary_check_bf_storage(block, "prj_boundary_write_bf_face");
    if (dir < 0 || dir >= 3 || !prj_boundary_storage_index_ok(i, j, k)) {
        fprintf(stderr, "prj_boundary_write_bf_face: invalid face index dir=%d i=%d j=%d k=%d\n",
            dir, i, j, k);
        exit(EXIT_FAILURE);
    }
    if (!isfinite(value) || fidelity < PRJ_MHD_FIDELITY_NONE ||
        fidelity > PRJ_MHD_FIDELITY_FINER) {
        fprintf(stderr, "prj_boundary_write_bf_face: invalid value or fidelity\n");
        exit(EXIT_FAILURE);
    }
    idx = IDX(i, j, k);
    if (fidelity < block->face_fidelity[dir][idx]) {
        return;
    }
    dst = prj_boundary_bf_array(block, dir, use_bf1);
    dst[idx] = value;
    if (fidelity > block->face_fidelity[dir][idx]) {
        block->face_fidelity[dir][idx] = fidelity;
    }
}

static double prj_boundary_restrict_bf_value(const prj_block *fine, int use_bf1,
    int dir, const double x[3])
{
    const double tol = 1.0e-8;
    const double *src;
    int normal;
    int base[3] = {0, 0, 0};
    int tan0 = (dir + 1) % 3;
    int tan1 = (dir + 2) % 3;
    double sum = 0.0;
    int a;
    int b;

    prj_boundary_check_bf_storage(fine, "prj_boundary_restrict_bf_value");
    src = prj_boundary_bf_array_const(fine, dir, use_bf1);
    for (a = 0; a < 3; ++a) {
        double q;

        if (fine->dx[a] <= 0.0) {
            fprintf(stderr, "prj_boundary_restrict_bf_value: invalid fine cell size\n");
            exit(EXIT_FAILURE);
        }
        q = (x[a] - fine->xmin[a]) / fine->dx[a];
        if (a == dir) {
            normal = prj_boundary_nearest_int(q);
            if (fabs(q - (double)normal) > tol) {
                fprintf(stderr, "prj_boundary_restrict_bf_value: normal face is not fine-aligned\n");
                exit(EXIT_FAILURE);
            }
            base[a] = normal;
        } else {
            int lower = prj_floor_to_int(q - 0.5);

            if (fabs((q - 0.5) - ((double)lower + 0.5)) > tol) {
                fprintf(stderr, "prj_boundary_restrict_bf_value: tangential face is not centered on two fine faces\n");
                exit(EXIT_FAILURE);
            }
            base[a] = lower;
        }
    }

    for (a = 0; a < 2; ++a) {
        for (b = 0; b < 2; ++b) {
            int idx[3];
            double value;

            idx[dir] = base[dir];
            idx[tan0] = base[tan0] + a;
            idx[tan1] = base[tan1] + b;
            if (!prj_boundary_storage_index_ok(idx[0], idx[1], idx[2])) {
                fprintf(stderr,
                    "prj_boundary_restrict_bf_value: fine index out of storage dir=%d i=%d j=%d k=%d\n",
                    dir, idx[0], idx[1], idx[2]);
                exit(EXIT_FAILURE);
            }
            value = src[IDX(idx[0], idx[1], idx[2])];
            if (!isfinite(value)) {
                fprintf(stderr, "prj_boundary_restrict_bf_value: non-finite fine face\n");
                exit(EXIT_FAILURE);
            }
            sum += value;
        }
    }
    return 0.25 * sum;
}

static void prj_boundary_copy_bf_same_level(const prj_block *src_block,
    prj_block *dst_block, int use_bf1)
{
    int dir;

    prj_boundary_check_bf_storage(src_block, "prj_boundary_copy_bf_same_level");
    prj_boundary_check_bf_storage(dst_block, "prj_boundary_copy_bf_same_level");
    for (dir = 0; dir < 3; ++dir) {
        const double *src = prj_boundary_bf_array_const(src_block, dir, use_bf1);
        int i;
        int j;
        int k;

        for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
            for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
                for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                    double x[3];
                    int sidx[3];
                    double value;

                    if (prj_boundary_bf_face_active(dir, i, j, k)) {
                        continue;
                    }
                    x[0] = prj_boundary_bf_face_coord(dst_block, dir, 0, i);
                    x[1] = prj_boundary_bf_face_coord(dst_block, dir, 1, j);
                    x[2] = prj_boundary_bf_face_coord(dst_block, dir, 2, k);
                    if (!prj_boundary_bf_point_inside(src_block, x)) {
                        continue;
                    }
                    if (!prj_boundary_bf_index_from_coord(src_block, dir, x, sidx)) {
                        fprintf(stderr, "prj_boundary_copy_bf_same_level: source face is not grid-aligned\n");
                        exit(EXIT_FAILURE);
                    }
                    value = src[IDX(sidx[0], sidx[1], sidx[2])];
                    prj_boundary_write_bf_face(dst_block, use_bf1, dir, i, j, k, value,
                        PRJ_MHD_FIDELITY_SAME);
                }
            }
        }
    }
}

static void prj_boundary_restrict_bf_to_coarse(const prj_block *fine,
    prj_block *coarse, int use_bf1)
{
    int dir;

    prj_boundary_check_bf_storage(fine, "prj_boundary_restrict_bf_to_coarse");
    prj_boundary_check_bf_storage(coarse, "prj_boundary_restrict_bf_to_coarse");
    for (dir = 0; dir < 3; ++dir) {
        int i;
        int j;
        int k;

        for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
            for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
                for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                    double x[3];
                    double value;

                    if (prj_boundary_bf_face_active(dir, i, j, k)) {
                        continue;
                    }
                    x[0] = prj_boundary_bf_face_coord(coarse, dir, 0, i);
                    x[1] = prj_boundary_bf_face_coord(coarse, dir, 1, j);
                    x[2] = prj_boundary_bf_face_coord(coarse, dir, 2, k);
                    if (!prj_boundary_bf_point_inside(fine, x)) {
                        continue;
                    }
                    value = prj_boundary_restrict_bf_value(fine, use_bf1, dir, x);
                    prj_boundary_write_bf_face(coarse, use_bf1, dir, i, j, k, value,
                        PRJ_MHD_FIDELITY_FINER);
                }
            }
        }
    }
}

static int prj_boundary_patch_fully_active(int fi, int fj, int fk)
{
    return fi >= 0 && fi + 1 < PRJ_BLOCK_SIZE &&
        fj >= 0 && fj + 1 < PRJ_BLOCK_SIZE &&
        fk >= 0 && fk + 1 < PRJ_BLOCK_SIZE;
}

static void prj_boundary_prolong_bf_to_fine(const prj_block *coarse,
    prj_block *fine, int use_bf1)
{
    int fi;
    int fj;
    int fk;

    prj_boundary_check_bf_storage(coarse, "prj_boundary_prolong_bf_to_fine");
    prj_boundary_check_bf_storage(fine, "prj_boundary_prolong_bf_to_fine");
    for (fi = -PRJ_NGHOST; fi + 2 < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++fi) {
        for (fj = -PRJ_NGHOST; fj + 2 < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++fj) {
            for (fk = -PRJ_NGHOST; fk + 2 < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++fk) {
                double lo[3];
                double hi[3];
                int ci;
                int cj;
                int ck;

                if (prj_boundary_patch_fully_active(fi, fj, fk)) {
                    continue;
                }
                lo[0] = fine->xmin[0] + (double)fi * fine->dx[0];
                lo[1] = fine->xmin[1] + (double)fj * fine->dx[1];
                lo[2] = fine->xmin[2] + (double)fk * fine->dx[2];
                hi[0] = lo[0] + 2.0 * fine->dx[0];
                hi[1] = lo[1] + 2.0 * fine->dx[1];
                hi[2] = lo[2] + 2.0 * fine->dx[2];
                if (lo[0] < coarse->xmin[0] - 1.0e-12 || hi[0] > coarse->xmax[0] + 1.0e-12 ||
                    lo[1] < coarse->xmin[1] - 1.0e-12 || hi[1] > coarse->xmax[1] + 1.0e-12 ||
                    lo[2] < coarse->xmin[2] - 1.0e-12 || hi[2] > coarse->xmax[2] + 1.0e-12) {
                    continue;
                }
                ci = prj_boundary_aligned_cell_index(coarse, 0, lo[0]);
                cj = prj_boundary_aligned_cell_index(coarse, 1, lo[1]);
                ck = prj_boundary_aligned_cell_index(coarse, 2, lo[2]);
                if (ci < -PRJ_NGHOST || ci >= PRJ_BLOCK_SIZE + PRJ_NGHOST ||
                    cj < -PRJ_NGHOST || cj >= PRJ_BLOCK_SIZE + PRJ_NGHOST ||
                    ck < -PRJ_NGHOST || ck >= PRJ_BLOCK_SIZE + PRJ_NGHOST) {
                    continue;
                }
                prj_mhd_bf_prolongate(coarse, fine, ci, cj, ck, fi, fj, fk, use_bf1);
            }
        }
    }
}

static void prj_boundary_init_face_fidelity(prj_mesh *mesh)
{
    int bidx;

    if (mesh == 0) {
        return;
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        int dir;

        if (!prj_boundary_active_block(block)) {
            continue;
        }
        prj_boundary_check_bf_storage(block, "prj_boundary_init_face_fidelity");
        for (dir = 0; dir < 3; ++dir) {
            int i;
            int j;
            int k;
            int n;

            for (n = 0; n < PRJ_BLOCK_NCELLS; ++n) {
                block->face_fidelity[dir][n] = PRJ_MHD_FIDELITY_NONE;
            }

            for (i = 0; i <= prj_boundary_bf_axis_active_max(dir, 0); ++i) {
                for (j = 0; j <= prj_boundary_bf_axis_active_max(dir, 1); ++j) {
                    for (k = 0; k <= prj_boundary_bf_axis_active_max(dir, 2); ++k) {
                        block->face_fidelity[dir][IDX(i, j, k)] = PRJ_MHD_FIDELITY_SAME;
                    }
                }
            }
        }
    }
}

static void prj_boundary_check_outflow_bf(int bc_type, const char *label)
{
    if (bc_type != PRJ_BC_OUTFLOW) {
        fprintf(stderr, "%s: only outflow magnetic-field boundary conditions are implemented\n", label);
        exit(EXIT_FAILURE);
    }
}

static void prj_boundary_apply_bf_axis(prj_block *block, int use_bf1,
    int axis, int side, int bc_type)
{
    int dir;

    prj_boundary_check_outflow_bf(bc_type, "prj_boundary_apply_bf_axis");
    prj_boundary_check_bf_storage(block, "prj_boundary_apply_bf_axis");
    for (dir = 0; dir < 3; ++dir) {
        int max_axis = prj_boundary_bf_axis_active_max(dir, axis);
        int i;
        int j;
        int k;

        for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
            for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
                for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                    int idx[3];
                    int src_idx[3];
                    int src_flat;
                    int fidelity;
                    double value;
                    const double *src;

                    idx[0] = i;
                    idx[1] = j;
                    idx[2] = k;
                    if (side == 0) {
                        if (idx[axis] >= 0) {
                            continue;
                        }
                        src_idx[axis] = 0;
                    } else {
                        if (idx[axis] <= max_axis) {
                            continue;
                        }
                        src_idx[axis] = max_axis;
                    }
                    src_idx[(axis + 1) % 3] = idx[(axis + 1) % 3];
                    src_idx[(axis + 2) % 3] = idx[(axis + 2) % 3];
                    if (!prj_boundary_storage_index_ok(src_idx[0], src_idx[1], src_idx[2])) {
                        continue;
                    }
                    src_flat = IDX(src_idx[0], src_idx[1], src_idx[2]);
                    fidelity = block->face_fidelity[dir][src_flat];
                    if (fidelity == PRJ_MHD_FIDELITY_NONE) {
                        continue;
                    }
                    src = prj_boundary_bf_array_const(block, dir, use_bf1);
                    value = src[src_flat];
                    prj_boundary_write_bf_face(block, use_bf1, dir, i, j, k, value, fidelity);
                }
            }
        }
    }
}

static void prj_boundary_physical_bf(const prj_bc *bc, prj_block *block, int use_bf1)
{
    const double tol = 1.0e-12;
    int pass;

    if (bc == 0 || prj_boundary_active_mesh == 0) {
        fprintf(stderr, "prj_boundary_physical_bf: missing boundary context\n");
        exit(EXIT_FAILURE);
    }
    for (pass = 0; pass < 3; ++pass) {
        if (prj_abs_double(block->xmin[0] - prj_boundary_active_mesh->coord.x1min) < tol) {
            prj_boundary_apply_bf_axis(block, use_bf1, 0, 0, bc->bc_x1_inner);
        }
        if (prj_abs_double(block->xmax[0] - prj_boundary_active_mesh->coord.x1max) < tol) {
            prj_boundary_apply_bf_axis(block, use_bf1, 0, 1, bc->bc_x1_outer);
        }
        if (prj_abs_double(block->xmin[1] - prj_boundary_active_mesh->coord.x2min) < tol) {
            prj_boundary_apply_bf_axis(block, use_bf1, 1, 0, bc->bc_x2_inner);
        }
        if (prj_abs_double(block->xmax[1] - prj_boundary_active_mesh->coord.x2max) < tol) {
            prj_boundary_apply_bf_axis(block, use_bf1, 1, 1, bc->bc_x2_outer);
        }
        if (prj_abs_double(block->xmin[2] - prj_boundary_active_mesh->coord.x3min) < tol) {
            prj_boundary_apply_bf_axis(block, use_bf1, 2, 0, bc->bc_x3_inner);
        }
        if (prj_abs_double(block->xmax[2] - prj_boundary_active_mesh->coord.x3max) < tol) {
            prj_boundary_apply_bf_axis(block, use_bf1, 2, 1, bc->bc_x3_outer);
        }
    }
}

void prj_boundary_send_bf(prj_block *block, int use_bf1, int fill_kind)
{
    int n;

    if (!prj_boundary_active_block(block)) {
        return;
    }
    prj_boundary_check_bf_storage(block, "prj_boundary_send_bf");
    for (n = 0; n < 56; ++n) {
        int id = block->slot[n].id;
        prj_block *neighbor;

        if (id < 0 || prj_boundary_active_mesh == 0 || id >= prj_boundary_active_mesh->nblocks) {
            continue;
        }
        neighbor = &prj_boundary_active_mesh->blocks[id];
        if (!prj_boundary_active_block(neighbor) || neighbor->rank != block->rank) {
            continue;
        }
        if (fill_kind == PRJ_BOUNDARY_FILL_NONRECON) {
            if (block->slot[n].rel_level == 0) {
                prj_boundary_copy_bf_same_level(block, neighbor, use_bf1);
            } else if (block->slot[n].rel_level < 0) {
                prj_boundary_restrict_bf_to_coarse(block, neighbor, use_bf1);
            }
        } else if (fill_kind == PRJ_BOUNDARY_FILL_RECON) {
            if (block->slot[n].rel_level > 0) {
                prj_boundary_prolong_bf_to_fine(block, neighbor, use_bf1);
            }
        } else if (fill_kind == PRJ_BOUNDARY_FILL_ALL) {
            if (block->slot[n].rel_level == 0) {
                prj_boundary_copy_bf_same_level(block, neighbor, use_bf1);
            } else if (block->slot[n].rel_level < 0) {
                prj_boundary_restrict_bf_to_coarse(block, neighbor, use_bf1);
            } else {
                prj_boundary_prolong_bf_to_fine(block, neighbor, use_bf1);
            }
        }
    }
}

void prj_boundary_fill_bf(prj_mesh *mesh, const prj_bc *bc, int use_bf1)
{
    prj_mpi *mpi = prj_mpi_current();
    int i;

    if (mesh == 0 || bc == 0) {
        fprintf(stderr, "prj_boundary_fill_bf: mesh or bc is null\n");
        exit(EXIT_FAILURE);
    }
    prj_boundary_active_mesh = mesh;
    prj_boundary_init_face_fidelity(mesh);
    for (i = 0; i < mesh->nblocks; ++i) {
        if (prj_boundary_active_block(&mesh->blocks[i])) {
            prj_boundary_physical_bf(bc, &mesh->blocks[i], use_bf1);
        }
    }
    for (i = 0; i < mesh->nblocks; ++i) {
        if (prj_boundary_active_block(&mesh->blocks[i])) {
            prj_boundary_send_bf(&mesh->blocks[i], use_bf1, PRJ_BOUNDARY_FILL_NONRECON);
        }
    }
    prj_mpi_exchange_bf(mesh, mpi, use_bf1, PRJ_BOUNDARY_FILL_NONRECON);
    for (i = 0; i < mesh->nblocks; ++i) {
        if (prj_boundary_active_block(&mesh->blocks[i])) {
            prj_boundary_send_bf(&mesh->blocks[i], use_bf1, PRJ_BOUNDARY_FILL_RECON);
        }
    }
    prj_mpi_exchange_bf(mesh, mpi, use_bf1, PRJ_BOUNDARY_FILL_RECON);
    for (i = 0; i < mesh->nblocks; ++i) {
        if (prj_boundary_active_block(&mesh->blocks[i])) {
            prj_boundary_physical_bf(bc, &mesh->blocks[i], use_bf1);
        }
    }
}
#endif

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
