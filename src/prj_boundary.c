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
    if ((cases[0] == 1 || cases[0] == 3) &&
        (cases[1] == 1 || cases[1] == 3) &&
        (cases[2] == 1 || cases[2] == 3)) {
        return PRJ_BOUNDARY_FILL_RECON;
    }
    return -1;
}

static double prj_boundary_read_prim(const double *src, int var, int i, int j, int k)
{
    if (i < -PRJ_NGHOST || i >= PRJ_BLOCK_SIZE + PRJ_NGHOST ||
        j < -PRJ_NGHOST || j >= PRJ_BLOCK_SIZE + PRJ_NGHOST ||
        k < -PRJ_NGHOST || k >= PRJ_BLOCK_SIZE + PRJ_NGHOST) {
        fprintf(stderr,
            "prj_boundary_read_prim: out-of-range access var=%d i=%d j=%d k=%d "
            "(valid [%d, %d])\n",
            var, i, j, k, -PRJ_NGHOST, PRJ_BLOCK_SIZE + PRJ_NGHOST - 1);
        exit(EXIT_FAILURE);
    }
    return src[VIDX(var, i, j, k)];
}

void prj_boundary_get_prim(const prj_block *block, int stage, double x1, double x2, double x3, double *w)
{
    const double *src = prj_boundary_stage_array_const(block, stage);
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
        fprintf(stderr, "prj_boundary_get_prim: unsupported sample location (%g, %g, %g)\n", x1, x2, x3);
        exit(EXIT_FAILURE);
    }

    for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
        if (cases[0] == 1 && cases[1] == 1 && cases[2] == 1) {
            int i = (int)(ox[0] >= 0.0 ? ox[0] + 0.5 : ox[0] - 0.5);
            int j = (int)(ox[1] >= 0.0 ? ox[1] + 0.5 : ox[1] - 0.5);
            int k = (int)(ox[2] >= 0.0 ? ox[2] + 0.5 : ox[2] - 0.5);

            w[v] = prj_boundary_read_prim(src, v, i, j, k);
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
                        sum += prj_boundary_read_prim(src, v, i + di, j + dj, k + dk);
                    }
                }
            }
            w[v] = 0.125 * sum;
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
            double base = prj_boundary_read_prim(src, v, i, j, k);

            stx[0] = prj_boundary_read_prim(src, v, i - 1, j, k);
            stx[1] = base;
            stx[2] = prj_boundary_read_prim(src, v, i + 1, j, k);
            sty[0] = prj_boundary_read_prim(src, v, i, j - 1, k);
            sty[1] = base;
            sty[2] = prj_boundary_read_prim(src, v, i, j + 1, k);
            stz[0] = prj_boundary_read_prim(src, v, i, j, k - 1);
            stz[1] = base;
            stz[2] = prj_boundary_read_prim(src, v, i, j, k + 1);
            sx = prj_reconstruct_slope(stx, block->dx[0]);
            sy = prj_reconstruct_slope(sty, block->dx[1]);
            sz = prj_reconstruct_slope(stz, block->dx[2]);
            w[v] = base +
                sx * (x1 - xcenter) +
                sy * (x2 - ycenter) +
                sz * (x3 - zcenter);
        } else {
            fprintf(stderr, "prj_boundary_get_prim: mixed unsupported sample location (%g, %g, %g)\n", x1, x2, x3);
            exit(EXIT_FAILURE);
        }
    }
}

void prj_boundary_send(prj_block *block, int stage, int fill_kind)
{
    int n;

    for (n = 0; n < 56; ++n) {
        int id = block->slot[n].id;

        if (id >= 0 && prj_boundary_active_mesh != 0 && id < prj_boundary_active_mesh->nblocks) {
            prj_block *neighbor = &prj_boundary_active_mesh->blocks[id];
            double *dst;
            int i;
            int j;
            int k;

            if (!prj_boundary_active_block(neighbor)) {
                continue;
            }
            dst = prj_boundary_stage_array(neighbor, stage);
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
                            if (sample_kind < 0 ||
                                (fill_kind != PRJ_BOUNDARY_FILL_ALL && sample_kind != fill_kind)) {
                                continue;
                            }
                            prj_boundary_get_prim(block, stage, x1, x2, x3, w);
                            if (neighbor->rank == block->rank) {
                                for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
                                    dst[VIDX(v, i, j, k)] = w[v];
                                }
                            } else {
                                /* MPI ghost exchange buffer path not implemented yet. */
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
