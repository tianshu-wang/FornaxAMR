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
    PRJ_BOUNDARY_PHYS_FACE_ONLY = 0,
    PRJ_BOUNDARY_PHYS_ALL = 1
};

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

void prj_boundary_get_slope(double *src, int v, int i, int j, int k, int is_eosvar, double *slope)
{
    double stx[3];
    double sty[3];
    double stz[3];
    double base = is_eosvar==1 ? src[EIDX(v, i, j, k)] : src[VIDX(v, i, j, k)];

    stx[0] = is_eosvar==1 ? src[EIDX(v, i-1, j, k)] : src[VIDX(v, i-1, j, k)];
    stx[1] = base;
    stx[2] = is_eosvar==1 ? src[EIDX(v, i+1, j, k)] : src[VIDX(v, i+1, j, k)];
    sty[0] = is_eosvar==1 ? src[EIDX(v, i, j-1, k)] : src[VIDX(v, i, j-1, k)];
    sty[1] = base;
    sty[2] = is_eosvar==1 ? src[EIDX(v, i, j+1, k)] : src[VIDX(v, i, j+1, k)];
    stz[0] = is_eosvar==1 ? src[EIDX(v, i, j, k-1)] : src[VIDX(v, i, j, k-1)];
    stz[1] = base;
    stz[2] = is_eosvar==1 ? src[EIDX(v, i, j, k+1)] : src[VIDX(v, i, j, k+1)];
    slope[0] = prj_reconstruct_slope(stx, 1.0);
    slope[1] = prj_reconstruct_slope(sty, 1.0);
    slope[2] = prj_reconstruct_slope(stz, 1.0);
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
    
    double *W_send = stage == 2 ? block->W1 : block->W;
    double *eos_send = block->eosvar;

    for (n = 0; n < 56; ++n) {
        int id = block->slot[n].id;

        if (id >= 0 && prj_boundary_active_mesh != 0 && id < prj_boundary_active_mesh->nblocks) {
            prj_block *neighbor = &prj_boundary_active_mesh->blocks[id];
            double *W_recv = stage == 2 ? neighbor->W1 : neighbor->W;
            double *eos_recv = neighbor->eosvar;
            int i;
            int j;
            int k;

            if (!prj_boundary_active_block(neighbor)) {
                continue;
            }
            if (neighbor->rank == block->rank) {
                const prj_neighbor *slot = &block->slot[n];
                int sample_kind;
                if(slot->rel_level<=0){sample_kind=PRJ_BOUNDARY_FILL_NONRECON;}
                else{sample_kind=PRJ_BOUNDARY_FILL_RECON;}
                if (fill_kind != PRJ_BOUNDARY_FILL_ALL && sample_kind != fill_kind) {
                    continue;
                }

                int eos_done = (slot->rel_level==0) ? 1 : 0;

                for (i = 0; i < slot->recv_loc_end[0]-slot->recv_loc_start[0]; ++i) {
                    for (j = 0; j < slot->recv_loc_end[1]-slot->recv_loc_start[1]; ++j) {
                        for (k = 0; k < slot->recv_loc_end[2]-slot->recv_loc_start[2]; ++k) {
                            int v;

                            if (slot->rel_level==0){
                                // Same level
                                for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
                                    W_recv[VIDX(v, i+slot->recv_loc_start[0], j+slot->recv_loc_start[1], k+slot->recv_loc_start[2])] = 
                                        W_send[VIDX(v, i+slot->send_loc_start[0], j+slot->send_loc_start[1], k+slot->send_loc_start[2])];
                                }
                                for (v = 0; v < PRJ_NVAR_EOSVAR; ++v) {
                                    eos_recv[EIDX(v, i+slot->recv_loc_start[0], j+slot->recv_loc_start[1], k+slot->recv_loc_start[2])] = 
                                        eos_send[EIDX(v, i+slot->send_loc_start[0], j+slot->send_loc_start[1], k+slot->send_loc_start[2])];
                                }
                            } else if (slot->rel_level==-1) {
                                // Neighbor is coarser, restriction
                                for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
                                    W_recv[VIDX(v, i+slot->recv_loc_start[0], j+slot->recv_loc_start[1], k+slot->recv_loc_start[2])] = 
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
                                    eos_recv[EIDX(v, i+slot->recv_loc_start[0], j+slot->recv_loc_start[1], k+slot->recv_loc_start[2])] = 
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
                                double slope[3]={0};
                                for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
                                    prj_boundary_get_slope(W_send, v, 
                                                           i/2+slot->send_loc_start[0], 
                                                           j/2+slot->send_loc_start[1], 
                                                           k/2+slot->send_loc_start[2], 
                                                           0, slope);
                                    W_recv[VIDX(v, i+slot->recv_loc_start[0], j+slot->recv_loc_start[1], k+slot->recv_loc_start[2])] =
                                        W_send[VIDX(v, i/2+slot->send_loc_start[0], j/2+slot->send_loc_start[1], k/2+slot->send_loc_start[2])]
                                       +((i%2==0) ? +0.25*slope[0] : -0.25*slope[0])
                                       +((j%2==0) ? +0.25*slope[1] : -0.25*slope[1])
                                       +((k%2==0) ? +0.25*slope[2] : -0.25*slope[2]);
                                }
                                for (v = 0; v < PRJ_NVAR_EOSVAR; ++v) {
                                    prj_boundary_get_slope(eos_send, v, 
                                                           i/2+slot->send_loc_start[0], 
                                                           j/2+slot->send_loc_start[1], 
                                                           k/2+slot->send_loc_start[2], 
                                                           1, slope);
                                    eos_recv[EIDX(v, i+slot->recv_loc_start[0], j+slot->recv_loc_start[1], k+slot->recv_loc_start[2])] =
                                        eos_send[EIDX(v, i/2+slot->send_loc_start[0], j/2+slot->send_loc_start[1], k/2+slot->send_loc_start[2])]
                                       +((i%2==0) ? +0.25*slope[0] : -0.25*slope[0])
                                       +((j%2==0) ? +0.25*slope[1] : -0.25*slope[1])
                                       +((k%2==0) ? +0.25*slope[2] : -0.25*slope[2]);
                                }
                            } else {
                              fprintf(stderr,"slot->rel_level unrecognized: %d\n", slot->rel_level);
                              exit(1);
                            }
                            neighbor->eos_done[IDX(i,j,k)] = eos_done;
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

static int prj_boundary_face_storage_index_ok(int dir, int i, int j, int k)
{
    int idx[3] = {i, j, k};
    int d;

    for (d = 0; d < 3; ++d) {
        int lo = -PRJ_NGHOST;
        int hi = dir == d ? PRJ_BLOCK_SIZE + PRJ_NGHOST : PRJ_BLOCK_SIZE + PRJ_NGHOST - 1;
        if (idx[d] < lo || idx[d] > hi) {
            return 0;
        }
    }
    return 1;
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

void prj_boundary_write_bf_face(prj_block *block, int use_bf1, int dir,
    int i, int j, int k, double value, int fidelity)
{
    int idx;
    double *dst;

    prj_boundary_check_bf_storage(block, "prj_boundary_write_bf_face");
    if (dir < 0 || dir >= 3 || !prj_boundary_face_storage_index_ok(dir, i, j, k)) {
        fprintf(stderr, "prj_boundary_write_bf_face: invalid face index dir=%d i=%d j=%d k=%d\n",
            dir, i, j, k);
        exit(EXIT_FAILURE);
    }
    if (!isfinite(value) || fidelity < PRJ_MHD_FIDELITY_NONE ||
        fidelity > PRJ_MHD_FIDELITY_FINER) {
        fprintf(stderr, "prj_boundary_write_bf_face: invalid value or fidelity\n");
        exit(EXIT_FAILURE);
    }
    idx = FACE_IDX(dir, i, j, k);
    if (fidelity < block->face_fidelity[dir][idx]) {
        return;
    }
    dst = prj_boundary_bf_array(block, dir, use_bf1);
    dst[idx] = value;
    if (fidelity > block->face_fidelity[dir][idx]) {
        block->face_fidelity[dir][idx] = fidelity;
    }
}

static void prj_boundary_copy_bf_same_level(const prj_block *src_block,
    prj_block *dst_block, int use_bf1, const prj_neighbor *slot)
{
    int dir;
    int max_count = 0;
    double *buffer;

    prj_boundary_check_bf_storage(src_block, "prj_boundary_copy_bf_same_level");
    prj_boundary_check_bf_storage(dst_block, "prj_boundary_copy_bf_same_level");

    for (dir = 0; dir < 3; ++dir) {
        int di = slot->recv_loc_end[0] - slot->recv_loc_start[0] + (dir == 0 ? 1 : 0);
        int dj = slot->recv_loc_end[1] - slot->recv_loc_start[1] + (dir == 1 ? 1 : 0);
        int dk = slot->recv_loc_end[2] - slot->recv_loc_start[2] + (dir == 2 ? 1 : 0);
        int count = di * dj * dk;
        if (count > max_count) max_count = count;
    }

    buffer = (double *)malloc(max_count * sizeof(double));
    if (buffer == 0) {
        fprintf(stderr, "prj_boundary_copy_bf_same_level: malloc failed\n");
        exit(EXIT_FAILURE);
    }

    for (dir = 0; dir < 3; ++dir) {
        const double *src = prj_boundary_bf_array_const(src_block, dir, use_bf1);
        int buf_idx = 0;
        int i;
        int j;
        int k;

        for (i = 0; i < slot->recv_loc_end[0]-slot->recv_loc_start[0]+((dir==0)?1:0); ++i) {
            for (j = 0; j < slot->recv_loc_end[1]-slot->recv_loc_start[1]+((dir==1)?1:0); ++j) {
                for (k = 0; k < slot->recv_loc_end[2]-slot->recv_loc_start[2]+((dir==2)?1:0); ++k) {
                    if (prj_boundary_bf_face_active(dir, i+slot->recv_loc_start[0], j+slot->recv_loc_start[1], k+slot->recv_loc_start[2])) {
                        continue;
                    }
                    buffer[buf_idx++] = src[FACE_IDX(dir, i+slot->send_loc_start[0], j+slot->send_loc_start[1], k+slot->send_loc_start[2])];
                }
            }
        }

        buf_idx = 0;
        for (i = 0; i < slot->recv_loc_end[0]-slot->recv_loc_start[0]+((dir==0)?1:0); ++i) {
            for (j = 0; j < slot->recv_loc_end[1]-slot->recv_loc_start[1]+((dir==1)?1:0); ++j) {
                for (k = 0; k < slot->recv_loc_end[2]-slot->recv_loc_start[2]+((dir==2)?1:0); ++k) {
                    if (prj_boundary_bf_face_active(dir, i+slot->recv_loc_start[0], j+slot->recv_loc_start[1], k+slot->recv_loc_start[2])) {
                        continue;
                    }
                    prj_boundary_write_bf_face(dst_block, use_bf1, dir,
                                               i+slot->recv_loc_start[0],
                                               j+slot->recv_loc_start[1],
                                               k+slot->recv_loc_start[2],
                                               buffer[buf_idx++], PRJ_MHD_FIDELITY_SAME);
                }
            }
        }
    }

    free(buffer);
}

static void prj_boundary_restrict_bf_to_coarse(const prj_block *fine,
    prj_block *coarse, int use_bf1, const prj_neighbor *slot)
{
    int dir;
    int max_count = 0;
    double *buffer;

    prj_boundary_check_bf_storage(fine, "prj_boundary_restrict_bf_to_coarse");
    prj_boundary_check_bf_storage(coarse, "prj_boundary_restrict_bf_to_coarse");

    for (dir = 0; dir < 3; ++dir) {
        int di = slot->send_loc_end[0] - slot->send_loc_start[0] + (dir == 0 ? 1 : 0) + 1;
        int dj = slot->send_loc_end[1] - slot->send_loc_start[1] + (dir == 1 ? 1 : 0) + 1;
        int dk = slot->send_loc_end[2] - slot->send_loc_start[2] + (dir == 2 ? 1 : 0) + 1;
        int count = di * dj * dk;
        if (count > max_count) max_count = count;
    }

    buffer = (double *)malloc(max_count * sizeof(double));
    if (buffer == 0) {
        fprintf(stderr, "prj_boundary_restrict_bf_to_coarse: malloc failed\n");
        exit(EXIT_FAILURE);
    }

    for (dir = 0; dir < 3; ++dir) {
        const double *src = prj_boundary_bf_array_const(fine, dir, use_bf1);
        int dj = slot->send_loc_end[1] - slot->send_loc_start[1] + (dir == 1 ? 1 : 0) + 1;
        int dk = slot->send_loc_end[2] - slot->send_loc_start[2] + (dir == 2 ? 1 : 0) + 1;
        int buf_idx = 0;
        int tan0 = (dir + 1) % 3;
        int tan1 = (dir + 2) % 3;
        int si;
        int sj;
        int sk;
        int i;
        int j;
        int k;

        for (si = slot->send_loc_start[0]; si <= slot->send_loc_end[0] + (dir == 0 ? 1 : 0); ++si) {
            for (sj = slot->send_loc_start[1]; sj <= slot->send_loc_end[1] + (dir == 1 ? 1 : 0); ++sj) {
                for (sk = slot->send_loc_start[2]; sk <= slot->send_loc_end[2] + (dir == 2 ? 1 : 0); ++sk) {
                    buffer[buf_idx++] = src[FACE_IDX(dir, si, sj, sk)];
                }
            }
        }

        buf_idx = 0;
        for (i = 0; i < slot->recv_loc_end[0]-slot->recv_loc_start[0]+((dir==0)?1:0); ++i) {
            for (j = 0; j < slot->recv_loc_end[1]-slot->recv_loc_start[1]+((dir==1)?1:0); ++j) {
                for (k = 0; k < slot->recv_loc_end[2]-slot->recv_loc_start[2]+((dir==2)?1:0); ++k) {
                    int it_send0[3] = {0, 0, 0};
                    int it_send1[3] = {0, 0, 0};
                    int it_send2[3] = {0, 0, 0};
                    int it_send3[3] = {0, 0, 0};
                    int b0, b1, b2, b3;
                    double value;

                    it_send0[0] = 2*i+slot->send_loc_start[0];
                    it_send0[1] = 2*j+slot->send_loc_start[1];
                    it_send0[2] = 2*k+slot->send_loc_start[2];

                    it_send1[dir]  = it_send0[dir];
                    it_send1[tan0] = it_send0[tan0];
                    it_send1[tan1] = it_send0[tan1]+1;

                    it_send2[dir]  = it_send0[dir];
                    it_send2[tan0] = it_send0[tan0]+1;
                    it_send2[tan1] = it_send0[tan1];

                    it_send3[dir]  = it_send0[dir];
                    it_send3[tan0] = it_send0[tan0]+1;
                    it_send3[tan1] = it_send0[tan1]+1;

                    b0 = ((it_send0[0]-slot->send_loc_start[0])*dj + (it_send0[1]-slot->send_loc_start[1]))*dk + (it_send0[2]-slot->send_loc_start[2]);
                    b1 = ((it_send1[0]-slot->send_loc_start[0])*dj + (it_send1[1]-slot->send_loc_start[1]))*dk + (it_send1[2]-slot->send_loc_start[2]);
                    b2 = ((it_send2[0]-slot->send_loc_start[0])*dj + (it_send2[1]-slot->send_loc_start[1]))*dk + (it_send2[2]-slot->send_loc_start[2]);
                    b3 = ((it_send3[0]-slot->send_loc_start[0])*dj + (it_send3[1]-slot->send_loc_start[1]))*dk + (it_send3[2]-slot->send_loc_start[2]);

                    value = 0.25*(buffer[b0] + buffer[b1] + buffer[b2] + buffer[b3]);
                    prj_boundary_write_bf_face(coarse, use_bf1, dir,
                                               i+slot->recv_loc_start[0],
                                               j+slot->recv_loc_start[1],
                                               k+slot->recv_loc_start[2],
                                               value, PRJ_MHD_FIDELITY_FINER);
                }
            }
        }
    }

    free(buffer);
}

static inline double prj_boundary_buf_read(const double *buf,
    const int buf_lo[3], const int buf_n[3], int i, int j, int k)
{
    return buf[((i - buf_lo[0]) * buf_n[1] + (j - buf_lo[1])) * buf_n[2] +
               (k - buf_lo[2])];
}

static inline double prj_boundary_buf_minmod(double left, double center,
    double right, double dx)
{
    double sl = (center - left) / dx;
    double sr = (right - center) / dx;

    if (sl * sr <= 0.0) return 0.0;
    return fabs(sl) < fabs(sr) ? sl : sr;
}

static inline double prj_boundary_buf_sign_half(int bit)
{
    return bit == 0 ? -1.0 : 1.0;
}

static inline double prj_boundary_interp_x1_buf(const double *buf,
    const int buf_lo[3], const int buf_n[3], const double dx[3],
    int i, int j, int k, int fine_j, int fine_k, double area)
{
    double base = prj_boundary_buf_read(buf, buf_lo, buf_n, i, j, k);
    double sy = prj_boundary_buf_minmod(
        prj_boundary_buf_read(buf, buf_lo, buf_n, i, j - 1, k),
        base,
        prj_boundary_buf_read(buf, buf_lo, buf_n, i, j + 1, k),
        dx[1]);
    double sz = prj_boundary_buf_minmod(
        prj_boundary_buf_read(buf, buf_lo, buf_n, i, j, k - 1),
        base,
        prj_boundary_buf_read(buf, buf_lo, buf_n, i, j, k + 1),
        dx[2]);
    double value = base +
        sy * (0.25 * prj_boundary_buf_sign_half(fine_j) * dx[1]) +
        sz * (0.25 * prj_boundary_buf_sign_half(fine_k) * dx[2]);

    return value * area;
}

static inline double prj_boundary_interp_x2_buf(const double *buf,
    const int buf_lo[3], const int buf_n[3], const double dx[3],
    int i, int j, int k, int fine_i, int fine_k, double area)
{
    double base = prj_boundary_buf_read(buf, buf_lo, buf_n, i, j, k);
    double sx = prj_boundary_buf_minmod(
        prj_boundary_buf_read(buf, buf_lo, buf_n, i - 1, j, k),
        base,
        prj_boundary_buf_read(buf, buf_lo, buf_n, i + 1, j, k),
        dx[0]);
    double sz = prj_boundary_buf_minmod(
        prj_boundary_buf_read(buf, buf_lo, buf_n, i, j, k - 1),
        base,
        prj_boundary_buf_read(buf, buf_lo, buf_n, i, j, k + 1),
        dx[2]);
    double value = base +
        sx * (0.25 * prj_boundary_buf_sign_half(fine_i) * dx[0]) +
        sz * (0.25 * prj_boundary_buf_sign_half(fine_k) * dx[2]);

    return value * area;
}

static inline double prj_boundary_interp_x3_buf(const double *buf,
    const int buf_lo[3], const int buf_n[3], const double dx[3],
    int i, int j, int k, int fine_i, int fine_j, double area)
{
    double base = prj_boundary_buf_read(buf, buf_lo, buf_n, i, j, k);
    double sx = prj_boundary_buf_minmod(
        prj_boundary_buf_read(buf, buf_lo, buf_n, i - 1, j, k),
        base,
        prj_boundary_buf_read(buf, buf_lo, buf_n, i + 1, j, k),
        dx[0]);
    double sy = prj_boundary_buf_minmod(
        prj_boundary_buf_read(buf, buf_lo, buf_n, i, j - 1, k),
        base,
        prj_boundary_buf_read(buf, buf_lo, buf_n, i, j + 1, k),
        dx[1]);
    double value = base +
        sx * (0.25 * prj_boundary_buf_sign_half(fine_i) * dx[0]) +
        sy * (0.25 * prj_boundary_buf_sign_half(fine_j) * dx[1]);

    return value * area;
}

static inline void prj_boundary_compute_inner_fluxes(double u[3][2][2],
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
        double si = prj_boundary_buf_sign_half(i);
        for (j = 0; j < 2; ++j) {
            double sj = prj_boundary_buf_sign_half(j);
            for (k = 0; k < 2; ++k) {
                double sk = prj_boundary_buf_sign_half(k);
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
        double sj = prj_boundary_buf_sign_half(j);
        for (k = 0; k < 2; ++k) {
            double sk = prj_boundary_buf_sign_half(k);
            u[1][j][k] = 0.5 * (u[2][j][k] + u[0][j][k]) +
                Uxx + sk * dx3s * Vxyz + sj * dx2s * Wxyz;
        }
    }
    for (i = 0; i < 2; ++i) {
        double si = prj_boundary_buf_sign_half(i);
        for (k = 0; k < 2; ++k) {
            double sk = prj_boundary_buf_sign_half(k);
            v[i][1][k] = 0.5 * (v[i][2][k] + v[i][0][k]) +
                Vyy + si * dx1s * Wxyz + sk * dx3s * Uxyz;
        }
    }
    for (i = 0; i < 2; ++i) {
        double si = prj_boundary_buf_sign_half(i);
        for (j = 0; j < 2; ++j) {
            double sj = prj_boundary_buf_sign_half(j);
            w[i][j][1] = 0.5 * (w[i][j][2] + w[i][j][0]) +
                Wzz + sj * dx2s * Uxyz + si * dx1s * Vxyz;
        }
    }
}

static void prj_boundary_prolong_bf_from_buffer(const double *buf[3],
    const int buf_lo[3][3], const int buf_n[3][3], const double coarse_dx[3],
    prj_block *fine, int ci, int cj, int ck, int fi, int fj, int fk,
    int use_bf1)
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

    prj_boundary_check_bf_storage(fine, "prj_boundary_prolong_bf_from_buffer");

    for (d = 0; d < 3; ++d) {
        dst[d] = prj_boundary_bf_array(fine, d, use_bf1);
    }

    area_u = fine->dx[1] * fine->dx[2];
    area_v = fine->dx[0] * fine->dx[2];
    area_w = fine->dx[0] * fine->dx[1];

    for (j = 0; j < 2; ++j) {
        for (k = 0; k < 2; ++k) {
            double flux;
            int idx;

            flux = prj_boundary_interp_x1_buf(buf[X1DIR], buf_lo[X1DIR],
                buf_n[X1DIR], coarse_dx, ci, cj, ck, j, k, area_u);
            idx = FACE_IDX(X1DIR, fi, fj + j, fk + k);
            if (fine->face_fidelity[X1DIR][idx] > PRJ_MHD_FIDELITY_COARSER) {
                u[0][j][k] = dst[X1DIR][idx] * area_u;
            } else {
                u[0][j][k] = flux;
                dst[X1DIR][idx] = flux / area_u;
                fine->face_fidelity[X1DIR][idx] = PRJ_MHD_FIDELITY_COARSER;
            }

            flux = prj_boundary_interp_x1_buf(buf[X1DIR], buf_lo[X1DIR],
                buf_n[X1DIR], coarse_dx, ci + 1, cj, ck, j, k, area_u);
            idx = FACE_IDX(X1DIR, fi + 2, fj + j, fk + k);
            if (fine->face_fidelity[X1DIR][idx] > PRJ_MHD_FIDELITY_COARSER) {
                u[2][j][k] = dst[X1DIR][idx] * area_u;
            } else {
                u[2][j][k] = flux;
                dst[X1DIR][idx] = flux / area_u;
                fine->face_fidelity[X1DIR][idx] = PRJ_MHD_FIDELITY_COARSER;
            }
        }
    }

    for (i = 0; i < 2; ++i) {
        for (k = 0; k < 2; ++k) {
            double flux;
            int idx;

            flux = prj_boundary_interp_x2_buf(buf[X2DIR], buf_lo[X2DIR],
                buf_n[X2DIR], coarse_dx, ci, cj, ck, i, k, area_v);
            idx = FACE_IDX(X2DIR, fi + i, fj, fk + k);
            if (fine->face_fidelity[X2DIR][idx] > PRJ_MHD_FIDELITY_COARSER) {
                v[i][0][k] = dst[X2DIR][idx] * area_v;
            } else {
                v[i][0][k] = flux;
                dst[X2DIR][idx] = flux / area_v;
                fine->face_fidelity[X2DIR][idx] = PRJ_MHD_FIDELITY_COARSER;
            }

            flux = prj_boundary_interp_x2_buf(buf[X2DIR], buf_lo[X2DIR],
                buf_n[X2DIR], coarse_dx, ci, cj + 1, ck, i, k, area_v);
            idx = FACE_IDX(X2DIR, fi + i, fj + 2, fk + k);
            if (fine->face_fidelity[X2DIR][idx] > PRJ_MHD_FIDELITY_COARSER) {
                v[i][2][k] = dst[X2DIR][idx] * area_v;
            } else {
                v[i][2][k] = flux;
                dst[X2DIR][idx] = flux / area_v;
                fine->face_fidelity[X2DIR][idx] = PRJ_MHD_FIDELITY_COARSER;
            }
        }
    }

    for (i = 0; i < 2; ++i) {
        for (j = 0; j < 2; ++j) {
            double flux;
            int idx;

            flux = prj_boundary_interp_x3_buf(buf[X3DIR], buf_lo[X3DIR],
                buf_n[X3DIR], coarse_dx, ci, cj, ck, i, j, area_w);
            idx = FACE_IDX(X3DIR, fi + i, fj + j, fk);
            if (fine->face_fidelity[X3DIR][idx] > PRJ_MHD_FIDELITY_COARSER) {
                w[i][j][0] = dst[X3DIR][idx] * area_w;
            } else {
                w[i][j][0] = flux;
                dst[X3DIR][idx] = flux / area_w;
                fine->face_fidelity[X3DIR][idx] = PRJ_MHD_FIDELITY_COARSER;
            }

            flux = prj_boundary_interp_x3_buf(buf[X3DIR], buf_lo[X3DIR],
                buf_n[X3DIR], coarse_dx, ci, cj, ck + 1, i, j, area_w);
            idx = FACE_IDX(X3DIR, fi + i, fj + j, fk + 2);
            if (fine->face_fidelity[X3DIR][idx] > PRJ_MHD_FIDELITY_COARSER) {
                w[i][j][2] = dst[X3DIR][idx] * area_w;
            } else {
                w[i][j][2] = flux;
                dst[X3DIR][idx] = flux / area_w;
                fine->face_fidelity[X3DIR][idx] = PRJ_MHD_FIDELITY_COARSER;
            }
        }
    }

    prj_boundary_compute_inner_fluxes(u, v, w,
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

static void prj_boundary_prolong_bf_to_fine(const prj_block *coarse,
    prj_block *fine, int use_bf1, const prj_neighbor *slot)
{
    int dir;
    int max_count = 0;
    double *buf[3];
    int buf_lo[3][3];
    int buf_n[3][3];
    int i;
    int j;
    int k;

    prj_boundary_check_bf_storage(coarse, "prj_boundary_prolong_bf_to_fine");
    prj_boundary_check_bf_storage(fine, "prj_boundary_prolong_bf_to_fine");

    for (dir = 0; dir < 3; ++dir) {
        int d;
        int count;
        for (d = 0; d < 3; ++d) {
            if (d == dir) {
                buf_lo[dir][d] = slot->send_loc_start[d];
                buf_n[dir][d] = slot->send_loc_end[d] - slot->send_loc_start[d] + 2;
            } else {
                buf_lo[dir][d] = slot->send_loc_start[d] - 1;
                buf_n[dir][d] = slot->send_loc_end[d] - slot->send_loc_start[d] + 3;
            }
        }
        count = buf_n[dir][0] * buf_n[dir][1] * buf_n[dir][2];
        if (count > max_count) max_count = count;
    }

    buf[0] = (double *)malloc(max_count * sizeof(double));
    buf[1] = (double *)malloc(max_count * sizeof(double));
    buf[2] = (double *)malloc(max_count * sizeof(double));
    if (buf[0] == 0 || buf[1] == 0 || buf[2] == 0) {
        fprintf(stderr, "prj_boundary_prolong_bf_to_fine: malloc failed\n");
        exit(EXIT_FAILURE);
    }

    for (dir = 0; dir < 3; ++dir) {
        const double *src = prj_boundary_bf_array_const(coarse, dir, use_bf1);
        int buf_idx = 0;
        int lo[3];
        int hi[3];
        int si;
        int sj;
        int sk;
        for (i = 0; i < 3; ++i) {
            lo[i] = buf_lo[dir][i];
            hi[i] = lo[i] + buf_n[dir][i] - 1;
        }
        for (si = lo[0]; si <= hi[0]; ++si) {
            for (sj = lo[1]; sj <= hi[1]; ++sj) {
                for (sk = lo[2]; sk <= hi[2]; ++sk) {
                    buf[dir][buf_idx++] = src[FACE_IDX(dir, si, sj, sk)];
                }
            }
        }
    }

    for (i = 0; i < slot->recv_loc_end[0]-slot->recv_loc_start[0]; i+=2) {
        for (j = 0; j < slot->recv_loc_end[1]-slot->recv_loc_start[1]; j+=2) {
            for (k = 0; k < slot->recv_loc_end[2]-slot->recv_loc_start[2]; k+=2) {
                const double *cbuf[3];
                int fi = i+slot->recv_loc_start[0];
                int fj = j+slot->recv_loc_start[1];
                int fk = k+slot->recv_loc_start[2];
                int ci = i/2+slot->send_loc_start[0];
                int cj = j/2+slot->send_loc_start[1];
                int ck = k/2+slot->send_loc_start[2];
                cbuf[0] = buf[0];
                cbuf[1] = buf[1];
                cbuf[2] = buf[2];
                prj_boundary_prolong_bf_from_buffer(cbuf, buf_lo, buf_n,
                    coarse->dx, fine, ci, cj, ck, fi, fj, fk, use_bf1);
            }
        }
    }

    free(buf[0]);
    free(buf[1]);
    free(buf[2]);
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

            for (n = 0; n < PRJ_BLOCK_NFACES; ++n) {
                block->face_fidelity[dir][n] = PRJ_MHD_FIDELITY_NONE;
            }

            for (i = 0; i <= prj_boundary_bf_axis_active_max(dir, 0); ++i) {
                for (j = 0; j <= prj_boundary_bf_axis_active_max(dir, 1); ++j) {
                    for (k = 0; k <= prj_boundary_bf_axis_active_max(dir, 2); ++k) {
                        block->face_fidelity[dir][FACE_IDX(dir, i, j, k)] = PRJ_MHD_FIDELITY_SAME;
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
                    if (!prj_boundary_face_storage_index_ok(dir, src_idx[0], src_idx[1], src_idx[2])) {
                        continue;
                    }
                    src_flat = FACE_IDX(dir, src_idx[0], src_idx[1], src_idx[2]);
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
        

        {
            const prj_neighbor *slot = &block->slot[n];

            switch (fill_kind) {
            case PRJ_BOUNDARY_FILL_SAME_LEVEL:
                if (slot->rel_level == 0) {
                    prj_boundary_copy_bf_same_level(block, neighbor, use_bf1, slot);
                }
                break;
            case PRJ_BOUNDARY_FILL_RESTRICTION:
                if (slot->rel_level < 0) {
                    prj_boundary_restrict_bf_to_coarse(block, neighbor, use_bf1, slot);
                }
                break;
            case PRJ_BOUNDARY_FILL_PROLONGATION:
                if (slot->rel_level > 0) {
                    prj_boundary_prolong_bf_to_fine(block, neighbor, use_bf1, slot);
                }
                break;
            case PRJ_BOUNDARY_FILL_NONRECON:
                if (slot->rel_level == 0) {
                    prj_boundary_copy_bf_same_level(block, neighbor, use_bf1, slot);
                } else if (slot->rel_level < 0) {
                    prj_boundary_restrict_bf_to_coarse(block, neighbor, use_bf1, slot);
                }
                break;
            case PRJ_BOUNDARY_FILL_RECON:
                if (slot->rel_level > 0) {
                    prj_boundary_prolong_bf_to_fine(block, neighbor, use_bf1, slot);
                }
                break;
            case PRJ_BOUNDARY_FILL_ALL:
                if (slot->rel_level == 0) {
                    prj_boundary_copy_bf_same_level(block, neighbor, use_bf1, slot);
                } else if (slot->rel_level < 0) {
                    prj_boundary_restrict_bf_to_coarse(block, neighbor, use_bf1, slot);
                } else {
                    prj_boundary_prolong_bf_to_fine(block, neighbor, use_bf1, slot);
                }
                break;
            default:
                fprintf(stderr, "prj_boundary_send_bf: invalid fill_kind=%d\n", fill_kind);
                exit(EXIT_FAILURE);
            }
        }
    }
}

void prj_boundary_fill_bf(prj_mesh *mesh, const prj_bc *bc, int use_bf1, prj_eos *eos)
{
    prj_mpi *mpi = prj_mpi_current();
    const int fill_passes[3] = {
        PRJ_BOUNDARY_FILL_SAME_LEVEL,
        PRJ_BOUNDARY_FILL_RESTRICTION,
        PRJ_BOUNDARY_FILL_PROLONGATION
    };
    int i;
    int pass;

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
    for (pass = 0; pass < 2; ++pass) {
        int fill_kind = fill_passes[pass];

        for (i = 0; i < mesh->nblocks; ++i) {
            if (prj_boundary_active_block(&mesh->blocks[i])) {
                prj_boundary_send_bf(&mesh->blocks[i], use_bf1, fill_kind);
            }
        }
    }
    prj_mpi_exchange_bf(mesh, mpi, use_bf1, PRJ_BOUNDARY_FILL_NONRECON);
    for (i = 0; i < mesh->nblocks; ++i) {
        if (prj_boundary_active_block(&mesh->blocks[i])) {
            prj_boundary_send_bf(&mesh->blocks[i], use_bf1, fill_passes[2]);
        }
    }
    prj_mpi_exchange_bf(mesh, mpi, use_bf1, PRJ_BOUNDARY_FILL_RECON);
    for (i = 0; i < mesh->nblocks; ++i) {
        if (prj_boundary_active_block(&mesh->blocks[i])) {
            prj_boundary_physical_bf(bc, &mesh->blocks[i], use_bf1);
        }
    }
    for (i = 0; i < mesh->nblocks; ++i) {
        if (prj_boundary_active_block(&mesh->blocks[i])) {
            prj_mhd_bf2bc(eos, &mesh->blocks[i], use_bf1);
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
