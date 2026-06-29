#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "prj.h"

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
    return prj_block_prim_stage(block, prj_stage_slot_from_stage_arg(stage));
}

static const double *prj_boundary_stage_array_const(const prj_block *block, int stage)
{
    return prj_block_prim_stage_const(block, prj_stage_slot_from_stage_arg(stage));
}

static int prj_boundary_active_block(const prj_mpi *mpi, const prj_block *block)
{
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
    PRJ_BOUNDARY_PHYS_FACE_ONLY = 0,  /* pure-face ghosts (1 transverse-interior slab) */
    PRJ_BOUNDARY_PHYS_ALL = 1,        /* faces + edges + corners */
    PRJ_BOUNDARY_PHYS_EDGE_CORNER = 2 /* edges + corners only (faces done separately) */
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
    return src[is_eosvar != 0 ? EIDX(var, i, j, k) : WIDX(var, i, j, k)];
}

static int prj_boundary_use_BJ(const prj_mesh *mesh)
{
    return mesh != 0 && mesh->use_BJ != 0;
}

static double prj_boundary_prolongate_value(const double *src, int var,
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
                    prj_boundary_read_value(src, var, i + di, j + dj, k + dk, is_eosvar);
            }
        }
    }
    return prj_reconstruct_cell_for_prolongate(stencil, target, use_BJ);
}

static void prj_boundary_get_values(const prj_block *block, const double *src, int nvar, int is_eosvar,
    int use_BJ, const char *label, double x1, double x2, double x3, double *dst)
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
            double xcenter = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
            double ycenter = block->xmin[1] + ((double)j + 0.5) * block->dx[1];
            double zcenter = block->xmin[2] + ((double)k + 0.5) * block->dx[2];
            double target[3];

            target[0] = (x1 - xcenter) / block->dx[0];
            target[1] = (x2 - ycenter) / block->dx[1];
            target[2] = (x3 - zcenter) / block->dx[2];
            dst[v] = prj_boundary_prolongate_value(src, v, i, j, k, is_eosvar, target, use_BJ);
        } else {
            fprintf(stderr, "%s: mixed unsupported sample location (%g, %g, %g)\n", label, x1, x2, x3);
            exit(EXIT_FAILURE);
        }
    }
}

void prj_boundary_get_prim(const prj_mesh *mesh, const prj_block *block, int stage,
    double x1, double x2, double x3, double *w)
{
    prj_boundary_get_values(block, prj_boundary_stage_array_const(block, stage),
        PRJ_NVAR_PRIM, 0, prj_boundary_use_BJ(mesh), "prj_boundary_get_prim", x1, x2, x3, w);
}

void prj_boundary_get_eosvar(const prj_mesh *mesh, const prj_block *block,
    double x1, double x2, double x3, double *eosv)
{
    prj_boundary_get_values(block, block != 0 ? block->eosvar : 0,
        PRJ_NVAR_EOSVAR, 1, prj_boundary_use_BJ(mesh), "prj_boundary_get_eosvar", x1, x2, x3, eosv);
}

void prj_boundary_send(prj_mesh *mesh, const prj_mpi *mpi, prj_block *block, int stage, int fill_kind)
{
    int n;
    int use_BJ = prj_boundary_use_BJ(mesh);
    
    double *W_send = prj_boundary_stage_array(block, stage);
    double *eos_send = block->eosvar;

    for (n = 0; n < 56; ++n) {
        int id = block->slot[n].id;

        if (id >= 0 && mesh != 0 && id < mesh->nblocks) {
            prj_block *neighbor = &mesh->blocks[id];
            double *W_recv = prj_boundary_stage_array(neighbor, stage);
            double *eos_recv = neighbor->eosvar;
            int i;
            int j;
            int k;

            if (!prj_boundary_active_block(mpi, neighbor)) {
                continue;
            }
            if (neighbor->rank == block->rank) {
                const prj_neighbor *slot = &block->slot[n];
                int sample_kind;
                if (slot->rel_level == 0) { sample_kind = PRJ_BOUNDARY_FILL_SAME_LEVEL; }
                else if (slot->rel_level < 0) { sample_kind = PRJ_BOUNDARY_FILL_RESTRICTION; }
                else { sample_kind = PRJ_BOUNDARY_FILL_PROLONGATION; }
                if (sample_kind != fill_kind) {
                    continue;
                }

                /* Hydro primitives + EOS vars: iterate the full recv box. */
                for (i = 0; i < slot->recv_loc_end[0]-slot->recv_loc_start[0]; ++i) {
                    for (j = 0; j < slot->recv_loc_end[1]-slot->recv_loc_start[1]; ++j) {
                        for (k = 0; k < slot->recv_loc_end[2]-slot->recv_loc_start[2]; ++k) {
                            int v;

                            if (slot->rel_level==0){
                                // Same level
                                for (v = 0; v < PRJ_NHYDRO; ++v) {
                                    W_recv[WIDX(v, i+slot->recv_loc_start[0], j+slot->recv_loc_start[1], k+slot->recv_loc_start[2])] =
                                        W_send[WIDX(v, i+slot->send_loc_start[0], j+slot->send_loc_start[1], k+slot->send_loc_start[2])];
                                }
                                for (v = 0; v < PRJ_NVAR_EOSVAR; ++v) {
                                    eos_recv[EIDX(v, i+slot->recv_loc_start[0], j+slot->recv_loc_start[1], k+slot->recv_loc_start[2])] =
                                        eos_send[EIDX(v, i+slot->send_loc_start[0], j+slot->send_loc_start[1], k+slot->send_loc_start[2])];
                                }
                            } else if (slot->rel_level==-1) {
                                // Neighbor is coarser, restriction
                                for (v = 0; v < PRJ_NHYDRO; ++v) {
                                    W_recv[WIDX(v, i+slot->recv_loc_start[0], j+slot->recv_loc_start[1], k+slot->recv_loc_start[2])] =
                                        0.125*
                                        (W_send[WIDX(v, 2*i+slot->send_loc_start[0],
                                                        2*j+slot->send_loc_start[1],
                                                        2*k+slot->send_loc_start[2])]
                                        +W_send[WIDX(v, 2*i+slot->send_loc_start[0],
                                                        2*j+slot->send_loc_start[1],
                                                        2*k+slot->send_loc_start[2]+1)]
                                        +W_send[WIDX(v, 2*i+slot->send_loc_start[0],
                                                        2*j+slot->send_loc_start[1]+1,
                                                        2*k+slot->send_loc_start[2])]
                                        +W_send[WIDX(v, 2*i+slot->send_loc_start[0],
                                                        2*j+slot->send_loc_start[1]+1,
                                                        2*k+slot->send_loc_start[2]+1)]
                                        +W_send[WIDX(v, 2*i+slot->send_loc_start[0]+1,
                                                        2*j+slot->send_loc_start[1],
                                                        2*k+slot->send_loc_start[2])]
                                        +W_send[WIDX(v, 2*i+slot->send_loc_start[0]+1,
                                                        2*j+slot->send_loc_start[1],
                                                        2*k+slot->send_loc_start[2]+1)]
                                        +W_send[WIDX(v, 2*i+slot->send_loc_start[0]+1,
                                                        2*j+slot->send_loc_start[1]+1,
                                                        2*k+slot->send_loc_start[2])]
                                        +W_send[WIDX(v, 2*i+slot->send_loc_start[0]+1,
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
                                double target[3];
                                int ai = i + slot->recv_loc_start[0];
                                int aj = j + slot->recv_loc_start[1];
                                int ak = k + slot->recv_loc_start[2];

                                target[0] = (ai % 2 == 0) ? -0.25 : 0.25;
                                target[1] = (aj % 2 == 0) ? -0.25 : 0.25;
                                target[2] = (ak % 2 == 0) ? -0.25 : 0.25;
                                for (v = 0; v < PRJ_NHYDRO; ++v) {
                                    W_recv[WIDX(v, ai, aj, ak)] =
                                        prj_boundary_prolongate_value(W_send, v,
                                            i/2+slot->send_loc_start[0],
                                            j/2+slot->send_loc_start[1],
                                            k/2+slot->send_loc_start[2],
                                            0, target, use_BJ);
                                }
                                for (v = 0; v < PRJ_NVAR_EOSVAR; ++v) {
                                    eos_recv[EIDX(v, ai, aj, ak)] =
                                        prj_boundary_prolongate_value(eos_send, v,
                                            i/2+slot->send_loc_start[0],
                                            j/2+slot->send_loc_start[1],
                                            k/2+slot->send_loc_start[2],
                                            1, target, use_BJ);
                                }
                            } else {
                              fprintf(stderr,"slot->rel_level unrecognized: %d\n", slot->rel_level);
                              exit(1);
                            }
                        }
                    }
                }

#if PRJ_NRAD > 0
                /* Radiation primitives: iterate the precomputed rad-clipped box. */
                for (i = 0; i < slot->recv_loc_end_rad[0]-slot->recv_loc_start_rad[0]; ++i) {
                    for (j = 0; j < slot->recv_loc_end_rad[1]-slot->recv_loc_start_rad[1]; ++j) {
                        for (k = 0; k < slot->recv_loc_end_rad[2]-slot->recv_loc_start_rad[2]; ++k) {
                            int v;

                            if (slot->rel_level==0){
                                for (v = PRJ_NHYDRO; v < PRJ_NVAR_PRIM; ++v) {
                                    W_recv[WIDX(v, i+slot->recv_loc_start_rad[0], j+slot->recv_loc_start_rad[1], k+slot->recv_loc_start_rad[2])] =
                                        W_send[WIDX(v, i+slot->send_loc_start_rad[0], j+slot->send_loc_start_rad[1], k+slot->send_loc_start_rad[2])];
                                }
                            } else if (slot->rel_level==-1) {
                                for (v = PRJ_NHYDRO; v < PRJ_NVAR_PRIM; ++v) {
                                    W_recv[WIDX(v, i+slot->recv_loc_start_rad[0], j+slot->recv_loc_start_rad[1], k+slot->recv_loc_start_rad[2])] =
                                        0.125*
                                        (W_send[WIDX(v, 2*i+slot->send_loc_start_rad[0],
                                                        2*j+slot->send_loc_start_rad[1],
                                                        2*k+slot->send_loc_start_rad[2])]
                                        +W_send[WIDX(v, 2*i+slot->send_loc_start_rad[0],
                                                        2*j+slot->send_loc_start_rad[1],
                                                        2*k+slot->send_loc_start_rad[2]+1)]
                                        +W_send[WIDX(v, 2*i+slot->send_loc_start_rad[0],
                                                        2*j+slot->send_loc_start_rad[1]+1,
                                                        2*k+slot->send_loc_start_rad[2])]
                                        +W_send[WIDX(v, 2*i+slot->send_loc_start_rad[0],
                                                        2*j+slot->send_loc_start_rad[1]+1,
                                                        2*k+slot->send_loc_start_rad[2]+1)]
                                        +W_send[WIDX(v, 2*i+slot->send_loc_start_rad[0]+1,
                                                        2*j+slot->send_loc_start_rad[1],
                                                        2*k+slot->send_loc_start_rad[2])]
                                        +W_send[WIDX(v, 2*i+slot->send_loc_start_rad[0]+1,
                                                        2*j+slot->send_loc_start_rad[1],
                                                        2*k+slot->send_loc_start_rad[2]+1)]
                                        +W_send[WIDX(v, 2*i+slot->send_loc_start_rad[0]+1,
                                                        2*j+slot->send_loc_start_rad[1]+1,
                                                        2*k+slot->send_loc_start_rad[2])]
                                        +W_send[WIDX(v, 2*i+slot->send_loc_start_rad[0]+1,
                                                        2*j+slot->send_loc_start_rad[1]+1,
                                                        2*k+slot->send_loc_start_rad[2]+1)]);
                                }
                            } else if (slot->rel_level==1) {
                                /* Fine-block offset parity is determined by the (rad-clipped) recv index;
                                 * recv_loc_start_rad differs from recv_loc_start by an even amount,
                                 * so absolute parity of (recv_loc_start_rad + i) matches the original. */
                                double target[3];
                                int ai = i + slot->recv_loc_start_rad[0];
                                int aj = j + slot->recv_loc_start_rad[1];
                                int ak = k + slot->recv_loc_start_rad[2];

                                target[0] = (ai % 2 == 0) ? -0.25 : 0.25;
                                target[1] = (aj % 2 == 0) ? -0.25 : 0.25;
                                target[2] = (ak % 2 == 0) ? -0.25 : 0.25;
                                for (v = PRJ_NHYDRO; v < PRJ_NVAR_PRIM; ++v) {
                                    W_recv[WIDX(v, ai, aj, ak)] =
                                        prj_boundary_prolongate_value(W_send, v,
                                            i/2+slot->send_loc_start_rad[0],
                                            j/2+slot->send_loc_start_rad[1],
                                            k/2+slot->send_loc_start_rad[2],
                                            0, target, use_BJ);
                                }
                            }
                        }
                    }
                }
#endif
            }
        }
    }
}

static int prj_boundary_idx_outside(int idx)
{
    return idx < 0 || idx >= PRJ_BLOCK_SIZE;
}

/* Copy one boundary band for variables [v_begin, v_end).  Iterates only the
 * NGHOST-thick slab normal to `axis` on `side` (no full-cube walk).  The
 * transverse extent is set by `region`:
 *   FACE_ONLY    -> transverse interior [0, BLOCK_SIZE)            (pure faces)
 *   EDGE_CORNER  -> full transverse band, skipping pure-face cells (edges/corners)
 *   ALL          -> full transverse band                          (faces+edges+corners)
 * `do_special` enables the velocity sign-flip / user override (hydro band only). */
static void prj_boundary_apply_axis_band(double *dst, int axis, int side, int bc_type,
    int region, int ng, int v_begin, int v_end, int do_special)
{
    const int axis1 = (axis + 1) % 3;
    const int axis2 = (axis + 2) % 3;
    const int normal_v = axis == 0 ? PRJ_PRIM_V1 : (axis == 1 ? PRJ_PRIM_V2 : PRJ_PRIM_V3);
    const int n_lo = side == 0 ? -ng : PRJ_BLOCK_SIZE;
    const int n_hi = side == 0 ? 0 : PRJ_BLOCK_SIZE + ng; /* exclusive */
    const int t_lo = region == PRJ_BOUNDARY_PHYS_FACE_ONLY ? 0 : -ng;
    const int t_hi = region == PRJ_BOUNDARY_PHYS_FACE_ONLY
        ? PRJ_BLOCK_SIZE : PRJ_BLOCK_SIZE + ng; /* exclusive */
    int a;
    int t1;
    int t2;

    for (a = n_lo; a < n_hi; ++a) {
        int src_a = side == 0 ? -a - 1 : 2 * PRJ_BLOCK_SIZE - 1 - a;

        for (t1 = t_lo; t1 < t_hi; ++t1) {
            for (t2 = t_lo; t2 < t_hi; ++t2) {
                int idx[3];
                int src_idx[3];
                int v;

                /* Edge/corner pass leaves the pure-face cells to the face pass. */
                if (region == PRJ_BOUNDARY_PHYS_EDGE_CORNER &&
                    !prj_boundary_idx_outside(t1) && !prj_boundary_idx_outside(t2)) {
                    continue;
                }
                idx[axis] = a;
                idx[axis1] = t1;
                idx[axis2] = t2;
                src_idx[axis] = src_a;
                src_idx[axis1] = t1;
                src_idx[axis2] = t2;
                for (v = v_begin; v < v_end; ++v) {
                    dst[WIDX(v, idx[0], idx[1], idx[2])] =
                        dst[WIDX(v, src_idx[0], src_idx[1], src_idx[2])];
                }
                if (do_special) {
                    if (bc_type == PRJ_BC_REFLECT) {
                        dst[WIDX(normal_v, idx[0], idx[1], idx[2])] =
                            -dst[WIDX(normal_v, src_idx[0], src_idx[1], src_idx[2])];
                    } else if (bc_type == PRJ_BC_USER) {
                        if (axis == 0 && side == 0) {
                            dst[WIDX(PRJ_PRIM_V1, idx[0], idx[1], idx[2])] = 5.0;
                        }
                    }
                }
            }
        }
    }
}

static void prj_boundary_apply_axis(double *dst, int axis, int side, int bc_type, int region)
{
    /* Hydro primitives use the full NGHOST band. */
    prj_boundary_apply_axis_band(dst, axis, side, bc_type, region, PRJ_NGHOST,
        0, PRJ_NHYDRO, 1);
#if PRJ_NRAD > 0
    /* Radiation primitives use the (possibly narrower) NGHOST_RAD band. */
    prj_boundary_apply_axis_band(dst, axis, side, bc_type, region, PRJ_NGHOST_RAD,
        PRJ_NHYDRO, PRJ_NVAR_PRIM, 0);
#endif
}

void prj_boundary_physical(const prj_mesh *mesh, const prj_bc *bc, prj_block *block, int stage, int mode)
{
    double *dst = prj_boundary_stage_array(block, stage);
    const double tol = 1.0e-12;
    int pass;
    /* Faces resolve in one pass; edges/corners need the per-axis copy to
     * propagate (faces->edges->corners), so those modes run three passes. */
    int npass = mode == PRJ_BOUNDARY_PHYS_FACE_ONLY ? 1 : 3;

    for (pass = 0; pass < npass; ++pass) {
        if (prj_abs_double(block->xmin[0] - mesh->coord.x1min) < tol) {
            prj_boundary_apply_axis(dst, 0, 0, bc->bc_x1_inner, mode);
        }
        if (prj_abs_double(block->xmax[0] - mesh->coord.x1max) < tol) {
            prj_boundary_apply_axis(dst, 0, 1, bc->bc_x1_outer, mode);
        }
        if (prj_abs_double(block->xmin[1] - mesh->coord.x2min) < tol) {
            prj_boundary_apply_axis(dst, 1, 0, bc->bc_x2_inner, mode);
        }
        if (prj_abs_double(block->xmax[1] - mesh->coord.x2max) < tol) {
            prj_boundary_apply_axis(dst, 1, 1, bc->bc_x2_outer, mode);
        }
        if (prj_abs_double(block->xmin[2] - mesh->coord.x3min) < tol) {
            prj_boundary_apply_axis(dst, 2, 0, bc->bc_x3_inner, mode);
        }
        if (prj_abs_double(block->xmax[2] - mesh->coord.x3max) < tol) {
            prj_boundary_apply_axis(dst, 2, 1, bc->bc_x3_outer, mode);
        }
    }
}

#if PRJ_MHD
static double *prj_boundary_bf_array(prj_block *block, int dir, int use_bf1)
{
    return prj_block_bf_stage(block, dir, prj_stage_slot_from_bf_arg(use_bf1));
}

static const double *prj_boundary_bf_array_const(const prj_block *block, int dir, int use_bf1)
{
    return prj_block_bf_stage_const(block, dir, prj_stage_slot_from_bf_arg(use_bf1));
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
        if (block->Bf[d] == 0) {
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

static void prj_boundary_prolong_bf_to_fine(const prj_block *coarse,
    prj_block *fine, int use_bf1, const prj_neighbor *slot, int use_BJ)
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
                prj_mhd_prolong_bf_from_buffer(cbuf, buf_lo, buf_n,
                    coarse->dx, fine, ci, cj, ck, fi, fj, fk, use_bf1,
                    use_BJ);
            }
        }
    }

    free(buf[0]);
    free(buf[1]);
    free(buf[2]);
}

static void prj_boundary_init_face_fidelity(prj_mesh *mesh, const prj_mpi *mpi)
{
    int bidx;

    if (mesh == 0) {
        return;
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        int dir;

        if (!prj_boundary_active_block(mpi, block)) {
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

static void prj_boundary_physical_bf(const prj_mesh *mesh, const prj_bc *bc, prj_block *block, int use_bf1)
{
    const double tol = 1.0e-12;
    int pass;

    if (mesh == 0 || bc == 0) {
        fprintf(stderr, "prj_boundary_physical_bf: missing boundary context\n");
        exit(EXIT_FAILURE);
    }
    for (pass = 0; pass < 3; ++pass) {
        if (prj_abs_double(block->xmin[0] - mesh->coord.x1min) < tol) {
            prj_boundary_apply_bf_axis(block, use_bf1, 0, 0, bc->bc_x1_inner);
        }
        if (prj_abs_double(block->xmax[0] - mesh->coord.x1max) < tol) {
            prj_boundary_apply_bf_axis(block, use_bf1, 0, 1, bc->bc_x1_outer);
        }
        if (prj_abs_double(block->xmin[1] - mesh->coord.x2min) < tol) {
            prj_boundary_apply_bf_axis(block, use_bf1, 1, 0, bc->bc_x2_inner);
        }
        if (prj_abs_double(block->xmax[1] - mesh->coord.x2max) < tol) {
            prj_boundary_apply_bf_axis(block, use_bf1, 1, 1, bc->bc_x2_outer);
        }
        if (prj_abs_double(block->xmin[2] - mesh->coord.x3min) < tol) {
            prj_boundary_apply_bf_axis(block, use_bf1, 2, 0, bc->bc_x3_inner);
        }
        if (prj_abs_double(block->xmax[2] - mesh->coord.x3max) < tol) {
            prj_boundary_apply_bf_axis(block, use_bf1, 2, 1, bc->bc_x3_outer);
        }
    }
}

void prj_boundary_send_bf(prj_mesh *mesh, const prj_mpi *mpi, prj_block *block, int use_bf1, int fill_kind)
{
    int n;
    int use_BJ = prj_boundary_use_BJ(mesh);

    if (!prj_boundary_active_block(mpi, block)) {
        return;
    }
    prj_boundary_check_bf_storage(block, "prj_boundary_send_bf");
    for (n = 0; n < 56; ++n) {
        int id = block->slot[n].id;
        prj_block *neighbor;

        if (id < 0 || mesh == 0 || id >= mesh->nblocks) {
            continue;
        }
        neighbor = &mesh->blocks[id];
        if (!prj_boundary_active_block(mpi, neighbor) || neighbor->rank != block->rank) {
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
                    prj_boundary_prolong_bf_to_fine(block, neighbor, use_bf1, slot, use_BJ);
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
                    prj_boundary_prolong_bf_to_fine(block, neighbor, use_bf1, slot, use_BJ);
                }
                break;
            case PRJ_BOUNDARY_FILL_ALL:
                if (slot->rel_level == 0) {
                    prj_boundary_copy_bf_same_level(block, neighbor, use_bf1, slot);
                } else if (slot->rel_level < 0) {
                    prj_boundary_restrict_bf_to_coarse(block, neighbor, use_bf1, slot);
                } else {
                    prj_boundary_prolong_bf_to_fine(block, neighbor, use_bf1, slot, use_BJ);
                }
                break;
            default:
                fprintf(stderr, "prj_boundary_send_bf: invalid fill_kind=%d\n", fill_kind);
                exit(EXIT_FAILURE);
            }
        }
    }
}

void prj_boundary_fill_bf(prj_mesh *mesh, prj_mpi *mpi, const prj_bc *bc, int use_bf1, prj_eos *eos)
{
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
    prj_boundary_init_face_fidelity(mesh, mpi);
    for (i = 0; i < mesh->nblocks; ++i) {
        if (prj_boundary_active_block(mpi, &mesh->blocks[i])) {
            prj_boundary_physical_bf(mesh, bc, &mesh->blocks[i], use_bf1);
        }
    }
    for (pass = 0; pass < 3; ++pass) {
        int fill_kind = fill_passes[pass];

        for (i = 0; i < mesh->nblocks; ++i) {
            if (prj_boundary_active_block(mpi, &mesh->blocks[i])) {
                prj_boundary_send_bf(mesh, mpi, &mesh->blocks[i], use_bf1, fill_kind);
            }
        }
        prj_mpi_exchange_bf(mesh, mpi, use_bf1, fill_kind);
    }
    for (i = 0; i < mesh->nblocks; ++i) {
        if (prj_boundary_active_block(mpi, &mesh->blocks[i])) {
            prj_boundary_physical_bf(mesh, bc, &mesh->blocks[i], use_bf1);
        }
    }
    for (i = 0; i < mesh->nblocks; ++i) {
        if (prj_boundary_active_block(mpi, &mesh->blocks[i])) {
            prj_mhd_bf2bc_all(eos, &mesh->blocks[i], use_bf1);
        }
    }
}
#endif

void prj_boundary_mpi_recv(prj_mesh *mesh, prj_mpi *mpi, int stage, int fill_kind)
{
    if (mpi == 0) {
        return;
    }
    prj_mpi_exchange_ghosts(mesh, mpi, stage, fill_kind);
}

void prj_boundary_fill_ghosts(prj_mesh *mesh, prj_mpi *mpi, const prj_bc *bc, int stage)
{
    const int fill_passes[3] = {
        PRJ_BOUNDARY_FILL_SAME_LEVEL,
        PRJ_BOUNDARY_FILL_RESTRICTION,
        PRJ_BOUNDARY_FILL_PROLONGATION
    };
    int i;
    int pass;

    for (i = 0; i < mesh->nblocks; ++i) {
        if (prj_boundary_active_block(mpi, &mesh->blocks[i])) {
            prj_boundary_physical(mesh, bc, &mesh->blocks[i], stage, PRJ_BOUNDARY_PHYS_FACE_ONLY);
        }
    }
    for (pass = 0; pass < 3; ++pass) {
        int fill_kind = fill_passes[pass];

        for (i = 0; i < mesh->nblocks; ++i) {
            if (prj_boundary_active_block(mpi, &mesh->blocks[i])) {
                prj_boundary_send(mesh, mpi, &mesh->blocks[i], stage, fill_kind);
            }
        }
        prj_boundary_mpi_recv(mesh, mpi, stage, fill_kind);
    }
    for (i = 0; i < mesh->nblocks; ++i) {
        if (prj_boundary_active_block(mpi, &mesh->blocks[i])) {
            prj_boundary_physical(mesh, bc, &mesh->blocks[i], stage, PRJ_BOUNDARY_PHYS_EDGE_CORNER);
        }
    }
}

void prj_boundary_fill_ghosts_and_bf(prj_mesh *mesh, prj_mpi *mpi, const prj_bc *bc,
    int stage, int use_bf1, prj_eos *eos, prj_grav *grav, prj_rad *rad, int timer_scope)
{
    const int fill_passes[3] = {
        PRJ_BOUNDARY_FILL_SAME_LEVEL,
        PRJ_BOUNDARY_FILL_RESTRICTION,
        PRJ_BOUNDARY_FILL_PROLONGATION
    };
    int i;
    int pass;

    (void)grav;
    (void)rad;
    (void)timer_scope;

    PRJ_SUBTIMER_START("sub_ghost_phys_face");
    for (i = 0; i < mesh->nblocks; ++i) {
        if (prj_boundary_active_block(mpi, &mesh->blocks[i])) {
            prj_boundary_physical(mesh, bc, &mesh->blocks[i], stage, PRJ_BOUNDARY_PHYS_FACE_ONLY);
        }
    }
    PRJ_SUBTIMER_STOP("sub_ghost_phys_face");
#if PRJ_MHD
    PRJ_SUBTIMER_START("sub_ghost_mhd_pre");
    prj_boundary_init_face_fidelity(mesh, mpi);
    for (i = 0; i < mesh->nblocks; ++i) {
        if (prj_boundary_active_block(mpi, &mesh->blocks[i])) {
            prj_boundary_physical_bf(mesh, bc, &mesh->blocks[i], use_bf1);
        }
    }
    PRJ_SUBTIMER_STOP("sub_ghost_mhd_pre");
#else
    (void)use_bf1;
    (void)eos;
#endif
    for (pass = 0; pass < 3; ++pass) {
        int fill_kind = fill_passes[pass];

        /* Post the inter-rank Isend/Irecv first, then do the intra-rank local
         * copies while the MPI messages are in flight. The pack inside the post
         * reads only active cells; the local copies write only ghost zones, so
         * the two never touch the same memory. The wait completes the exchange
         * after the local work is done. */
        PRJ_SUBTIMER_START("sub_ghost_post");
        prj_mpi_post_ghosts_and_bf(mesh, mpi, stage, fill_kind, use_bf1);
        PRJ_SUBTIMER_STOP("sub_ghost_post");
        PRJ_SUBTIMER_START("sub_ghost_send");
        for (i = 0; i < mesh->nblocks; ++i) {
            if (prj_boundary_active_block(mpi, &mesh->blocks[i])) {
                prj_boundary_send(mesh, mpi, &mesh->blocks[i], stage, fill_kind);
            }
        }
        PRJ_SUBTIMER_STOP("sub_ghost_send");
#if PRJ_MHD
        PRJ_SUBTIMER_START("sub_ghost_send_bf");
        for (i = 0; i < mesh->nblocks; ++i) {
            if (prj_boundary_active_block(mpi, &mesh->blocks[i])) {
                prj_boundary_send_bf(mesh, mpi, &mesh->blocks[i], use_bf1, fill_kind);
            }
        }
        PRJ_SUBTIMER_STOP("sub_ghost_send_bf");
#endif
#if PRJ_USE_GRAVITY
        /* Overlap the gravity radial reduce/integrate with the in-flight
         * same-level exchange. Both read only active cells (the reduce's
         * stage matches this fill's stage); the reduce writes radial profiles
         * and the integrate writes per-cell lapse/grav from geometry, so
         * neither touches the W ghost zones the exchange is filling. The two
         * allreduces inside the reduce progress the posted Isend/Irecv. */
        if (grav != 0 && fill_kind == PRJ_BOUNDARY_FILL_SAME_LEVEL) {
            PRJ_SUBTIMER_START("sub_ghost_grav");
            prj_gravity_monopole_reduce(mesh, grav, mpi, stage);
            prj_gravity_monopole_integrate(mesh, grav, mpi);
            PRJ_SUBTIMER_STOP("sub_ghost_grav");
        }
#endif
        /* Transport opacity for active cells overlaps the same-level exchange
         * too: it reads only active W/eosvar (stable after eos_fill_active) and
         * writes the active region of kappa_cell/sigma_cell. The 1-ghost halo
         * is filled by the caller after eos_fill_mesh. */
        if (rad != 0 && fill_kind == PRJ_BOUNDARY_FILL_SAME_LEVEL) {
            PRJ_SUBTIMER_START("sub_ghost_opac");
            prj_flux_fill_transport_opacity_active(mesh, rad, mpi, stage);
            PRJ_SUBTIMER_STOP("sub_ghost_opac");
        }
        PRJ_SUBTIMER_START("sub_ghost_wait");
        prj_mpi_wait_ghosts_and_bf(mesh, mpi, stage, fill_kind, use_bf1);
        PRJ_SUBTIMER_STOP("sub_ghost_wait");
    }
    PRJ_SUBTIMER_START("sub_ghost_phys_all");
    for (i = 0; i < mesh->nblocks; ++i) {
        if (prj_boundary_active_block(mpi, &mesh->blocks[i])) {
            prj_boundary_physical(mesh, bc, &mesh->blocks[i], stage, PRJ_BOUNDARY_PHYS_EDGE_CORNER);
        }
    }
    PRJ_SUBTIMER_STOP("sub_ghost_phys_all");
#if PRJ_MHD
    PRJ_SUBTIMER_START("sub_ghost_mhd_post");
    for (i = 0; i < mesh->nblocks; ++i) {
        if (prj_boundary_active_block(mpi, &mesh->blocks[i])) {
            prj_boundary_physical_bf(mesh, bc, &mesh->blocks[i], use_bf1);
        }
    }
    for (i = 0; i < mesh->nblocks; ++i) {
        if (prj_boundary_active_block(mpi, &mesh->blocks[i])) {
            prj_mhd_bf2bc_all(eos, &mesh->blocks[i], use_bf1);
        }
    }
    PRJ_SUBTIMER_STOP("sub_ghost_mhd_post");
#endif
}
