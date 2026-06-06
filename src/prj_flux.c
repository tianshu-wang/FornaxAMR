#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "prj.h"

static inline void prj_flux_face_cells(int dir, int i, int j, int k,
    int *il, int *jl, int *kl, int *ir, int *jr, int *kr)
{
    *il = i;
    *jl = j;
    *kl = k;
    *ir = i;
    *jr = j;
    *kr = k;
    if (dir == X1DIR) {
        *il = i - 1;
        *ir = i;
    } else if (dir == X2DIR) {
        *jl = j - 1;
        *jr = j;
    } else {
        *kl = k - 1;
        *kr = k;
    }
}

static inline int prj_flux_prim_var_radiation(int v);

/* Reconstruct a single primitive face value.  The variable-type decisions
 * (radiation? which reconstruction scheme?) are loop-invariant for a fixed
 * variable, so callers precompute them once per variable via
 * prj_flux_prim_var_kind() and pass them in here, keeping the per-face inner loop
 * branch-free on variable type.  The reconstruction-scheme value is returned
 * as-is (no positivity guard). */
static inline double prj_flux_prim_face_value(const double *W, int v, int dir,
    int i, int j, int k, double target, int is_rad, int recon_scheme)
{
    double q[PRJ_RECON_NCELLS];
    int n;

#if !(PRJ_NRAD > 0 && PRJ_MIXED_PRECISION)
    (void)is_rad;
#endif

    for (n = 0; n < PRJ_RECON_NCELLS; ++n) {
        int offset = n - (PRJ_RECON_NCELLS / 2);

        if (dir == X1DIR) {
            q[n] = W[VIDX(v, i + offset, j, k)];
        } else if (dir == X2DIR) {
            q[n] = W[VIDX(v, i, j + offset, k)];
        } else {
            q[n] = W[VIDX(v, i, j, k + offset)];
        }
#if PRJ_NRAD > 0 && PRJ_MIXED_PRECISION
        if (is_rad) {
            q[n] = prj_rad_mixed_round(q[n]);
        }
#endif
    }

    return prj_reconstruct_face_value_for_scheme(q, target, recon_scheme);
}

/* Resolve the (radiation, scheme) reconstruction kind for a variable.
 * Mirrors the dispatch that prj_flux_prim_face_value() previously did per face. */
static inline void prj_flux_prim_var_kind(int v, int *is_rad, int *recon_scheme)
{
    *is_rad = prj_flux_prim_var_radiation(v);
    *recon_scheme = *is_rad ? PRJ_RECON_RADIATION : PRJ_RECON_HYDRO;
}

static inline int prj_flux_prim_var_radiation(int v)
{
#if PRJ_NRAD > 0
    return v >= PRJ_NHYDRO && v < PRJ_NVAR_PRIM;
#else
    (void)v;
    return 0;
#endif
}

static inline double prj_flux_eos_face_value(const double *eosvar, int v, int dir,
    int i, int j, int k, double target)
{
    double q[PRJ_RECON_NCELLS];
    const int half = PRJ_RECON_NCELLS / 2;
    int n;

    for (n = 0; n < PRJ_RECON_NCELLS; ++n) {
        int offset = n - half;

        if (dir == X1DIR) {
            q[n] = eosvar[EIDX(v, i + offset, j, k)];
        } else if (dir == X2DIR) {
            q[n] = eosvar[EIDX(v, i, j + offset, k)];
        } else {
            q[n] = eosvar[EIDX(v, i, j, k + offset)];
        }
    }

    return prj_reconstruct_face_value_for_scheme(q, target, PRJ_RECON_HYDRO);
}

/* Build the local->global primitive index permutation for direction `dir`.
 * The sweep-normal velocity/B/flux component is rotated into the local X1 slot
 * so the Riemann solver always sees the sweep direction as its normal. The map
 * is a permutation over all PRJ_NVAR_PRIM indices. */
static void prj_flux_build_var_map(int dir, int *gmap)
{
    int field = 0;
    int group = 0;

    gmap[PRJ_PRIM_RHO] = PRJ_PRIM_RHO;
    gmap[PRJ_PRIM_EINT] = PRJ_PRIM_EINT;
    gmap[PRJ_PRIM_YE] = PRJ_PRIM_YE;

    if (dir == X1DIR) {
        gmap[PRJ_PRIM_V1] = PRJ_PRIM_V1;
        gmap[PRJ_PRIM_V2] = PRJ_PRIM_V2;
        gmap[PRJ_PRIM_V3] = PRJ_PRIM_V3;
#if PRJ_MHD
        gmap[PRJ_PRIM_B1] = PRJ_PRIM_B1;
        gmap[PRJ_PRIM_B2] = PRJ_PRIM_B2;
        gmap[PRJ_PRIM_B3] = PRJ_PRIM_B3;
#endif
    } else if (dir == X2DIR) {
        gmap[PRJ_PRIM_V1] = PRJ_PRIM_V2;
        gmap[PRJ_PRIM_V2] = PRJ_PRIM_V3;
        gmap[PRJ_PRIM_V3] = PRJ_PRIM_V1;
#if PRJ_MHD
        gmap[PRJ_PRIM_B1] = PRJ_PRIM_B2;
        gmap[PRJ_PRIM_B2] = PRJ_PRIM_B3;
        gmap[PRJ_PRIM_B3] = PRJ_PRIM_B1;
#endif
    } else {
        gmap[PRJ_PRIM_V1] = PRJ_PRIM_V3;
        gmap[PRJ_PRIM_V2] = PRJ_PRIM_V1;
        gmap[PRJ_PRIM_V3] = PRJ_PRIM_V2;
#if PRJ_MHD
        gmap[PRJ_PRIM_B1] = PRJ_PRIM_B3;
        gmap[PRJ_PRIM_B2] = PRJ_PRIM_B1;
        gmap[PRJ_PRIM_B3] = PRJ_PRIM_B2;
#endif
    }

#if PRJ_NRAD > 0
    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            gmap[PRJ_PRIM_RAD_E(field, group)] = PRJ_PRIM_RAD_E(field, group);
            if (dir == X1DIR) {
                gmap[PRJ_PRIM_RAD_F1(field, group)] = PRJ_PRIM_RAD_F1(field, group);
                gmap[PRJ_PRIM_RAD_F2(field, group)] = PRJ_PRIM_RAD_F2(field, group);
                gmap[PRJ_PRIM_RAD_F3(field, group)] = PRJ_PRIM_RAD_F3(field, group);
            } else if (dir == X2DIR) {
                gmap[PRJ_PRIM_RAD_F1(field, group)] = PRJ_PRIM_RAD_F2(field, group);
                gmap[PRJ_PRIM_RAD_F2(field, group)] = PRJ_PRIM_RAD_F3(field, group);
                gmap[PRJ_PRIM_RAD_F3(field, group)] = PRJ_PRIM_RAD_F1(field, group);
            } else {
                gmap[PRJ_PRIM_RAD_F1(field, group)] = PRJ_PRIM_RAD_F3(field, group);
                gmap[PRJ_PRIM_RAD_F2(field, group)] = PRJ_PRIM_RAD_F1(field, group);
                gmap[PRJ_PRIM_RAD_F3(field, group)] = PRJ_PRIM_RAD_F2(field, group);
            }
        }
    }
#else
    (void)field;
    (void)group;
#endif
}

/* Reconstruct every primitive (and the pressure/gamma EOS vars) onto the dir
 * faces of a block, with the *variable* loop outermost.  Results are written to
 * block-wide buffers: WL_block/WR_block use layout [local_var * PRJ_BLOCK_NCELLS
 * + IDX(i, j, k)]; the pL/pR/gL/gR buffers are indexed by IDX(i, j, k).  The
 * Riemann sweep in prj_flux_update consumes these per face. */
static void prj_flux_reconstruct_block(double *W, const double *eosvar, int dir,
    int istart, int iend, int jstart, int jend, int kstart, int kend,
    double *WL_block, double *WR_block,
    double *pL_block, double *pR_block, double *gL_block, double *gR_block)
{
    int gmap[PRJ_NVAR_PRIM];
    int lv;
    int i;
    int j;
    int k;

    prj_flux_build_var_map(dir, gmap);

    for (lv = 0; lv < PRJ_NVAR_PRIM; ++lv) {
        int gv = gmap[lv];
        double *wl = &WL_block[(size_t)lv * PRJ_BLOCK_NCELLS];
        double *wr = &WR_block[(size_t)lv * PRJ_BLOCK_NCELLS];
        int is_rad;
        int recon_scheme;

        /* Variable-type dispatch hoisted out of the per-face loop below. */
        prj_flux_prim_var_kind(gv, &is_rad, &recon_scheme);

        for (i = istart; i <= iend; ++i) {
            for (j = jstart; j <= jend; ++j) {
                for (k = kstart; k <= kend; ++k) {
                    int il;
                    int jl;
                    int kl;
                    int ir;
                    int jr;
                    int kr;

                    prj_flux_face_cells(dir, i, j, k, &il, &jl, &kl, &ir, &jr, &kr);
                    wl[IDX(i, j, k)] = prj_flux_prim_face_value(W, gv, dir, il, jl, kl, 0.5,
                        is_rad, recon_scheme);
                    wr[IDX(i, j, k)] = prj_flux_prim_face_value(W, gv, dir, ir, jr, kr, -0.5,
                        is_rad, recon_scheme);
                }
            }
        }
    }

    if (eosvar == 0) {
        return;
    }

    {
        const int eos_var[2] = { PRJ_EOSVAR_PRESSURE, PRJ_EOSVAR_GAMMA };
        double *eosL[2];
        double *eosR[2];
        int e;

        eosL[0] = pL_block;
        eosL[1] = gL_block;
        eosR[0] = pR_block;
        eosR[1] = gR_block;

        for (e = 0; e < 2; ++e) {
            int ev = eos_var[e];
            double *el = eosL[e];
            double *er = eosR[e];

            for (i = istart; i <= iend; ++i) {
                for (j = jstart; j <= jend; ++j) {
                    for (k = kstart; k <= kend; ++k) {
                        int il;
                        int jl;
                        int kl;
                        int ir;
                        int jr;
                        int kr;

                        prj_flux_face_cells(dir, i, j, k, &il, &jl, &kl, &ir, &jr, &kr);
                        el[IDX(i, j, k)] = prj_flux_eos_face_value(eosvar, ev, dir, il, jl, kl, 0.5);
                        er[IDX(i, j, k)] = prj_flux_eos_face_value(eosvar, ev, dir, ir, jr, kr, -0.5);
                    }
                }
            }
        }
    }
}

static void prj_flux_store_face_velocity(prj_block *block, int dir, int i, int j, int k,
    const double v_face_loc[3])
{
    double *dst;

    if (block->v_riemann[dir] == 0) {
        return;
    }
    dst = block->v_riemann[dir];
    if (dir == X1DIR) {
        dst[0 * PRJ_BLOCK_NCELLS + IDX(i, j, k)] = v_face_loc[0];
        dst[1 * PRJ_BLOCK_NCELLS + IDX(i, j, k)] = v_face_loc[1];
        dst[2 * PRJ_BLOCK_NCELLS + IDX(i, j, k)] = v_face_loc[2];
    } else if (dir == X2DIR) {
        dst[1 * PRJ_BLOCK_NCELLS + IDX(i, j, k)] = v_face_loc[0];
        dst[2 * PRJ_BLOCK_NCELLS + IDX(i, j, k)] = v_face_loc[1];
        dst[0 * PRJ_BLOCK_NCELLS + IDX(i, j, k)] = v_face_loc[2];
    } else {
        dst[2 * PRJ_BLOCK_NCELLS + IDX(i, j, k)] = v_face_loc[0];
        dst[0 * PRJ_BLOCK_NCELLS + IDX(i, j, k)] = v_face_loc[1];
        dst[1 * PRJ_BLOCK_NCELLS + IDX(i, j, k)] = v_face_loc[2];
    }
}

static void prj_flux_store_local_flux(double *dst, int dir, int i, int j, int k,
    const double *Fl)
{
    int field=0;
    int group=0;

    dst[VIDX(PRJ_CONS_RHO, i, j, k)] = Fl[PRJ_CONS_RHO];
    dst[VIDX(PRJ_CONS_ETOT, i, j, k)] = Fl[PRJ_CONS_ETOT];
    dst[VIDX(PRJ_CONS_YE, i, j, k)] = Fl[PRJ_CONS_YE];

    if (dir == X1DIR) {
        dst[VIDX(PRJ_CONS_MOM1, i, j, k)] = Fl[PRJ_CONS_MOM1];
        dst[VIDX(PRJ_CONS_MOM2, i, j, k)] = Fl[PRJ_CONS_MOM2];
        dst[VIDX(PRJ_CONS_MOM3, i, j, k)] = Fl[PRJ_CONS_MOM3];
#if PRJ_MHD
        dst[VIDX(PRJ_CONS_B1, i, j, k)] = Fl[PRJ_CONS_B1];
        dst[VIDX(PRJ_CONS_B2, i, j, k)] = Fl[PRJ_CONS_B2];
        dst[VIDX(PRJ_CONS_B3, i, j, k)] = Fl[PRJ_CONS_B3];
#endif
    } else if (dir == X2DIR) {
        dst[VIDX(PRJ_CONS_MOM2, i, j, k)] = Fl[PRJ_CONS_MOM1];
        dst[VIDX(PRJ_CONS_MOM3, i, j, k)] = Fl[PRJ_CONS_MOM2];
        dst[VIDX(PRJ_CONS_MOM1, i, j, k)] = Fl[PRJ_CONS_MOM3];
#if PRJ_MHD
        dst[VIDX(PRJ_CONS_B2, i, j, k)] = Fl[PRJ_CONS_B1];
        dst[VIDX(PRJ_CONS_B3, i, j, k)] = Fl[PRJ_CONS_B2];
        dst[VIDX(PRJ_CONS_B1, i, j, k)] = Fl[PRJ_CONS_B3];
#endif
    } else {
        dst[VIDX(PRJ_CONS_MOM3, i, j, k)] = Fl[PRJ_CONS_MOM1];
        dst[VIDX(PRJ_CONS_MOM1, i, j, k)] = Fl[PRJ_CONS_MOM2];
        dst[VIDX(PRJ_CONS_MOM2, i, j, k)] = Fl[PRJ_CONS_MOM3];
#if PRJ_MHD
        dst[VIDX(PRJ_CONS_B3, i, j, k)] = Fl[PRJ_CONS_B1];
        dst[VIDX(PRJ_CONS_B1, i, j, k)] = Fl[PRJ_CONS_B2];
        dst[VIDX(PRJ_CONS_B2, i, j, k)] = Fl[PRJ_CONS_B3];
#endif
    }

#if PRJ_NRAD > 0
    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            dst[VIDX(PRJ_CONS_RAD_E(field, group), i, j, k)] = Fl[PRJ_CONS_RAD_E(field, group)];
            if (dir == X1DIR) {
                dst[VIDX(PRJ_CONS_RAD_F1(field, group), i, j, k)] = Fl[PRJ_CONS_RAD_F1(field, group)];
                dst[VIDX(PRJ_CONS_RAD_F2(field, group), i, j, k)] = Fl[PRJ_CONS_RAD_F2(field, group)];
                dst[VIDX(PRJ_CONS_RAD_F3(field, group), i, j, k)] = Fl[PRJ_CONS_RAD_F3(field, group)];
            } else if (dir == X2DIR) {
                dst[VIDX(PRJ_CONS_RAD_F2(field, group), i, j, k)] = Fl[PRJ_CONS_RAD_F1(field, group)];
                dst[VIDX(PRJ_CONS_RAD_F3(field, group), i, j, k)] = Fl[PRJ_CONS_RAD_F2(field, group)];
                dst[VIDX(PRJ_CONS_RAD_F1(field, group), i, j, k)] = Fl[PRJ_CONS_RAD_F3(field, group)];
            } else {
                dst[VIDX(PRJ_CONS_RAD_F3(field, group), i, j, k)] = Fl[PRJ_CONS_RAD_F1(field, group)];
                dst[VIDX(PRJ_CONS_RAD_F1(field, group), i, j, k)] = Fl[PRJ_CONS_RAD_F2(field, group)];
                dst[VIDX(PRJ_CONS_RAD_F2(field, group), i, j, k)] = Fl[PRJ_CONS_RAD_F3(field, group)];
            }
        }
    }
#else
    (void)field;
    (void)group;
#endif
}

static void prj_flux_velocity_deltas(double *W, int dir,
    int il, int jl, int kl, int ir, int jr, int kr,
    double *deltau, double *deltav, double *deltaw)
{
    if (dir == X1DIR) {
        *deltau = W[VIDX(PRJ_PRIM_V1, ir, jr, kr)] -
            W[VIDX(PRJ_PRIM_V1, il, jl, kl)];
        *deltav = PRJ_MIN(
            PRJ_MIN(W[VIDX(PRJ_PRIM_V2, il, jl, kl)] -
                    W[VIDX(PRJ_PRIM_V2, il, jl - 1, kl)],
                    W[VIDX(PRJ_PRIM_V2, il, jl + 1, kl)] -
                    W[VIDX(PRJ_PRIM_V2, il, jl, kl)]),
            PRJ_MIN(W[VIDX(PRJ_PRIM_V2, ir, jr, kr)] -
                    W[VIDX(PRJ_PRIM_V2, ir, jr - 1, kr)],
                    W[VIDX(PRJ_PRIM_V2, ir, jr + 1, kr)] -
                    W[VIDX(PRJ_PRIM_V2, ir, jr, kr)]));
        *deltaw = PRJ_MIN(
            PRJ_MIN(W[VIDX(PRJ_PRIM_V3, il, jl, kl)] -
                    W[VIDX(PRJ_PRIM_V3, il, jl, kl - 1)],
                    W[VIDX(PRJ_PRIM_V3, il, jl, kl + 1)] -
                    W[VIDX(PRJ_PRIM_V3, il, jl, kl)]),
            PRJ_MIN(W[VIDX(PRJ_PRIM_V3, ir, jr, kr)] -
                    W[VIDX(PRJ_PRIM_V3, ir, jr, kr - 1)],
                    W[VIDX(PRJ_PRIM_V3, ir, jr, kr + 1)] -
                    W[VIDX(PRJ_PRIM_V3, ir, jr, kr)]));
    } else if (dir == X2DIR) {
        *deltau = W[VIDX(PRJ_PRIM_V2, ir, jr, kr)] -
            W[VIDX(PRJ_PRIM_V2, il, jl, kl)];
        *deltav = PRJ_MIN(
            PRJ_MIN(W[VIDX(PRJ_PRIM_V3, il, jl, kl)] -
                    W[VIDX(PRJ_PRIM_V3, il, jl, kl - 1)],
                    W[VIDX(PRJ_PRIM_V3, il, jl, kl + 1)] -
                    W[VIDX(PRJ_PRIM_V3, il, jl, kl)]),
            PRJ_MIN(W[VIDX(PRJ_PRIM_V3, ir, jr, kr)] -
                    W[VIDX(PRJ_PRIM_V3, ir, jr, kr - 1)],
                    W[VIDX(PRJ_PRIM_V3, ir, jr, kr + 1)] -
                    W[VIDX(PRJ_PRIM_V3, ir, jr, kr)]));
        *deltaw = PRJ_MIN(
            PRJ_MIN(W[VIDX(PRJ_PRIM_V1, il, jl, kl)] -
                    W[VIDX(PRJ_PRIM_V1, il - 1, jl, kl)],
                    W[VIDX(PRJ_PRIM_V1, il + 1, jl, kl)] -
                    W[VIDX(PRJ_PRIM_V1, il, jl, kl)]),
            PRJ_MIN(W[VIDX(PRJ_PRIM_V1, ir, jr, kr)] -
                    W[VIDX(PRJ_PRIM_V1, ir - 1, jr, kr)],
                    W[VIDX(PRJ_PRIM_V1, ir + 1, jr, kr)] -
                    W[VIDX(PRJ_PRIM_V1, ir, jr, kr)]));
    } else {
        *deltau = W[VIDX(PRJ_PRIM_V3, ir, jr, kr)] -
            W[VIDX(PRJ_PRIM_V3, il, jl, kl)];
        *deltav = PRJ_MIN(
            PRJ_MIN(W[VIDX(PRJ_PRIM_V1, il, jl, kl)] -
                    W[VIDX(PRJ_PRIM_V1, il - 1, jl, kl)],
                    W[VIDX(PRJ_PRIM_V1, il + 1, jl, kl)] -
                    W[VIDX(PRJ_PRIM_V1, il, jl, kl)]),
            PRJ_MIN(W[VIDX(PRJ_PRIM_V1, ir, jr, kr)] -
                    W[VIDX(PRJ_PRIM_V1, ir - 1, jr, kr)],
                    W[VIDX(PRJ_PRIM_V1, ir + 1, jr, kr)] -
                    W[VIDX(PRJ_PRIM_V1, ir, jr, kr)]));
        *deltaw = PRJ_MIN(
            PRJ_MIN(W[VIDX(PRJ_PRIM_V2, il, jl, kl)] -
                    W[VIDX(PRJ_PRIM_V2, il, jl - 1, kl)],
                    W[VIDX(PRJ_PRIM_V2, il, jl + 1, kl)] -
                    W[VIDX(PRJ_PRIM_V2, il, jl, kl)]),
            PRJ_MIN(W[VIDX(PRJ_PRIM_V2, ir, jr, kr)] -
                    W[VIDX(PRJ_PRIM_V2, ir, jr - 1, kr)],
                    W[VIDX(PRJ_PRIM_V2, ir, jr + 1, kr)] -
                    W[VIDX(PRJ_PRIM_V2, ir, jr, kr)]));
    }
}

#if PRJ_NRAD > 0
/* Transport absorption/scattering opacity for one cell, stored per group in
 * block->kappa_cell / block->sigma_cell. Used by the radiation flux. */
static void prj_flux_opacity_cell(prj_rad *rad, prj_block *block, const double *W,
    const double *eosvar, int ii, int jj, int kk)
{
    const size_t stride = (size_t)PRJ_NRAD * (size_t)PRJ_NEGROUP;
    double rho = W[VIDX(PRJ_PRIM_RHO, ii, jj, kk)];
    double ye = W[VIDX(PRJ_PRIM_YE, ii, jj, kk)];
    double T = (eosvar != 0) ? eosvar[EIDX(PRJ_EOSVAR_TEMPERATURE, ii, jj, kk)] : 0.0;
    size_t off = (size_t)IDX(ii, jj, kk) * stride;

    prj_rad3_opac_lookup(rad, rho, T, ye,
        &block->kappa_cell[off], &block->sigma_cell[off], 0, 0);
}

static int prj_flux_opacity_block_ready(const prj_block *block)
{
    return block != 0 && block->kappa_cell != 0 && block->sigma_cell != 0;
}
#endif

/* Fill transport opacity over the active cells only (0..PRJ_BLOCK_SIZE-1).
 * Reads active W/eosvar, which are stable once eos_fill_active has run, so this
 * can be overlapped with the same-level ghost exchange. */
void prj_flux_fill_transport_opacity_active(prj_mesh *mesh, prj_rad *rad,
    const prj_mpi *mpi, int stage)
{
#if PRJ_NRAD > 0
    int bidx;

    if (mesh == 0 || rad == 0) {
        return;
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        double *W = stage == 2 ? block->W1 : block->W;
        int ii;
        int jj;
        int kk;

        if (block->id < 0 || block->active != 1 || W == 0 ||
            !prj_flux_opacity_block_ready(block)) {
            continue;
        }
        if (mpi != 0 && block->rank != mpi->rank) {
            continue;
        }
        for (ii = 0; ii < PRJ_BLOCK_SIZE; ++ii) {
            for (jj = 0; jj < PRJ_BLOCK_SIZE; ++jj) {
                for (kk = 0; kk < PRJ_BLOCK_SIZE; ++kk) {
                    prj_flux_opacity_cell(rad, block, W, block->eosvar, ii, jj, kk);
                }
            }
        }
    }
#else
    (void)mesh;
    (void)rad;
    (void)mpi;
    (void)stage;
#endif
}

/* Fill transport opacity over the 1-cell ghost halo (the shell of the
 * [-1, PRJ_BLOCK_SIZE] box). Needs ghost W/eosvar, so it must run after all
 * ghost filling and eos_fill_mesh are done. */
void prj_flux_fill_transport_opacity_halo(prj_mesh *mesh, prj_rad *rad,
    const prj_mpi *mpi, int stage)
{
#if PRJ_NRAD > 0
    int bidx;

    if (mesh == 0 || rad == 0) {
        return;
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        double *W = stage == 2 ? block->W1 : block->W;
        int ii;
        int jj;
        int kk;

        if (block->id < 0 || block->active != 1 || W == 0 ||
            !prj_flux_opacity_block_ready(block)) {
            continue;
        }
        if (mpi != 0 && block->rank != mpi->rank) {
            continue;
        }
        for (ii = -1; ii <= PRJ_BLOCK_SIZE; ++ii) {
            for (jj = -1; jj <= PRJ_BLOCK_SIZE; ++jj) {
                for (kk = -1; kk <= PRJ_BLOCK_SIZE; ++kk) {
                    if (ii >= 0 && ii < PRJ_BLOCK_SIZE &&
                        jj >= 0 && jj < PRJ_BLOCK_SIZE &&
                        kk >= 0 && kk < PRJ_BLOCK_SIZE) {
                        continue;
                    }
                    prj_flux_opacity_cell(rad, block, W, block->eosvar, ii, jj, kk);
                }
            }
        }
    }
#else
    (void)mesh;
    (void)rad;
    (void)mpi;
    (void)stage;
#endif
}

void prj_flux_update(prj_eos *eos, prj_rad *rad, prj_block *block, double *W,
    double *eosvar, double *flux[3], int use_bf1)
{
    /* Block-wide face-state scratch for the variable-outermost reconstruction
     * pass.  Reused across calls (this code path is serial); WL/WR use layout
     * [local_var * PRJ_BLOCK_NCELLS + IDX(i, j, k)]. */
    static double WL_block[PRJ_NVAR_PRIM * PRJ_BLOCK_NCELLS];
    static double WR_block[PRJ_NVAR_PRIM * PRJ_BLOCK_NCELLS];
    static double pL_block[PRJ_BLOCK_NCELLS];
    static double pR_block[PRJ_BLOCK_NCELLS];
    static double gL_block[PRJ_BLOCK_NCELLS];
    static double gR_block[PRJ_BLOCK_NCELLS];
    int dir;
#if PRJ_NRAD == 0
    (void)rad;
#endif
#if !PRJ_MHD
    (void)use_bf1;
#endif

    /* Transport opacity (kappa_cell/sigma_cell) is now filled ahead of the flux:
     * active cells in the same-level ghost shadow and the 1-ghost halo after
     * eos_fill_mesh. See prj_flux_fill_transport_opacity_active/halo. */

    for (dir = 0; dir < 3; ++dir) {
        int i;
        int j;
        int k;
        int istart = 0;
        int iend = PRJ_BLOCK_SIZE;
        int jstart = 0;
        int jend = PRJ_BLOCK_SIZE;
        int kstart = 0;
        int kend = PRJ_BLOCK_SIZE;

#if PRJ_MHD
        /* MHD needs one extra transverse face layer (CT/EMF reconstruction). */
        if (dir != X1DIR) {
            istart = -1;
        }
        if (dir != X2DIR) {
            jstart = -1;
        }
        if (dir != X3DIR) {
            kstart = -1;
        }
#else
        /* Transverse directions carry one face per cell ([0, PRJ_BLOCK_SIZE - 1]);
         * only the sweep-normal direction has PRJ_BLOCK_SIZE + 1 faces. */
        if (dir != X1DIR) {
            iend = PRJ_BLOCK_SIZE - 1;
        }
        if (dir != X2DIR) {
            jend = PRJ_BLOCK_SIZE - 1;
        }
        if (dir != X3DIR) {
            kend = PRJ_BLOCK_SIZE - 1;
        }
#endif

        /* Variable-outermost reconstruction pass over the block. */
        prj_flux_reconstruct_block(W, eosvar, dir,
            istart, iend, jstart, jend, kstart, kend,
            WL_block, WR_block, pL_block, pR_block, gL_block, gR_block);

        for (i = istart; i <= iend; ++i) {
            for (j = jstart; j <= jend; ++j) {
                for (k = kstart; k <= kend; ++k) {
                    double WL[PRJ_NVAR_PRIM];
                    double WR[PRJ_NVAR_PRIM];
                    double pL;
                    double pR;
                    double gL;
                    double gR;
                    double Fl[PRJ_NVAR_CONS];
                    int il;
                    int jl;
                    int kl;
                    int ir;
                    int jr;
                    int kr;
                    int v;
                    size_t fidx;

                    double v_face_loc[3] = {0.0, 0.0, 0.0};
                    double deltau;
                    double deltav;
                    double deltaw;

                    /* Gather the reconstructed left/right states for this face. */
                    fidx = (size_t)IDX(i, j, k);
                    for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
                        WL[v] = WL_block[(size_t)v * PRJ_BLOCK_NCELLS + fidx];
                        WR[v] = WR_block[(size_t)v * PRJ_BLOCK_NCELLS + fidx];
                    }
                    if (eosvar != 0) {
                        pL = pL_block[fidx];
                        pR = pR_block[fidx];
                        gL = gL_block[fidx];
                        gR = gR_block[fidx];
                    } else {
                        pL = 0.0;
                        pR = 0.0;
                        gL = 0.0;
                        gR = 0.0;
                    }
                    prj_flux_face_cells(dir, i, j, k, &il, &jl, &kl, &ir, &jr, &kr);
                    prj_flux_velocity_deltas(W, dir, il, jl, kl, ir, jr, kr,
                        &deltau, &deltav, &deltaw);
#if PRJ_MHD
                    {
                        double *bf_dir = use_bf1 != 0 ? block->Bf1[dir] : block->Bf[dir];
                        double bv1 = 0.0;
                        double bv2 = 0.0;
                        double bn;

                        if (bf_dir == 0 || block->Bv1[dir] == 0 || block->Bv2[dir] == 0) {
                            fprintf(stderr, "prj_flux_update: missing MHD face storage for dir=%d\n", dir);
                            exit(1);
                        }
                        bn = bf_dir[FACE_IDX(dir, i, j, k)];
                        WL[PRJ_PRIM_B1] = bn;
                        WR[PRJ_PRIM_B1] = bn;

#if PRJ_LHLLD_RIEMANN
                        /* PRJ_LHLLD_RIEMANN=1 is the default Minoshima &
                         * Miyoshi LHLLD path.  Set LHLLD_RIEMANN=0 at build
                         * time for the legacy current-HLLD comparison path. */
                        prj_riemann_lhlld(WL, WR, pL, pR, gL, gR, eos, bn, Fl,
                            v_face_loc, &bv1, &bv2, deltau, deltav, deltaw);
#else
                        /* Legacy current-HLLD comparison path. */
                        prj_riemann_hlld(WL, WR, pL, pR, gL, gR, eos, bn, Fl,
                            v_face_loc, &bv1, &bv2, deltau, deltav, deltaw);
#endif
                        block->Bv1[dir][IDX(i, j, k)] = bv1;
                        block->Bv2[dir][IDX(i, j, k)] = bv2;
                    }
#else
                    prj_riemann_hllc(WL, WR, pL, pR, gL, gR, eos, Fl,
                        v_face_loc, deltau, deltav, deltaw);
#endif
                    prj_flux_store_face_velocity(block, dir, i, j, k, v_face_loc);
#if PRJ_NRAD > 0
                    {
                        double dx_dir;
                        double lapse_face = 1.0;
                        double chi_face[PRJ_NRAD * PRJ_NEGROUP];
                        const size_t stride = (size_t)PRJ_NRAD * (size_t)PRJ_NEGROUP;
                        size_t off_L = (size_t)IDX(il, jl, kl) * stride;
                        size_t off_R = (size_t)IDX(ir, jr, kr) * stride;
                        const double *kappa_L = &block->kappa_cell[off_L];
                        const double *sigma_L = &block->sigma_cell[off_L];
                        const double *kappa_R = &block->kappa_cell[off_R];
                        const double *sigma_R = &block->sigma_cell[off_R];
                        int idx;

                        if (block->lapse != 0) {
                            lapse_face = 0.5 *
                                (block->lapse[IDX(il, jl, kl)] + block->lapse[IDX(ir, jr, kr)]);
                        }
                        dx_dir = block->dx[dir];

                        for (idx = 0; idx < PRJ_NRAD * PRJ_NEGROUP; ++idx) {
                            double kL = kappa_L[idx];
                            double kR = kappa_R[idx];
                            double sL_o = sigma_L[idx];
                            double sR_o = sigma_R[idx];
                            double k_sum = kL + kR;
                            double s_sum = sL_o + sR_o;
                            double k_face = 2.0 * kL * kR / k_sum;
                            double s_face = 2.0 * sL_o * sR_o / s_sum;
                            chi_face[idx] = k_face + s_face;
                        }

                        prj_rad_flux(rad, WL, WR, lapse_face, chi_face, dx_dir, v_face_loc[0], Fl);
                    }
#endif
                    prj_flux_store_local_flux(flux[dir], dir, i, j, k, Fl);
                }
            }
        }
    }
}

void prj_flux_div(double *flux[3], double area[3], double vol, int i, int j, int k, double *fluxdiv)
{
    int v;

    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        fluxdiv[v] =
            -(area[X1DIR] * (flux[X1DIR][VIDX(v, i + 1, j, k)] - flux[X1DIR][VIDX(v, i, j, k)]) +
                area[X2DIR] * (flux[X2DIR][VIDX(v, i, j + 1, k)] - flux[X2DIR][VIDX(v, i, j, k)]) +
                area[X3DIR] * (flux[X3DIR][VIDX(v, i, j, k + 1)] - flux[X3DIR][VIDX(v, i, j, k)])) / vol;
    }
}
