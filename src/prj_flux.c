#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "prj.h"

#if defined(PRJ_ENABLE_MPI)
#include <mpi.h>
#endif

/* The per-cell MC slope fast path (prj_flux_recon_var_mc) replaces the per-face
 * reconstruction when both physics use the MC scheme. The per-face face-value
 * helpers and macros below are only needed for the WENO fallback. */
#if PRJ_RECON_HYDRO == PRJ_RECON_MC && PRJ_RECON_RADIATION == PRJ_RECON_MC
#define PRJ_FLUX_RECON_MC_FASTPATH 1
#else
#define PRJ_FLUX_RECON_MC_FASTPATH 0
#endif

static inline void prj_flux_face_cells_x1(int i, int j, int k,
    int *il, int *jl, int *kl, int *ir, int *jr, int *kr)
{
    *il = i - 1;
    *jl = j;
    *kl = k;
    *ir = i;
    *jr = j;
    *kr = k;
}

static inline void prj_flux_face_cells_x2(int i, int j, int k,
    int *il, int *jl, int *kl, int *ir, int *jr, int *kr)
{
    *il = i;
    *jl = j - 1;
    *kl = k;
    *ir = i;
    *jr = j;
    *kr = k;
}

static inline void prj_flux_face_cells_x3(int i, int j, int k,
    int *il, int *jl, int *kl, int *ir, int *jr, int *kr)
{
    *il = i;
    *jl = j;
    *kl = k - 1;
    *ir = i;
    *jr = j;
    *kr = k;
}

static inline void prj_flux_face_cells(int dir, int i, int j, int k,
    int *il, int *jl, int *kl, int *ir, int *jr, int *kr)
{
    if (dir == X1DIR) {
        prj_flux_face_cells_x1(i, j, k, il, jl, kl, ir, jr, kr);
    } else if (dir == X2DIR) {
        prj_flux_face_cells_x2(i, j, k, il, jl, kl, ir, jr, kr);
    } else {
        prj_flux_face_cells_x3(i, j, k, il, jl, kl, ir, jr, kr);
    }
}

#if PRJ_NRAD > 0
/* The radiation flux at a face is consumed by the flux divergence only when the
 * face's two transverse indices lie in the active range [0, PRJ_BLOCK_SIZE - 1].
 * For MHD the face loop carries an extra transverse ghost layer (CT/EMF
 * reconstruction); radiation never reads those faces, so the M1 flux there is
 * pure waste. Non-MHD loop bounds already exclude them, making this a no-op. */
static inline int prj_flux_rad_face_needed(int dir, int i, int j, int k)
{
    if (dir == X1DIR) {
        return j >= 0 && j < PRJ_BLOCK_SIZE && k >= 0 && k < PRJ_BLOCK_SIZE;
    }
    if (dir == X2DIR) {
        return i >= 0 && i < PRJ_BLOCK_SIZE && k >= 0 && k < PRJ_BLOCK_SIZE;
    }
    return i >= 0 && i < PRJ_BLOCK_SIZE && j >= 0 && j < PRJ_BLOCK_SIZE;
}
#endif

#if !PRJ_FLUX_RECON_MC_FASTPATH
static inline double prj_flux_radiation_stencil_value(double q)
{
#if PRJ_NRAD > 0 && PRJ_MIXED_PRECISION_FLUX
    return prj_rad_mixed_round(q);
#else
    return q;
#endif
}

static inline double prj_flux_hydro_prim_face_value_x1(const double *W, int v,
    int i, int j, int k, double target)
{
    double q[PRJ_RECON_HYDRO_NCELLS];
    const int half = PRJ_RECON_HYDRO_NCELLS / 2;
    int n;

    for (n = 0; n < PRJ_RECON_HYDRO_NCELLS; ++n) {
        int offset = n - half;
        q[n] = W[WIDX(v, i + offset, j, k)];
    }

    return prj_reconstruct_hydro_face_value(q, target);
}

static inline double prj_flux_hydro_prim_face_value_x2(const double *W, int v,
    int i, int j, int k, double target)
{
    double q[PRJ_RECON_HYDRO_NCELLS];
    const int half = PRJ_RECON_HYDRO_NCELLS / 2;
    int n;

    for (n = 0; n < PRJ_RECON_HYDRO_NCELLS; ++n) {
        int offset = n - half;
        q[n] = W[WIDX(v, i, j + offset, k)];
    }

    return prj_reconstruct_hydro_face_value(q, target);
}

static inline double prj_flux_hydro_prim_face_value_x3(const double *W, int v,
    int i, int j, int k, double target)
{
    double q[PRJ_RECON_HYDRO_NCELLS];
    const int half = PRJ_RECON_HYDRO_NCELLS / 2;
    int n;

    for (n = 0; n < PRJ_RECON_HYDRO_NCELLS; ++n) {
        int offset = n - half;
        q[n] = W[WIDX(v, i, j, k + offset)];
    }

    return prj_reconstruct_hydro_face_value(q, target);
}

static inline double prj_flux_radiation_prim_face_value_x1(const double *W, int v,
    int i, int j, int k, double target)
{
    double q[PRJ_RECON_RADIATION_NCELLS];
    const int half = PRJ_RECON_RADIATION_NCELLS / 2;
    int n;

    for (n = 0; n < PRJ_RECON_RADIATION_NCELLS; ++n) {
        int offset = n - half;
        q[n] = prj_flux_radiation_stencil_value(W[WIDX(v, i + offset, j, k)]);
    }

    return prj_reconstruct_radiation_face_value(q, target);
}

static inline double prj_flux_radiation_prim_face_value_x2(const double *W, int v,
    int i, int j, int k, double target)
{
    double q[PRJ_RECON_RADIATION_NCELLS];
    const int half = PRJ_RECON_RADIATION_NCELLS / 2;
    int n;

    for (n = 0; n < PRJ_RECON_RADIATION_NCELLS; ++n) {
        int offset = n - half;
        q[n] = prj_flux_radiation_stencil_value(W[WIDX(v, i, j + offset, k)]);
    }

    return prj_reconstruct_radiation_face_value(q, target);
}

static inline double prj_flux_radiation_prim_face_value_x3(const double *W, int v,
    int i, int j, int k, double target)
{
    double q[PRJ_RECON_RADIATION_NCELLS];
    const int half = PRJ_RECON_RADIATION_NCELLS / 2;
    int n;

    for (n = 0; n < PRJ_RECON_RADIATION_NCELLS; ++n) {
        int offset = n - half;
        q[n] = prj_flux_radiation_stencil_value(W[WIDX(v, i, j, k + offset)]);
    }

    return prj_reconstruct_radiation_face_value(q, target);
}

static inline double prj_flux_eos_face_value_x1(const double *eosvar, int v,
    int i, int j, int k, double target)
{
    double q[PRJ_RECON_HYDRO_NCELLS];
    const int half = PRJ_RECON_HYDRO_NCELLS / 2;
    int n;

    for (n = 0; n < PRJ_RECON_HYDRO_NCELLS; ++n) {
        int offset = n - half;
        q[n] = eosvar[EIDX(v, i + offset, j, k)];
    }

    return prj_reconstruct_hydro_face_value(q, target);
}

static inline double prj_flux_eos_face_value_x2(const double *eosvar, int v,
    int i, int j, int k, double target)
{
    double q[PRJ_RECON_HYDRO_NCELLS];
    const int half = PRJ_RECON_HYDRO_NCELLS / 2;
    int n;

    for (n = 0; n < PRJ_RECON_HYDRO_NCELLS; ++n) {
        int offset = n - half;
        q[n] = eosvar[EIDX(v, i, j + offset, k)];
    }

    return prj_reconstruct_hydro_face_value(q, target);
}

static inline double prj_flux_eos_face_value_x3(const double *eosvar, int v,
    int i, int j, int k, double target)
{
    double q[PRJ_RECON_HYDRO_NCELLS];
    const int half = PRJ_RECON_HYDRO_NCELLS / 2;
    int n;

    for (n = 0; n < PRJ_RECON_HYDRO_NCELLS; ++n) {
        int offset = n - half;
        q[n] = eosvar[EIDX(v, i, j, k + offset)];
    }

    return prj_reconstruct_hydro_face_value(q, target);
}
#endif /* !PRJ_FLUX_RECON_MC_FASTPATH */

/* Build the local->global primitive index permutation for direction `dir`.
 * The sweep-normal velocity/B/flux component is rotated into the local X1 slot
 * so the Riemann solver always sees the sweep direction as its normal. The map
 * is a permutation over all PRJ_NVAR_PRIM indices. */
static void prj_flux_build_var_map(int dir, int *gmap)
{
    int field = 0;
    int group = 0;
#if PRJ_USE_RADIATION_FSA
    int angle = 0;
#endif

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
#if PRJ_USE_RADIATION_FSA
    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            for (angle = 0; angle < PRJ_NANGLE; ++angle) {
                gmap[PRJ_PRIM_RAD_I(field, group, angle)] =
                    PRJ_PRIM_RAD_I(field, group, angle);
            }
        }
    }
#else
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
#endif
#else
    (void)field;
    (void)group;
#endif
}

#if PRJ_FLUX_RECON_MC_FASTPATH
enum prj_flux_recon_src_kind {
    PRJ_FLUX_RECON_SRC_BLOCK = 0,
    PRJ_FLUX_RECON_SRC_PRIM_STAGE = 1
};

static inline double prj_flux_recon_stencil_val(double q, int is_rad)
{
#if PRJ_NRAD > 0 && PRJ_MIXED_PRECISION_FLUX
    return is_rad ? prj_rad_mixed_round(q) : q;
#else
    (void)is_rad;
    return q;
#endif
}

/* MC fast path for one variable: compute each cell's limited slope ONCE along the
 * sweep, then emit both faces it touches (q0 - 0.5*slope as the right state of the
 * lower face, q0 + 0.5*slope as the left state of the upper face). This halves the
 * limiter evaluations vs reconstructing each face/side independently, and is
 * bit-identical because both faces share the same slope. `slab` is the source
 * variable plane (W or eosvar), `wl`/`wr` the left/right output planes. */
static void prj_flux_recon_var_mc(const double *src, int nc, int gv,
    double *wl, double *wr,
    int dir, int istart, int iend, int jstart, int jend, int kstart, int kend,
    int is_rad, enum prj_flux_recon_src_kind src_kind)
{
    /* Component gv is first gathered from the block buffer `src` (nc components
     * per cell, SoA or AoSoA) into a contiguous untiled SoA plane `slab`.  The
     * MC sweep then streams `slab` with stride 1/BS/BS^2 and full cache-line
     * use, instead of the AoSoA per-variable interleave that wastes ~half of
     * every line and would read each cell three times in the stencil.  Reading
     * the AoSoA data once (the gather) is cheaper than the 3x scattered stencil
     * reads, and the slab sweep is the original bit-identical SoA kernel. */
    static double slab[PRJ_BLOCK_NCELLS];
    int i;
    int j;
    int k;

    for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
        for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
            for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                slab[LIDX(i, j, k)] = src_kind == PRJ_FLUX_RECON_SRC_PRIM_STAGE ?
                    src[WIDX(gv, i, j, k)] : src[BIDX(nc, gv, i, j, k)];
            }
        }
    }

    if (dir == X1DIR) {
        const size_t S = (size_t)PRJ_BS * PRJ_BS;
        for (i = istart - 1; i <= iend; ++i) {
            for (j = jstart; j <= jend; ++j) {
                for (k = kstart; k <= kend; ++k) {
                    size_t idx = (size_t)LIDX(i, j, k);
                    double qm = prj_flux_recon_stencil_val(slab[idx - S], is_rad);
                    double q0 = prj_flux_recon_stencil_val(slab[idx], is_rad);
                    double qp = prj_flux_recon_stencil_val(slab[idx + S], is_rad);
                    double h = 0.5 * prj_reconstruct_mc_face_slope(qm, q0, qp);
                    if (i >= istart) wr[idx] = q0 - h;
                    if (i < iend) wl[idx + S] = q0 + h;
                }
            }
        }
    } else if (dir == X2DIR) {
        const size_t S = (size_t)PRJ_BS;
        for (i = istart; i <= iend; ++i) {
            for (j = jstart - 1; j <= jend; ++j) {
                for (k = kstart; k <= kend; ++k) {
                    size_t idx = (size_t)LIDX(i, j, k);
                    double qm = prj_flux_recon_stencil_val(slab[idx - S], is_rad);
                    double q0 = prj_flux_recon_stencil_val(slab[idx], is_rad);
                    double qp = prj_flux_recon_stencil_val(slab[idx + S], is_rad);
                    double h = 0.5 * prj_reconstruct_mc_face_slope(qm, q0, qp);
                    if (j >= jstart) wr[idx] = q0 - h;
                    if (j < jend) wl[idx + S] = q0 + h;
                }
            }
        }
    } else {
        const size_t S = 1;
        for (i = istart; i <= iend; ++i) {
            for (j = jstart; j <= jend; ++j) {
                for (k = kstart - 1; k <= kend; ++k) {
                    size_t idx = (size_t)LIDX(i, j, k);
                    double qm = prj_flux_recon_stencil_val(slab[idx - S], is_rad);
                    double q0 = prj_flux_recon_stencil_val(slab[idx], is_rad);
                    double qp = prj_flux_recon_stencil_val(slab[idx + S], is_rad);
                    double h = 0.5 * prj_reconstruct_mc_face_slope(qm, q0, qp);
                    if (k >= kstart) wr[idx] = q0 - h;
                    if (k < kend) wl[idx + S] = q0 + h;
                }
            }
        }
    }
}
#endif

/* Reconstruct every primitive (and the pressure/gamma EOS vars) onto the dir
 * faces of a block, with the *variable* loop outermost.  Results are written to
 * block-wide buffers: WL_block/WR_block use layout [local_var * PRJ_BLOCK_NCELLS
 * + IDX(i, j, k)]; the pL/pR/gL/gR buffers are indexed by IDX(i, j, k).  The
 * Riemann sweep in prj_flux_update consumes these per face. */
#define PRJ_FLUX_RECONSTRUCT_PRIM_RANGE(FACE_CELLS, PRIM_FACE_VALUE, LV_BEGIN, LV_END) do { \
        int lv; \
        for (lv = (LV_BEGIN); lv < (LV_END); ++lv) { \
            int gv = gmap[lv]; \
            double *wl = &WL_block[VPLANE((size_t)lv)]; \
            double *wr = &WR_block[VPLANE((size_t)lv)]; \
            int i; \
            for (i = istart; i <= iend; ++i) { \
                int j; \
                for (j = jstart; j <= jend; ++j) { \
                    int k; \
                    for (k = kstart; k <= kend; ++k) { \
                        int il; \
                        int jl; \
                        int kl; \
                        int ir; \
                        int jr; \
                        int kr; \
                        FACE_CELLS(i, j, k, &il, &jl, &kl, &ir, &jr, &kr); \
                        wl[LIDX(i, j, k)] = PRIM_FACE_VALUE(W, gv, il, jl, kl, 0.5); \
                        wr[LIDX(i, j, k)] = PRIM_FACE_VALUE(W, gv, ir, jr, kr, -0.5); \
                    } \
                } \
            } \
        } \
    } while (0)

#define PRJ_FLUX_RECONSTRUCT_PRIM_DIR(FACE_CELLS, HYDRO_FACE_VALUE, RADIATION_FACE_VALUE) do { \
        PRJ_FLUX_RECONSTRUCT_PRIM_RANGE(FACE_CELLS, HYDRO_FACE_VALUE, 0, PRJ_NHYDRO); \
        PRJ_FLUX_RECONSTRUCT_PRIM_RANGE(FACE_CELLS, RADIATION_FACE_VALUE, \
            PRJ_NHYDRO, PRJ_NVAR_PRIM); \
    } while (0)

#define PRJ_FLUX_RECONSTRUCT_EOS_DIR(FACE_CELLS, EOS_FACE_VALUE) do { \
        if (eosvar != 0) { \
            const int eos_var[2] = { PRJ_EOSVAR_PRESSURE, PRJ_EOSVAR_GAMMA }; \
            double *eosL[2]; \
            double *eosR[2]; \
            int e; \
            eosL[0] = pL_block; \
            eosL[1] = gL_block; \
            eosR[0] = pR_block; \
            eosR[1] = gR_block; \
            for (e = 0; e < 2; ++e) { \
                int ev = eos_var[e]; \
                double *el = eosL[e]; \
                double *er = eosR[e]; \
                int i; \
                for (i = istart; i <= iend; ++i) { \
                    int j; \
                    for (j = jstart; j <= jend; ++j) { \
                        int k; \
                        for (k = kstart; k <= kend; ++k) { \
                            int il; \
                            int jl; \
                            int kl; \
                            int ir; \
                            int jr; \
                            int kr; \
                            FACE_CELLS(i, j, k, &il, &jl, &kl, &ir, &jr, &kr); \
                            el[LIDX(i, j, k)] = EOS_FACE_VALUE(eosvar, ev, il, jl, kl, 0.5); \
                            er[LIDX(i, j, k)] = EOS_FACE_VALUE(eosvar, ev, ir, jr, kr, -0.5); \
                        } \
                    } \
                } \
            } \
        } \
    } while (0)

static void prj_flux_reconstruct_block(double *W, const double *eosvar, int dir,
    int istart, int iend, int jstart, int jend, int kstart, int kend,
    double *WL_block, double *WR_block,
    double *pL_block, double *pR_block, double *gL_block, double *gR_block)
{
    int gmap[PRJ_NVAR_PRIM];

    prj_flux_build_var_map(dir, gmap);

#if PRJ_FLUX_RECON_MC_FASTPATH
    {
        int lv;

        for (lv = 0; lv < PRJ_NHYDRO; ++lv) {
            int gv = gmap[lv];
            prj_flux_recon_var_mc(W, PRJ_NVAR_CONS, gv,
                &WL_block[VPLANE((size_t)lv)],
                &WR_block[VPLANE((size_t)lv)],
                dir, istart, iend, jstart, jend, kstart, kend, 0,
                PRJ_FLUX_RECON_SRC_PRIM_STAGE);
        }
        for (lv = PRJ_NHYDRO; lv < PRJ_NVAR_PRIM; ++lv) {
            int gv = gmap[lv];
            prj_flux_recon_var_mc(W, PRJ_NVAR_CONS, gv,
                &WL_block[VPLANE((size_t)lv)],
                &WR_block[VPLANE((size_t)lv)],
                dir, istart, iend, jstart, jend, kstart, kend, 1,
                PRJ_FLUX_RECON_SRC_PRIM_STAGE);
        }
        if (eosvar != 0) {
            prj_flux_recon_var_mc(eosvar, PRJ_NVAR_EOSVAR, PRJ_EOSVAR_PRESSURE,
                pL_block, pR_block, dir, istart, iend, jstart, jend, kstart, kend, 0,
                PRJ_FLUX_RECON_SRC_BLOCK);
            prj_flux_recon_var_mc(eosvar, PRJ_NVAR_EOSVAR, PRJ_EOSVAR_GAMMA,
                gL_block, gR_block, dir, istart, iend, jstart, jend, kstart, kend, 0,
                PRJ_FLUX_RECON_SRC_BLOCK);
        }
    }
#else
    if (dir == X1DIR) {
        PRJ_FLUX_RECONSTRUCT_PRIM_DIR(prj_flux_face_cells_x1,
            prj_flux_hydro_prim_face_value_x1, prj_flux_radiation_prim_face_value_x1);
        PRJ_FLUX_RECONSTRUCT_EOS_DIR(prj_flux_face_cells_x1, prj_flux_eos_face_value_x1);
    } else if (dir == X2DIR) {
        PRJ_FLUX_RECONSTRUCT_PRIM_DIR(prj_flux_face_cells_x2,
            prj_flux_hydro_prim_face_value_x2, prj_flux_radiation_prim_face_value_x2);
        PRJ_FLUX_RECONSTRUCT_EOS_DIR(prj_flux_face_cells_x2, prj_flux_eos_face_value_x2);
    } else {
        PRJ_FLUX_RECONSTRUCT_PRIM_DIR(prj_flux_face_cells_x3,
            prj_flux_hydro_prim_face_value_x3, prj_flux_radiation_prim_face_value_x3);
        PRJ_FLUX_RECONSTRUCT_EOS_DIR(prj_flux_face_cells_x3, prj_flux_eos_face_value_x3);
    }
#endif
}

#undef PRJ_FLUX_RECONSTRUCT_PRIM_RANGE
#undef PRJ_FLUX_RECONSTRUCT_PRIM_DIR
#undef PRJ_FLUX_RECONSTRUCT_EOS_DIR

static void prj_flux_store_face_velocity(prj_block *block, int dir, int i, int j, int k,
    const double v_face_loc[3])
{
    double *dst;

    if (block->v_riemann[dir] == 0) {
        return;
    }
    dst = block->v_riemann[dir];
    if (dir == X1DIR) {
        dst[VRIDX(0, i, j, k)] = v_face_loc[0];
        dst[VRIDX(1, i, j, k)] = v_face_loc[1];
        dst[VRIDX(2, i, j, k)] = v_face_loc[2];
    } else if (dir == X2DIR) {
        dst[VRIDX(1, i, j, k)] = v_face_loc[0];
        dst[VRIDX(2, i, j, k)] = v_face_loc[1];
        dst[VRIDX(0, i, j, k)] = v_face_loc[2];
    } else {
        dst[VRIDX(2, i, j, k)] = v_face_loc[0];
        dst[VRIDX(0, i, j, k)] = v_face_loc[1];
        dst[VRIDX(1, i, j, k)] = v_face_loc[2];
    }
}

static void prj_flux_store_local_flux(double *dst, int dir, int i, int j, int k,
    const double *Fl)
{
    int field=0;
    int group=0;
#if PRJ_USE_RADIATION_FSA
    int angle=0;
#endif

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
    /* The radiation flux divergence only reads faces with transverse indices in
     * the active range; on MHD transverse ghost faces these stores are never
     * consumed, so skip the radiation-variable copy there. */
    if (prj_flux_rad_face_needed(dir, i, j, k)) {
#if PRJ_USE_RADIATION_FSA
    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            for (angle = 0; angle < PRJ_NANGLE; ++angle) {
                dst[VIDX(PRJ_CONS_RAD_I(field, group, angle), i, j, k)] =
                    Fl[PRJ_CONS_RAD_I(field, group, angle)];
            }
        }
    }
#else
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
#endif
    }
#else
    (void)field;
    (void)group;
#endif
}

#if PRJ_DYNAMIC_GR
static const char *prj_flux_gr_status_name(int status)
{
    switch (status) {
    case PRJ_EOS_GR_OK:
        return "ok";
    case PRJ_EOS_GR_NULL_ARG:
        return "null_arg";
    case PRJ_EOS_GR_BAD_METRIC:
        return "bad_metric";
    case PRJ_EOS_GR_BAD_STATE:
        return "bad_state";
    case PRJ_EOS_GR_NO_CONVERGE:
        return "no_converge";
    default:
        return "missing_geometry";
    }
}

static const char *prj_flux_prim_name(int var)
{
    static const char *const names[6] = {
        "RHO", "V1", "V2", "V3", "EINT", "YE"
    };

    return var >= 0 && var < 6 ? names[var] : "RAD";
}

static const char *prj_flux_cons_name(int var)
{
    static const char *const names[6] = {
        "RHO", "MOM1", "MOM2", "MOM3", "ETOT", "YE"
    };

    return var >= 0 && var < 6 ? names[var] : "RAD";
}

static const char *prj_flux_eosvar_name(int var)
{
    static const char *const names[3] = {
        "pressure", "temperature", "gamma"
    };

    return var >= 0 && var < 3 ? names[var] : "unknown";
}

static void prj_flux_gr_abort(void)
{
#if defined(PRJ_ENABLE_MPI)
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#endif
    exit(EXIT_FAILURE);
}

static void prj_flux_gr_fail(const char *op, int status, int dir, int i, int j, int k)
{
    fprintf(stderr,
        "dynamic GR flux %s failed at dir=%d face (%d,%d,%d): status=%d (%s)\n",
        op, dir, i, j, k, status, prj_flux_gr_status_name(status));
    fflush(stderr);
    prj_flux_gr_abort();
}

static void prj_flux_local_axes(int dir, int axis[3])
{
    if (dir == X1DIR) {
        axis[0] = 0;
        axis[1] = 1;
        axis[2] = 2;
    } else if (dir == X2DIR) {
        axis[0] = 1;
        axis[1] = 2;
        axis[2] = 0;
    } else {
        axis[0] = 2;
        axis[1] = 0;
        axis[2] = 1;
    }
}

static double prj_flux_det3(const double g[3][3])
{
    return g[0][0] * (g[1][1] * g[2][2] - g[1][2] * g[2][1])
        - g[0][1] * (g[1][0] * g[2][2] - g[1][2] * g[2][0])
        + g[0][2] * (g[1][0] * g[2][1] - g[1][1] * g[2][0]);
}

/* Cofactor inverse of a 3-metric; returns det(g). On a non-finite or
 * non-positive determinant the inverse is left untouched (the caller keeps its
 * fallback and the degenerate metric is caught downstream). Keeping gamma_inv
 * the exact inverse of gamma is what makes face norms, raising/lowering and
 * densitization mutually consistent. */
static double prj_flux_inv3(const double g[3][3], double ginv[3][3])
{
    double det = prj_flux_det3(g);
    double inv_det;

    if (!isfinite(det) || det <= 0.0) {
        return det;
    }
    inv_det = 1.0 / det;
    ginv[0][0] = (g[1][1] * g[2][2] - g[1][2] * g[2][1]) * inv_det;
    ginv[0][1] = (g[0][2] * g[2][1] - g[0][1] * g[2][2]) * inv_det;
    ginv[0][2] = (g[0][1] * g[1][2] - g[0][2] * g[1][1]) * inv_det;
    ginv[1][0] = (g[1][2] * g[2][0] - g[1][0] * g[2][2]) * inv_det;
    ginv[1][1] = (g[0][0] * g[2][2] - g[0][2] * g[2][0]) * inv_det;
    ginv[1][2] = (g[0][2] * g[1][0] - g[0][0] * g[1][2]) * inv_det;
    ginv[2][0] = (g[1][0] * g[2][1] - g[1][1] * g[2][0]) * inv_det;
    ginv[2][1] = (g[0][1] * g[2][0] - g[0][0] * g[2][1]) * inv_det;
    ginv[2][2] = (g[0][0] * g[1][1] - g[0][1] * g[1][0]) * inv_det;
    return det;
}

static void prj_flux_gr_print_face_location(const prj_block *block,
    int dir, int i, int j, int k)
{
    double x[3];

    if (block == 0) {
        return;
    }
    x[0] = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
    x[1] = block->xmin[1] + ((double)j + 0.5) * block->dx[1];
    x[2] = block->xmin[2] + ((double)k + 0.5) * block->dx[2];
    x[dir] = block->xmin[dir] + (double)(dir == X1DIR ? i : dir == X2DIR ? j : k) *
        block->dx[dir];
    fprintf(stderr, "  face x=(%.17e, %.17e, %.17e) face_index=%zu\n",
        x[0], x[1], x[2], (size_t)FACE_IDX(dir, i, j, k));
}

static void prj_flux_gr_print_face_geom(const prj_z4c_hydro_geom *geom)
{
    double det;
    int a;

    if (geom == 0) {
        fprintf(stderr, "  face geometry is NULL\n");
        return;
    }
    det = prj_flux_det3(geom->gamma);
    fprintf(stderr,
        "  face geometry: alpha=%.17e beta=(%.17e, %.17e, %.17e) "
        "sqrt_gamma=%.17e det(gamma)=%.17e\n",
        geom->alpha, geom->beta[0], geom->beta[1], geom->beta[2],
        geom->sqrt_gamma, det);
    fprintf(stderr, "  face gamma_ij in local flux axes:\n");
    for (a = 0; a < 3; ++a) {
        fprintf(stderr, "    [%.17e, %.17e, %.17e]\n",
            geom->gamma[a][0], geom->gamma[a][1], geom->gamma[a][2]);
    }
    fprintf(stderr, "  face gamma^ij in local flux axes:\n");
    for (a = 0; a < 3; ++a) {
        fprintf(stderr, "    [%.17e, %.17e, %.17e]\n",
            geom->gamma_inv[a][0], geom->gamma_inv[a][1], geom->gamma_inv[a][2]);
    }
}

static void prj_flux_gr_print_prim_vector(const char *title, const double *W)
{
    int v;

    if (W == 0) {
        fprintf(stderr, "  %s is NULL\n", title);
        return;
    }
    fprintf(stderr, "  %s:\n", title);
    for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
        fprintf(stderr, "    W[%2d] %-5s = %.17e\n", v,
            prj_flux_prim_name(v), W[v]);
    }
}

static void prj_flux_gr_print_cons_vector(const char *title, const double *U)
{
    int v;

    if (U == 0) {
        fprintf(stderr, "  %s is NULL\n", title);
        return;
    }
    fprintf(stderr, "  %s:\n", title);
    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        fprintf(stderr, "    U[%2d] %-5s = %.17e\n", v,
            prj_flux_cons_name(v), U[v]);
    }
}

static void prj_flux_gr_print_velocity_diag(const prj_z4c_hydro_geom *geom,
    const double *W)
{
    double beta[3];
    double beta2 = 0.0;
    int a;
    int b;

    if (geom == 0 || W == 0) {
        fprintf(stderr, "  face velocity diagnostics are unavailable\n");
        return;
    }
    beta[0] = W[PRJ_PRIM_V1] / PRJ_CLIGHT;
    beta[1] = W[PRJ_PRIM_V2] / PRJ_CLIGHT;
    beta[2] = W[PRJ_PRIM_V3] / PRJ_CLIGHT;
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            beta2 += geom->gamma[a][b] * beta[a] * beta[b];
        }
    }
    fprintf(stderr,
        "  face velocity diagnostics: beta=(%.17e, %.17e, %.17e) "
        "gamma_ij beta^i beta^j=%.17e lorentz=% .17e\n",
        beta[0], beta[1], beta[2], beta2,
        beta2 >= 0.0 && beta2 < 1.0 ? 1.0 / sqrt(1.0 - beta2) : NAN);
}

static void prj_flux_gr_print_source_cell(const char *side, const prj_block *block,
    const double *W_block, const double *eosvar, int z4c_stage,
    int i, int j, int k)
{
    double x[3];
    int v;

    if (block == 0) {
        fprintf(stderr, "  %s source cell: block is NULL\n", side);
        return;
    }

    x[0] = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
    x[1] = block->xmin[1] + ((double)j + 0.5) * block->dx[1];
    x[2] = block->xmin[2] + ((double)k + 0.5) * block->dx[2];
    fprintf(stderr,
        "  %s source cell: cell=(%d,%d,%d) cell_lidx=%zu storage_idx=%zu "
        "z4c_lidx=%zu x=(%.17e, %.17e, %.17e)\n",
        side, i, j, k, (size_t)LIDX(i, j, k), (size_t)IDX(i, j, k),
        (size_t)Z4CLIDX(i, j, k), x[0], x[1], x[2]);

    if (W_block != 0) {
        fprintf(stderr, "  %s stored primitive cell state:\n", side);
        for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
            fprintf(stderr, "    Wcell[%2d] %-5s = %.17e\n", v,
                prj_flux_prim_name(v), W_block[WIDX(v, i, j, k)]);
        }
    } else {
        fprintf(stderr, "  %s stored primitive cell state is unavailable\n", side);
    }

    if (eosvar != 0) {
        fprintf(stderr, "  %s stored EOS cell variables:\n", side);
        for (v = 0; v < PRJ_NVAR_EOSVAR; ++v) {
            fprintf(stderr, "    eosvar[%d] %-11s = %.17e\n", v,
                prj_flux_eosvar_name(v), eosvar[EIDX(v, i, j, k)]);
        }
    } else {
        fprintf(stderr, "  %s stored EOS cell variables are unavailable\n", side);
    }

    if (prj_block_has_cons_storage(block)) {
        fprintf(stderr, "  %s stored conserved cell state:\n", side);
        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
            fprintf(stderr, "    Ucell[%2d] %-5s = %.17e\n", v,
                prj_flux_cons_name(v), prj_block_cons_value_const(block, v, i, j, k));
        }
    } else {
        fprintf(stderr, "  %s stored conserved cell state is unavailable\n", side);
    }

    {
        const double *z = prj_block_z4c_stage_const(block, z4c_stage);

        if (z != 0) {
            fprintf(stderr, "  %s Z4c variables (stage %d) at cell (%d,%d,%d):\n",
                side, z4c_stage, i, j, k);
            for (v = 0; v < PRJ_NZ4C; ++v) {
                fprintf(stderr, "    z4c[%2d] %-16s = %.17e\n", v,
                    prj_z4c_var_name(v), z[Z4CIDX(v, i, j, k)]);
            }
        } else {
            fprintf(stderr, "  %s Z4c array is unavailable for stage %d\n",
                side, z4c_stage);
        }
    }
}

static void prj_flux_gr_state_fail(const char *op, int status,
    const prj_block *block, const double *W_block, const double *eosvar,
    int z4c_stage, const prj_z4c_hydro_geom *geom, const double *W_face,
    double pressure, double gas_gamma, const double *U_face, int dir,
    int face_i, int face_j, int face_k, const char *side,
    int cell_i, int cell_j, int cell_k)
{
    fprintf(stderr,
        "dynamic GR flux %s failed at dir=%d face (%d,%d,%d): "
        "status=%d (%s) side=%s source_cell=(%d,%d,%d)\n",
        op, dir, face_i, face_j, face_k, status,
        prj_flux_gr_status_name(status), side, cell_i, cell_j, cell_k);
    fprintf(stderr,
        "  block id=%d rank=%d level=%d active=%d parent=%d z4c_stage=%d "
        "xmin=(%.17e, %.17e, %.17e) xmax=(%.17e, %.17e, %.17e) "
        "dx=(%.17e, %.17e, %.17e)\n",
        block != 0 ? block->id : -1, block != 0 ? block->rank : -1,
        block != 0 ? block->level : -1, block != 0 ? block->active : -1,
        block != 0 ? block->parent : -1, z4c_stage,
        block != 0 ? block->xmin[0] : 0.0, block != 0 ? block->xmin[1] : 0.0,
        block != 0 ? block->xmin[2] : 0.0, block != 0 ? block->xmax[0] : 0.0,
        block != 0 ? block->xmax[1] : 0.0, block != 0 ? block->xmax[2] : 0.0,
        block != 0 ? block->dx[0] : 0.0, block != 0 ? block->dx[1] : 0.0,
        block != 0 ? block->dx[2] : 0.0);
    prj_flux_gr_print_face_location(block, dir, face_i, face_j, face_k);
    fprintf(stderr, "  reconstructed face EOS: pressure=%.17e gamma=%.17e\n",
        pressure, gas_gamma);
    prj_flux_gr_print_prim_vector("reconstructed face primitive state", W_face);
    prj_flux_gr_print_velocity_diag(geom, W_face);
    prj_flux_gr_print_cons_vector("prim2cons output U (partial/invalid on failure)",
        U_face);
    prj_flux_gr_print_face_geom(geom);
    prj_flux_gr_print_source_cell(side, block, W_block, eosvar, z4c_stage,
        cell_i, cell_j, cell_k);
    fflush(stderr);
    prj_flux_gr_abort();
}

static void prj_flux_gr_face_geom_from_cells(const prj_z4c_hydro_geom *gl,
    const prj_z4c_hydro_geom *gr, int dir, prj_z4c_hydro_geom *geom)
{
    int axis[3];
    int a;
    int b;
    int d;

    prj_flux_local_axes(dir, axis);
    geom->alpha = 0.5 * (gl->alpha + gr->alpha);
    geom->sqrt_gamma = 0.5 * (gl->sqrt_gamma + gr->sqrt_gamma);
    for (a = 0; a < 3; ++a) {
        geom->beta[a] = 0.5 * (gl->beta[axis[a]] + gr->beta[axis[a]]);
        geom->dalpha[a] = 0.5 * (gl->dalpha[axis[a]] + gr->dalpha[axis[a]]);
        for (b = 0; b < 3; ++b) {
            geom->gamma[a][b] =
                0.5 * (gl->gamma[axis[a]][axis[b]] + gr->gamma[axis[a]][axis[b]]);
            geom->gamma_inv[a][b] =
                0.5 * (gl->gamma_inv[axis[a]][axis[b]] + gr->gamma_inv[axis[a]][axis[b]]);
            geom->K_dd[a][b] =
                0.5 * (gl->K_dd[axis[a]][axis[b]] + gr->K_dd[axis[a]][axis[b]]);
            geom->dbeta[a][b] =
                0.5 * (gl->dbeta[axis[a]][axis[b]] + gr->dbeta[axis[a]][axis[b]]);
            for (d = 0; d < 3; ++d) {
                geom->dgamma[d][a][b] =
                    0.5 * (gl->dgamma[axis[d]][axis[a]][axis[b]] +
                        gr->dgamma[axis[d]][axis[a]][axis[b]]);
            }
        }
    }
    /* Derive the face inverse-metric and sqrt(det) from the averaged gamma so
     * they stay mutually consistent (gamma_inv == inverse(gamma),
     * sqrt_gamma == sqrt(det gamma)). The independent averages above are kept
     * only as a fallback for a degenerate metric; using them directly leaves
     * the fields inconsistent and corrupts norms, raising/lowering, flux
     * densitization and the GRMHD HLLD input on curved/off-diagonal metrics. */
    {
        double det = prj_flux_inv3((const double (*)[3])geom->gamma, geom->gamma_inv);

        if (isfinite(det) && det > 0.0) {
            geom->sqrt_gamma = sqrt(det);
        }
    }
}

static int prj_flux_gr_face_geom(const prj_mesh *mesh, const prj_block *block,
    int stage, int dir, int il, int jl, int kl, int ir, int jr, int kr,
    prj_z4c_hydro_geom *geom)
{
    prj_z4c_hydro_geom gl;
    prj_z4c_hydro_geom gr;

    if (!prj_z4c_load_hydro_geom(mesh, block, stage, il, jl, kl, &gl) ||
        !prj_z4c_load_hydro_geom(mesh, block, stage, ir, jr, kr, &gr)) {
        return 0;
    }
    prj_flux_gr_face_geom_from_cells(&gl, &gr, dir, geom);
    return 1;
}

static int prj_flux_gr_face_geom_cached(const prj_z4c_hydro_geom *cell_geom,
    int dir, int il, int jl, int kl, int ir, int jr, int kr,
    prj_z4c_hydro_geom *geom)
{
    if (cell_geom == 0 || geom == 0 ||
        il < -PRJ_NGHOST || il >= PRJ_BLOCK_SIZE + PRJ_NGHOST ||
        jl < -PRJ_NGHOST || jl >= PRJ_BLOCK_SIZE + PRJ_NGHOST ||
        kl < -PRJ_NGHOST || kl >= PRJ_BLOCK_SIZE + PRJ_NGHOST ||
        ir < -PRJ_NGHOST || ir >= PRJ_BLOCK_SIZE + PRJ_NGHOST ||
        jr < -PRJ_NGHOST || jr >= PRJ_BLOCK_SIZE + PRJ_NGHOST ||
        kr < -PRJ_NGHOST || kr >= PRJ_BLOCK_SIZE + PRJ_NGHOST) {
        return 0;
    }
    prj_flux_gr_face_geom_from_cells(
        &cell_geom[LIDX(il, jl, kl)], &cell_geom[LIDX(ir, jr, kr)], dir, geom);
    return 1;
}

static void prj_flux_gr_hydro_state_flux(prj_eos *eos, const prj_block *block,
    const double *W_block, const double *eosvar, int z4c_stage,
    const prj_z4c_hydro_geom *geom, const double *W, double pressure,
    double gas_gamma, double *U, double *F, double *speed, int dir, int i, int j, int k,
    const char *side, int cell_i, int cell_j, int cell_k)
{
    prj_eos_gr_geom egeom;
    double Uloc[PRJ_NVAR_CONS];
    double det;
    double sqrtg;
    double vn;
    double cs2;
    double cs;
    double gnn_inv;
    int status;
    int a;
    int b;
    int v;

    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        U[v] = NAN;
    }
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            egeom.gamma[a][b] = geom->gamma[a][b];
        }
    }
    status = prj_eos_gr_prim2cons(eos, &egeom, W, U, PRJ_EOS_CTX_MAIN);
    if (status != PRJ_EOS_GR_OK) {
        prj_flux_gr_state_fail("prim2cons", status, block, W_block, eosvar,
            z4c_stage, geom, W, pressure, gas_gamma, U, dir, i, j, k, side,
            cell_i, cell_j, cell_k);
    }
    det = prj_flux_det3(geom->gamma);
    if (!isfinite(det) || det <= 0.0) {
        prj_flux_gr_state_fail("metric", PRJ_EOS_GR_BAD_METRIC, block, W_block,
            eosvar, z4c_stage, geom, W, pressure, gas_gamma, U, dir, i, j, k,
            side, cell_i, cell_j, cell_k);
    }
    sqrtg = sqrt(det);
    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        Uloc[v] = U[v] / sqrtg;
        F[v] = 0.0;
    }
    vn = geom->alpha * W[PRJ_PRIM_V1] - PRJ_CLIGHT * geom->beta[0];
    F[PRJ_CONS_RHO] = sqrtg * Uloc[PRJ_CONS_RHO] * vn;
    F[PRJ_CONS_MOM1] = sqrtg * (Uloc[PRJ_CONS_MOM1] * vn + geom->alpha * pressure);
    F[PRJ_CONS_MOM2] = sqrtg * Uloc[PRJ_CONS_MOM2] * vn;
    F[PRJ_CONS_MOM3] = sqrtg * Uloc[PRJ_CONS_MOM3] * vn;
    F[PRJ_CONS_ETOT] =
        sqrtg * (Uloc[PRJ_CONS_ETOT] * vn + geom->alpha * pressure * W[PRJ_PRIM_V1]);
    F[PRJ_CONS_YE] = sqrtg * Uloc[PRJ_CONS_YE] * vn;

    cs2 = gas_gamma * pressure / W[PRJ_PRIM_RHO];
    if (!isfinite(cs2) || cs2 < 0.0) {
        cs2 = 0.0;
    }
    if (cs2 > PRJ_CLIGHT * PRJ_CLIGHT) {
        cs2 = PRJ_CLIGHT * PRJ_CLIGHT;
    }
    cs = sqrt(cs2);
    /* The normal characteristic speed scales with the norm of the face-normal
     * one-form, sqrt(gamma^nn) = sqrt(gamma_inv[0][0]). 1/sqrt(gamma[0][0])
     * only coincides with it for a diagonal face metric. */
    gnn_inv = geom->gamma_inv[0][0];
    if (!isfinite(gnn_inv) || gnn_inv <= 0.0) {
        gnn_inv = 1.0;
    }
    *speed = fabs(vn) + geom->alpha * cs * sqrt(gnn_inv);
}

static void prj_flux_gr_hydro_hll(prj_eos *eos, const prj_mesh *mesh,
    const prj_block *block, const double *W_block, const double *eosvar,
    int z4c_stage, int dir, int i, int j, int k, int il, int jl, int kl,
    int ir, int jr, int kr, const double *WL, const double *WR,
    double pL, double pR, double gL, double gR, double *Fl,
    double v_face_loc[3])
{
    prj_z4c_hydro_geom geom;
    double UL[PRJ_NVAR_CONS];
    double UR[PRJ_NVAR_CONS];
    double FL[PRJ_NVAR_CONS];
    double FR[PRJ_NVAR_CONS];
    double sL;
    double sR;
    double smax;
    int v;

    if (!prj_flux_gr_face_geom(mesh, block, z4c_stage, dir,
            il, jl, kl, ir, jr, kr, &geom)) {
        prj_flux_gr_fail("geometry load", -1, dir, i, j, k);
    }
    prj_flux_gr_hydro_state_flux(eos, block, W_block, eosvar, z4c_stage,
        &geom, WL, pL, gL, UL, FL, &sL, dir, i, j, k, "left", il, jl, kl);
    prj_flux_gr_hydro_state_flux(eos, block, W_block, eosvar, z4c_stage,
        &geom, WR, pR, gR, UR, FR, &sR, dir, i, j, k, "right", ir, jr, kr);
    smax = fmax(sL, sR);
    if (!isfinite(smax)) {
        prj_flux_gr_fail("wavespeed", PRJ_EOS_GR_BAD_STATE, dir, i, j, k);
    }
    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        Fl[v] = 0.5 * (FL[v] + FR[v]) - 0.5 * smax * (UR[v] - UL[v]);
    }
    if (Fl[PRJ_CONS_RHO] >= 0.0) {
        v_face_loc[0] = WL[PRJ_PRIM_V1];
        v_face_loc[1] = WL[PRJ_PRIM_V2];
        v_face_loc[2] = WL[PRJ_PRIM_V3];
    } else {
        v_face_loc[0] = WR[PRJ_PRIM_V1];
        v_face_loc[1] = WR[PRJ_PRIM_V2];
        v_face_loc[2] = WR[PRJ_PRIM_V3];
    }
}

#if PRJ_USE_RADIATION_M1
static void prj_flux_gr_m1_raise_vec(const prj_z4c_hydro_geom *geom,
    const double vcov[3], double vcon[3])
{
    int a;
    int b;

    for (a = 0; a < 3; ++a) {
        vcon[a] = 0.0;
        for (b = 0; b < 3; ++b) {
            vcon[a] += geom->gamma_inv[a][b] * vcov[b];
        }
    }
}

static double prj_flux_gr_m1_chi_exact(double f)
{
    if (f <= 0.0) {
        return 1.0 / 3.0;
    }
    if (f >= 1.0) {
        return 1.0;
    }
    return (3.0 + 4.0 * f * f) / (5.0 + 2.0 * sqrt(4.0 - 3.0 * f * f));
}

static double prj_flux_gr_m1_chi(const prj_rad *rad, double f)
{
    double scaled;
    double w;
    double chi;
    int idx;

    if (!isfinite(f) || f <= 0.0) {
        f = 0.0;
    } else if (f > 1.0) {
        f = 1.0;
    }
    if (rad == 0 || !isfinite(rad->chi[0]) || !isfinite(rad->chi[NCLOSURE]) ||
        rad->chi[0] <= 0.0 || rad->chi[NCLOSURE] <= 0.0) {
        return prj_flux_gr_m1_chi_exact(f);
    }
    if (f <= 0.0) {
        return rad->chi[0];
    }
    if (f >= 1.0) {
        return rad->chi[NCLOSURE];
    }

    scaled = f * (double)NCLOSURE;
    idx = (int)scaled;
    w = scaled - (double)idx;
    chi = rad->chi[idx] + w * (rad->chi[idx + 1] - rad->chi[idx]);
    return isfinite(chi) ? chi : prj_flux_gr_m1_chi_exact(f);
}

static void prj_flux_gr_m1_limit_state(const prj_z4c_hydro_geom *geom,
    double *E, double Fcov[3], double Fcon[3], double *Fmag_out)
{
    double F2;
    double Fmag;
    double cE;
    double scale;
    int d;

    if (!isfinite(*E) || *E < 0.0) {
        *E = 0.0;
    }
    prj_flux_gr_m1_raise_vec(geom, Fcov, Fcon);
    F2 = 0.0;
    for (d = 0; d < 3; ++d) {
        F2 += Fcov[d] * Fcon[d];
    }
    if (!isfinite(F2) || F2 < 0.0) {
        F2 = 0.0;
    }
    Fmag = sqrt(F2);
    cE = PRJ_CLIGHT * (*E);
    if (Fmag > cE && Fmag > 0.0) {
        scale = cE / Fmag;
        for (d = 0; d < 3; ++d) {
            Fcov[d] *= scale;
            Fcon[d] *= scale;
        }
        Fmag = cE;
    }
    *Fmag_out = Fmag;
}

static void prj_flux_gr_m1_wavespeeds(const prj_rad *rad,
    const prj_z4c_hydro_geom *geom, const prj_rad_gr_m1_side_data *side,
    const double Fcon[3], double Fmag, double zeta,
    double *smin, double *smax)
{
    double alpha = geom->alpha;
    double beta_d = geom->beta[0];
    double gdd = geom->gamma_inv[0][0];
    double wlor = side->wlor;
    double wlor2;
    double p;
    double r2;
    double r;
    double den;
    double lambda_thin_l;
    double lambda_thin_r;
    double lambda_thick_l_a;
    double lambda_thick_r_a;
    double lambda_fluid;
    double lambda_thick_l;
    double lambda_thick_r;
    double chi;
    double thin_w;
    double thick_w;
    double lambda_l;
    double lambda_r;

    if (!isfinite(alpha) || alpha <= 0.0) {
        alpha = 1.0;
    }
    if (!isfinite(beta_d)) {
        beta_d = 0.0;
    }
    if (!isfinite(gdd) || gdd <= 0.0) {
        gdd = 1.0;
    }
    if (!isfinite(wlor) || wlor <= 0.0) {
        wlor = 1.0;
    }

    if (Fmag > 0.0 && isfinite(Fmag)) {
        double thin_speed = alpha * fabs(Fcon[0]) / Fmag;

        lambda_thin_l = -beta_d - thin_speed;
        lambda_thin_r = -beta_d + thin_speed;
    } else {
        lambda_thin_l = -beta_d;
        lambda_thin_r = -beta_d;
    }

    wlor2 = wlor * wlor;
    p = alpha * side->vcon[0] / wlor;
    if (!isfinite(p)) {
        p = 0.0;
    }
    r2 = alpha * alpha * gdd * (2.0 * wlor2 + 1.0) -
        2.0 * wlor2 * p * p;
    if (!isfinite(r2) || r2 < 0.0) {
        r2 = 0.0;
    }
    r = sqrt(r2);
    den = 2.0 * wlor2 + 1.0;
    lambda_fluid = -beta_d + p;
    lambda_thick_l_a = -beta_d + (2.0 * p * wlor2 - r) / den;
    lambda_thick_r_a = -beta_d + (2.0 * p * wlor2 + r) / den;
    lambda_thick_l = lambda_thick_l_a < lambda_fluid ?
        lambda_thick_l_a : lambda_fluid;
    lambda_thick_r = lambda_thick_r_a > lambda_fluid ?
        lambda_thick_r_a : lambda_fluid;

    chi = prj_flux_gr_m1_chi(rad, zeta);
    thin_w = 0.5 * (3.0 * chi - 1.0);
    thick_w = 1.5 * (1.0 - chi);
    lambda_l = thin_w * lambda_thin_l + thick_w * lambda_thick_l;
    lambda_r = thin_w * lambda_thin_r + thick_w * lambda_thick_r;

    if (!isfinite(lambda_l) || !isfinite(lambda_r) || lambda_l >= lambda_r) {
        *smin = -PRJ_CLIGHT;
        *smax = PRJ_CLIGHT;
        return;
    }
    *smin = PRJ_CLIGHT * lambda_l;
    *smax = PRJ_CLIGHT * lambda_r;
}

static int prj_flux_gr_m1_cell_in_storage(int i, int j, int k)
{
    return i >= -PRJ_NGHOST && i < PRJ_BLOCK_SIZE + PRJ_NGHOST &&
        j >= -PRJ_NGHOST && j < PRJ_BLOCK_SIZE + PRJ_NGHOST &&
        k >= -PRJ_NGHOST && k < PRJ_BLOCK_SIZE + PRJ_NGHOST;
}

/* Central-difference d(v_comp/c)/dx^gt at cell (ci,cj,ck), scaled by the
 * caller-supplied inv2dx = 1/(2 dx^gt c). Returns 0 (leaving *out untouched) if
 * the 3-point stencil leaves storage or hits a non-finite velocity. */
static int prj_flux_gr_m1_cell_transverse_dv(const double *W_block, int comp,
    int gt, int ci, int cj, int ck, double inv2dx, double *out)
{
    int oi = gt == X1DIR ? 1 : 0;
    int oj = gt == X2DIR ? 1 : 0;
    int ok = gt == X3DIR ? 1 : 0;
    double vp;
    double vm;

    if (!prj_flux_gr_m1_cell_in_storage(ci + oi, cj + oj, ck + ok) ||
        !prj_flux_gr_m1_cell_in_storage(ci - oi, cj - oj, ck - ok)) {
        return 0;
    }
    vp = W_block[WIDX(PRJ_PRIM_V1 + comp, ci + oi, cj + oj, ck + ok)];
    vm = W_block[WIDX(PRJ_PRIM_V1 + comp, ci - oi, cj - oj, ck - ok)];
    if (!isfinite(vp) || !isfinite(vm)) {
        return 0;
    }
    *out = (vp - vm) * inv2dx;
    return 1;
}

/* Face-centred velocity-gradient tensor dvdx[d][a] = d(v_{axis[a]}/c)/dx^{axis[d]}
 * in the local flux axes (axis[0] = normal `dir`). The normal row (d = 0) is the
 * two-cell difference across the face; the transverse rows (d = 1,2) average the
 * cell-centred central differences of the two cells sharing the face, each side
 * contributing only where its stencil stays in storage. This mirrors the
 * coordinate-basis full gradient the source-term closure uses
 * (prj_src_gr_m1_dvdx), rather than leaving the tangential shear rows zero. */
static int prj_flux_gr_m1_face_dvdx(const double *W_block, int dir,
    int il, int jl, int kl, int ir, int jr, int kr, const double *dx,
    double dvdx[3][3])
{
    int axis[3];
    int a;
    int d;
    double dx_dir;

    for (d = 0; d < 3; ++d) {
        for (a = 0; a < 3; ++a) {
            dvdx[d][a] = 0.0;
        }
    }
    if (W_block == 0 || dx == 0) {
        return 0;
    }
    dx_dir = dx[dir];
    if (!isfinite(dx_dir) || dx_dir <= 0.0) {
        return 0;
    }

    prj_flux_local_axes(dir, axis);
    for (a = 0; a < 3; ++a) {
        int v = PRJ_PRIM_V1 + axis[a];
        double vL = W_block[WIDX(v, il, jl, kl)];
        double vR = W_block[WIDX(v, ir, jr, kr)];

        if (!isfinite(vL) || !isfinite(vR)) {
            for (d = 0; d < 3; ++d) {
                int q;

                for (q = 0; q < 3; ++q) {
                    dvdx[d][q] = 0.0;
                }
            }
            return 0;
        }
        dvdx[0][a] = (vR - vL) / (dx_dir * PRJ_CLIGHT);
    }

    for (d = 1; d < 3; ++d) {
        int gt = axis[d];
        double dxt = dx[gt];
        double inv2dx;

        if (!isfinite(dxt) || dxt <= 0.0) {
            continue;
        }
        inv2dx = 1.0 / (2.0 * dxt * PRJ_CLIGHT);
        for (a = 0; a < 3; ++a) {
            int comp = axis[a];
            double acc = 0.0;
            double contrib;
            int count = 0;

            if (prj_flux_gr_m1_cell_transverse_dv(W_block, comp, gt,
                    il, jl, kl, inv2dx, &contrib)) {
                acc += contrib;
                ++count;
            }
            if (prj_flux_gr_m1_cell_transverse_dv(W_block, comp, gt,
                    ir, jr, kr, inv2dx, &contrib)) {
                acc += contrib;
                ++count;
            }
            if (count > 0) {
                dvdx[d][a] = acc / (double)count;
            }
        }
    }
    return 1;
}

static void prj_flux_gr_m1_closure_ctx(const prj_z4c_hydro_geom *geom,
    const double *W, double opacity, const double dvdx[3][3], int have_dvdx,
    prj_rad_gr_m1_closure_ctx *ctx)
{
    int a;
    int b;
    int d;

    memset(ctx, 0, sizeof(*ctx));
    for (a = 0; a < 3; ++a) {
        ctx->vcon[a] = W != 0 ? W[PRJ_PRIM_V1 + a] / PRJ_CLIGHT : 0.0;
        for (b = 0; b < 3; ++b) {
            ctx->gamma[a][b] = geom->gamma[a][b];
            ctx->gamma_inv[a][b] = geom->gamma_inv[a][b];
            ctx->K_dd[a][b] = geom->K_dd[a][b];
            for (d = 0; d < 3; ++d) {
                ctx->dgamma[d][a][b] = geom->dgamma[d][a][b];
            }
        }
    }
    if (have_dvdx && dvdx != 0) {
        for (d = 0; d < 3; ++d) {
            for (a = 0; a < 3; ++a) {
                ctx->dvdx[d][a] = dvdx[d][a];
            }
        }
    }
    ctx->opacity = opacity;
    /* Face transport carries the full velocity-gradient tensor: normal row from
     * the across-face difference and tangential rows from the averaged
     * cell-centred transverse stencils (prj_flux_gr_m1_face_dvdx). */
    ctx->have_shear = have_dvdx;
}

/* Per-group state flux for one side of a face. The closure context (metric,
 * dgamma, dvdx, vcon) is invariant across all (field, group) pairs of this
 * side, so it is built once by the caller and passed in with only ctx->opacity
 * refreshed per group; this hoists a memset + ~54-double geometry copy out of
 * the NRAD*NEGROUP inner loop (bit-identical result). */
static void prj_flux_gr_m1_state_flux(const prj_rad *rad,
    const prj_z4c_hydro_geom *geom, const double *W, int field, int group,
    const prj_rad_gr_m1_closure_ctx *closure,
    const prj_rad_gr_m1_side_data *side, double U[4], double F[4],
    double *smin, double *smax)
{
    const double c = PRJ_CLIGHT;
    const double c2 = PRJ_CLIGHT * PRJ_CLIGHT;
    double E;
    double Fcov[3];
    double Fcon[3];
    double Pcon[3][3];
    double Fmag;
    double fbar = 0.0;
    int i;

    E = W[PRJ_PRIM_RAD_E(field, group)];
    Fcov[0] = W[PRJ_PRIM_RAD_F1(field, group)];
    Fcov[1] = W[PRJ_PRIM_RAD_F2(field, group)];
    Fcov[2] = W[PRJ_PRIM_RAD_F3(field, group)];
    for (i = 0; i < 3; ++i) {
        if (!isfinite(Fcov[i])) {
            Fcov[i] = 0.0;
        }
    }
    prj_flux_gr_m1_limit_state(geom, &E, Fcov, Fcon, &Fmag);

    prj_rad_gr_m1_pressure_fbar_cached(rad, closure, side, E, Fcov, Pcon,
        &fbar);

    U[0] = geom->sqrt_gamma * E;
    F[0] = geom->sqrt_gamma * (geom->alpha * Fcon[0] -
        c * geom->beta[0] * E);
    for (i = 0; i < 3; ++i) {
        double Pn_i = 0.0;
        int a;

        U[1 + i] = geom->sqrt_gamma * Fcov[i];
        for (a = 0; a < 3; ++a) {
            Pn_i += geom->gamma[i][a] * Pcon[0][a];
        }
        F[1 + i] = geom->sqrt_gamma *
            (geom->alpha * c2 * Pn_i - c * geom->beta[0] * Fcov[i]);
    }

    prj_flux_gr_m1_wavespeeds(rad, geom, side, Fcon, Fmag, fbar, smin, smax);
}

static void prj_flux_gr_m1(const prj_rad *rad, const double *WL,
    const double *WR, const double *W_block, int dir,
    int il, int jl, int kl, int ir, int jr, int kr,
    const prj_z4c_hydro_geom *geom, const double *chi_face, const double *dx,
    double *flux)
{
    double dvdx_face[3][3];
    double dx_dir = dx != 0 ? dx[dir] : 0.0;
    int have_dvdx;
    int field;
    int group;
    /* Closure context is invariant across (field, group) for each side of the
     * face; build it once here (opacity refreshed per group below) instead of
     * rebuilding it in every prj_flux_gr_m1_state_flux call. */
    prj_rad_gr_m1_closure_ctx closureL;
    prj_rad_gr_m1_closure_ctx closureR;
    /* Per-side velocity/geometry kinematics (shear, Christoffel, ...) are
     * invariant across all NRAD*NEGROUP groups, so build them once per side
     * here and reuse in every prj_rad_gr_m1_pressure_fbar_cached call below. */
    prj_rad_gr_m1_side_data sideL;
    prj_rad_gr_m1_side_data sideR;

    if (rad == 0 || WL == 0 || WR == 0 || geom == 0 || flux == 0) {
        return;
    }
    have_dvdx = prj_flux_gr_m1_face_dvdx(W_block, dir, il, jl, kl, ir, jr, kr,
        dx, dvdx_face);
    prj_flux_gr_m1_closure_ctx(geom, WL, 0.0, dvdx_face, have_dvdx, &closureL);
    prj_flux_gr_m1_closure_ctx(geom, WR, 0.0, dvdx_face, have_dvdx, &closureR);
    prj_rad_gr_m1_prepare_side(&closureL, &sideL);
    prj_rad_gr_m1_prepare_side(&closureR, &sideR);
    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            int idx = field * PRJ_NEGROUP + group;
            double UL[4];
            double UR[4];
            double FphysL[4];
            double FphysR[4];
            double lamL_min;
            double lamL_max;
            double lamR_min;
            double lamR_max;
            double sL;
            double sR;
            double denom;
            double inv_denom;
            double chi_ext;
            double tau;
            double eps;
            double eps2;
            int q;

            chi_ext = chi_face != 0 ? chi_face[idx] : 0.0;
            if (!isfinite(chi_ext) || chi_ext < 0.0) {
                chi_ext = 0.0;
            }
            closureL.opacity = chi_ext;
            closureR.opacity = chi_ext;

            prj_flux_gr_m1_state_flux(rad, geom, WL, field, group, &closureL,
                &sideL, UL, FphysL, &lamL_min, &lamL_max);
            prj_flux_gr_m1_state_flux(rad, geom, WR, field, group, &closureR,
                &sideR, UR, FphysR, &lamR_min, &lamR_max);

            sL = lamL_min < lamR_min ? lamL_min : lamR_min;
            sR = lamL_max > lamR_max ? lamL_max : lamR_max;
            if (sL > 0.0) {
                sL = 0.0;
            }
            if (sR < 0.0) {
                sR = 0.0;
            }
            if (sR - sL < 1.0e-30) {
                sL = -PRJ_CLIGHT;
                sR = PRJ_CLIGHT;
            }
            denom = sR - sL;
            inv_denom = 1.0 / denom;

            tau = chi_ext * dx_dir;
            eps = 3.0 / (5.0 * tau + 1.0e-10);
            if (eps > 1.0) {
                eps = 1.0;
            }
            eps2 = eps * eps;

            flux[PRJ_CONS_RAD_E(field, group)] =
                (sR * FphysL[0] - sL * FphysR[0] +
                    eps * sL * sR * (UR[0] - UL[0])) * inv_denom;
            for (q = 0; q < 3; ++q) {
                int v = PRJ_CONS_RAD_F1(field, group) + q;

                flux[v] =
                    (eps2 * (sR * FphysL[1 + q] - sL * FphysR[1 + q]) +
                        eps * sL * sR * (UR[1 + q] - UL[1 + q])) *
                    inv_denom +
                    (1.0 - eps2) * 0.5 * (FphysL[1 + q] + FphysR[1 + q]);
            }
        }
    }
}
#endif
#endif

#if !PRJ_MHD
static void prj_flux_velocity_deltas(double *W, int dir,
    int il, int jl, int kl, int ir, int jr, int kr,
    double *deltau, double *deltav, double *deltaw)
{
    if (dir == X1DIR) {
        *deltau = W[WIDX(PRJ_PRIM_V1, ir, jr, kr)] -
            W[WIDX(PRJ_PRIM_V1, il, jl, kl)];
        *deltav = PRJ_MIN(
            PRJ_MIN(W[WIDX(PRJ_PRIM_V2, il, jl, kl)] -
                    W[WIDX(PRJ_PRIM_V2, il, jl - 1, kl)],
                    W[WIDX(PRJ_PRIM_V2, il, jl + 1, kl)] -
                    W[WIDX(PRJ_PRIM_V2, il, jl, kl)]),
            PRJ_MIN(W[WIDX(PRJ_PRIM_V2, ir, jr, kr)] -
                    W[WIDX(PRJ_PRIM_V2, ir, jr - 1, kr)],
                    W[WIDX(PRJ_PRIM_V2, ir, jr + 1, kr)] -
                    W[WIDX(PRJ_PRIM_V2, ir, jr, kr)]));
        *deltaw = PRJ_MIN(
            PRJ_MIN(W[WIDX(PRJ_PRIM_V3, il, jl, kl)] -
                    W[WIDX(PRJ_PRIM_V3, il, jl, kl - 1)],
                    W[WIDX(PRJ_PRIM_V3, il, jl, kl + 1)] -
                    W[WIDX(PRJ_PRIM_V3, il, jl, kl)]),
            PRJ_MIN(W[WIDX(PRJ_PRIM_V3, ir, jr, kr)] -
                    W[WIDX(PRJ_PRIM_V3, ir, jr, kr - 1)],
                    W[WIDX(PRJ_PRIM_V3, ir, jr, kr + 1)] -
                    W[WIDX(PRJ_PRIM_V3, ir, jr, kr)]));
    } else if (dir == X2DIR) {
        *deltau = W[WIDX(PRJ_PRIM_V2, ir, jr, kr)] -
            W[WIDX(PRJ_PRIM_V2, il, jl, kl)];
        *deltav = PRJ_MIN(
            PRJ_MIN(W[WIDX(PRJ_PRIM_V3, il, jl, kl)] -
                    W[WIDX(PRJ_PRIM_V3, il, jl, kl - 1)],
                    W[WIDX(PRJ_PRIM_V3, il, jl, kl + 1)] -
                    W[WIDX(PRJ_PRIM_V3, il, jl, kl)]),
            PRJ_MIN(W[WIDX(PRJ_PRIM_V3, ir, jr, kr)] -
                    W[WIDX(PRJ_PRIM_V3, ir, jr, kr - 1)],
                    W[WIDX(PRJ_PRIM_V3, ir, jr, kr + 1)] -
                    W[WIDX(PRJ_PRIM_V3, ir, jr, kr)]));
        *deltaw = PRJ_MIN(
            PRJ_MIN(W[WIDX(PRJ_PRIM_V1, il, jl, kl)] -
                    W[WIDX(PRJ_PRIM_V1, il - 1, jl, kl)],
                    W[WIDX(PRJ_PRIM_V1, il + 1, jl, kl)] -
                    W[WIDX(PRJ_PRIM_V1, il, jl, kl)]),
            PRJ_MIN(W[WIDX(PRJ_PRIM_V1, ir, jr, kr)] -
                    W[WIDX(PRJ_PRIM_V1, ir - 1, jr, kr)],
                    W[WIDX(PRJ_PRIM_V1, ir + 1, jr, kr)] -
                    W[WIDX(PRJ_PRIM_V1, ir, jr, kr)]));
    } else {
        *deltau = W[WIDX(PRJ_PRIM_V3, ir, jr, kr)] -
            W[WIDX(PRJ_PRIM_V3, il, jl, kl)];
        *deltav = PRJ_MIN(
            PRJ_MIN(W[WIDX(PRJ_PRIM_V1, il, jl, kl)] -
                    W[WIDX(PRJ_PRIM_V1, il - 1, jl, kl)],
                    W[WIDX(PRJ_PRIM_V1, il + 1, jl, kl)] -
                    W[WIDX(PRJ_PRIM_V1, il, jl, kl)]),
            PRJ_MIN(W[WIDX(PRJ_PRIM_V1, ir, jr, kr)] -
                    W[WIDX(PRJ_PRIM_V1, ir - 1, jr, kr)],
                    W[WIDX(PRJ_PRIM_V1, ir + 1, jr, kr)] -
                    W[WIDX(PRJ_PRIM_V1, ir, jr, kr)]));
        *deltaw = PRJ_MIN(
            PRJ_MIN(W[WIDX(PRJ_PRIM_V2, il, jl, kl)] -
                    W[WIDX(PRJ_PRIM_V2, il, jl - 1, kl)],
                    W[WIDX(PRJ_PRIM_V2, il, jl + 1, kl)] -
                    W[WIDX(PRJ_PRIM_V2, il, jl, kl)]),
            PRJ_MIN(W[WIDX(PRJ_PRIM_V2, ir, jr, kr)] -
                    W[WIDX(PRJ_PRIM_V2, ir, jr - 1, kr)],
                    W[WIDX(PRJ_PRIM_V2, ir, jr + 1, kr)] -
                    W[WIDX(PRJ_PRIM_V2, ir, jr, kr)]));
    }
}
#endif

#if PRJ_NRAD > 0
/* Transport absorption/scattering opacity for one cell, stored per group in
 * block->kappa_cell / block->sigma_cell. Used by the radiation flux. */
static void prj_flux_opacity_cell(prj_rad *rad, prj_block *block, const double *W,
    const double *eosvar, int ii, int jj, int kk)
{
    const size_t stride = (size_t)PRJ_NRAD * (size_t)PRJ_NEGROUP;
    double rho = W[WIDX(PRJ_PRIM_RHO, ii, jj, kk)];
    double ye = W[WIDX(PRJ_PRIM_YE, ii, jj, kk)];
    double T = eosvar[EIDX(PRJ_EOSVAR_TEMPERATURE, ii, jj, kk)];
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
        double *W = prj_block_prim_stage(block, prj_stage_slot_from_stage_arg(stage));
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
        double *W = prj_block_prim_stage(block, prj_stage_slot_from_stage_arg(stage));
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

void prj_flux_update(prj_eos *eos, prj_rad *rad, const prj_mesh *mesh,
    prj_block *block, double *W, double *eosvar, double *flux[3], int use_bf1)
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
    /* Per-pencil (fixed i, j) face-major transpose of the variable-major block
     * buffers above. Filled once per (i, j) strip as contiguous k-streams, then
     * consumed face-by-face zero-copy: WLp[kk] holds all PRJ_NVAR_PRIM vars of
     * face kk = k - kstart contiguously, the layout the Riemann/rad solvers
     * want. Strip length kend - kstart + 1 <= PRJ_BS. */
    static double WLp[PRJ_BS * PRJ_NVAR_PRIM];
    static double WRp[PRJ_BS * PRJ_NVAR_PRIM];
    static double pLp[PRJ_BS];
    static double pRp[PRJ_BS];
    static double gLp[PRJ_BS];
    static double gRp[PRJ_BS];
#if PRJ_DYNAMIC_GR
    static prj_z4c_hydro_geom gr_cell_geom[PRJ_BLOCK_NCELLS];
    int full_dynamic_gr = prj_eos_full_dynamic_gr_enabled(mesh);
    int z4c_stage = prj_stage_slot_from_bf_arg(use_bf1);
    int gr_cell_geom_ready = 0;
#endif
    int dir;
#if PRJ_NRAD == 0
    (void)rad;
#endif
#if !PRJ_DYNAMIC_GR
    (void)mesh;
#endif
#if !PRJ_MHD && !PRJ_DYNAMIC_GR
    (void)use_bf1;
#endif

    /* Transport opacity (kappa_cell/sigma_cell) is now filled ahead of the flux:
     * active cells in the same-level ghost shadow and the 1-ghost halo after
     * eos_fill_mesh. See prj_flux_fill_transport_opacity_active/halo. */

#if PRJ_DYNAMIC_GR
    if (full_dynamic_gr) {
        int gi;
        int gj;
        int gk;

        for (gi = -PRJ_NGHOST; gi < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++gi) {
            for (gj = -PRJ_NGHOST; gj < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++gj) {
                for (gk = -PRJ_NGHOST; gk < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++gk) {
                    if (!prj_z4c_load_hydro_geom(mesh, block, z4c_stage,
                            gi, gj, gk, &gr_cell_geom[LIDX(gi, gj, gk)])) {
                        prj_flux_gr_fail("geometry cache load", -1, 0, gi, gj, gk);
                    }
                }
            }
        }
        gr_cell_geom_ready = 1;
    }
#endif

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
        PRJ_SUBTIMER_START("sub_flux_reconstruct");
        prj_flux_reconstruct_block(W, eosvar, dir,
            istart, iend, jstart, jend, kstart, kend,
            WL_block, WR_block, pL_block, pR_block, gL_block, gR_block);
        PRJ_SUBTIMER_STOP("sub_flux_reconstruct");

        for (i = istart; i <= iend; ++i) {
            for (j = jstart; j <= jend; ++j) {
                int strip = kend - kstart + 1;

                /* Transpose this (i, j) pencil from variable-major block buffers
                 * into face-major strip buffers. The inner k-run is contiguous
                 * in the block buffers (k is the fastest IDX axis), so each
                 * variable streams sequentially instead of striding 13.8 KB per
                 * face as the old per-face gather did. */
                PRJ_SUBTIMER_START("sub_flux_gather");
                {
                    size_t base = (size_t)LIDX(i, j, kstart);
                    int v;
                    int kk;
                    for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
                        const double *src_l =
                            &WL_block[VPLANE((size_t)v) + base];
                        const double *src_r =
                            &WR_block[VPLANE((size_t)v) + base];
                        for (kk = 0; kk < strip; ++kk) {
                            WLp[(size_t)kk * PRJ_NVAR_PRIM + v] = src_l[kk];
                            WRp[(size_t)kk * PRJ_NVAR_PRIM + v] = src_r[kk];
                        }
                    }
                    for (kk = 0; kk < strip; ++kk) {
                        pLp[kk] = pL_block[base + kk];
                        pRp[kk] = pR_block[base + kk];
                        gLp[kk] = gL_block[base + kk];
                        gRp[kk] = gR_block[base + kk];
                    }
                }
                PRJ_SUBTIMER_STOP("sub_flux_gather");

                for (k = kstart; k <= kend; ++k) {
                    int kk = k - kstart;
                    double *WL = &WLp[(size_t)kk * PRJ_NVAR_PRIM];
                    double *WR = &WRp[(size_t)kk * PRJ_NVAR_PRIM];
                    double pL = pLp[kk];
                    double pR = pRp[kk];
                    double gL = gLp[kk];
                    double gR = gRp[kk];
                    double Fl[PRJ_NVAR_CONS];
                    int il;
                    int jl;
                    int kl;
                    int ir;
                    int jr;
                    int kr;

                    double v_face_loc[3] = {0.0, 0.0, 0.0};
                    double deltau;
                    double deltav;
                    double deltaw;

                    prj_flux_face_cells(dir, i, j, k, &il, &jl, &kl, &ir, &jr, &kr);

#if PRJ_MHD
                    deltau = 0.0;
                    deltav = 0.0;
                    deltaw = 0.0;
#else
                    PRJ_SUBTIMER_START("sub_flux_veldelta");
                    prj_flux_velocity_deltas(W, dir, il, jl, kl, ir, jr, kr,
                        &deltau, &deltav, &deltaw);
                    PRJ_SUBTIMER_STOP("sub_flux_veldelta");
#endif
                    PRJ_SUBTIMER_START("sub_flux_riemann");
#if PRJ_MHD
                    {
                        double *bf_dir = prj_block_bf_stage(block, dir,
                            prj_stage_slot_from_bf_arg(use_bf1));
                        double bv1 = 0.0;
                        double bv2 = 0.0;
                        double bn;

                        bn = bf_dir[FACE_IDX(dir, i, j, k)];
#if PRJ_DYNAMIC_GR
                        if (full_dynamic_gr) {
                            prj_z4c_hydro_geom geom;

                            if (!gr_cell_geom_ready ||
                                !prj_flux_gr_face_geom_cached(gr_cell_geom, dir,
                                    il, jl, kl, ir, jr, kr, &geom)) {
                                prj_flux_gr_fail("geometry load", -1, dir, i, j, k);
                            }
                            prj_riemann_gr_hlld(WL, WR, pL, pR, gL, gR,
                                eos, (const double (*)[3])geom.gamma,
                                geom.sqrt_gamma, geom.alpha,
                                geom.beta, bn, Fl, v_face_loc, &bv1, &bv2,
                                deltau, deltav, deltaw);
                        } else
#endif
                        {
                            WL[PRJ_PRIM_B1] = bn;
                            WR[PRJ_PRIM_B1] = bn;
                            prj_riemann_hlld(WL, WR, pL, pR, gL, gR, eos, bn, Fl,
                                v_face_loc, &bv1, &bv2, deltau, deltav, deltaw);
                        }
                        block->Bv1[dir][IDX(i, j, k)] = bv1;
                        block->Bv2[dir][IDX(i, j, k)] = bv2;
                    }
#else
#if PRJ_DYNAMIC_GR
                    if (prj_eos_full_dynamic_gr_enabled(mesh)) {
                        prj_flux_gr_hydro_hll(eos, mesh, block,
                            W, eosvar, prj_stage_slot_from_bf_arg(use_bf1),
                            dir, i, j, k, il, jl, kl, ir, jr, kr, WL, WR,
                            pL, pR, gL, gR, Fl, v_face_loc);
                    } else
#endif
                    prj_riemann_hllc(WL, WR, pL, pR, gL, gR, eos, Fl,
                        v_face_loc, deltau, deltav, deltaw);
#endif
                    PRJ_SUBTIMER_STOP("sub_flux_riemann");
                    prj_flux_store_face_velocity(block, dir, i, j, k, v_face_loc);
#if PRJ_NRAD > 0
                    PRJ_SUBTIMER_START("sub_flux_rad");
                    if (prj_flux_rad_face_needed(dir, i, j, k)) {
#if PRJ_USE_RADIATION_FSA
                        double lapse_face = 0.5 *
                            (block->lapse[IDX(il, jl, kl)] + block->lapse[IDX(ir, jr, kr)]);

                        prj_rad_flux_fsa(rad, block, WL, WR, lapse_face, dir,
                            v_face_loc[0], il, jl, kl, ir, jr, kr, Fl);
#else
                        double dx_dir;
                        double chi_face[PRJ_NRAD * PRJ_NEGROUP];
                        const size_t stride = (size_t)PRJ_NRAD * (size_t)PRJ_NEGROUP;
                        size_t off_L = (size_t)IDX(il, jl, kl) * stride;
                        size_t off_R = (size_t)IDX(ir, jr, kr) * stride;
                        const double *kappa_L = &block->kappa_cell[off_L];
                        const double *sigma_L = &block->sigma_cell[off_L];
                        const double *kappa_R = &block->kappa_cell[off_R];
                        const double *sigma_R = &block->sigma_cell[off_R];
                        int idx;
                        dx_dir = block->dx[dir];

                        for (idx = 0; idx < PRJ_NRAD * PRJ_NEGROUP; ++idx) {
                            double kL = kappa_L[idx];
                            double kR = kappa_R[idx];
                            double sL_o = sigma_L[idx];
                            double sR_o = sigma_R[idx];
                            double k_sum = kL + kR;
                            double s_sum = sL_o + sR_o;
                            double k_face = k_sum > 0.0 ? 2.0 * kL * kR / k_sum : 0.0;
                            double s_face = s_sum > 0.0 ? 2.0 * sL_o * sR_o / s_sum : 0.0;
                            chi_face[idx] = k_face + s_face;
                        }

#if PRJ_DYNAMIC_GR
                        if (full_dynamic_gr) {
                            prj_z4c_hydro_geom geom;

                            if (!gr_cell_geom_ready ||
                                !prj_flux_gr_face_geom_cached(gr_cell_geom, dir,
                                    il, jl, kl, ir, jr, kr, &geom)) {
                                prj_flux_gr_fail("radiation geometry load", -1,
                                    dir, i, j, k);
                            }
                            prj_flux_gr_m1(rad, WL, WR, W, dir,
                                il, jl, kl, ir, jr, kr, &geom, chi_face,
                                block->dx, Fl);
                        } else
#endif
                        {
                            double lapse_face = 0.5 *
                                (block->lapse[IDX(il, jl, kl)] +
                                    block->lapse[IDX(ir, jr, kr)]);

                            prj_rad_flux(rad, WL, WR, lapse_face, chi_face,
                                dx_dir, v_face_loc[0], Fl);
                        }
#endif
                    }
                    /* Transverse ghost faces (MHD CT/EMF layer) are skipped: the
                     * radiation flux divergence never reads them, and the store
                     * (prj_flux_store_local_flux) gates the radiation copy on the
                     * same prj_flux_rad_face_needed() predicate. */
                    PRJ_SUBTIMER_STOP("sub_flux_rad");
#endif
                    PRJ_SUBTIMER_START("sub_flux_store");
                    prj_flux_store_local_flux(flux[dir], dir, i, j, k, Fl);
                    PRJ_SUBTIMER_STOP("sub_flux_store");
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
