#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "prj.h"
#include "prj_rad_inel.h"

#if PRJ_NRAD > 0
static double prj_timeint_cell_lapse(const prj_block *block, int i, int j, int k)
{
    if (block != 0 && block->lapse != 0) {
        return block->lapse[IDX(i, j, k)];
    }
    return 1.0;
}
#endif

static void prj_timeint_unavailable(const char *name)
{
    fprintf(stderr, "%s is unavailable for the selected TIME_INTEGRATION\n", name);
    exit(EXIT_FAILURE);
}

static int prj_timeint_local_block(const prj_mpi *mpi, const prj_block *block)
{
    return block != 0 && block->id >= 0 && block->active == 1 &&
        (mpi == 0 || block->rank == mpi->rank);
}

static void prj_timeint_cell_prim(const double *src, int i, int j, int k, double *w)
{
    int v;

    for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
        w[v] = src[WIDX(v, i, j, k)];
    }
}

static void prj_timeint_cell_prim_store(double *dst, int i, int j, int k, const double *w)
{
    int v;

    for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
        dst[WIDX(v, i, j, k)] = w[v];
    }
}

#if PRJ_MHD
static void prj_timeint_mhd_fail(const char *message)
{
    fprintf(stderr, "%s\n", message);
    exit(EXIT_FAILURE);
}


static int prj_timeint_mhd_face_axis_max(int dir, int axis)
{
    return dir == axis ? PRJ_BLOCK_SIZE : PRJ_BLOCK_SIZE - 1;
}

static int prj_timeint_mhd_edge_axis_max(int dir, int axis)
{
    return dir == axis ? PRJ_BLOCK_SIZE - 1 : PRJ_BLOCK_SIZE;
}

static double prj_timeint_mhd_cell_emf(const prj_mesh *mesh, const prj_block *block,
    int stage, const double *W, int dir, int i, int j, int k)
{
    double b1;
    double b2;
    double b3;
    double v1;
    double v2;
    double v3;
    double emf;

    if (W == 0 || dir < 0 || dir >= 3) {
        prj_timeint_mhd_fail("prj_timeint_mhd_cell_emf: invalid input");
    }
    b1 = W[WIDX(PRJ_PRIM_B1, i, j, k)];
    b2 = W[WIDX(PRJ_PRIM_B2, i, j, k)];
    b3 = W[WIDX(PRJ_PRIM_B3, i, j, k)];
    v1 = W[WIDX(PRJ_PRIM_V1, i, j, k)];
    v2 = W[WIDX(PRJ_PRIM_V2, i, j, k)];
    v3 = W[WIDX(PRJ_PRIM_V3, i, j, k)];
#if PRJ_DYNAMIC_GR
    if (prj_eos_full_dynamic_gr_enabled(mesh)) {
        prj_z4c_hydro_geom geom;

        if (!prj_z4c_load_hydro_geom(mesh, block, stage, i, j, k, &geom)) {
            prj_timeint_mhd_fail("prj_timeint_mhd_cell_emf: failed to load full-GR geometry");
        }
        b1 *= geom.sqrt_gamma;
        b2 *= geom.sqrt_gamma;
        b3 *= geom.sqrt_gamma;
        v1 = geom.alpha * v1 - PRJ_CLIGHT * geom.beta[0];
        v2 = geom.alpha * v2 - PRJ_CLIGHT * geom.beta[1];
        v3 = geom.alpha * v3 - PRJ_CLIGHT * geom.beta[2];
    }
#else
    (void)mesh;
    (void)block;
    (void)stage;
#endif
    if (dir == X1DIR) {
        emf = b2 * v3 - b3 * v2;
    } else if (dir == X2DIR) {
        emf = b3 * v1 - b1 * v3;
    } else {
        emf = b1 * v2 - b2 * v1;
    }
    if (!isfinite(emf)) {
        prj_timeint_mhd_fail("prj_timeint_mhd_cell_emf: non-finite cell emf");
    }
    return emf;
}

static double prj_timeint_mhd_face_emf(const prj_block *block, int face_dir, int emf_dir,
    int i, int j, int k)
{
    double value = 0.0;

    if (block == 0 || face_dir < 0 || face_dir >= 3 || emf_dir < 0 || emf_dir >= 3 ||
        block->Bv1[face_dir] == 0 || block->Bv2[face_dir] == 0) {
        prj_timeint_mhd_fail("prj_timeint_mhd_face_emf: invalid input");
    }
    if (emf_dir == (face_dir + 1) % 3) {
        value = block->Bv1[face_dir][IDX(i, j, k)];
    } else if (emf_dir == (face_dir + 2) % 3) {
        value = block->Bv2[face_dir][IDX(i, j, k)];
    } else {
        prj_timeint_mhd_fail("prj_timeint_mhd_face_emf: requested emf is normal to face");
    }
    if (!isfinite(value)) {
        prj_timeint_mhd_fail("prj_timeint_mhd_face_emf: non-finite face emf");
    }
    return value;
}

static double prj_timeint_mhd_face_vnorm(const prj_block *block, int face_dir,
    int i, int j, int k)
{
    double value;

    if (block == 0 || face_dir < 0 || face_dir >= 3 || block->v_riemann[face_dir] == 0) {
        prj_timeint_mhd_fail("prj_timeint_mhd_face_vnorm: invalid input");
    }
    value = block->v_riemann[face_dir][VRIDX(face_dir, i, j, k)];
    if (!isfinite(value)) {
        prj_timeint_mhd_fail("prj_timeint_mhd_face_vnorm: non-finite face velocity");
    }
    return value;
}

/* Anomalous CT resistivity applied in dense, proton-rich EOS-transition material.
 * The per-cell artificial eta grows smoothly with both rho/RHO_FLOOR and
 * Ye/YE_THR, and the EMF edge uses the maximum eta from its incident cells.
 *
 * Applying the value on all incident edges (rather than only those that straddle
 * a feature) makes the damping orientation-independent. */
#define PRJ_MHD_ETA_BOOST 1.0e12
#define PRJ_MHD_ETA_RHO_FLOOR 1.0e13
#define PRJ_MHD_ETA_YE_THR 0.25

/* Fill every cell over the active region plus one ghost layer -- the range the
 * EMF edge loop reads.  The current stage primitives already include ghost zones
 * copied from neighbour blocks, so adjacent blocks agree on shared edges with no
 * extra exchange. */
static void prj_timeint_mhd_fill_eta_mask(prj_block *block, const double *W)
{
    int i;
    int j;
    int k;

    if (block == 0 || block->eta_mask == 0 || W == 0) {
        prj_timeint_mhd_fail("prj_timeint_mhd_fill_eta_mask: invalid input");
    }
    for (i = -1; i <= PRJ_BLOCK_SIZE; ++i) {
        for (j = -1; j <= PRJ_BLOCK_SIZE; ++j) {
            for (k = -1; k <= PRJ_BLOCK_SIZE; ++k) {
                double rho = W[WIDX(PRJ_PRIM_RHO, i, j, k)];
                double ye = W[WIDX(PRJ_PRIM_YE, i, j, k)];
                double rho_factor = 0.0;
                double ye_factor = 0.0;
                double eta_art = 0.0;

                if (isfinite(rho) && rho > 0.0) {
                    rho_factor = PRJ_MAX(0.0, log(rho / PRJ_MHD_ETA_RHO_FLOOR));
                }
                if (isfinite(ye) && ye > 0.0) {
                    ye_factor = PRJ_MAX(0.0, log(ye / PRJ_MHD_ETA_YE_THR));
                }
                eta_art = PRJ_MHD_ETA_BOOST * rho_factor * ye_factor;
                if (!isfinite(eta_art)) {
                    prj_timeint_mhd_fail("prj_timeint_mhd_fill_eta_mask: non-finite artificial eta");
                }
                block->eta_mask[IDX(i, j, k)] = eta_art;
            }
        }
    }
}

/* Maximum artificial eta from the four cells incident to an EMF edge. */
static double prj_timeint_mhd_edge_eta_mask(const prj_block *block,
    const int cell0[3], const int cell1[3], const int cell2[3], const int cell3[3])
{
    const int *cells[4];
    int n;
    double eta = 0.0;

    if (block == 0 || block->eta_mask == 0) {
        return 0.0;
    }
    cells[0] = cell0;
    cells[1] = cell1;
    cells[2] = cell2;
    cells[3] = cell3;
    for (n = 0; n < 4; ++n) {
        double cell_eta = block->eta_mask[IDX(cells[n][0], cells[n][1], cells[n][2])];

        if (!isfinite(cell_eta) || cell_eta < 0.0) {
            prj_timeint_mhd_fail("prj_timeint_mhd_edge_eta_mask: invalid artificial eta");
        }
        if (cell_eta > eta) {
            eta = cell_eta;
        }
    }
    return eta;
}

static void prj_timeint_mhd_update_emf(const prj_mesh *mesh, prj_block *block,
    int stage, double eta)
{
    double *W;
    const double *bf[3];
    int dir;

    if (block == 0) {
        prj_timeint_mhd_fail("prj_timeint_mhd_update_emf: invalid input");
    }
    W = prj_block_prim_stage(block, stage);
    if (W == 0) {
        prj_timeint_mhd_fail("prj_timeint_mhd_update_emf: missing prim storage");
    }
    for (dir = 0; dir < 3; ++dir) {
        bf[dir] = prj_block_bf_stage(block, dir, stage);
        if (bf[dir] == 0) {
            prj_timeint_mhd_fail("prj_timeint_mhd_update_emf: missing Bf storage");
        }
    }
    prj_timeint_mhd_fill_eta_mask(block, W);
    for (dir = 0; dir < 3; ++dir) {
        int up = (dir + 1) % 3;
        int right = (dir + 2) % 3;
        int e[3];

        for (e[0] = 0; e[0] <= prj_timeint_mhd_edge_axis_max(dir, 0); ++e[0]) {
            for (e[1] = 0; e[1] <= prj_timeint_mhd_edge_axis_max(dir, 1); ++e[1]) {
                for (e[2] = 0; e[2] <= prj_timeint_mhd_edge_axis_max(dir, 2); ++e[2]) {
                    int cell0[3];
                    int cell1[3];
                    int cell2[3];
                    int cell3[3];
                    int face0[3];
                    int face1[3];
                    int face2[3];
                    int face3[3];
                    double emf_face[4];
                    double emf_cell[4];
                    double v_norm[4];
                    int d;

                    for (d = 0; d < 3; ++d) {
                        cell2[d] = e[d];
                        face0[d] = e[d];
                        face1[d] = e[d];
                        face2[d] = e[d];
                        face3[d] = e[d];
                    }
                    cell2[up] -= 1;
                    cell2[right] -= 1;
                    for (d = 0; d < 3; ++d) {
                        cell1[d] = cell2[d];
                        cell3[d] = cell2[d];
                        cell0[d] = cell2[d];
                    }
                    cell1[right] += 1;
                    cell3[up] += 1;
                    cell0[up] += 1;
                    cell0[right] += 1;

                    face2[up] -= 1;
                    face3[right] -= 1;

                    emf_face[0] = prj_timeint_mhd_face_emf(block, right, dir,
                        face0[0], face0[1], face0[2]);
                    emf_face[1] = prj_timeint_mhd_face_emf(block, up, dir,
                        face1[0], face1[1], face1[2]);
                    emf_face[2] = prj_timeint_mhd_face_emf(block, right, dir,
                        face2[0], face2[1], face2[2]);
                    emf_face[3] = prj_timeint_mhd_face_emf(block, up, dir,
                        face3[0], face3[1], face3[2]);
                    emf_cell[0] = prj_timeint_mhd_cell_emf(mesh, block, stage, W, dir,
                        cell0[0], cell0[1], cell0[2]);
                    emf_cell[1] = prj_timeint_mhd_cell_emf(mesh, block, stage, W, dir,
                        cell1[0], cell1[1], cell1[2]);
                    emf_cell[2] = prj_timeint_mhd_cell_emf(mesh, block, stage, W, dir,
                        cell2[0], cell2[1], cell2[2]);
                    emf_cell[3] = prj_timeint_mhd_cell_emf(mesh, block, stage, W, dir,
                        cell3[0], cell3[1], cell3[2]);
                    v_norm[0] = prj_timeint_mhd_face_vnorm(block, right,
                        face0[0], face0[1], face0[2]);
                    v_norm[1] = prj_timeint_mhd_face_vnorm(block, up,
                        face1[0], face1[1], face1[2]);
                    v_norm[2] = prj_timeint_mhd_face_vnorm(block, right,
                        face2[0], face2[1], face2[2]);
                    v_norm[3] = prj_timeint_mhd_face_vnorm(block, up,
                        face3[0], face3[1], face3[2]);
                    double edge_eta = eta + prj_timeint_mhd_edge_eta_mask(block,
                        cell0, cell1, cell2, cell3);

                    prj_mhd_emf_upwind(block, dir, e[0], e[1], e[2],
                        emf_face, emf_cell, v_norm, bf, edge_eta);
                }
            }
        }
    }
}

static double prj_timeint_mhd_ct_value(const prj_block *block, double *src[3],
    int dir, int i, int j, int k, double dt)
{
    double value;

    if (block == 0 || src == 0 || src[dir] == 0 || block->emf[0] == 0 ||
        block->emf[1] == 0 || block->emf[2] == 0) {
        prj_timeint_mhd_fail("prj_timeint_mhd_ct_value: invalid input");
    }
    value = src[dir][FACE_IDX(dir, i, j, k)];
    if (dir == X1DIR) {
        value += -dt * (block->emf[X3DIR][EDGE_IDX(X3DIR, i, j + 1, k)] -
                        block->emf[X3DIR][EDGE_IDX(X3DIR, i, j, k)]) / block->dx[1];
        value +=  dt * (block->emf[X2DIR][EDGE_IDX(X2DIR, i, j, k + 1)] -
                        block->emf[X2DIR][EDGE_IDX(X2DIR, i, j, k)]) / block->dx[2];
    } else if (dir == X2DIR) {
        value += -dt * (block->emf[X1DIR][EDGE_IDX(X1DIR, i, j, k + 1)] -
                        block->emf[X1DIR][EDGE_IDX(X1DIR, i, j, k)]) / block->dx[2];
        value +=  dt * (block->emf[X3DIR][EDGE_IDX(X3DIR, i + 1, j, k)] -
                        block->emf[X3DIR][EDGE_IDX(X3DIR, i, j, k)]) / block->dx[0];
    } else {
        value += -dt * (block->emf[X2DIR][EDGE_IDX(X2DIR, i + 1, j, k)] -
                        block->emf[X2DIR][EDGE_IDX(X2DIR, i, j, k)]) / block->dx[0];
        value +=  dt * (block->emf[X1DIR][EDGE_IDX(X1DIR, i, j + 1, k)] -
                        block->emf[X1DIR][EDGE_IDX(X1DIR, i, j, k)]) / block->dx[1];
    }
    if (!isfinite(value)) {
        prj_timeint_mhd_fail("prj_timeint_mhd_ct_value: non-finite updated magnetic field");
    }
    return value;
}

#if TIME_INTEGRATION != PRJ_TIMEINT_IMEX
static void prj_timeint_mhd_update_bf(prj_block *block, int stage, double dt)
{
    double *src[3];
    double *dst[3];
    double *old[3];
    int dir;

    if (block == 0 || (stage != 0 && stage != 1 && stage != 2)) {
        prj_timeint_mhd_fail("prj_timeint_mhd_update_bf: invalid input");
    }
    for (dir = 0; dir < 3; ++dir) {
        if (block->Bf[dir] == 0) {
            prj_timeint_mhd_fail("prj_timeint_mhd_update_bf: missing Bf storage");
        }
        old[dir] = prj_block_bf_stage(block, dir, 0);
        if (stage == 0) {
            src[dir] = prj_block_bf_stage(block, dir, 0);
            dst[dir] = prj_block_bf_stage(block, dir, 0);
        } else if (stage == 1) {
            src[dir] = prj_block_bf_stage(block, dir, 0);
            dst[dir] = prj_block_bf_stage(block, dir, 1);
        } else {
            src[dir] = prj_block_bf_stage(block, dir, 1);
            dst[dir] = prj_block_bf_stage(block, dir, 0);
        }
    }
    for (dir = 0; dir < 3; ++dir) {
        int i;
        int j;
        int k;

        for (i = 0; i <= prj_timeint_mhd_face_axis_max(dir, 0); ++i) {
            for (j = 0; j <= prj_timeint_mhd_face_axis_max(dir, 1); ++j) {
                for (k = 0; k <= prj_timeint_mhd_face_axis_max(dir, 2); ++k) {
                    double value = prj_timeint_mhd_ct_value(block, src, dir, i, j, k, dt);

                    if (stage == 2) {
                        value = 0.5 * old[dir][FACE_IDX(dir, i, j, k)] + 0.5 * value;
                    }
                    dst[dir][FACE_IDX(dir, i, j, k)] = value;
                }
            }
        }
    }
}
#endif

static void prj_timeint_mhd_set_cons_b_from_bf(const prj_block *block,
    double *bf[3], int i, int j, int k, double *u)
{
    if (block == 0 || bf == 0 || u == 0 || bf[0] == 0 || bf[1] == 0 || bf[2] == 0) {
        prj_timeint_mhd_fail("prj_timeint_mhd_set_cons_b_from_bf: invalid input");
    }
    u[PRJ_CONS_B1] = 0.5 * (bf[X1DIR][FACE_IDX(X1DIR, i, j, k)] + bf[X1DIR][FACE_IDX(X1DIR, i + 1, j, k)]);
    u[PRJ_CONS_B2] = 0.5 * (bf[X2DIR][FACE_IDX(X2DIR, i, j, k)] + bf[X2DIR][FACE_IDX(X2DIR, i, j + 1, k)]);
    u[PRJ_CONS_B3] = 0.5 * (bf[X3DIR][FACE_IDX(X3DIR, i, j, k)] + bf[X3DIR][FACE_IDX(X3DIR, i, j, k + 1)]);
    if (!isfinite(u[PRJ_CONS_B1]) || !isfinite(u[PRJ_CONS_B2]) || !isfinite(u[PRJ_CONS_B3])) {
        prj_timeint_mhd_fail("prj_timeint_mhd_set_cons_b_from_bf: non-finite cell-centered magnetic field");
    }
}

#if TIME_INTEGRATION != PRJ_TIMEINT_IMEX
static void prj_timeint_mhd_update_mesh_emf(prj_mesh *mesh, const prj_mpi *mpi,
    int stage)
{
    int bidx;

    if (mesh == 0) {
        prj_timeint_mhd_fail("prj_timeint_mhd_update_mesh_emf: invalid input");
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (prj_timeint_local_block(mpi, block)) {
            prj_timeint_mhd_update_emf(mesh, block, stage, mesh->mhd_eta);
        }
    }
    prj_mhd_emf_send(mesh, mpi);
}
#endif
#endif

static void prj_timeint_update_dt_src_values(const prj_mesh *mesh, const prj_grav *grav,
    const prj_block *block, double rho, double mom1, double mom2, double mom3,
    double etot, int i, int j, int k, const prj_mpi *mpi, double *dt_src)
{
    double dt_src_local;

    if (dt_src == 0) {
        return;
    }

    dt_src_local = 0.02 * etot /
        (fabs(prj_block_rhs_value_const(block, PRJ_CONS_ETOT, i, j, k)) + 1.0e-50);
    if (dt_src_local < *dt_src) {
        *dt_src = dt_src_local;
    }
#if !PRJ_GRAV_DEBUG
    (void)mesh;
    (void)grav;
    (void)mpi;
    (void)rho;
    (void)mom1;
    (void)mom2;
    (void)mom3;
#endif
#if PRJ_GRAV_DEBUG
    if (dt_src_local<1e-7) {
        int rank = mpi != 0 ? mpi->rank : 0;
        double inv_rho = rho != 0.0 ? 1.0 / rho : 0.0;
        double x1 = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
        double x2 = block->xmin[1] + ((double)j + 0.5) * block->dx[1];
        double x3 = block->xmin[2] + ((double)k + 0.5) * block->dx[2];
        int gidx = prj_block_cache_index(i, j, k);
        double dx1 = x1 - mesh->x_com[0];
        double dx2 = x2 - mesh->x_com[1];
        double dx3 = x3 - mesh->x_com[2];
        double rcom = block->r_com != 0 ? block->r_com[gidx] :
            sqrt(dx1 * dx1 + dx2 * dx2 + dx3 * dx3);
        double accel = (block->grav[0] != 0 && rcom > 0.0) ?
            (block->grav[0][gidx] * dx1 + block->grav[1][gidx] * dx2 +
             block->grav[2][gidx] * dx3) / rcom : 0.0;

        fprintf(stderr,
            "[grav debug] shortest dt_src rank=%d block=%d cell=(%d,%d,%d) "
            "x=(%.17e,%.17e,%.17e) rho=%.17e etot=%.17e "
            "v=(%.17e,%.17e,%.17e) accel=%.17e dt_src=%.17e\n",
            rank, block->id, i, j, k, x1, x2, x3, rho, etot,
            mom1 * inv_rho, mom2 * inv_rho, mom3 * inv_rho, accel, dt_src_local);
    }
#endif
}

#if !(PRJ_MHD && PRJ_NRAD > 0) && TIME_INTEGRATION != PRJ_TIMEINT_IMEX
static void prj_timeint_update_dt_src(const prj_mesh *mesh, const prj_grav *grav,
    const prj_block *block, const double *u, int i, int j, int k,
    const prj_mpi *mpi, double *dt_src)
{
    if (u == 0) {
        return;
    }
    prj_timeint_update_dt_src_values(mesh, grav, block, u[PRJ_CONS_RHO], u[PRJ_CONS_MOM1],
        u[PRJ_CONS_MOM2], u[PRJ_CONS_MOM3], u[PRJ_CONS_ETOT], i, j, k, mpi, dt_src);
}
#endif

#if PRJ_MHD && PRJ_NRAD > 0 && TIME_INTEGRATION != PRJ_TIMEINT_IMEX
static inline double prj_timeint_flux_div_var(const prj_block *block,
    int v, int i, int j, int k)
{
    return -(block->area[X1DIR] *
            (block->flux[X1DIR][VIDX(v, i + 1, j, k)] -
                block->flux[X1DIR][VIDX(v, i, j, k)]) +
        block->area[X2DIR] *
            (block->flux[X2DIR][VIDX(v, i, j + 1, k)] -
                block->flux[X2DIR][VIDX(v, i, j, k)]) +
        block->area[X3DIR] *
            (block->flux[X3DIR][VIDX(v, i, j, k + 1)] -
                block->flux[X3DIR][VIDX(v, i, j, k)])) / block->vol;
}

static inline void prj_timeint_mhd_hydro_cons_from_prim(const double *W,
    int i, int j, int k, double *rho, double *mom1, double *mom2, double *mom3,
    double *etot, double *ye, double *b1, double *b2, double *b3)
{
    double r = W[WIDX(PRJ_PRIM_RHO, i, j, k)];
    double v1 = W[WIDX(PRJ_PRIM_V1, i, j, k)];
    double v2 = W[WIDX(PRJ_PRIM_V2, i, j, k)];
    double v3 = W[WIDX(PRJ_PRIM_V3, i, j, k)];
    double eint = W[WIDX(PRJ_PRIM_EINT, i, j, k)];
    double cb1 = W[WIDX(PRJ_PRIM_B1, i, j, k)];
    double cb2 = W[WIDX(PRJ_PRIM_B2, i, j, k)];
    double cb3 = W[WIDX(PRJ_PRIM_B3, i, j, k)];

    *rho = r;
    *mom1 = r * v1;
    *mom2 = r * v2;
    *mom3 = r * v3;
    *etot = r * eint + 0.5 * r * (v1 * v1 + v2 * v2 + v3 * v3) +
        0.5 * (cb1 * cb1 + cb2 * cb2 + cb3 * cb3);
    *ye = r * W[WIDX(PRJ_PRIM_YE, i, j, k)];
    *b1 = cb1;
    *b2 = cb2;
    *b3 = cb3;
}

static void prj_timeint_cell_cons_from_prim_mhd_rad(const double *W_mhd, const double *W_rad,
    int i, int j, int k, double *u)
{
    int field;
    int group;
#if PRJ_USE_RADIATION_FSA
    int angle;
#endif

    prj_timeint_mhd_hydro_cons_from_prim(W_mhd, i, j, k,
        &u[PRJ_CONS_RHO], &u[PRJ_CONS_MOM1], &u[PRJ_CONS_MOM2],
        &u[PRJ_CONS_MOM3], &u[PRJ_CONS_ETOT], &u[PRJ_CONS_YE],
        &u[PRJ_CONS_B1], &u[PRJ_CONS_B2], &u[PRJ_CONS_B3]);

    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
#if PRJ_USE_RADIATION_FSA
            for (angle = 0; angle < PRJ_NANGLE; ++angle) {
                int iv = PRJ_CONS_RAD_I(field, group, angle);

                u[iv] = W_rad[WIDX(PRJ_RAD_PRIM_I(field, group, angle), i, j, k)];
            }
#else
            u[PRJ_CONS_RAD_E(field, group)] =
                W_rad[WIDX(PRJ_RAD_PRIM_E(field, group), i, j, k)];
            u[PRJ_CONS_RAD_F1(field, group)] =
                W_rad[WIDX(PRJ_RAD_PRIM_F1(field, group), i, j, k)];
            u[PRJ_CONS_RAD_F2(field, group)] =
                W_rad[WIDX(PRJ_RAD_PRIM_F2(field, group), i, j, k)];
            u[PRJ_CONS_RAD_F3(field, group)] =
                W_rad[WIDX(PRJ_RAD_PRIM_F3(field, group), i, j, k)];
#endif
        }
    }
}

static void prj_timeint_store_mhd_rad_cell(double *W_mhd, double *W_rad,
    int i, int j, int k, double *u)
{
    double rho = u[PRJ_CONS_RHO];
    int field;
    int group;
#if PRJ_USE_RADIATION_FSA
    int angle;
#endif

    if (rho == 0.0) {
        W_mhd[WIDX(PRJ_PRIM_RHO, i, j, k)] = 0.0;
        W_mhd[WIDX(PRJ_PRIM_V1, i, j, k)] = 0.0;
        W_mhd[WIDX(PRJ_PRIM_V2, i, j, k)] = 0.0;
        W_mhd[WIDX(PRJ_PRIM_V3, i, j, k)] = 0.0;
        W_mhd[WIDX(PRJ_PRIM_EINT, i, j, k)] = 0.0;
        W_mhd[WIDX(PRJ_PRIM_YE, i, j, k)] = 0.0;
        W_mhd[WIDX(PRJ_PRIM_B1, i, j, k)] = 0.0;
        W_mhd[WIDX(PRJ_PRIM_B2, i, j, k)] = 0.0;
        W_mhd[WIDX(PRJ_PRIM_B3, i, j, k)] = 0.0;
    } else {
        double v1 = u[PRJ_CONS_MOM1] / rho;
        double v2 = u[PRJ_CONS_MOM2] / rho;
        double v3 = u[PRJ_CONS_MOM3] / rho;
        double kinetic = 0.5 * (v1 * v1 + v2 * v2 + v3 * v3);
        double magnetic_density = 0.5 * (u[PRJ_CONS_B1] * u[PRJ_CONS_B1] +
            u[PRJ_CONS_B2] * u[PRJ_CONS_B2] + u[PRJ_CONS_B3] * u[PRJ_CONS_B3]);
        double magnetic = magnetic_density / rho;
        double eint = u[PRJ_CONS_ETOT] / rho - kinetic - magnetic;

        W_mhd[WIDX(PRJ_PRIM_RHO, i, j, k)] = rho;
        W_mhd[WIDX(PRJ_PRIM_V1, i, j, k)] = v1;
        W_mhd[WIDX(PRJ_PRIM_V2, i, j, k)] = v2;
        W_mhd[WIDX(PRJ_PRIM_V3, i, j, k)] = v3;
        W_mhd[WIDX(PRJ_PRIM_EINT, i, j, k)] = eint;
        W_mhd[WIDX(PRJ_PRIM_YE, i, j, k)] = u[PRJ_CONS_YE] / rho;
        W_mhd[WIDX(PRJ_PRIM_B1, i, j, k)] = u[PRJ_CONS_B1];
        W_mhd[WIDX(PRJ_PRIM_B2, i, j, k)] = u[PRJ_CONS_B2];
        W_mhd[WIDX(PRJ_PRIM_B3, i, j, k)] = u[PRJ_CONS_B3];
    }

    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
#if PRJ_USE_RADIATION_FSA
            for (angle = 0; angle < PRJ_NANGLE; ++angle) {
                int iv = PRJ_CONS_RAD_I(field, group, angle);

                W_rad[WIDX(PRJ_RAD_PRIM_I(field, group, angle), i, j, k)] = u[iv];
            }
#else
            int e = PRJ_CONS_RAD_E(field, group);
            int f1 = PRJ_CONS_RAD_F1(field, group);
            int f2 = PRJ_CONS_RAD_F2(field, group);
            int f3 = PRJ_CONS_RAD_F3(field, group);

            W_rad[WIDX(PRJ_RAD_PRIM_E(field, group), i, j, k)] = u[e];
            W_rad[WIDX(PRJ_RAD_PRIM_F1(field, group), i, j, k)] = u[f1];
            W_rad[WIDX(PRJ_RAD_PRIM_F2(field, group), i, j, k)] = u[f2];
            W_rad[WIDX(PRJ_RAD_PRIM_F3(field, group), i, j, k)] = u[f3];
#endif
        }
    }
}
#endif

#if TIME_INTEGRATION != PRJ_TIMEINT_IMEX
static void prj_timeint_update_cell_stage1_mhd_rad(const prj_mesh *mesh, prj_rad *rad, prj_eos *eos,
    prj_block *block, int i, int j, int k,
    double dt, const prj_grav *grav, const prj_mpi *mpi, double *dt_src, double *Wdst,
    double *Wdst_rad, int use_bf1, double dtau_cm)
{
    (void)Wdst_rad;
    (void)dtau_cm;
#if PRJ_MHD && PRJ_NRAD > 0
    double *bf_dst[3];
    double u[PRJ_NVAR_CONS];
    int d;
    int field;
    int group;

    for (d = 0; d < 3; ++d) {
        bf_dst[d] = prj_block_bf_stage(block, d, use_bf1 != 0 ? 1 : 0);
    }

    PRJ_SUBTIMER_START("sub_cell_cons_from_prim");
    prj_timeint_cell_cons_from_prim_mhd_rad(block->W_mhd, block->W_rad, i, j, k, u);
    PRJ_SUBTIMER_STOP("sub_cell_cons_from_prim");
    PRJ_SUBTIMER_START("sub_cell_dt_src");
    prj_timeint_update_dt_src_values(mesh, grav, block, u[PRJ_CONS_RHO], u[PRJ_CONS_MOM1],
        u[PRJ_CONS_MOM2], u[PRJ_CONS_MOM3], u[PRJ_CONS_ETOT], i, j, k, mpi, dt_src);
    PRJ_SUBTIMER_STOP("sub_cell_dt_src");

#define PRJ_STAGE1_UPDATE(v) do { \
        int vv = (v); \
        u[vv] += dt * (prj_block_rhs_value_const(block, vv, i, j, k) + \
            prj_timeint_flux_div_var(block, vv, i, j, k)); \
    } while (0)

    PRJ_SUBTIMER_START("sub_cell_flux_src_update");
    PRJ_STAGE1_UPDATE(PRJ_CONS_RHO);
    PRJ_STAGE1_UPDATE(PRJ_CONS_MOM1);
    PRJ_STAGE1_UPDATE(PRJ_CONS_MOM2);
    PRJ_STAGE1_UPDATE(PRJ_CONS_MOM3);
    PRJ_STAGE1_UPDATE(PRJ_CONS_ETOT);
    PRJ_STAGE1_UPDATE(PRJ_CONS_YE);

    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
#if PRJ_USE_RADIATION_FSA
            int angle;

            for (angle = 0; angle < PRJ_NANGLE; ++angle) {
                PRJ_STAGE1_UPDATE(PRJ_CONS_RAD_I(field, group, angle));
            }
#else
            PRJ_STAGE1_UPDATE(PRJ_CONS_RAD_E(field, group));
            PRJ_STAGE1_UPDATE(PRJ_CONS_RAD_F1(field, group));
            PRJ_STAGE1_UPDATE(PRJ_CONS_RAD_F2(field, group));
            PRJ_STAGE1_UPDATE(PRJ_CONS_RAD_F3(field, group));
#endif
        }
    }
    PRJ_SUBTIMER_STOP("sub_cell_flux_src_update");

#undef PRJ_STAGE1_UPDATE

    PRJ_SUBTIMER_START("sub_cell_mhd_set_b");
    prj_timeint_mhd_set_cons_b_from_bf(block, bf_dst, i, j, k, u);
    PRJ_SUBTIMER_STOP("sub_cell_mhd_set_b");
    {
        double T_cell = block->eosvar[EIDX(PRJ_EOSVAR_TEMPERATURE, i, j, k)];
        double lapse_cell = prj_timeint_cell_lapse(block, i, j, k);
#if PRJ_USE_RADIATION_M1
        double kappa[PRJ_NRAD * PRJ_NEGROUP];
#endif

        PRJ_SUBTIMER_START("sub_cell_radiation");
#if PRJ_USE_RADIATION_FSA
        /* Discard unphysical negative-J undershoots before the radiation update;
           no matter back-reaction (see prj_rad_fsa_clamp_intensities). */
        prj_rad_fsa_clamp_intensities(u);
#endif
        PRJ_SUBTIMER_START("sub_rad_freq_flux");
        prj_rad_freq_flux_apply(rad, block, block->W_rad, u, i, j, k, lapse_cell, dt);
        PRJ_SUBTIMER_STOP("sub_rad_freq_flux");
        PRJ_SUBTIMER_START("sub_rad_ang_flux");
        prj_rad_ang_flux_apply(rad, block, block->W_rad, u, i, j, k, lapse_cell, dt);
        PRJ_SUBTIMER_STOP("sub_rad_ang_flux");
#if PRJ_USE_RADIATION_FSA
        PRJ_SUBTIMER_START("sub_rad_inel");
        prj_rad_inel_fsa(rad, block, i, j, k, eos, u, dt, T_cell);
        PRJ_SUBTIMER_STOP("sub_rad_inel");
        PRJ_SUBTIMER_START("sub_rad_energy_momentum");
        prj_rad_energy_momentum_update_fsa(rad, block, i, j, k, eos, u, dt, lapse_cell);
        PRJ_SUBTIMER_STOP("sub_rad_energy_momentum");
#if DO_FFC
        PRJ_SUBTIMER_START("sub_rad_ffc");
        prj_rad_ffc_fsa(rad, u, dt);
        PRJ_SUBTIMER_STOP("sub_rad_ffc");
#endif
#else
        PRJ_SUBTIMER_START("sub_rad_eleinel");
        prj_rad_eleinel_step(rad, eos, u, dt, T_cell);
        PRJ_SUBTIMER_STOP("sub_rad_eleinel");
        PRJ_SUBTIMER_START("sub_rad_nucinel");
        prj_rad_nucinel_step(rad, eos, u, dt, T_cell);
        PRJ_SUBTIMER_STOP("sub_rad_nucinel");
        PRJ_SUBTIMER_START("sub_rad_energy");
        prj_rad_energy_update(rad, eos, u, dt, lapse_cell, &T_cell, kappa);
        PRJ_SUBTIMER_STOP("sub_rad_energy");
        PRJ_SUBTIMER_START("sub_rad_momentum");
        prj_rad_momentum_update(rad, eos, u, dt, lapse_cell, T_cell, kappa);
        PRJ_SUBTIMER_STOP("sub_rad_momentum");
#endif
        PRJ_SUBTIMER_STOP("sub_cell_radiation");
    }
    PRJ_SUBTIMER_START("sub_cell_store");
    prj_timeint_store_mhd_rad_cell(Wdst, Wdst_rad, i, j, k, u);
    PRJ_SUBTIMER_STOP("sub_cell_store");
#else
    double w[PRJ_NVAR_PRIM];
    double u[PRJ_NVAR_CONS];
    double u1[PRJ_NVAR_CONS];
    double fluxdiv[PRJ_NVAR_CONS];
    int v;

    PRJ_SUBTIMER_START("sub_cu_fluxdiv");
    prj_flux_div(block->flux, block->area, block->vol, i, j, k, fluxdiv);
    PRJ_SUBTIMER_STOP("sub_cu_fluxdiv");
    prj_timeint_cell_prim(block->W_mhd, i, j, k, w);
    PRJ_SUBTIMER_START("sub_cu_prim2cons");
    prj_eos_cell_prim2cons(eos, mesh, block, 0, i, j, k, w, u, PRJ_EOS_CTX_MAIN);
    PRJ_SUBTIMER_STOP("sub_cu_prim2cons");
    PRJ_SUBTIMER_START("sub_cu_dt_src");
    prj_timeint_update_dt_src(mesh, grav, block, u, i, j, k, mpi, dt_src);
    PRJ_SUBTIMER_STOP("sub_cu_dt_src");
    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        u1[v] = u[v] + dt * (prj_block_rhs_value_const(block, v, i, j, k) + fluxdiv[v]);
    }
#if PRJ_MHD
    {
        double *bf_dst[3];
        int d;

        for (d = 0; d < 3; ++d) {
            bf_dst[d] = prj_block_bf_stage(block, d, use_bf1 != 0 ? 1 : 0);
        }
        prj_timeint_mhd_set_cons_b_from_bf(block, bf_dst, i, j, k, u1);
    }
#else
    (void)use_bf1;
#endif
#if PRJ_NRAD > 0
    {
        double T_cell = block->eosvar[EIDX(PRJ_EOSVAR_TEMPERATURE, i, j, k)];
        double lapse_cell = prj_timeint_cell_lapse(block, i, j, k);
#if PRJ_USE_RADIATION_M1
        double kappa[PRJ_NRAD * PRJ_NEGROUP];
#endif

#if PRJ_USE_RADIATION_FSA
        /* Discard unphysical negative-J undershoots before the radiation update;
           no matter back-reaction (see prj_rad_fsa_clamp_intensities). */
        prj_rad_fsa_clamp_intensities(u1);
#endif
        prj_rad_freq_flux_apply(rad, block, block->W_rad, u1, i, j, k, lapse_cell, dt);
        prj_rad_ang_flux_apply(rad, block, block->W_rad, u1, i, j, k, lapse_cell, dt);
#if PRJ_USE_RADIATION_FSA
        prj_rad_inel_fsa(rad, block, i, j, k, eos, u1, dt, T_cell);
        prj_rad_energy_momentum_update_fsa(rad, block, i, j, k, eos, u1, dt, lapse_cell);
#if DO_FFC
        prj_rad_ffc_fsa(rad, u1, dt);
#endif
#else
        prj_rad_eleinel_step(rad, eos, u1, dt, T_cell);
        prj_rad_nucinel_step(rad, eos, u1, dt, T_cell);
        prj_rad_energy_update(rad, eos, u1, dt, lapse_cell, &T_cell, kappa);
        prj_rad_momentum_update(rad, eos, u1, dt, lapse_cell, T_cell, kappa);
#endif
    }
#else
    (void)rad;
#endif
    /* Fuse the stage geometry advance into the per-cell update: densitization
     * (prim2cons above) used slot 0 = g^n; recover against the advanced metric
     * g_new = g^n + dtau*rhs0 written in place to the recovery slot. The Z4c
     * point-update is pointwise, so mutating this cell's slot never disturbs
     * another cell's still-g^n densitization. See prj_z4c_update_linear_cell. */
    if (prj_z4c_runtime_enabled(mesh)) {
        PRJ_SUBTIMER_START("sub_cu_z4c_linear");
        prj_z4c_update_linear_cell(block, prj_stage_slot_from_bf_arg(use_bf1),
            0, 1.0, 0, 0.0, 0, dtau_cm, i, j, k);
        PRJ_SUBTIMER_STOP("sub_cu_z4c_linear");
    }
    PRJ_SUBTIMER_START("sub_cu_cons2prim");
    prj_eos_cell_cons2prim(eos, mesh, block,
        prj_stage_slot_from_bf_arg(use_bf1), i, j, k, u1, w, PRJ_EOS_CTX_MAIN);
    PRJ_SUBTIMER_STOP("sub_cu_cons2prim");
    prj_timeint_cell_prim_store(Wdst, i, j, k, w);
#endif
}
#endif

#if TIME_INTEGRATION == RK2
static void prj_timeint_update_cell_stage2_mhd_rad(const prj_mesh *mesh, prj_rad *rad, prj_eos *eos,
    prj_block *block, int i, int j, int k,
    double dt, const prj_grav *grav, const prj_mpi *mpi, double *dt_src, double dtau_cm)
{
#if PRJ_MHD && PRJ_NRAD > 0
    double u[PRJ_NVAR_CONS];
    double rho0;
    double mom10;
    double mom20;
    double mom30;
    double etot0;
    double ye0;
    double b10;
    double b20;
    double b30;
    double rho1;
    double mom11;
    double mom21;
    double mom31;
    double etot1;
    double ye1;
    double b11;
    double b21;
    double b31;
    double *W_stage1 = prj_block_mhd_stage(block, 1);
    double *W_rad_stage1 = prj_block_rad_stage(block, 1);
    int field;
    int group;

    (void)dtau_cm;
    PRJ_SUBTIMER_START("sub_cell_cons_from_prim");
    prj_timeint_mhd_hydro_cons_from_prim(block->W_mhd, i, j, k,
        &rho0, &mom10, &mom20, &mom30, &etot0, &ye0, &b10, &b20, &b30);
    prj_timeint_mhd_hydro_cons_from_prim(W_stage1, i, j, k,
        &rho1, &mom11, &mom21, &mom31, &etot1, &ye1, &b11, &b21, &b31);
    PRJ_SUBTIMER_STOP("sub_cell_cons_from_prim");
    PRJ_SUBTIMER_START("sub_cell_dt_src");
    prj_timeint_update_dt_src_values(mesh, grav, block, rho0, mom10, mom20, mom30, etot0,
        i, j, k, mpi, dt_src);
    PRJ_SUBTIMER_STOP("sub_cell_dt_src");

#define PRJ_STAGE2_UPDATE(v, u0v, u1v) do { \
        int vv = (v); \
        u[vv] = 0.5 * (u0v) + 0.5 * ((u1v) + dt * \
            (prj_block_rhs_value_const(block, vv, i, j, k) + \
                prj_timeint_flux_div_var(block, vv, i, j, k))); \
    } while (0)

    PRJ_SUBTIMER_START("sub_cell_flux_src_update");
    PRJ_STAGE2_UPDATE(PRJ_CONS_RHO, rho0, rho1);
    PRJ_STAGE2_UPDATE(PRJ_CONS_MOM1, mom10, mom11);
    PRJ_STAGE2_UPDATE(PRJ_CONS_MOM2, mom20, mom21);
    PRJ_STAGE2_UPDATE(PRJ_CONS_MOM3, mom30, mom31);
    PRJ_STAGE2_UPDATE(PRJ_CONS_ETOT, etot0, etot1);
    PRJ_STAGE2_UPDATE(PRJ_CONS_YE, ye0, ye1);

    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
#if PRJ_USE_RADIATION_FSA
            int angle;

            for (angle = 0; angle < PRJ_NANGLE; ++angle) {
                int iv = PRJ_CONS_RAD_I(field, group, angle);

                PRJ_STAGE2_UPDATE(iv, block->W_rad[WIDX(PRJ_RAD_PRIM_I(field, group, angle), i, j, k)],
                    W_rad_stage1[WIDX(PRJ_RAD_PRIM_I(field, group, angle), i, j, k)]);
            }
#else
            int e = PRJ_CONS_RAD_E(field, group);
            int f1 = PRJ_CONS_RAD_F1(field, group);
            int f2 = PRJ_CONS_RAD_F2(field, group);
            int f3 = PRJ_CONS_RAD_F3(field, group);

            PRJ_STAGE2_UPDATE(e, block->W_rad[WIDX(PRJ_RAD_PRIM_E(field, group), i, j, k)], W_rad_stage1[WIDX(PRJ_RAD_PRIM_E(field, group), i, j, k)]);
            PRJ_STAGE2_UPDATE(f1, block->W_rad[WIDX(PRJ_RAD_PRIM_F1(field, group), i, j, k)], W_rad_stage1[WIDX(PRJ_RAD_PRIM_F1(field, group), i, j, k)]);
            PRJ_STAGE2_UPDATE(f2, block->W_rad[WIDX(PRJ_RAD_PRIM_F2(field, group), i, j, k)], W_rad_stage1[WIDX(PRJ_RAD_PRIM_F2(field, group), i, j, k)]);
            PRJ_STAGE2_UPDATE(f3, block->W_rad[WIDX(PRJ_RAD_PRIM_F3(field, group), i, j, k)], W_rad_stage1[WIDX(PRJ_RAD_PRIM_F3(field, group), i, j, k)]);
#endif
        }
    }
    PRJ_SUBTIMER_STOP("sub_cell_flux_src_update");

#undef PRJ_STAGE2_UPDATE

    (void)b10;
    (void)b20;
    (void)b30;
    (void)b11;
    (void)b21;
    (void)b31;
    PRJ_SUBTIMER_START("sub_cell_mhd_set_b");
    prj_timeint_mhd_set_cons_b_from_bf(block, block->Bf, i, j, k, u);
    PRJ_SUBTIMER_STOP("sub_cell_mhd_set_b");
    {
        double T_cell = block->eosvar[EIDX(PRJ_EOSVAR_TEMPERATURE, i, j, k)];
        double lapse_cell = prj_timeint_cell_lapse(block, i, j, k);
#if PRJ_USE_RADIATION_M1
        double kappa[PRJ_NRAD * PRJ_NEGROUP];
#endif

        PRJ_SUBTIMER_START("sub_cell_radiation");
#if PRJ_USE_RADIATION_FSA
        /* Discard unphysical negative-J undershoots before the radiation update;
           no matter back-reaction (see prj_rad_fsa_clamp_intensities). */
        prj_rad_fsa_clamp_intensities(u);
#endif
        PRJ_SUBTIMER_START("sub_rad_freq_flux");
        prj_rad_freq_flux_apply(rad, block, W_rad_stage1, u, i, j, k, lapse_cell, 0.5 * dt);
        PRJ_SUBTIMER_STOP("sub_rad_freq_flux");
        PRJ_SUBTIMER_START("sub_rad_ang_flux");
        prj_rad_ang_flux_apply(rad, block, W_rad_stage1, u, i, j, k, lapse_cell, 0.5 * dt);
        PRJ_SUBTIMER_STOP("sub_rad_ang_flux");
#if PRJ_USE_RADIATION_FSA
        PRJ_SUBTIMER_START("sub_rad_inel");
        prj_rad_inel_fsa(rad, block, i, j, k, eos, u, 0.5 * dt, T_cell);
        PRJ_SUBTIMER_STOP("sub_rad_inel");
        PRJ_SUBTIMER_START("sub_rad_energy_momentum");
        prj_rad_energy_momentum_update_fsa(rad, block, i, j, k, eos, u, 0.5 * dt, lapse_cell);
        PRJ_SUBTIMER_STOP("sub_rad_energy_momentum");
#if DO_FFC
        PRJ_SUBTIMER_START("sub_rad_ffc");
        prj_rad_ffc_fsa(rad, u, 0.5 * dt);
        PRJ_SUBTIMER_STOP("sub_rad_ffc");
#endif
#else
        PRJ_SUBTIMER_START("sub_rad_eleinel");
        prj_rad_eleinel_step(rad, eos, u, 0.5 * dt, T_cell);
        PRJ_SUBTIMER_STOP("sub_rad_eleinel");
        PRJ_SUBTIMER_START("sub_rad_nucinel");
        prj_rad_nucinel_step(rad, eos, u, 0.5 * dt, T_cell);
        PRJ_SUBTIMER_STOP("sub_rad_nucinel");
        PRJ_SUBTIMER_START("sub_rad_energy");
        prj_rad_energy_update(rad, eos, u, 0.5 * dt, lapse_cell, &T_cell, kappa);
        PRJ_SUBTIMER_STOP("sub_rad_energy");
        PRJ_SUBTIMER_START("sub_rad_momentum");
        prj_rad_momentum_update(rad, eos, u, 0.5 * dt, lapse_cell, T_cell, kappa);
        PRJ_SUBTIMER_STOP("sub_rad_momentum");
#endif
        PRJ_SUBTIMER_STOP("sub_cell_radiation");
    }
    PRJ_SUBTIMER_START("sub_cell_store");
    prj_timeint_store_mhd_rad_cell(block->W_mhd, block->W_rad, i, j, k, u);
    PRJ_SUBTIMER_STOP("sub_cell_store");
#else
    double w[PRJ_NVAR_PRIM];
    double u[PRJ_NVAR_CONS];
    double u1[PRJ_NVAR_CONS];
    double fluxdiv[PRJ_NVAR_CONS];
    double *W_stage1 = prj_block_prim_stage(block, 1);
    int v;

    prj_flux_div(block->flux, block->area, block->vol, i, j, k, fluxdiv);
    prj_timeint_cell_prim(block->W_mhd, i, j, k, w);
    prj_eos_cell_prim2cons(eos, mesh, block, 0, i, j, k, w, u, PRJ_EOS_CTX_MAIN);
    prj_timeint_cell_prim(W_stage1, i, j, k, w);
    prj_eos_cell_prim2cons(eos, mesh, block, 1, i, j, k, w, u1, PRJ_EOS_CTX_MAIN);
    prj_timeint_update_dt_src(mesh, grav, block, u, i, j, k, mpi, dt_src);
    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        u[v] = 0.5 * u[v] + 0.5 *
            (u1[v] + dt * (prj_block_rhs_value_const(block, v, i, j, k) + fluxdiv[v]));
    }
#if PRJ_MHD
    prj_timeint_mhd_set_cons_b_from_bf(block, block->Bf, i, j, k, u);
#endif
#if PRJ_NRAD > 0
    {
        double T_cell = block->eosvar[EIDX(PRJ_EOSVAR_TEMPERATURE, i, j, k)];
        double lapse_cell = prj_timeint_cell_lapse(block, i, j, k);
#if PRJ_USE_RADIATION_M1
        double kappa[PRJ_NRAD * PRJ_NEGROUP];
#endif

#if PRJ_USE_RADIATION_FSA
        /* Discard unphysical negative-J undershoots before the radiation update;
           no matter back-reaction (see prj_rad_fsa_clamp_intensities). */
        prj_rad_fsa_clamp_intensities(u);
#endif
        prj_rad_freq_flux_apply(rad, block, W_rad_stage1, u, i, j, k, lapse_cell, 0.5 * dt);
        prj_rad_ang_flux_apply(rad, block, W_rad_stage1, u, i, j, k, lapse_cell, 0.5 * dt);
#if PRJ_USE_RADIATION_FSA
        prj_rad_inel_fsa(rad, block, i, j, k, eos, u, 0.5 * dt, T_cell);
        prj_rad_energy_momentum_update_fsa(rad, block, i, j, k, eos, u, 0.5 * dt, lapse_cell);
#if DO_FFC
        prj_rad_ffc_fsa(rad, u, 0.5 * dt);
#endif
#else
        prj_rad_eleinel_step(rad, eos, u, 0.5 * dt, T_cell);
        prj_rad_nucinel_step(rad, eos, u, 0.5 * dt, T_cell);
        prj_rad_energy_update(rad, eos, u, 0.5 * dt, lapse_cell, &T_cell, kappa);
        prj_rad_momentum_update(rad, eos, u, 0.5 * dt, lapse_cell, T_cell, kappa);
#endif
    }
#else
    (void)rad;
#endif
    /* Fuse the final geometry advance: the two prim2cons above densitized with
     * slot 0 = g^n (u0) and slot 1 = g1 (u1). Recover the RK2-combined state
     * against the final metric g^{n+1} = 1/2 g^n + 1/2 g1 + 1/2 dtau*rhs1
     * (rhs1 lives in rhs slot 0 here), written in place to slot 0. Pointwise, so
     * later cells still read their own g^n/g1 for their prim2cons. */
    if (prj_z4c_runtime_enabled(mesh)) {
        prj_z4c_update_linear_cell(block, 0, 0, 0.5, 1, 0.5, 0, 0.5 * dtau_cm, i, j, k);
    }
    prj_eos_cell_cons2prim(eos, mesh, block, 0, i, j, k, u, w, PRJ_EOS_CTX_MAIN);
    prj_timeint_cell_prim_store(block->W_mhd, i, j, k, w);
#endif
}
#endif

#if TIME_INTEGRATION == PRJ_TIMEINT_IMEX
#if defined(PRJ_DEBUG)
static int prj_timeint_imex_debug_tracking_active = 0;
static int prj_timeint_imex_debug_nstages = 0;
static int prj_timeint_imex_debug_stage_solved[PRJ_BLOCK_NSTAGES];
#endif

static void prj_timeint_imex_fail(const char *message)
{
    fprintf(stderr, "%s\n", message);
    exit(EXIT_FAILURE);
}

static void prj_timeint_imex_validate(const prj_timeint_imex_tableau *tableau)
{
    int i;
    int j;

    if (tableau == 0) {
        prj_timeint_imex_fail("prj_timeint_imex: missing tableau");
    }
    if (tableau->nstages < 1 || tableau->nstages > PRJ_BLOCK_NSTAGES) {
        prj_timeint_imex_fail("prj_timeint_imex: tableau stage count exceeds block storage");
    }
    if (tableau->a_ex == 0 || tableau->a_im == 0 ||
        tableau->b_ex == 0 || tableau->b_im == 0) {
        prj_timeint_imex_fail("prj_timeint_imex: incomplete tableau");
    }
    if (!(tableau->r > 0.0) || !isfinite(tableau->r)) {
        prj_timeint_imex_fail("prj_timeint_imex: invalid SSP coefficient");
    }
    for (i = 0; i < tableau->nstages; ++i) {
        if (!isfinite(tableau->b_ex[i]) || !isfinite(tableau->b_im[i])) {
            prj_timeint_imex_fail("prj_timeint_imex: non-finite final coefficient");
        }
        /* A zero implicit diagonal is allowed: step_im treats dt_implicit == 0
         * as an explicit (no-op) stage, leaving u unchanged and deriv = 0. */
        for (j = 0; j < tableau->nstages; ++j) {
            double a_ex = tableau->a_ex[(size_t)i * (size_t)tableau->nstages + (size_t)j];
            double a_im = tableau->a_im[(size_t)i * (size_t)tableau->nstages + (size_t)j];

            if (!isfinite(a_ex) || !isfinite(a_im)) {
                prj_timeint_imex_fail("prj_timeint_imex: non-finite stage coefficient");
            }
            if (j >= i && a_ex != 0.0) {
                prj_timeint_imex_fail("prj_timeint_imex: explicit tableau is not strictly lower triangular");
            }
            if (j > i && a_im != 0.0) {
                prj_timeint_imex_fail("prj_timeint_imex: implicit tableau is not lower triangular");
            }
        }
    }
}

static inline double prj_timeint_imex_a_ex(const prj_timeint_imex_tableau *tableau,
    int row, int col)
{
    return tableau->a_ex[(size_t)row * (size_t)tableau->nstages + (size_t)col];
}

static inline double prj_timeint_imex_a_im(const prj_timeint_imex_tableau *tableau,
    int row, int col)
{
    return tableau->a_im[(size_t)row * (size_t)tableau->nstages + (size_t)col];
}

static void prj_timeint_imex_cons_from_stage(prj_eos *eos, const prj_mesh *mesh,
    const prj_block *block, int stage, int i, int j, int k, double *u)
{
    double w[PRJ_NVAR_PRIM];
    double *W_stage;

    if (block == 0 || u == 0) {
        return;
    }
    W_stage = prj_block_prim_stage((prj_block *)block, stage);
    prj_timeint_cell_prim(W_stage, i, j, k, w);
    prj_eos_cell_prim2cons(eos, mesh, block, stage, i, j, k, w, u,
        PRJ_EOS_CTX_MAIN);
#if PRJ_MHD
    {
        double *bf[3];
        int d;

        for (d = 0; d < 3; ++d) {
            bf[d] = prj_block_bf_stage((prj_block *)block, d, stage);
        }
        prj_timeint_mhd_set_cons_b_from_bf(block, bf, i, j, k, u);
    }
#endif
}

static void prj_timeint_imex_store_cons_to_stage(prj_eos *eos, const prj_mesh *mesh,
    prj_block *block, int stage, int i, int j, int k, double *u)
{
    double w[PRJ_NVAR_PRIM];
    double *W_stage;

    if (block == 0 || u == 0) {
        return;
    }
    W_stage = prj_block_prim_stage(block, stage);
#if PRJ_MHD
    {
        double *bf[3];
        int d;

        for (d = 0; d < 3; ++d) {
            bf[d] = prj_block_bf_stage(block, d, stage);
        }
        prj_timeint_mhd_set_cons_b_from_bf(block, bf, i, j, k, u);
    }
#endif
    prj_eos_cell_cons2prim(eos, mesh, block, stage, i, j, k, u, w,
        PRJ_EOS_CTX_MAIN);
    prj_timeint_cell_prim_store(W_stage, i, j, k, w);
}

static void prj_timeint_imex_save_base_cons(prj_mesh *mesh, prj_eos *eos, const prj_mpi *mpi)
{
    int bidx;

    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        int i;
        int j;
        int k;

        if (!prj_timeint_local_block(mpi, block) || !prj_block_has_cons_storage(block)) {
            continue;
        }
        for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
            for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    double u[PRJ_NVAR_CONS];

                    prj_timeint_imex_cons_from_stage(eos, mesh, block, 0, i, j, k, u);
                    prj_block_store_cons_cell(block, i, j, k, u);
                }
            }
        }
    }
}

#if PRJ_MHD
static void prj_timeint_mhd_update_mesh_emf_stage(prj_mesh *mesh, const prj_mpi *mpi,
    int stage)
{
    int bidx;

    if (mesh == 0) {
        prj_timeint_mhd_fail("prj_timeint_mhd_update_mesh_emf_stage: invalid input");
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (prj_timeint_local_block(mpi, block)) {
            prj_timeint_mhd_update_emf(mesh, block, stage, mesh->mhd_eta);
        }
    }
    prj_mhd_emf_send(mesh, mpi);
}

static void prj_timeint_mhd_save_deriv_bf(prj_block *block, int stage)
{
    double *src[3];
    int dir;

    if (block == 0) {
        prj_timeint_mhd_fail("prj_timeint_mhd_save_deriv_bf: invalid input");
    }
    for (dir = 0; dir < 3; ++dir) {
        src[dir] = prj_block_bf_stage(block, dir, stage);
        if (src[dir] == 0 || prj_block_deriv_bf_stage(block, dir, stage) == 0) {
            prj_timeint_mhd_fail("prj_timeint_mhd_save_deriv_bf: missing Bf derivative storage");
        }
    }
    for (dir = 0; dir < 3; ++dir) {
        double *dst = prj_block_deriv_bf_stage(block, dir, stage);
        int i;
        int j;
        int k;

        for (i = 0; i <= prj_timeint_mhd_face_axis_max(dir, 0); ++i) {
            for (j = 0; j <= prj_timeint_mhd_face_axis_max(dir, 1); ++j) {
                for (k = 0; k <= prj_timeint_mhd_face_axis_max(dir, 2); ++k) {
                    int idx = FACE_IDX(dir, i, j, k);

                    dst[idx] = prj_timeint_mhd_ct_value(block, src, dir, i, j, k, 1.0) -
                        src[dir][idx];
                }
            }
        }
    }
}
#endif

static void prj_timeint_imex_update_bf_stage(prj_block *block,
    const prj_timeint_imex_tableau *tableau, int dst_stage, int coeff_row, int final_stage,
    int nterms, double dt)
{
#if PRJ_MHD
    int dir;

    if (block == 0) {
        return;
    }
    for (dir = 0; dir < 3; ++dir) {
        double *dst = prj_block_bf_stage(block, dir, dst_stage);
        double *base = prj_block_bf_stage(block, dir, 0);
        int i;
        int j;
        int k;

        if (dst == 0 || base == 0) {
            prj_timeint_mhd_fail("prj_timeint_imex_update_bf_stage: missing Bf storage");
        }
        for (i = 0; i <= prj_timeint_mhd_face_axis_max(dir, 0); ++i) {
            for (j = 0; j <= prj_timeint_mhd_face_axis_max(dir, 1); ++j) {
                for (k = 0; k <= prj_timeint_mhd_face_axis_max(dir, 2); ++k) {
                    int idx = FACE_IDX(dir, i, j, k);
                    double value = base[idx];
                    int s;

                    for (s = 0; s < nterms; ++s) {
                        double coeff = final_stage != 0 ? tableau->b_ex[s] :
                            prj_timeint_imex_a_ex(tableau, coeff_row, s);
                        const double *deriv = prj_block_deriv_bf_stage_const(block, dir, s);

                        if (coeff == 0.0) {
                            continue;
                        }
                        value += dt * coeff * deriv[idx];
                    }
                    dst[idx] = value;
                }
            }
        }
    }
#else
    (void)block;
    (void)tableau;
    (void)dst_stage;
    (void)coeff_row;
    (void)final_stage;
    (void)nterms;
    (void)dt;
#endif
}

static void prj_timeint_imex_assemble_stage(prj_mesh *mesh, prj_eos *eos,
    const prj_mpi *mpi, const prj_timeint_imex_tableau *tableau,
    int dst_stage, int coeff_row, int final_stage, int nterms, double dt)
{
    int bidx;

    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        int i;
        int j;
        int k;

        if (!prj_timeint_local_block(mpi, block)) {
            continue;
        }
        prj_timeint_imex_update_bf_stage(block, tableau, dst_stage, coeff_row, final_stage,
            nterms, dt);
        for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
            for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    double u[PRJ_NVAR_CONS];
                    int v;
                    int s;

                    prj_block_load_cons_cell_const(block, i, j, k, u);
                    for (s = 0; s < nterms; ++s) {
                        double coeff_ex = final_stage != 0 ? tableau->b_ex[s] :
                            prj_timeint_imex_a_ex(tableau, coeff_row, s);
                        double coeff_im = final_stage != 0 ? tableau->b_im[s] :
                            prj_timeint_imex_a_im(tableau, coeff_row, s);
                        const double *deriv_ex = prj_block_deriv_stage_const(block->deriv_ex, s);
                        const double *deriv_im = prj_block_deriv_stage_const(block->deriv_im, s);

                        if (coeff_ex == 0.0 && coeff_im == 0.0) {
                            continue;
                        }
                        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
                            size_t idx = (size_t)VIDX(v, i, j, k);

                            u[v] += dt * (coeff_ex * deriv_ex[idx] + coeff_im * deriv_im[idx]);
                        }
                    }
                    prj_timeint_imex_store_cons_to_stage(eos, mesh, block, dst_stage,
                        i, j, k, u);
                }
            }
        }
    }
}

#if PRJ_NRAD > 0
static void prj_timeint_imex_add_explicit_rad_deriv(prj_eos *eos, prj_rad *rad,
    const prj_mesh *mesh, prj_block *block, int stage, int i, int j, int k,
    double dt, double *deriv)
{
    double u0[PRJ_NVAR_CONS];
    double u1[PRJ_NVAR_CONS];
    double *W_stage;
    double *W_rad_stage;
#if PRJ_USE_RADIATION_M1 || PRJ_USE_RADIATION_FSA
    double T_cell;
#endif
    double lapse_cell;
    int v;

    if (dt == 0.0 || deriv == 0) {
        return;
    }
    W_stage = prj_block_mhd_stage(block, stage);
    W_rad_stage = prj_block_rad_stage(block, stage);
    if (W_stage == 0) {
        prj_timeint_imex_fail("prj_timeint_step_ex: missing stage storage");
    }
    prj_timeint_imex_cons_from_stage(eos, mesh, block, stage, i, j, k, u0);
    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        u1[v] = u0[v];
    }
    lapse_cell = prj_timeint_cell_lapse(block, i, j, k);

#if PRJ_USE_RADIATION_FSA
    /* Discard unphysical negative-J undershoots before the radiation update;
       no matter back-reaction (see prj_rad_fsa_clamp_intensities). */
    prj_rad_fsa_clamp_intensities(u1);
#endif
    PRJ_SUBTIMER_START("sub_rad_freq_flux");
    prj_rad_freq_flux_apply(rad, block, W_rad_stage, u1, i, j, k, lapse_cell, dt);
    PRJ_SUBTIMER_STOP("sub_rad_freq_flux");
    PRJ_SUBTIMER_START("sub_rad_ang_flux");
    prj_rad_ang_flux_apply(rad, block, W_rad_stage, u1, i, j, k, lapse_cell, dt);
    PRJ_SUBTIMER_STOP("sub_rad_ang_flux");
#if PRJ_USE_RADIATION_M1 || PRJ_USE_RADIATION_FSA
    T_cell = block->eosvar[EIDX(PRJ_EOSVAR_TEMPERATURE, i, j, k)];
#if PRJ_USE_RADIATION_FSA
    PRJ_SUBTIMER_START("sub_rad_inel");
    prj_rad_inel_fsa(rad, block, i, j, k, eos, u1, dt, T_cell);
    PRJ_SUBTIMER_STOP("sub_rad_inel");
#else
    PRJ_SUBTIMER_START("sub_rad_eleinel");
    prj_rad_eleinel_step(rad, eos, u1, dt, T_cell);
    PRJ_SUBTIMER_STOP("sub_rad_eleinel");
    PRJ_SUBTIMER_START("sub_rad_nucinel");
    prj_rad_nucinel_step(rad, eos, u1, dt, T_cell);
    PRJ_SUBTIMER_STOP("sub_rad_nucinel");
#endif
#endif

    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        size_t idx = (size_t)VIDX(v, i, j, k);

        deriv[idx] += (u1[v] - u0[v]) / dt;
    }
}
#endif

static void prj_timeint_imex_fill_updated_stage(prj_mesh *mesh, const prj_bc *bc,
    prj_eos *eos, prj_rad *rad, prj_grav *grav, prj_mpi *mpi, prj_timer *timer,
    int stage, int fill_opacity_halo)
{
    PRJ_TIMER_BARRIER_START(timer, mpi, "imex_eos_fill_active");
    prj_eos_fill_active_cells(mesh, eos, mpi, stage + 1, PRJ_EOS_CTX_MAIN);
    PRJ_TIMER_BARRIER_STOP(timer, mpi, "imex_eos_fill_active");

    PRJ_TIMER_BARRIER_START(timer, mpi, "imex_ghost_grav_opac");
    prj_boundary_fill_ghosts_and_bf(mesh, mpi, bc, stage + 1, stage, eos, grav, rad,
        PRJ_BOUNDARY_TIMER_SCOPE_NONE);
    PRJ_TIMER_BARRIER_STOP(timer, mpi, "imex_ghost_grav_opac");

    PRJ_TIMER_BARRIER_START(timer, mpi, "imex_eos_fill_mesh");
    prj_eos_fill_mesh(mesh, eos, mpi, stage + 1, PRJ_EOS_CTX_MAIN);
    PRJ_TIMER_BARRIER_STOP(timer, mpi, "imex_eos_fill_mesh");

    if (fill_opacity_halo) {
        prj_flux_fill_transport_opacity_halo(mesh, rad, mpi, stage + 1);
    }
}
#endif

void prj_timeint_init(const prj_timeint_imex_tableau *tableau)
{
#if TIME_INTEGRATION == PRJ_TIMEINT_IMEX
    prj_timeint_imex_validate(tableau);
#else
    (void)tableau;
#endif
}

#if PRJ_USE_RADIATION_M1
static double prj_timeint_cell_rad_denom(const double *w, const double dx[3])
{
    double max_denom = 0.0;
    int field;
    int group;

    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            double E = w[PRJ_PRIM_RAD_E(field, group)];
            double F1 = w[PRJ_PRIM_RAD_F1(field, group)];
            double F2 = w[PRJ_PRIM_RAD_F2(field, group)];
            double F3 = w[PRJ_PRIM_RAD_F3(field, group)];
            double denom = 0.0;
            double lam_min;
            double lam_max;
            double c_abs;
            double cE;
            double Fmag;
            double inv_Fmag;
            double f;

            /* Fmag, inv_Fmag and f are rotation-invariant, so compute them once
             * and reuse across the three direction-normal wavespeed evaluations. */
            cE = PRJ_CLIGHT * (E > 0.0 ? E : 0.0);
            Fmag = sqrt(F1 * F1 + F2 * F2 + F3 * F3);
            inv_Fmag = (Fmag > 0.0) ? (1.0 / Fmag) : 0.0;
            if (cE > 0.0) {
                f = Fmag / cE;
                if (f > 1.0) {
                    f = 1.0;
                }
            } else {
                f = 0.0;
            }

            /* The first flux argument is the direction-normal component. */
            prj_rad_m1_wavespeeds_with_fluxmag(E, F1, Fmag, inv_Fmag, f, &lam_min, &lam_max);
            c_abs = fabs(lam_min);
            if (fabs(lam_max) > c_abs) {
                c_abs = fabs(lam_max);
            }
            c_abs *= PRJ_CLIGHT;
            denom += c_abs / dx[0];

            prj_rad_m1_wavespeeds_with_fluxmag(E, F2, Fmag, inv_Fmag, f, &lam_min, &lam_max);
            c_abs = fabs(lam_min);
            if (fabs(lam_max) > c_abs) {
                c_abs = fabs(lam_max);
            }
            c_abs *= PRJ_CLIGHT;
            denom += c_abs / dx[1];

            prj_rad_m1_wavespeeds_with_fluxmag(E, F3, Fmag, inv_Fmag, f, &lam_min, &lam_max);
            c_abs = fabs(lam_min);
            if (fabs(lam_max) > c_abs) {
                c_abs = fabs(lam_max);
            }
            c_abs *= PRJ_CLIGHT;
            denom += c_abs / dx[2];

            /* Preserve the correlation between directional wavespeeds belonging
             * to one M1 state. Different groups are independent systems, so the
             * cell CFL denominator is the largest complete group denominator. */
            if (denom > max_denom) {
                max_denom = denom;
            }
        }
    }

    return max_denom;
}
#endif

double prj_timeint_calc_dt(const prj_mesh *mesh, prj_eos *eos, const prj_rad *rad,
    const prj_mpi *mpi, double cfl, const prj_timeint_imex_tableau *tableau)
{
    double dt_min = 1.0e99;
    double dt_global;
    int bidx;

    (void)eos;
    (void)rad;
#if TIME_INTEGRATION != PRJ_TIMEINT_IMEX
    (void)tableau;
#endif

    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        const prj_block *block = &mesh->blocks[bidx];
        int i;
        int j;
        int k;

        if (!prj_timeint_local_block(mpi, block)) {
            continue;
        }
        for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
            for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    double q[PRJ_EOS_NQUANT];
                    double w[PRJ_NVAR_PRIM];
                    double denom;
                    double cs;
                    double dt_cell;

                    prj_timeint_cell_prim(block->W_mhd, i, j, k, w);
                    q[PRJ_EOS_PRESSURE] = block->eosvar[EIDX(PRJ_EOSVAR_PRESSURE, i, j, k)];
                    q[PRJ_EOS_GAMMA] = block->eosvar[EIDX(PRJ_EOSVAR_GAMMA, i, j, k)];
                    cs = sqrt(q[PRJ_EOS_GAMMA] * q[PRJ_EOS_PRESSURE] / w[PRJ_PRIM_RHO]);
                    denom =
                        (fabs(w[PRJ_PRIM_V1]) + cs) / block->dx[0] +
                        (fabs(w[PRJ_PRIM_V2]) + cs) / block->dx[1] +
                        (fabs(w[PRJ_PRIM_V3]) + cs) / block->dx[2];
                    dt_cell = cfl / denom;
#if PRJ_USE_RADIATION_M1
                    {
                        double rad_denom = prj_timeint_cell_rad_denom(w, block->dx);

                        if (rad_denom > 0.0) {
                            double dt_rad = cfl / rad_denom;

                            if (dt_rad < dt_cell) {
                                dt_cell = dt_rad;
                            }
                        }
                    }
#endif
#if PRJ_USE_RADIATION_FSA
                    {
                        double lapse_c = prj_timeint_cell_lapse(block, i, j, k) * PRJ_CLIGHT;
                        double rad_denom = 0.0;
                        int a;

                        /* The angular-cell directions n0 are shared across all
                         * radiation species and energy groups, so iterate over
                         * angles only and take the fastest projected signal. */
                        for (a = 0; a < PRJ_NANGLE; ++a) {
                            double n[3];
                            double denom_a =
                                0.0;

                            prj_rad_fsa_rotated_angle_dir(rad, block, a, i, j, k, n);
                            denom_a =
                                (fabs(w[PRJ_PRIM_V1]) + lapse_c * n[0]) / block->dx[0] +
                                (fabs(w[PRJ_PRIM_V2]) + lapse_c * n[1]) / block->dx[1] +
                                (fabs(w[PRJ_PRIM_V3]) + lapse_c * n[2]) / block->dx[2];

                            if (a == 0 || denom_a > rad_denom) {
                                rad_denom = denom_a;
                            }
                        }
                        if (rad_denom > 0.0) {
                            double dt_rad = cfl / rad_denom;

                            if (dt_rad < dt_cell) {
                                dt_cell = dt_rad;
                            }
                        }
                    }
#endif
                    if (dt_cell < dt_min) {
                        dt_min = dt_cell;
                    }
                }
            }
        }
    }

    dt_global = prj_mpi_min_dt(mpi, dt_min);
    if (prj_z4c_runtime_enabled(mesh)) {
        double dt_z4c = prj_z4c_calc_dt_seconds(mesh, mpi, cfl);

        if (dt_z4c < dt_global) {
            dt_global = dt_z4c;
        }
    }
#if TIME_INTEGRATION == PRJ_TIMEINT_ESSPRK
    dt_global *= (double)(PRJ_TIMEINT_ESSPRK_N - 1);
#elif TIME_INTEGRATION == PRJ_TIMEINT_ESSPRK9_3
    dt_global *= 6.0;
#elif TIME_INTEGRATION == PRJ_TIMEINT_IMEX
    dt_global *= tableau->r;
#endif
    return dt_global;
}

#if PRJ_TIMEINT_USES_ESSPRK_STEP
static void prj_timeint_eSSPRK_bad_saved_state(const char *name)
{
    fprintf(stderr, "%s: invalid saved state\n", name);
    exit(EXIT_FAILURE);
}

static double *prj_timeint_eSSPRK_saved_W(prj_block *block, int saved_state)
{
    if (block == 0) {
        prj_timeint_eSSPRK_bad_saved_state("prj_timeint_eSSPRK_saved_W");
    }
    if (saved_state == 1) {
        return prj_block_prim_stage(block, 1);
    }
#if PRJ_TIMEINT_EXTRA_SAVED_STATES
    if (saved_state == 2) {
        return prj_block_prim_stage(block, 2);
    }
    if (saved_state == 3) {
        return prj_block_prim_stage(block, 3);
    }
#endif
    prj_timeint_eSSPRK_bad_saved_state("prj_timeint_eSSPRK_saved_W");
    return 0;
}

#if PRJ_MHD
static double *prj_timeint_eSSPRK_saved_bf_dir(prj_block *block, int saved_state, int dir)
{
    if (block == 0 || dir < 0 || dir >= 3) {
        prj_timeint_eSSPRK_bad_saved_state("prj_timeint_eSSPRK_saved_bf");
    }
    if (saved_state == 1) {
        return prj_block_bf_stage(block, dir, 1);
    }
#if PRJ_TIMEINT_EXTRA_SAVED_STATES
    if (saved_state == 2) {
        return prj_block_bf_stage(block, dir, 2);
    }
    if (saved_state == 3) {
        return prj_block_bf_stage(block, dir, 3);
    }
#endif
    prj_timeint_eSSPRK_bad_saved_state("prj_timeint_eSSPRK_saved_bf");
    return 0;
}
#endif

static void prj_timeint_eSSPRK_save_state(prj_mesh *mesh, const prj_mpi *mpi, int saved_state)
{
    int bidx;

    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        double *W_saved;
        int i;
        int j;
        int k;
        int v;

        if (!prj_timeint_local_block(mpi, block)) {
            continue;
        }
        W_saved = prj_timeint_eSSPRK_saved_W(block, saved_state);
        for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
            for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
                        W_saved[WIDX(v, i, j, k)] = block->W_mhd[WIDX(v, i, j, k)];
                    }
                }
            }
        }
#if PRJ_MHD
        {
            int d;
            int n;

            for (d = 0; d < 3; ++d) {
                double *bf_saved = prj_timeint_eSSPRK_saved_bf_dir(block, saved_state, d);

                for (n = 0; n < PRJ_BLOCK_NFACES; ++n) {
                    bf_saved[n] = block->Bf[d][n];
                }
            }
        }
#endif
    }
    prj_z4c_save_stage(mesh, mpi, saved_state, 0);
}

#if PRJ_MHD
static void prj_timeint_eSSPRK_blend_bf(prj_block *block, int saved_state, double saved_weight)
{
    double current_weight = 1.0 - saved_weight;
    int dir;

    for (dir = 0; dir < 3; ++dir) {
        double *bf_saved = prj_timeint_eSSPRK_saved_bf_dir(block, saved_state, dir);
        int i;
        int j;
        int k;

        for (i = 0; i <= prj_timeint_mhd_face_axis_max(dir, 0); ++i) {
            for (j = 0; j <= prj_timeint_mhd_face_axis_max(dir, 1); ++j) {
                for (k = 0; k <= prj_timeint_mhd_face_axis_max(dir, 2); ++k) {
                    int idx = FACE_IDX(dir, i, j, k);

                    block->Bf[dir][idx] =
                        saved_weight * bf_saved[idx] + current_weight * block->Bf[dir][idx];
                }
            }
        }
    }
}
#endif

static void prj_timeint_eSSPRK_final_cell(prj_eos *eos, const prj_mesh *mesh,
    prj_block *block, int i, int j, int k, int saved_state, double saved_weight)
{
    double *W_saved = prj_timeint_eSSPRK_saved_W(block, saved_state);
    double w[PRJ_NVAR_PRIM];
    double u[PRJ_NVAR_CONS];
    double u_saved[PRJ_NVAR_CONS];
    double current_weight = 1.0 - saved_weight;
    int v;

    prj_timeint_cell_prim(block->W_mhd, i, j, k, w);
    prj_eos_cell_prim2cons(eos, mesh, block, 0, i, j, k, w, u, PRJ_EOS_CTX_MAIN);
    prj_timeint_cell_prim(W_saved, i, j, k, w);
    prj_eos_cell_prim2cons(eos, mesh, block, saved_state, i, j, k, w,
        u_saved, PRJ_EOS_CTX_MAIN);
    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        u[v] = saved_weight * u_saved[v] + current_weight * u[v];
    }
#if PRJ_MHD
    prj_timeint_mhd_set_cons_b_from_bf(block, block->Bf, i, j, k, u);
#endif
#if PRJ_MHD && PRJ_NRAD > 0
    prj_timeint_store_mhd_rad_cell(block->W_mhd, block->W_rad, i, j, k, u);
#else
    prj_eos_cell_cons2prim(eos, mesh, block, 0, i, j, k, u, w, PRJ_EOS_CTX_MAIN);
    prj_timeint_cell_prim_store(block->W_mhd, i, j, k, w);
#endif
}

static void prj_timeint_eSSPRK_blend_with_saved(prj_mesh *mesh, const prj_bc *bc,
    prj_eos *eos, prj_rad *rad, prj_grav *grav, prj_mpi *mpi, prj_timer *timer,
    int saved_state, double saved_weight, const char *cell_timer, const char *active_timer,
    const char *ghost_timer, const char *mesh_timer)
{
    int bidx;

    PRJ_TIMER_BARRIER_START(timer, mpi, cell_timer);
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        int i;
        int j;
        int k;

        if (!prj_timeint_local_block(mpi, block)) {
            continue;
        }
#if PRJ_MHD
        prj_timeint_eSSPRK_blend_bf(block, saved_state, saved_weight);
#endif
        for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
            for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    prj_timeint_eSSPRK_final_cell(eos, mesh, block, i, j, k,
                        saved_state, saved_weight);
                }
            }
        }
    }
    PRJ_TIMER_BARRIER_STOP(timer, mpi, cell_timer);

    prj_z4c_blend_with_saved(mesh, mpi, bc, saved_state, saved_weight);

    PRJ_TIMER_BARRIER_START(timer, mpi, active_timer);
    prj_eos_fill_active_cells(mesh, eos, mpi, 1, PRJ_EOS_CTX_MAIN);
    PRJ_TIMER_BARRIER_STOP(timer, mpi, active_timer);

    PRJ_TIMER_BARRIER_START(timer, mpi, ghost_timer);
    prj_boundary_fill_ghosts_and_bf(mesh, mpi, bc, 1, 0, eos, grav, rad,
        PRJ_BOUNDARY_TIMER_SCOPE_NONE);
    PRJ_TIMER_BARRIER_STOP(timer, mpi, ghost_timer);

    PRJ_TIMER_BARRIER_START(timer, mpi, mesh_timer);
    prj_eos_fill_mesh(mesh, eos, mpi, 1, PRJ_EOS_CTX_MAIN);
    PRJ_TIMER_BARRIER_STOP(timer, mpi, mesh_timer);

    prj_flux_fill_transport_opacity_halo(mesh, rad, mpi, 1);
}
#endif

void prj_timeint_eSSPRK_step(prj_mesh *mesh, const prj_coord *coord, const prj_bc *bc, prj_eos *eos,
    prj_rad *rad, prj_grav *grav, prj_mpi *mpi, double dt, double *dt_src, prj_timer *timer)
{
#if PRJ_TIMEINT_USES_ESSPRK_STEP
    int bidx;

    (void)coord;

    PRJ_TIMER_BARRIER_START(timer, mpi, "essprk_step_flux_send");
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (prj_timeint_local_block(mpi, block)) {
            prj_flux_update(eos, rad, mesh, block, block->W_mhd, block->eosvar,
                block->flux, 0);
        }
    }
    PRJ_SUBTIMER_START("sub_riemann_flux_send");
    prj_riemann_flux_send(mesh, mpi);
    PRJ_SUBTIMER_STOP("sub_riemann_flux_send");
    PRJ_TIMER_BARRIER_STOP(timer, mpi, "essprk_step_flux_send");

#if PRJ_MHD
    PRJ_TIMER_BARRIER_START(timer, mpi, "essprk_step_mhd_emf");
    prj_timeint_mhd_update_mesh_emf(mesh, mpi, 0);
    PRJ_TIMER_BARRIER_STOP(timer, mpi, "essprk_step_mhd_emf");
#endif

    PRJ_TIMER_BARRIER_START(timer, mpi, "essprk_step_flux_emf_exchange");
    prj_mpi_exchange_fluxes_and_emf(mesh, mpi);
    PRJ_TIMER_BARRIER_STOP(timer, mpi, "essprk_step_flux_emf_exchange");

#if PRJ_MHD
    PRJ_TIMER_BARRIER_START(timer, mpi, "essprk_step_mhd_update_bf");
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (prj_timeint_local_block(mpi, block)) {
            prj_timeint_mhd_update_bf(block, 0, dt);
        }
    }
    PRJ_TIMER_BARRIER_STOP(timer, mpi, "essprk_step_mhd_update_bf");
#endif

    if (prj_z4c_runtime_enabled(mesh)) {
        double tau_cm = PRJ_CLIGHT * mesh->time_seconds;

        prj_z4c_compute_rhs(mesh, mpi, rad, 0, 0, tau_cm);
        prj_z4c_apply_sommerfeld_rhs(mesh, mpi, bc, 0, 0);
    }

    PRJ_TIMER_BARRIER_START(timer, mpi, "essprk_step_src_cell_update");
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (prj_timeint_local_block(mpi, block)) {
            int i;
            int j;
            int k;

            PRJ_SUBTIMER_START("sub_cell_src_update");
            prj_src_update(eos, rad, grav, mesh, block, 0, block->W_mhd, block->W_rad,
                block->mhd_rhs, block->rad_rhs);
            PRJ_SUBTIMER_STOP("sub_cell_src_update");
            PRJ_SUBTIMER_START("sub_cell_essprk_step_update");
            for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
                for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                    for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                        prj_timeint_update_cell_stage1_mhd_rad(mesh, rad, eos, block,
                            i, j, k, dt, grav, mpi, dt_src, block->W_mhd, block->W_rad, 0,
                            PRJ_CLIGHT * dt);
                    }
                }
            }
            PRJ_SUBTIMER_STOP("sub_cell_essprk_step_update");
        }
    }
    PRJ_TIMER_BARRIER_STOP(timer, mpi, "essprk_step_src_cell_update");

    if (prj_z4c_runtime_enabled(mesh)) {
        /* Geometry slot 0 was advanced in place, per cell, inside the update
         * loop above (fused with the hydro recovery); only enforce + fill
         * ghosts here for the downstream flux/reconstruction. */
        prj_z4c_finalize_stage(mesh, mpi, bc, 0);
    }

    PRJ_TIMER_BARRIER_START(timer, mpi, "essprk_step_eos_fill_active");
    prj_eos_fill_active_cells(mesh, eos, mpi, 1, PRJ_EOS_CTX_MAIN);
    PRJ_TIMER_BARRIER_STOP(timer, mpi, "essprk_step_eos_fill_active");

    PRJ_TIMER_BARRIER_START(timer, mpi, "essprk_step_ghost_grav_opac");
    prj_boundary_fill_ghosts_and_bf(mesh, mpi, bc, 1, 0, eos, grav, rad,
        PRJ_BOUNDARY_TIMER_SCOPE_NONE);
    PRJ_TIMER_BARRIER_STOP(timer, mpi, "essprk_step_ghost_grav_opac");

    PRJ_TIMER_BARRIER_START(timer, mpi, "essprk_step_eos_fill_mesh");
    prj_eos_fill_mesh(mesh, eos, mpi, 1, PRJ_EOS_CTX_MAIN);
    PRJ_TIMER_BARRIER_STOP(timer, mpi, "essprk_step_eos_fill_mesh");

    prj_flux_fill_transport_opacity_halo(mesh, rad, mpi, 1);
#else
    (void)mesh;
    (void)coord;
    (void)bc;
    (void)eos;
    (void)rad;
    (void)grav;
    (void)mpi;
    (void)dt;
    (void)dt_src;
    (void)timer;
    prj_timeint_unavailable("prj_timeint_eSSPRK_step");
#endif
}

void prj_timeint_eSSPRK_final(prj_mesh *mesh, const prj_coord *coord, const prj_bc *bc,
    prj_eos *eos, prj_rad *rad, prj_grav *grav, prj_mpi *mpi, prj_timer *timer)
{
#if TIME_INTEGRATION == PRJ_TIMEINT_ESSPRK
    double inv_n = 1.0 / (double)PRJ_TIMEINT_ESSPRK_N;

    (void)coord;

    prj_timeint_eSSPRK_blend_with_saved(mesh, bc, eos, rad, grav, mpi, timer,
        1, inv_n,
        "essprk_final_cell_update",
        "essprk_final_eos_fill_active",
        "essprk_final_ghost_grav_opac",
        "essprk_final_eos_fill_mesh");
#else
    (void)mesh;
    (void)coord;
    (void)bc;
    (void)eos;
    (void)rad;
    (void)grav;
    (void)mpi;
    (void)timer;
    prj_timeint_unavailable("prj_timeint_eSSPRK_final");
#endif
}

#if TIME_INTEGRATION == PRJ_TIMEINT_ESSPRK9_3
static void prj_timeint_eSSPRK9_3_advance(prj_mesh *mesh, const prj_coord *coord,
    const prj_bc *bc, prj_eos *eos, prj_rad *rad, prj_grav *grav, prj_mpi *mpi,
    double dt, double *dt_src, prj_timer *timer)
{
    double dt_sub = dt / 6.0;
    int n;

    prj_timeint_eSSPRK_save_state(mesh, mpi, 1);
    prj_timeint_eSSPRK_step(mesh, coord, bc, eos, rad, grav, mpi, dt_sub, dt_src, timer);

    prj_timeint_eSSPRK_save_state(mesh, mpi, 2);
    prj_timeint_eSSPRK_step(mesh, coord, bc, eos, rad, grav, mpi, dt_sub, dt_src, timer);

    prj_timeint_eSSPRK_save_state(mesh, mpi, 3);
    for (n = 0; n < 3; ++n) {
        prj_timeint_eSSPRK_step(mesh, coord, bc, eos, rad, grav, mpi, dt_sub, dt_src, timer);
    }

    prj_timeint_eSSPRK_blend_with_saved(mesh, bc, eos, rad, grav, mpi, timer,
        1, 1.0 / 5.0,
        "essprk9_3_blend1_cell_update",
        "essprk9_3_blend1_eos_fill_active",
        "essprk9_3_blend1_ghost_grav_opac",
        "essprk9_3_blend1_eos_fill_mesh");

    prj_timeint_eSSPRK_step(mesh, coord, bc, eos, rad, grav, mpi, dt_sub, dt_src, timer);

    prj_timeint_eSSPRK_blend_with_saved(mesh, bc, eos, rad, grav, mpi, timer,
        2, 1.0 / 4.0,
        "essprk9_3_blend2_cell_update",
        "essprk9_3_blend2_eos_fill_active",
        "essprk9_3_blend2_ghost_grav_opac",
        "essprk9_3_blend2_eos_fill_mesh");

    prj_timeint_eSSPRK_step(mesh, coord, bc, eos, rad, grav, mpi, dt_sub, dt_src, timer);

    prj_timeint_eSSPRK_blend_with_saved(mesh, bc, eos, rad, grav, mpi, timer,
        3, 1.0 / 3.0,
        "essprk9_3_blend3_cell_update",
        "essprk9_3_blend3_eos_fill_active",
        "essprk9_3_blend3_ghost_grav_opac",
        "essprk9_3_blend3_eos_fill_mesh");

    for (n = 0; n < 2; ++n) {
        prj_timeint_eSSPRK_step(mesh, coord, bc, eos, rad, grav, mpi, dt_sub, dt_src, timer);
    }
}
#endif

void prj_timeint_stage1(prj_mesh *mesh, const prj_coord *coord, const prj_bc *bc, prj_eos *eos,
    prj_rad *rad, prj_grav *grav, prj_mpi *mpi, double dt, double *dt_src, prj_timer *timer)
{
#if TIME_INTEGRATION == RK2
    int bidx;

    (void)coord;
    (void)bc;

    PRJ_TIMER_BARRIER_START(timer, mpi, "stage1_flux_send");
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (prj_timeint_local_block(mpi, block)) {
            prj_flux_update(eos, rad, mesh, block, block->W_mhd, block->eosvar,
                block->flux, 0);
        }
    }
    PRJ_SUBTIMER_START("sub_riemann_flux_send");
    prj_riemann_flux_send(mesh, mpi);
    PRJ_SUBTIMER_STOP("sub_riemann_flux_send");
    PRJ_TIMER_BARRIER_STOP(timer, mpi, "stage1_flux_send");

#if PRJ_MHD
    PRJ_TIMER_BARRIER_START(timer, mpi, "stage1_mhd_emf");
    prj_timeint_mhd_update_mesh_emf(mesh, mpi, 0);
    PRJ_TIMER_BARRIER_STOP(timer, mpi, "stage1_mhd_emf");
#endif

    PRJ_TIMER_BARRIER_START(timer, mpi, "stage1_flux_emf_exchange");
    prj_mpi_exchange_fluxes_and_emf(mesh, mpi);
    PRJ_TIMER_BARRIER_STOP(timer, mpi, "stage1_flux_emf_exchange");

#if PRJ_MHD
    PRJ_TIMER_BARRIER_START(timer, mpi, "stage1_mhd_update_bf");
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (prj_timeint_local_block(mpi, block)) {
            prj_timeint_mhd_update_bf(block, 1, dt);
        }
    }
    PRJ_TIMER_BARRIER_STOP(timer, mpi, "stage1_mhd_update_bf");
#endif

    if (prj_z4c_runtime_enabled(mesh)) {
        double tau_cm = PRJ_CLIGHT * mesh->time_seconds;

        PRJ_SUBTIMER_START("sub_z4c_rhs");
        prj_z4c_compute_rhs(mesh, mpi, rad, 0, 0, tau_cm);
        PRJ_SUBTIMER_STOP("sub_z4c_rhs");
        PRJ_SUBTIMER_START("sub_z4c_sommerfeld");
        prj_z4c_apply_sommerfeld_rhs(mesh, mpi, bc, 0, 0);
        PRJ_SUBTIMER_STOP("sub_z4c_sommerfeld");
    }

    PRJ_TIMER_BARRIER_START(timer, mpi, "stage1_src_cell_update");
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (prj_timeint_local_block(mpi, block)) {
            int i;
            int j;
            int k;

            PRJ_SUBTIMER_START("sub_cell_src_update");
            prj_src_update(eos, rad, grav, mesh, block, 0, block->W_mhd, block->W_rad,
                block->mhd_rhs, block->rad_rhs);
            PRJ_SUBTIMER_STOP("sub_cell_src_update");
            PRJ_SUBTIMER_START("sub_cell_stage1_update");
            for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
                for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                    for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                        prj_timeint_update_cell_stage1_mhd_rad(mesh, rad, eos, block,
                            i, j, k, dt, grav, mpi, dt_src,
                            prj_block_mhd_stage(block, 1), prj_block_rad_stage(block, 1), 1,
                            PRJ_CLIGHT * dt);
                    }
                }
            }
            PRJ_SUBTIMER_STOP("sub_cell_stage1_update");
        }
    }
    PRJ_TIMER_BARRIER_STOP(timer, mpi, "stage1_src_cell_update");

    if (prj_z4c_runtime_enabled(mesh)) {
        /* Z4c slot 1 was advanced in place, per cell, inside the update loop
         * above (fused with the hydro recovery); only enforce + fill ghosts. */
        PRJ_SUBTIMER_START("sub_z4c_finalize");
        prj_z4c_finalize_stage(mesh, mpi, bc, 1);
        PRJ_SUBTIMER_STOP("sub_z4c_finalize");
    }

    PRJ_TIMER_BARRIER_START(timer, mpi, "stage1_eos_fill_active");
    prj_eos_fill_active_cells(mesh, eos, mpi, 2, PRJ_EOS_CTX_MAIN);
    PRJ_TIMER_BARRIER_STOP(timer, mpi, "stage1_eos_fill_active");

    /* The gravity radial reduce/integrate and the active-cell transport opacity
     * run inside the same-level pass of this ghost fill, overlapped with the
     * in-flight exchange. */
    PRJ_TIMER_BARRIER_START(timer, mpi, "stage1_ghost_grav_opac");
    prj_boundary_fill_ghosts_and_bf(mesh, mpi, bc, 2, 1, eos, grav, rad,
        PRJ_BOUNDARY_TIMER_SCOPE_STAGE1);
    PRJ_TIMER_BARRIER_STOP(timer, mpi, "stage1_ghost_grav_opac");

    PRJ_TIMER_BARRIER_START(timer, mpi, "stage1_eos_fill_mesh");
    prj_eos_fill_mesh(mesh, eos, mpi, 2, PRJ_EOS_CTX_MAIN);
    PRJ_TIMER_BARRIER_STOP(timer, mpi, "stage1_eos_fill_mesh");

    /* The 1-ghost halo opacity needs ghost W/eosvar, so fill it now that all
     * ghost filling and eos_fill_mesh are done. */
    prj_flux_fill_transport_opacity_halo(mesh, rad, mpi, 2);
#else
    (void)mesh;
    (void)coord;
    (void)bc;
    (void)eos;
    (void)rad;
    (void)grav;
    (void)mpi;
    (void)dt;
    (void)dt_src;
    (void)timer;
    prj_timeint_unavailable("prj_timeint_stage1");
#endif
}

void prj_timeint_stage2(prj_mesh *mesh, const prj_coord *coord, const prj_bc *bc, prj_eos *eos,
    prj_rad *rad, prj_grav *grav, prj_mpi *mpi, double dt, double *dt_src, prj_timer *timer)
{
#if TIME_INTEGRATION == RK2
    int bidx;

    (void)coord;

    PRJ_TIMER_BARRIER_START(timer, mpi, "stage2_flux_send");
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (prj_timeint_local_block(mpi, block)) {
            prj_flux_update(eos, rad, mesh, block, prj_block_prim_stage(block, 1),
                block->eosvar, block->flux, 1);
        }
    }
    PRJ_SUBTIMER_START("sub_riemann_flux_send");
    prj_riemann_flux_send(mesh, mpi);
    PRJ_SUBTIMER_STOP("sub_riemann_flux_send");
    PRJ_TIMER_BARRIER_STOP(timer, mpi, "stage2_flux_send");

#if PRJ_MHD
    PRJ_TIMER_BARRIER_START(timer, mpi, "stage2_mhd_emf");
    prj_timeint_mhd_update_mesh_emf(mesh, mpi, 1);
    PRJ_TIMER_BARRIER_STOP(timer, mpi, "stage2_mhd_emf");
#endif

    PRJ_TIMER_BARRIER_START(timer, mpi, "stage2_flux_emf_exchange");
    prj_mpi_exchange_fluxes_and_emf(mesh, mpi);
    PRJ_TIMER_BARRIER_STOP(timer, mpi, "stage2_flux_emf_exchange");

#if PRJ_MHD
    PRJ_TIMER_BARRIER_START(timer, mpi, "stage2_mhd_update_bf");
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (prj_timeint_local_block(mpi, block)) {
            prj_timeint_mhd_update_bf(block, 2, dt);
        }
    }
    PRJ_TIMER_BARRIER_STOP(timer, mpi, "stage2_mhd_update_bf");
#endif

    if (prj_z4c_runtime_enabled(mesh)) {
        double tau_cm = PRJ_CLIGHT * mesh->time_seconds;

        PRJ_SUBTIMER_START("sub_z4c_rhs");
        prj_z4c_compute_rhs(mesh, mpi, rad, 1, 0, tau_cm);
        PRJ_SUBTIMER_STOP("sub_z4c_rhs");
        PRJ_SUBTIMER_START("sub_z4c_sommerfeld");
        prj_z4c_apply_sommerfeld_rhs(mesh, mpi, bc, 1, 0);
        PRJ_SUBTIMER_STOP("sub_z4c_sommerfeld");
    }

    PRJ_TIMER_BARRIER_START(timer, mpi, "stage2_src_cell_update");
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (prj_timeint_local_block(mpi, block)) {
            int i;
            int j;
            int k;

            PRJ_SUBTIMER_START("sub_cell_src_update");
            prj_src_update(eos, rad, grav, mesh, block, 1, prj_block_mhd_stage(block, 1),
                prj_block_rad_stage(block, 1), block->mhd_rhs, block->rad_rhs);
            PRJ_SUBTIMER_STOP("sub_cell_src_update");
            PRJ_SUBTIMER_START("sub_cell_stage2_update");
            for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
                for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                    for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                        prj_timeint_update_cell_stage2_mhd_rad(mesh, rad, eos, block,
                            i, j, k, dt, grav, mpi, dt_src, PRJ_CLIGHT * dt);
                    }
                }
            }
            PRJ_SUBTIMER_STOP("sub_cell_stage2_update");
        }
    }
    PRJ_TIMER_BARRIER_STOP(timer, mpi, "stage2_src_cell_update");

    if (prj_z4c_runtime_enabled(mesh)) {
        /* Final geometry slot 0 (g^{n+1}) was advanced in place, per cell,
         * inside the update loop above; only enforce + fill ghosts here. */
        PRJ_SUBTIMER_START("sub_z4c_finalize");
        prj_z4c_finalize_stage(mesh, mpi, bc, 0);
        PRJ_SUBTIMER_STOP("sub_z4c_finalize");
    }

    PRJ_TIMER_BARRIER_START(timer, mpi, "stage2_eos_fill_active");
    prj_eos_fill_active_cells(mesh, eos, mpi, 1, PRJ_EOS_CTX_MAIN);
    PRJ_TIMER_BARRIER_STOP(timer, mpi, "stage2_eos_fill_active");

    /* The gravity radial reduce/integrate and the active-cell transport opacity
     * run inside the same-level pass of this ghost fill, overlapped with the
     * in-flight exchange. */
    PRJ_TIMER_BARRIER_START(timer, mpi, "stage2_ghost_grav_opac");
    prj_boundary_fill_ghosts_and_bf(mesh, mpi, bc, 1, 0, eos, grav, rad,
        PRJ_BOUNDARY_TIMER_SCOPE_STAGE2);
    PRJ_TIMER_BARRIER_STOP(timer, mpi, "stage2_ghost_grav_opac");

    PRJ_TIMER_BARRIER_START(timer, mpi, "stage2_eos_fill_mesh");
    prj_eos_fill_mesh(mesh, eos, mpi, 1, PRJ_EOS_CTX_MAIN);
    PRJ_TIMER_BARRIER_STOP(timer, mpi, "stage2_eos_fill_mesh");

    /* The 1-ghost halo opacity needs ghost W/eosvar, so fill it now that all
     * ghost filling and eos_fill_mesh are done. */
    prj_flux_fill_transport_opacity_halo(mesh, rad, mpi, 1);
#else
    (void)mesh;
    (void)coord;
    (void)bc;
    (void)eos;
    (void)rad;
    (void)grav;
    (void)mpi;
    (void)dt;
    (void)dt_src;
    (void)timer;
    prj_timeint_unavailable("prj_timeint_stage2");
#endif
}

void prj_timeint_step_ex(prj_mesh *mesh, const prj_coord *coord, const prj_bc *bc, prj_eos *eos,
    prj_rad *rad, prj_grav *grav, prj_mpi *mpi, const prj_timeint_imex_tableau *tableau,
    int stage, double dt, double *dt_src, prj_timer *timer)
{
#if TIME_INTEGRATION == PRJ_TIMEINT_IMEX
    int bidx;

    (void)coord;
    (void)bc;
#if PRJ_NRAD == 0
    (void)dt;
#endif

    if (tableau == 0 || stage < 0 || stage >= tableau->nstages) {
        prj_timeint_imex_fail("prj_timeint_step_ex: invalid stage");
    }
#if defined(PRJ_DEBUG)
    if (prj_timeint_imex_debug_tracking_active &&
        (stage >= prj_timeint_imex_debug_nstages ||
         !prj_timeint_imex_debug_stage_solved[stage])) {
        prj_timeint_imex_fail("prj_timeint_step_ex: stage was not solved implicitly");
    }
#endif

    PRJ_TIMER_BARRIER_START(timer, mpi, "imex_step_ex_flux_send");
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (prj_timeint_local_block(mpi, block)) {
            prj_flux_update(eos, rad, mesh, block, prj_block_prim_stage(block, stage),
                block->eosvar, block->flux, stage);
        }
    }
    PRJ_SUBTIMER_START("sub_riemann_flux_send");
    prj_riemann_flux_send(mesh, mpi);
    PRJ_SUBTIMER_STOP("sub_riemann_flux_send");
    PRJ_TIMER_BARRIER_STOP(timer, mpi, "imex_step_ex_flux_send");

#if PRJ_MHD
    PRJ_TIMER_BARRIER_START(timer, mpi, "imex_step_ex_mhd_emf");
    prj_timeint_mhd_update_mesh_emf_stage(mesh, mpi, stage);
    PRJ_TIMER_BARRIER_STOP(timer, mpi, "imex_step_ex_mhd_emf");
#endif

    PRJ_TIMER_BARRIER_START(timer, mpi, "imex_step_ex_flux_emf_exchange");
    prj_mpi_exchange_fluxes_and_emf(mesh, mpi);
    PRJ_TIMER_BARRIER_STOP(timer, mpi, "imex_step_ex_flux_emf_exchange");

#if PRJ_MHD
    PRJ_TIMER_BARRIER_START(timer, mpi, "imex_step_ex_bf_deriv");
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (prj_timeint_local_block(mpi, block)) {
            prj_timeint_mhd_save_deriv_bf(block, stage);
        }
    }
    PRJ_TIMER_BARRIER_STOP(timer, mpi, "imex_step_ex_bf_deriv");
#endif

    PRJ_TIMER_BARRIER_START(timer, mpi, "imex_step_ex_src_deriv");
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        double *W_stage;
        double *W_rad_stage;
        double *deriv;
        int i;
        int j;
        int k;

        if (!prj_timeint_local_block(mpi, block)) {
            continue;
        }
        W_stage = prj_block_prim_stage(block, stage);
        W_rad_stage = prj_block_rad_stage(block, stage);
        deriv = prj_block_deriv_stage(block->deriv_ex, stage);
        if (W_stage == 0 || deriv == 0) {
            prj_timeint_imex_fail("prj_timeint_step_ex: missing derivative storage");
        }
        PRJ_SUBTIMER_START("sub_cell_src_update");
        prj_src_update(eos, rad, grav, mesh, block, stage, W_stage, W_rad_stage,
            block->mhd_rhs, block->rad_rhs);
        PRJ_SUBTIMER_STOP("sub_cell_src_update");
        PRJ_SUBTIMER_START("sub_cell_ex_deriv");
        for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
            for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    double fluxdiv[PRJ_NVAR_CONS];
                    double u[PRJ_NVAR_CONS];
                    int v;

                    prj_flux_div(block->flux, block->area, block->vol, i, j, k, fluxdiv);
                    prj_timeint_imex_cons_from_stage(eos, mesh, block, stage, i, j, k, u);
                    prj_timeint_update_dt_src_values(mesh, grav, block,
                        u[PRJ_CONS_RHO], u[PRJ_CONS_MOM1], u[PRJ_CONS_MOM2],
                        u[PRJ_CONS_MOM3], u[PRJ_CONS_ETOT], i, j, k, mpi, dt_src);
                    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
                        deriv[VIDX(v, i, j, k)] =
                            prj_block_rhs_value_const(block, v, i, j, k) + fluxdiv[v];
                    }
#if PRJ_NRAD > 0
                    prj_timeint_imex_add_explicit_rad_deriv(eos, rad, mesh, block,
                        stage, i, j, k, dt, deriv);
#else
                    (void)rad;
#endif
                }
            }
        }
        PRJ_SUBTIMER_STOP("sub_cell_ex_deriv");
    }
    PRJ_TIMER_BARRIER_STOP(timer, mpi, "imex_step_ex_src_deriv");
#else
    (void)mesh;
    (void)coord;
    (void)bc;
    (void)eos;
    (void)rad;
    (void)grav;
    (void)mpi;
    (void)tableau;
    (void)stage;
    (void)dt;
    (void)dt_src;
    (void)timer;
    prj_timeint_unavailable("prj_timeint_step_ex");
#endif
}

void prj_timeint_step_im(prj_mesh *mesh, const prj_coord *coord, const prj_bc *bc, prj_eos *eos,
    prj_rad *rad, prj_grav *grav, prj_mpi *mpi, const prj_timeint_imex_tableau *tableau,
    int stage, double dt, double dt_implicit, double *dt_src, prj_timer *timer)
{
#if TIME_INTEGRATION == PRJ_TIMEINT_IMEX
    int bidx;

    (void)coord;
    (void)bc;
    (void)grav;
    (void)dt;
    (void)dt_src;

    if (tableau == 0 || stage < 0 || stage >= tableau->nstages) {
        prj_timeint_imex_fail("prj_timeint_step_im: invalid stage");
    }

    PRJ_TIMER_BARRIER_START(timer, mpi, "imex_step_im_deriv");
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        double *deriv;
        int i;
        int j;
        int k;

        if (!prj_timeint_local_block(mpi, block)) {
            continue;
        }
        deriv = prj_block_deriv_stage(block->deriv_im, stage);
        if (prj_block_prim_stage(block, stage) == 0 || deriv == 0) {
            prj_timeint_imex_fail("prj_timeint_step_im: missing derivative storage");
        }
        for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
            for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    double u0[PRJ_NVAR_CONS];
                    double u1[PRJ_NVAR_CONS];
                    int v;

                    prj_timeint_imex_cons_from_stage(eos, mesh, block, stage, i, j, k, u0);
                    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
                        u1[v] = u0[v];
                        deriv[VIDX(v, i, j, k)] = 0.0;
                    }
#if PRJ_NRAD > 0
                    if (dt_implicit != 0.0) {
                        double lapse_cell = prj_timeint_cell_lapse(block, i, j, k);
#if PRJ_USE_RADIATION_FSA
                        PRJ_SUBTIMER_START("sub_rad_energy_momentum");
                        prj_rad_energy_momentum_update_fsa(rad, block, i, j, k, eos, u1,
                            dt_implicit, lapse_cell);
                        PRJ_SUBTIMER_STOP("sub_rad_energy_momentum");
#if DO_FFC
                        PRJ_SUBTIMER_START("sub_rad_ffc");
                        prj_rad_ffc_fsa(rad, u1, dt_implicit);
                        PRJ_SUBTIMER_STOP("sub_rad_ffc");
#endif
#else
                        double T_cell = 0.0;
                        double kappa[PRJ_NRAD * PRJ_NEGROUP];

                        PRJ_SUBTIMER_START("sub_rad_energy");
                        prj_rad_energy_update(rad, eos, u1, dt_implicit, lapse_cell,
                            &T_cell, kappa);
                        PRJ_SUBTIMER_STOP("sub_rad_energy");
                        PRJ_SUBTIMER_START("sub_rad_momentum");
                        prj_rad_momentum_update(rad, eos, u1, dt_implicit, lapse_cell,
                            T_cell, kappa);
                        PRJ_SUBTIMER_STOP("sub_rad_momentum");
#endif
                        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
                            deriv[VIDX(v, i, j, k)] = (u1[v] - u0[v]) / dt_implicit;
                        }
                    }
#else
                    (void)rad;
                    (void)dt_implicit;
#endif
                    prj_timeint_imex_store_cons_to_stage(eos, mesh, block, stage,
                        i, j, k, u1);
                }
            }
        }
    }
    PRJ_TIMER_BARRIER_STOP(timer, mpi, "imex_step_im_deriv");
#if defined(PRJ_DEBUG)
    if (prj_timeint_imex_debug_tracking_active &&
        stage >= 0 && stage < PRJ_BLOCK_NSTAGES) {
        prj_timeint_imex_debug_stage_solved[stage] = 1;
    }
#endif
#else
    (void)mesh;
    (void)coord;
    (void)bc;
    (void)eos;
    (void)rad;
    (void)grav;
    (void)mpi;
    (void)tableau;
    (void)stage;
    (void)dt;
    (void)dt_implicit;
    (void)dt_src;
    (void)timer;
    prj_timeint_unavailable("prj_timeint_step_im");
#endif
}

void prj_timeint_step(prj_mesh *mesh, const prj_coord *coord, const prj_bc *bc, prj_eos *eos,
    prj_rad *rad, prj_grav *grav, prj_mpi *mpi, const prj_timeint_imex_tableau *tableau,
    double dt, double *dt_src, prj_timer *timer)
{
#if TIME_INTEGRATION == RK2
    (void)tableau;
    prj_timeint_stage1(mesh, coord, bc, eos, rad, grav, mpi, dt, dt_src, timer);
    prj_timeint_stage2(mesh, coord, bc, eos, rad, grav, mpi, dt, dt_src, timer);
#elif TIME_INTEGRATION == PRJ_TIMEINT_ESSPRK
    int n;
    double dt_sub = dt / (double)(PRJ_TIMEINT_ESSPRK_N - 1);

    (void)tableau;
    prj_timeint_eSSPRK_save_state(mesh, mpi, 1);
    for (n = 0; n < PRJ_TIMEINT_ESSPRK_N; ++n) {
        prj_timeint_eSSPRK_step(mesh, coord, bc, eos, rad, grav, mpi, dt_sub, dt_src, timer);
    }
    prj_timeint_eSSPRK_final(mesh, coord, bc, eos, rad, grav, mpi, timer);
#elif TIME_INTEGRATION == PRJ_TIMEINT_ESSPRK9_3
    (void)tableau;
    prj_timeint_eSSPRK9_3_advance(mesh, coord, bc, eos, rad, grav, mpi, dt, dt_src, timer);
#elif TIME_INTEGRATION == PRJ_TIMEINT_IMEX
    int stage;

    if (tableau == 0) {
        prj_timeint_imex_fail("prj_timeint_step: missing tableau");
    }
    prj_timeint_imex_save_base_cons(mesh, eos, mpi);
#if defined(PRJ_DEBUG)
    prj_timeint_imex_debug_tracking_active = 1;
    prj_timeint_imex_debug_nstages = tableau->nstages;
    for (stage = 0; stage < PRJ_BLOCK_NSTAGES; ++stage) {
        prj_timeint_imex_debug_stage_solved[stage] = 0;
    }
#endif
    for (stage = 0; stage < tableau->nstages; ++stage) {
        double diag = prj_timeint_imex_a_im(tableau, stage, stage);
        double dt_implicit = dt * diag;

        prj_timeint_imex_assemble_stage(mesh, eos, mpi, tableau,
            stage, stage, 0, stage, dt);
        prj_gravity_monopole_update_lapse_active(mesh, eos, rad, grav, mpi, stage);
        prj_timeint_step_im(mesh, coord, bc, eos, rad, grav, mpi, tableau,
            stage, dt, dt_implicit, dt_src, timer);
        prj_timeint_imex_fill_updated_stage(mesh, bc, eos, rad, grav, mpi, timer,
            stage, 1);
        prj_timeint_step_ex(mesh, coord, bc, eos, rad, grav, mpi, tableau,
            stage, dt, dt_src, timer);
    }
#if defined(PRJ_DEBUG)
    prj_timeint_imex_debug_tracking_active = 0;
#endif
    prj_timeint_imex_assemble_stage(mesh, eos, mpi, tableau,
        0, 0, 1, tableau->nstages, dt);
    prj_timeint_imex_fill_updated_stage(mesh, bc, eos, rad, grav, mpi, timer,
        0, 0);
#else
#error "Unsupported TIME_INTEGRATION value"
#endif
}
