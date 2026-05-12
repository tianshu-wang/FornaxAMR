#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "prj.h"
#include "prj_rad_inel.h"

#ifndef PRJ_FAST_MHD_RAD_CELL_UPDATE
#define PRJ_FAST_MHD_RAD_CELL_UPDATE 1
#endif

#if PRJ_NRAD > 0
static double prj_timeint_cell_lapse(const prj_block *block, int i, int j, int k)
{
    return prj_gravity_block_lapse_at(block, i, j, k);
}
#endif

static int prj_timeint_local_block(const prj_block *block)
{
    prj_mpi *mpi = prj_mpi_current();

    return block != 0 && block->id >= 0 && block->active == 1 &&
        (mpi == 0 || block->rank == mpi->rank);
}

static void prj_timeint_cell_prim(const double *src, int i, int j, int k, double *w)
{
    int v;

    for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
        w[v] = src[VIDX(v, i, j, k)];
    }
}

#if !(PRJ_MHD && PRJ_NRAD > 0 && PRJ_FAST_MHD_RAD_CELL_UPDATE)
static void prj_timeint_cell_cons_store(double *dst, int i, int j, int k, const double *u)
{
    int v;

    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        dst[VIDX(v, i, j, k)] = u[v];
    }
}

static void prj_timeint_cell_prim_store(double *dst, int i, int j, int k, const double *w)
{
    int v;

    for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
        dst[VIDX(v, i, j, k)] = w[v];
    }
}

static void prj_timeint_apply_eint_floor(const prj_mesh *mesh, double *u, double *w)
{
    double rho;
    double kinetic;
#if PRJ_MHD
    double magnetic;
#endif

    if (mesh == 0 || u == 0 || w == 0 || mesh->E_floor <= 0.0) {
        return;
    }

    rho = w[PRJ_PRIM_RHO];
    if (rho <= 0.0 || w[PRJ_PRIM_EINT] >= mesh->E_floor) {
        return;
    }

    kinetic = 0.5 * (w[PRJ_PRIM_V1] * w[PRJ_PRIM_V1] +
        w[PRJ_PRIM_V2] * w[PRJ_PRIM_V2] +
        w[PRJ_PRIM_V3] * w[PRJ_PRIM_V3]);
#if PRJ_MHD
    magnetic = 0.5 * (w[PRJ_PRIM_B1] * w[PRJ_PRIM_B1] +
        w[PRJ_PRIM_B2] * w[PRJ_PRIM_B2] +
        w[PRJ_PRIM_B3] * w[PRJ_PRIM_B3]);
#endif
    w[PRJ_PRIM_EINT] = mesh->E_floor;
    u[PRJ_CONS_ETOT] = rho * (mesh->E_floor + kinetic)
#if PRJ_MHD
        + magnetic
#endif
        ;
}
#endif

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

static double prj_timeint_mhd_cell_emf(const double *W, int dir, int i, int j, int k)
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
    b1 = W[VIDX(PRJ_PRIM_B1, i, j, k)];
    b2 = W[VIDX(PRJ_PRIM_B2, i, j, k)];
    b3 = W[VIDX(PRJ_PRIM_B3, i, j, k)];
    v1 = W[VIDX(PRJ_PRIM_V1, i, j, k)];
    v2 = W[VIDX(PRJ_PRIM_V2, i, j, k)];
    v3 = W[VIDX(PRJ_PRIM_V3, i, j, k)];
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
    value = block->v_riemann[face_dir][face_dir * PRJ_BLOCK_NCELLS + IDX(i, j, k)];
    if (!isfinite(value)) {
        prj_timeint_mhd_fail("prj_timeint_mhd_face_vnorm: non-finite face velocity");
    }
    return value;
}

static void prj_timeint_mhd_update_emf(prj_block *block, double *W)
{
    int dir;

    if (block == 0 || W == 0) {
        prj_timeint_mhd_fail("prj_timeint_mhd_update_emf: invalid input");
    }
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
                    emf_cell[0] = prj_timeint_mhd_cell_emf(W, dir,
                        cell0[0], cell0[1], cell0[2]);
                    emf_cell[1] = prj_timeint_mhd_cell_emf(W, dir,
                        cell1[0], cell1[1], cell1[2]);
                    emf_cell[2] = prj_timeint_mhd_cell_emf(W, dir,
                        cell2[0], cell2[1], cell2[2]);
                    emf_cell[3] = prj_timeint_mhd_cell_emf(W, dir,
                        cell3[0], cell3[1], cell3[2]);
                    v_norm[0] = prj_timeint_mhd_face_vnorm(block, right,
                        face0[0], face0[1], face0[2]);
                    v_norm[1] = prj_timeint_mhd_face_vnorm(block, up,
                        face1[0], face1[1], face1[2]);
                    v_norm[2] = prj_timeint_mhd_face_vnorm(block, right,
                        face2[0], face2[1], face2[2]);
                    v_norm[3] = prj_timeint_mhd_face_vnorm(block, up,
                        face3[0], face3[1], face3[2]);
                    prj_mhd_emf_upwind(block, dir, e[0], e[1], e[2],
                        emf_face, emf_cell, v_norm);
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

static void prj_timeint_mhd_update_bf(prj_block *block, int stage, double dt)
{
    double *src[3];
    double *dst[3];
    double *old[3];
    int dir;

    if (block == 0 || (stage != 1 && stage != 2)) {
        prj_timeint_mhd_fail("prj_timeint_mhd_update_bf: invalid input");
    }
    for (dir = 0; dir < 3; ++dir) {
        if (block->Bf[dir] == 0 || block->Bf1[dir] == 0) {
            prj_timeint_mhd_fail("prj_timeint_mhd_update_bf: missing Bf storage");
        }
        old[dir] = block->Bf[dir];
        if (stage == 1) {
            src[dir] = block->Bf[dir];
            dst[dir] = block->Bf1[dir];
        } else {
            src[dir] = block->Bf1[dir];
            dst[dir] = block->Bf[dir];
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

static void prj_timeint_mhd_update_mesh_emf(prj_mesh *mesh, double *(*stage_array)(prj_block *),
    prj_timer *timer)
{
    int bidx;

    if (mesh == 0 || stage_array == 0) {
        prj_timeint_mhd_fail("prj_timeint_mhd_update_mesh_emf: invalid input");
    }
    PRJ_TIMER_START(timer, "mhd_emf_local");
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (prj_timeint_local_block(block)) {
            prj_timeint_mhd_update_emf(block, stage_array(block));
        }
    }
    PRJ_TIMER_STOP(timer, "mhd_emf_local");
    PRJ_TIMER_START(timer, "mhd_emf_send");
    prj_mhd_emf_send(mesh);
    PRJ_TIMER_STOP(timer, "mhd_emf_send");
}

static double *prj_timeint_stage1_array(prj_block *block)
{
    return block != 0 ? block->W : 0;
}

static double *prj_timeint_stage2_array(prj_block *block)
{
    return block != 0 ? block->W1 : 0;
}
#endif

static void prj_timeint_update_dt_src_values(const prj_block *block,
    double rho, double mom1, double mom2, double mom3, double etot,
    int i, int j, int k, double *dt_src)
{
    double dt_src_local;

    if (dt_src == 0) {
        return;
    }

    dt_src_local = 0.02 * etot /
        (fabs(block->dUdt[VIDX(PRJ_CONS_ETOT, i, j, k)]) + 1.0e-50);
    if (dt_src_local < *dt_src) {
        *dt_src = dt_src_local;
    }
#if !PRJ_GRAV_DEBUG
    (void)rho;
    (void)mom1;
    (void)mom2;
    (void)mom3;
#endif
#if PRJ_GRAV_DEBUG
    if (dt_src_local<1e-7) {
        prj_mpi *mpi = prj_mpi_current();
        int rank = mpi != 0 ? mpi->rank : 0;
        double inv_rho = rho != 0.0 ? 1.0 / rho : 0.0;
        double x1 = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
        double x2 = block->xmin[1] + ((double)j + 0.5) * block->dx[1];
        double x3 = block->xmin[2] + ((double)k + 0.5) * block->dx[2];
        double accel = prj_gravity_block_accel_at(block, i, j, k);

        fprintf(stderr,
            "[grav debug] shortest dt_src rank=%d block=%d cell=(%d,%d,%d) "
            "x=(%.17e,%.17e,%.17e) rho=%.17e etot=%.17e "
            "v=(%.17e,%.17e,%.17e) accel=%.17e dt_src=%.17e\n",
            rank, block->id, i, j, k, x1, x2, x3, rho, etot,
            mom1 * inv_rho, mom2 * inv_rho, mom3 * inv_rho, accel, dt_src_local);
    }
#endif
}

#if !(PRJ_MHD && PRJ_NRAD > 0 && PRJ_FAST_MHD_RAD_CELL_UPDATE)
static void prj_timeint_update_dt_src(const prj_block *block, const double *u,
    int i, int j, int k, double *dt_src)
{
    if (u == 0) {
        return;
    }
    prj_timeint_update_dt_src_values(block, u[PRJ_CONS_RHO], u[PRJ_CONS_MOM1],
        u[PRJ_CONS_MOM2], u[PRJ_CONS_MOM3], u[PRJ_CONS_ETOT], i, j, k, dt_src);
}
#endif

#if PRJ_MHD && PRJ_NRAD > 0 && PRJ_FAST_MHD_RAD_CELL_UPDATE
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
    double r = W[VIDX(PRJ_PRIM_RHO, i, j, k)];
    double v1 = W[VIDX(PRJ_PRIM_V1, i, j, k)];
    double v2 = W[VIDX(PRJ_PRIM_V2, i, j, k)];
    double v3 = W[VIDX(PRJ_PRIM_V3, i, j, k)];
    double eint = W[VIDX(PRJ_PRIM_EINT, i, j, k)];
    double cb1 = W[VIDX(PRJ_PRIM_B1, i, j, k)];
    double cb2 = W[VIDX(PRJ_PRIM_B2, i, j, k)];
    double cb3 = W[VIDX(PRJ_PRIM_B3, i, j, k)];

    *rho = r;
    *mom1 = r * v1;
    *mom2 = r * v2;
    *mom3 = r * v3;
    *etot = r * eint + 0.5 * r * (v1 * v1 + v2 * v2 + v3 * v3) +
        0.5 * (cb1 * cb1 + cb2 * cb2 + cb3 * cb3);
    *ye = r * W[VIDX(PRJ_PRIM_YE, i, j, k)];
    *b1 = cb1;
    *b2 = cb2;
    *b3 = cb3;
}

static void prj_timeint_cell_cons_from_prim_mhd_rad(const double *W,
    int i, int j, int k, double *u)
{
    int field;
    int group;

    prj_timeint_mhd_hydro_cons_from_prim(W, i, j, k,
        &u[PRJ_CONS_RHO], &u[PRJ_CONS_MOM1], &u[PRJ_CONS_MOM2],
        &u[PRJ_CONS_MOM3], &u[PRJ_CONS_ETOT], &u[PRJ_CONS_YE],
        &u[PRJ_CONS_B1], &u[PRJ_CONS_B2], &u[PRJ_CONS_B3]);

    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            u[PRJ_CONS_RAD_E(field, group)] =
                W[VIDX(PRJ_PRIM_RAD_E(field, group), i, j, k)];
            u[PRJ_CONS_RAD_F1(field, group)] =
                W[VIDX(PRJ_PRIM_RAD_F1(field, group), i, j, k)];
            u[PRJ_CONS_RAD_F2(field, group)] =
                W[VIDX(PRJ_PRIM_RAD_F2(field, group), i, j, k)];
            u[PRJ_CONS_RAD_F3(field, group)] =
                W[VIDX(PRJ_PRIM_RAD_F3(field, group), i, j, k)];
        }
    }
}

static void prj_timeint_store_mhd_rad_cell(const prj_mesh *mesh, prj_block *block,
    double *Wdst, int i, int j, int k, double *u)
{
    double rho = u[PRJ_CONS_RHO];
    int field;
    int group;

    if (rho == 0.0) {
        Wdst[VIDX(PRJ_PRIM_RHO, i, j, k)] = 0.0;
        Wdst[VIDX(PRJ_PRIM_V1, i, j, k)] = 0.0;
        Wdst[VIDX(PRJ_PRIM_V2, i, j, k)] = 0.0;
        Wdst[VIDX(PRJ_PRIM_V3, i, j, k)] = 0.0;
        Wdst[VIDX(PRJ_PRIM_EINT, i, j, k)] = 0.0;
        Wdst[VIDX(PRJ_PRIM_YE, i, j, k)] = 0.0;
        Wdst[VIDX(PRJ_PRIM_B1, i, j, k)] = 0.0;
        Wdst[VIDX(PRJ_PRIM_B2, i, j, k)] = 0.0;
        Wdst[VIDX(PRJ_PRIM_B3, i, j, k)] = 0.0;
    } else {
        double v1 = u[PRJ_CONS_MOM1] / rho;
        double v2 = u[PRJ_CONS_MOM2] / rho;
        double v3 = u[PRJ_CONS_MOM3] / rho;
        double kinetic = 0.5 * (v1 * v1 + v2 * v2 + v3 * v3);
        double magnetic_density = 0.5 * (u[PRJ_CONS_B1] * u[PRJ_CONS_B1] +
            u[PRJ_CONS_B2] * u[PRJ_CONS_B2] + u[PRJ_CONS_B3] * u[PRJ_CONS_B3]);
        double magnetic = magnetic_density / rho;
        double eint = u[PRJ_CONS_ETOT] / rho - kinetic - magnetic;

        if (mesh != 0 && mesh->E_floor > 0.0 && eint < mesh->E_floor) {
            eint = mesh->E_floor;
            u[PRJ_CONS_ETOT] = rho * (mesh->E_floor + kinetic) + magnetic_density;
        }

        Wdst[VIDX(PRJ_PRIM_RHO, i, j, k)] = rho;
        Wdst[VIDX(PRJ_PRIM_V1, i, j, k)] = v1;
        Wdst[VIDX(PRJ_PRIM_V2, i, j, k)] = v2;
        Wdst[VIDX(PRJ_PRIM_V3, i, j, k)] = v3;
        Wdst[VIDX(PRJ_PRIM_EINT, i, j, k)] = eint;
        Wdst[VIDX(PRJ_PRIM_YE, i, j, k)] = u[PRJ_CONS_YE] / rho;
        Wdst[VIDX(PRJ_PRIM_B1, i, j, k)] = u[PRJ_CONS_B1];
        Wdst[VIDX(PRJ_PRIM_B2, i, j, k)] = u[PRJ_CONS_B2];
        Wdst[VIDX(PRJ_PRIM_B3, i, j, k)] = u[PRJ_CONS_B3];
    }

    block->U[VIDX(PRJ_CONS_RHO, i, j, k)] = u[PRJ_CONS_RHO];
    block->U[VIDX(PRJ_CONS_MOM1, i, j, k)] = u[PRJ_CONS_MOM1];
    block->U[VIDX(PRJ_CONS_MOM2, i, j, k)] = u[PRJ_CONS_MOM2];
    block->U[VIDX(PRJ_CONS_MOM3, i, j, k)] = u[PRJ_CONS_MOM3];
    block->U[VIDX(PRJ_CONS_ETOT, i, j, k)] = u[PRJ_CONS_ETOT];
    block->U[VIDX(PRJ_CONS_YE, i, j, k)] = u[PRJ_CONS_YE];
    block->U[VIDX(PRJ_CONS_B1, i, j, k)] = u[PRJ_CONS_B1];
    block->U[VIDX(PRJ_CONS_B2, i, j, k)] = u[PRJ_CONS_B2];
    block->U[VIDX(PRJ_CONS_B3, i, j, k)] = u[PRJ_CONS_B3];

    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            int e = PRJ_CONS_RAD_E(field, group);
            int f1 = PRJ_CONS_RAD_F1(field, group);
            int f2 = PRJ_CONS_RAD_F2(field, group);
            int f3 = PRJ_CONS_RAD_F3(field, group);

            block->U[VIDX(e, i, j, k)] = u[e];
            block->U[VIDX(f1, i, j, k)] = u[f1];
            block->U[VIDX(f2, i, j, k)] = u[f2];
            block->U[VIDX(f3, i, j, k)] = u[f3];
            Wdst[VIDX(e, i, j, k)] = u[e];
            Wdst[VIDX(f1, i, j, k)] = u[f1];
            Wdst[VIDX(f2, i, j, k)] = u[f2];
            Wdst[VIDX(f3, i, j, k)] = u[f3];
        }
    }
}

static void prj_timeint_update_cell_stage1_mhd_rad(const prj_mesh *mesh,
    prj_rad *rad, prj_eos *eos, prj_block *block, int i, int j, int k,
    double dt, double *dt_src)
{
    double u[PRJ_NVAR_CONS];
    int field;
    int group;

    prj_timeint_cell_cons_from_prim_mhd_rad(block->W, i, j, k, u);
    prj_timeint_update_dt_src_values(block, u[PRJ_CONS_RHO], u[PRJ_CONS_MOM1],
        u[PRJ_CONS_MOM2], u[PRJ_CONS_MOM3], u[PRJ_CONS_ETOT], i, j, k, dt_src);

#define PRJ_STAGE1_UPDATE(v) do { \
        int vv = (v); \
        u[vv] += dt * (block->dUdt[VIDX(vv, i, j, k)] + \
            prj_timeint_flux_div_var(block, vv, i, j, k)); \
    } while (0)

    PRJ_STAGE1_UPDATE(PRJ_CONS_RHO);
    PRJ_STAGE1_UPDATE(PRJ_CONS_MOM1);
    PRJ_STAGE1_UPDATE(PRJ_CONS_MOM2);
    PRJ_STAGE1_UPDATE(PRJ_CONS_MOM3);
    PRJ_STAGE1_UPDATE(PRJ_CONS_ETOT);
    PRJ_STAGE1_UPDATE(PRJ_CONS_YE);

    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            PRJ_STAGE1_UPDATE(PRJ_CONS_RAD_E(field, group));
            PRJ_STAGE1_UPDATE(PRJ_CONS_RAD_F1(field, group));
            PRJ_STAGE1_UPDATE(PRJ_CONS_RAD_F2(field, group));
            PRJ_STAGE1_UPDATE(PRJ_CONS_RAD_F3(field, group));
        }
    }

#undef PRJ_STAGE1_UPDATE

    prj_timeint_mhd_set_cons_b_from_bf(block, block->Bf1, i, j, k, u);
    {
        double T_cell = block->eosvar[EIDX(PRJ_EOSVAR_TEMPERATURE, i, j, k)];
        double lapse_cell = prj_timeint_cell_lapse(block, i, j, k);

        prj_rad_freq_flux_apply(rad, block, block->W, u, i, j, k, lapse_cell, dt);
        PRJ_TIMER_START(prj_timer_current(), "rad_nucinel");
        prj_rad_nucinel_step(rad, eos, u, dt, T_cell);
        PRJ_TIMER_STOP(prj_timer_current(), "rad_nucinel");
        PRJ_TIMER_START(prj_timer_current(), "rad_eleinel");
        prj_rad_eleinel_step(rad, eos, u, dt, T_cell);
        PRJ_TIMER_STOP(prj_timer_current(), "rad_eleinel");
        prj_rad_energy_update(rad, eos, u, dt, lapse_cell, &T_cell);
        prj_rad_momentum_update(rad, eos, u, dt, lapse_cell, T_cell);
    }
    prj_timeint_store_mhd_rad_cell(mesh, block, block->W1, i, j, k, u);
}

static void prj_timeint_update_cell_stage2_mhd_rad(const prj_mesh *mesh,
    prj_rad *rad, prj_eos *eos, prj_block *block, int i, int j, int k,
    double dt, double *dt_src)
{
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
    int field;
    int group;

    prj_timeint_mhd_hydro_cons_from_prim(block->W, i, j, k,
        &rho0, &mom10, &mom20, &mom30, &etot0, &ye0, &b10, &b20, &b30);
    prj_timeint_mhd_hydro_cons_from_prim(block->W1, i, j, k,
        &rho1, &mom11, &mom21, &mom31, &etot1, &ye1, &b11, &b21, &b31);
    prj_timeint_update_dt_src_values(block, rho0, mom10, mom20, mom30, etot0,
        i, j, k, dt_src);

#define PRJ_STAGE2_UPDATE(v, u0v, u1v) do { \
        int vv = (v); \
        u[vv] = 0.5 * (u0v) + 0.5 * ((u1v) + dt * \
            (block->dUdt[VIDX(vv, i, j, k)] + \
                prj_timeint_flux_div_var(block, vv, i, j, k))); \
    } while (0)

    PRJ_STAGE2_UPDATE(PRJ_CONS_RHO, rho0, rho1);
    PRJ_STAGE2_UPDATE(PRJ_CONS_MOM1, mom10, mom11);
    PRJ_STAGE2_UPDATE(PRJ_CONS_MOM2, mom20, mom21);
    PRJ_STAGE2_UPDATE(PRJ_CONS_MOM3, mom30, mom31);
    PRJ_STAGE2_UPDATE(PRJ_CONS_ETOT, etot0, etot1);
    PRJ_STAGE2_UPDATE(PRJ_CONS_YE, ye0, ye1);

    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            int e = PRJ_CONS_RAD_E(field, group);
            int f1 = PRJ_CONS_RAD_F1(field, group);
            int f2 = PRJ_CONS_RAD_F2(field, group);
            int f3 = PRJ_CONS_RAD_F3(field, group);

            PRJ_STAGE2_UPDATE(e, block->W[VIDX(e, i, j, k)], block->W1[VIDX(e, i, j, k)]);
            PRJ_STAGE2_UPDATE(f1, block->W[VIDX(f1, i, j, k)], block->W1[VIDX(f1, i, j, k)]);
            PRJ_STAGE2_UPDATE(f2, block->W[VIDX(f2, i, j, k)], block->W1[VIDX(f2, i, j, k)]);
            PRJ_STAGE2_UPDATE(f3, block->W[VIDX(f3, i, j, k)], block->W1[VIDX(f3, i, j, k)]);
        }
    }

#undef PRJ_STAGE2_UPDATE

    (void)b10;
    (void)b20;
    (void)b30;
    (void)b11;
    (void)b21;
    (void)b31;
    prj_timeint_mhd_set_cons_b_from_bf(block, block->Bf, i, j, k, u);
    {
        double T_cell = block->eosvar[EIDX(PRJ_EOSVAR_TEMPERATURE, i, j, k)];
        double lapse_cell = prj_timeint_cell_lapse(block, i, j, k);

        prj_rad_freq_flux_apply(rad, block, block->W1, u, i, j, k, lapse_cell, 0.5 * dt);
        PRJ_TIMER_START(prj_timer_current(), "rad_nucinel");
        prj_rad_nucinel_step(rad, eos, u, 0.5 * dt, T_cell);
        PRJ_TIMER_STOP(prj_timer_current(), "rad_nucinel");
        PRJ_TIMER_START(prj_timer_current(), "rad_eleinel");
        prj_rad_eleinel_step(rad, eos, u, 0.5 * dt, T_cell);
        PRJ_TIMER_STOP(prj_timer_current(), "rad_eleinel");
        prj_rad_energy_update(rad, eos, u, 0.5 * dt, lapse_cell, &T_cell);
        prj_rad_momentum_update(rad, eos, u, 0.5 * dt, lapse_cell, T_cell);
    }
    prj_timeint_store_mhd_rad_cell(mesh, block, block->W, i, j, k, u);
}
#endif

#if PRJ_NRAD > 0
static double prj_timeint_cell_rad_denom(const double *w, const double dx[3])
{
    double cdir[3] = {0.0, 0.0, 0.0};
    int field;
    int group;

    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            double E = w[PRJ_PRIM_RAD_E(field, group)];
            double F1 = w[PRJ_PRIM_RAD_F1(field, group)];
            double F2 = w[PRJ_PRIM_RAD_F2(field, group)];
            double F3 = w[PRJ_PRIM_RAD_F3(field, group)];
            double lam_min;
            double lam_max;
            double c_abs;

            /* prj_rad_m1_wavespeeds treats F1 as the direction-normal flux. */
            prj_rad_m1_wavespeeds(E, F1, F2, F3, &lam_min, &lam_max);
            c_abs = fabs(lam_min);
            if (fabs(lam_max) > c_abs) {
                c_abs = fabs(lam_max);
            }
            c_abs *= PRJ_CLIGHT;
            if (c_abs > cdir[0]) {
                cdir[0] = c_abs;
            }

            prj_rad_m1_wavespeeds(E, F2, F3, F1, &lam_min, &lam_max);
            c_abs = fabs(lam_min);
            if (fabs(lam_max) > c_abs) {
                c_abs = fabs(lam_max);
            }
            c_abs *= PRJ_CLIGHT;
            if (c_abs > cdir[1]) {
                cdir[1] = c_abs;
            }

            prj_rad_m1_wavespeeds(E, F3, F1, F2, &lam_min, &lam_max);
            c_abs = fabs(lam_min);
            if (fabs(lam_max) > c_abs) {
                c_abs = fabs(lam_max);
            }
            c_abs *= PRJ_CLIGHT;
            if (c_abs > cdir[2]) {
                cdir[2] = c_abs;
            }
        }
    }

    return cdir[0] / dx[0] + cdir[1] / dx[1] + cdir[2] / dx[2];
}
#endif

double prj_timeint_calc_dt(const prj_mesh *mesh, prj_eos *eos, double cfl)
{
    double dt_min = 1.0e99;
    int bidx;

    (void)eos;

    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        const prj_block *block = &mesh->blocks[bidx];
        int i;
        int j;
        int k;

        if (!prj_timeint_local_block(block)) {
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

                    prj_timeint_cell_prim(block->W, i, j, k, w);
                    q[PRJ_EOS_PRESSURE] = block->eosvar[EIDX(PRJ_EOSVAR_PRESSURE, i, j, k)];
                    q[PRJ_EOS_GAMMA] = block->eosvar[EIDX(PRJ_EOSVAR_GAMMA, i, j, k)];
                    cs = sqrt(q[PRJ_EOS_GAMMA] * q[PRJ_EOS_PRESSURE] / w[PRJ_PRIM_RHO]);
                    denom =
                        (fabs(w[PRJ_PRIM_V1]) + cs) / block->dx[0] +
                        (fabs(w[PRJ_PRIM_V2]) + cs) / block->dx[1] +
                        (fabs(w[PRJ_PRIM_V3]) + cs) / block->dx[2];
                    dt_cell = cfl / denom;
#if PRJ_NRAD > 0
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
                    if (dt_cell < dt_min) {
                        dt_min = dt_cell;
                    }
                }
            }
        }
    }

    return prj_mpi_min_dt(dt_min);
}

void prj_timeint_stage1(prj_mesh *mesh, const prj_coord *coord, const prj_bc *bc, prj_eos *eos,
    prj_rad *rad, double dt, double *dt_src, prj_timer *timer)
{
    int bidx;

    (void)coord;
    (void)bc;
    PRJ_TIMER_START(timer, "stage1");
    PRJ_TIMER_START(timer, "flux_stage1");
    PRJ_TIMER_START(timer, "flux_local_stage1");
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (prj_timeint_local_block(block)) {
            prj_flux_update(eos, rad, block, block->W, block->eosvar, block->flux, 0);
        }
    }
    PRJ_TIMER_STOP(timer, "flux_local_stage1");
    PRJ_TIMER_START(timer, "flux_exchange_stage1");
    prj_riemann_flux_send(mesh);
    PRJ_TIMER_STOP(timer, "flux_exchange_stage1");
    PRJ_TIMER_STOP(timer, "flux_stage1");
#if PRJ_MHD
    PRJ_TIMER_START(timer, "mhd_emf_stage1");
    prj_timeint_mhd_update_mesh_emf(mesh, prj_timeint_stage1_array, timer);
    PRJ_TIMER_STOP(timer, "mhd_emf_stage1");
#endif
    PRJ_TIMER_START(timer, "mpi_exchange_flux_emf_stage1");
    prj_mpi_exchange_fluxes_and_emf(mesh, prj_mpi_current());
    PRJ_TIMER_STOP(timer, "mpi_exchange_flux_emf_stage1");
#if PRJ_MHD
    PRJ_TIMER_START(timer, "mhd_update_bf_stage1");
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (prj_timeint_local_block(block)) {
            prj_timeint_mhd_update_bf(block, 1, dt);
        }
    }
    PRJ_TIMER_STOP(timer, "mhd_update_bf_stage1");
#endif
    PRJ_TIMER_START(timer, "update_stage1");
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (prj_timeint_local_block(block)) {
            int i;
            int j;
            int k;

            PRJ_TIMER_START(timer, "src_update_stage1");
            prj_src_update(eos, block, block->W, block->dUdt);
            PRJ_TIMER_STOP(timer, "src_update_stage1");
            PRJ_TIMER_START(timer, "cell_update_stage1");
            for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
                for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                    for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
#if PRJ_MHD && PRJ_NRAD > 0 && PRJ_FAST_MHD_RAD_CELL_UPDATE
                        prj_timeint_update_cell_stage1_mhd_rad(mesh, rad, eos, block,
                            i, j, k, dt, dt_src);
#else
                        double w[PRJ_NVAR_PRIM];
                        double u[PRJ_NVAR_CONS];
                        double u1[PRJ_NVAR_CONS];
                        double fluxdiv[PRJ_NVAR_CONS];
                        int v;

                        prj_flux_div(block->flux, block->area, block->vol, i, j, k, fluxdiv);
                        prj_timeint_cell_prim(block->W, i, j, k, w);
                        prj_eos_prim2cons(eos, w, u);
                        prj_timeint_update_dt_src(block, u, i, j, k, dt_src);
                        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
                            u1[v] = u[v] + dt * (block->dUdt[VIDX(v, i, j, k)] + fluxdiv[v]);
                        }
#if PRJ_MHD
                        prj_timeint_mhd_set_cons_b_from_bf(block, block->Bf1, i, j, k, u1);
#endif
#if PRJ_NRAD > 0
                        {
                            double T_cell = block->eosvar[EIDX(PRJ_EOSVAR_TEMPERATURE, i, j, k)];
                            double lapse_cell = prj_timeint_cell_lapse(block, i, j, k);
                             prj_rad_freq_flux_apply(rad, block, block->W, u1, i, j, k, lapse_cell, dt);
                             PRJ_TIMER_START(prj_timer_current(), "rad_nucinel");
                             prj_rad_nucinel_step(rad, eos, u1, dt, T_cell);
                             PRJ_TIMER_STOP(prj_timer_current(), "rad_nucinel");
                             PRJ_TIMER_START(prj_timer_current(), "rad_eleinel");
                             prj_rad_eleinel_step(rad, eos, u1, dt, T_cell);
                             PRJ_TIMER_STOP(prj_timer_current(), "rad_eleinel");
                             prj_rad_energy_update(rad, eos, u1, dt, lapse_cell, &T_cell);
                             prj_rad_momentum_update(rad, eos, u1, dt, lapse_cell, T_cell);
                        }
#endif
                        prj_eos_cons2prim(eos, u1, w);
                        prj_timeint_apply_eint_floor(mesh, u1, w);
                        prj_timeint_cell_cons_store(block->U, i, j, k, u1);
                        prj_timeint_cell_prim_store(block->W1, i, j, k, w);
#endif
                    }
                }
            }
            PRJ_TIMER_STOP(timer, "cell_update_stage1");
        }
    }
    PRJ_TIMER_STOP(timer, "update_stage1");
    PRJ_TIMER_START(timer, "ghost_fill_stage1");
    PRJ_TIMER_START(timer, "eos_active_stage1");
    prj_eos_fill_active_cells(mesh, eos, 2);
    PRJ_TIMER_STOP(timer, "eos_active_stage1");
    PRJ_TIMER_START(timer, "boundary_ghost_bf_stage1");
    prj_boundary_fill_ghosts_and_bf(mesh, bc, 2, 1, eos);
    PRJ_TIMER_STOP(timer, "boundary_ghost_bf_stage1");
    PRJ_TIMER_START(timer, "eos_mesh_stage1");
    prj_eos_fill_mesh(mesh, eos, 2);
    PRJ_TIMER_STOP(timer, "eos_mesh_stage1");
    PRJ_TIMER_STOP(timer, "ghost_fill_stage1");
#if PRJ_USE_GRAVITY
    PRJ_TIMER_START(timer, "gravity_reduce_stage1");
    prj_gravity_monopole_reduce(mesh, 2);
    PRJ_TIMER_STOP(timer, "gravity_reduce_stage1");
    PRJ_TIMER_START(timer, "gravity_integrate_stage1");
    prj_gravity_monopole_integrate(mesh);
    PRJ_TIMER_STOP(timer, "gravity_integrate_stage1");
#endif
    PRJ_TIMER_STOP(timer, "stage1");
}

void prj_timeint_stage2(prj_mesh *mesh, const prj_coord *coord, const prj_bc *bc, prj_eos *eos,
    prj_rad *rad, double dt, double *dt_src, prj_timer *timer)
{
    int bidx;

    (void)coord;
    PRJ_TIMER_START(timer, "stage2");
    PRJ_TIMER_START(timer, "flux_stage2");
    PRJ_TIMER_START(timer, "flux_local_stage2");
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (prj_timeint_local_block(block)) {
            prj_flux_update(eos, rad, block, block->W1, block->eosvar, block->flux, 1);
        }
    }
    PRJ_TIMER_STOP(timer, "flux_local_stage2");
    PRJ_TIMER_START(timer, "flux_exchange_stage2");
    prj_riemann_flux_send(mesh);
    PRJ_TIMER_STOP(timer, "flux_exchange_stage2");
    PRJ_TIMER_STOP(timer, "flux_stage2");
#if PRJ_MHD
    PRJ_TIMER_START(timer, "mhd_emf_stage2");
    prj_timeint_mhd_update_mesh_emf(mesh, prj_timeint_stage2_array, timer);
    PRJ_TIMER_STOP(timer, "mhd_emf_stage2");
#endif
    PRJ_TIMER_START(timer, "mpi_exchange_flux_emf_stage2");
    prj_mpi_exchange_fluxes_and_emf(mesh, prj_mpi_current());
    PRJ_TIMER_STOP(timer, "mpi_exchange_flux_emf_stage2");
#if PRJ_MHD
    PRJ_TIMER_START(timer, "mhd_update_bf_stage2");
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (prj_timeint_local_block(block)) {
            prj_timeint_mhd_update_bf(block, 2, dt);
        }
    }
    PRJ_TIMER_STOP(timer, "mhd_update_bf_stage2");
#endif
    PRJ_TIMER_START(timer, "update_stage2");
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (prj_timeint_local_block(block)) {
            int i;
            int j;
            int k;

            PRJ_TIMER_START(timer, "src_update_stage2");
            prj_src_update(eos, block, block->W1, block->dUdt);
            PRJ_TIMER_STOP(timer, "src_update_stage2");
            PRJ_TIMER_START(timer, "cell_update_stage2");
            for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
                for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                    for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
#if PRJ_MHD && PRJ_NRAD > 0 && PRJ_FAST_MHD_RAD_CELL_UPDATE
                        prj_timeint_update_cell_stage2_mhd_rad(mesh, rad, eos, block,
                            i, j, k, dt, dt_src);
#else
                        double w[PRJ_NVAR_PRIM];
                        double u[PRJ_NVAR_CONS];
                        double u1[PRJ_NVAR_CONS];
                        double fluxdiv[PRJ_NVAR_CONS];
                        int v;

                        prj_flux_div(block->flux, block->area, block->vol, i, j, k, fluxdiv);
                        prj_timeint_cell_prim(block->W, i, j, k, w);
                        prj_eos_prim2cons(eos, w, u);
                        prj_timeint_cell_prim(block->W1, i, j, k, w);
                        prj_eos_prim2cons(eos, w, u1);
                        prj_timeint_update_dt_src(block, u, i, j, k, dt_src);
                        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
                            u[v] = 0.5 * u[v] + 0.5 * (u1[v] + dt * (block->dUdt[VIDX(v, i, j, k)] + fluxdiv[v]));
                        }
#if PRJ_MHD
                        prj_timeint_mhd_set_cons_b_from_bf(block, block->Bf, i, j, k, u);
#endif
#if PRJ_NRAD > 0
                        {
                            double T_cell = block->eosvar[EIDX(PRJ_EOSVAR_TEMPERATURE, i, j, k)];
                            double lapse_cell = prj_timeint_cell_lapse(block, i, j, k);
                            prj_rad_freq_flux_apply(rad, block, block->W1, u, i, j, k, lapse_cell, 0.5 * dt);
                            PRJ_TIMER_START(prj_timer_current(), "rad_nucinel");
                            prj_rad_nucinel_step(rad, eos, u, 0.5 * dt, T_cell);
                            PRJ_TIMER_STOP(prj_timer_current(), "rad_nucinel");
                            PRJ_TIMER_START(prj_timer_current(), "rad_eleinel");
                            prj_rad_eleinel_step(rad, eos, u, 0.5 * dt, T_cell);
                            PRJ_TIMER_STOP(prj_timer_current(), "rad_eleinel");
                            prj_rad_energy_update(rad, eos, u, 0.5 * dt, lapse_cell, &T_cell);
                            prj_rad_momentum_update(rad, eos, u, 0.5 * dt, lapse_cell, T_cell);
                        }
#endif
                        prj_eos_cons2prim(eos, u, w);
                        prj_timeint_apply_eint_floor(mesh, u, w);
                        prj_timeint_cell_cons_store(block->U, i, j, k, u);
                        prj_timeint_cell_prim_store(block->W, i, j, k, w);
#endif
                    }
                }
            }
            PRJ_TIMER_STOP(timer, "cell_update_stage2");
        }
    }
    PRJ_TIMER_STOP(timer, "update_stage2");
    PRJ_TIMER_START(timer, "ghost_fill_stage2");
    PRJ_TIMER_START(timer, "eos_active_stage2");
    prj_eos_fill_active_cells(mesh, eos, 1);
    PRJ_TIMER_STOP(timer, "eos_active_stage2");
    PRJ_TIMER_START(timer, "boundary_ghost_bf_stage2");
    prj_boundary_fill_ghosts_and_bf(mesh, bc, 1, 0, eos);
    PRJ_TIMER_STOP(timer, "boundary_ghost_bf_stage2");
    PRJ_TIMER_START(timer, "eos_mesh_stage2");
    prj_eos_fill_mesh(mesh, eos, 1);
    PRJ_TIMER_STOP(timer, "eos_mesh_stage2");
    PRJ_TIMER_STOP(timer, "ghost_fill_stage2");
#if PRJ_USE_GRAVITY
    PRJ_TIMER_START(timer, "gravity_reduce_stage2");
    prj_gravity_monopole_reduce(mesh, 1);
    PRJ_TIMER_STOP(timer, "gravity_reduce_stage2");
    PRJ_TIMER_START(timer, "gravity_integrate_stage2");
    prj_gravity_monopole_integrate(mesh);
    PRJ_TIMER_STOP(timer, "gravity_integrate_stage2");
#endif
    PRJ_TIMER_STOP(timer, "stage2");
}

void prj_timeint_step(prj_mesh *mesh, const prj_coord *coord, const prj_bc *bc, prj_eos *eos,
    prj_rad *rad, double dt, double *dt_src, prj_timer *timer)
{
    prj_timeint_stage1(mesh, coord, bc, eos, rad, dt, dt_src, timer);
    prj_timeint_stage2(mesh, coord, bc, eos, rad, dt, dt_src, timer);
}
