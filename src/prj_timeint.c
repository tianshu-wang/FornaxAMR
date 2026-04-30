#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "prj.h"

#if PRJ_TIMER
#define PRJ_TIMER_START(timer, name) prj_timer_start((timer), (name))
#define PRJ_TIMER_STOP(timer, name) prj_timer_stop((timer), (name))
#else
#define PRJ_TIMER_START(timer, name) ((void)(timer), (void)(name))
#define PRJ_TIMER_STOP(timer, name) ((void)(timer), (void)(name))
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
    value = src[dir][IDX(i, j, k)];
    if (dir == X1DIR) {
        value += -dt * (block->emf[X3DIR][IDX(i, j + 1, k)] -
            block->emf[X3DIR][IDX(i, j, k)]) / block->dx[1];
        value += dt * (block->emf[X2DIR][IDX(i, j, k + 1)] -
            block->emf[X2DIR][IDX(i, j, k)]) / block->dx[2];
    } else if (dir == X2DIR) {
        value += -dt * (block->emf[X1DIR][IDX(i, j, k + 1)] -
            block->emf[X1DIR][IDX(i, j, k)]) / block->dx[2];
        value += dt * (block->emf[X3DIR][IDX(i + 1, j, k)] -
            block->emf[X3DIR][IDX(i, j, k)]) / block->dx[0];
    } else {
        value += -dt * (block->emf[X2DIR][IDX(i + 1, j, k)] -
            block->emf[X2DIR][IDX(i, j, k)]) / block->dx[0];
        value += dt * (block->emf[X1DIR][IDX(i, j + 1, k)] -
            block->emf[X1DIR][IDX(i, j, k)]) / block->dx[1];
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
                        value = 0.5 * old[dir][IDX(i, j, k)] + 0.5 * value;
                    }
                    dst[dir][IDX(i, j, k)] = value;
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
    u[PRJ_CONS_B1] = 0.5 * (bf[X1DIR][IDX(i, j, k)] + bf[X1DIR][IDX(i + 1, j, k)]);
    u[PRJ_CONS_B2] = 0.5 * (bf[X2DIR][IDX(i, j, k)] + bf[X2DIR][IDX(i, j + 1, k)]);
    u[PRJ_CONS_B3] = 0.5 * (bf[X3DIR][IDX(i, j, k)] + bf[X3DIR][IDX(i, j, k + 1)]);
    if (!isfinite(u[PRJ_CONS_B1]) || !isfinite(u[PRJ_CONS_B2]) || !isfinite(u[PRJ_CONS_B3])) {
        prj_timeint_mhd_fail("prj_timeint_mhd_set_cons_b_from_bf: non-finite cell-centered magnetic field");
    }
}

static void prj_timeint_mhd_update_mesh_emf(prj_mesh *mesh, double *(*stage_array)(prj_block *))
{
    int bidx;

    if (mesh == 0 || stage_array == 0) {
        prj_timeint_mhd_fail("prj_timeint_mhd_update_mesh_emf: invalid input");
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (prj_timeint_local_block(block)) {
            prj_timeint_mhd_update_emf(block, stage_array(block));
        }
    }
    prj_mhd_emf_send(mesh);
    prj_mpi_exchange_emf(mesh, prj_mpi_current());
#if PRJ_MHD_DEBUG
    prj_mhd_debug_check_emf(mesh);
#endif
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

static void prj_timeint_update_dt_src(const prj_block *block, const double *u, int i, int j, int k, double *dt_src)
{
    double dt_src_local;

    if (dt_src == 0) {
        return;
    }

    dt_src_local = 0.02 * u[PRJ_CONS_ETOT] /
        (fabs(block->dUdt[VIDX(PRJ_CONS_ETOT, i, j, k)]) + 1.0e-50);
    if (dt_src_local < *dt_src) {
#if PRJ_GRAV_DEBUG
        prj_mpi *mpi = prj_mpi_current();
        int rank = mpi != 0 ? mpi->rank : 0;
        double rho = u[PRJ_CONS_RHO];
        double inv_rho = rho != 0.0 ? 1.0 / rho : 0.0;
        double x1 = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
        double x2 = block->xmin[1] + ((double)j + 0.5) * block->dx[1];
        double x3 = block->xmin[2] + ((double)k + 0.5) * block->dx[2];
        double accel = prj_gravity_block_accel_at(block, i, j, k);

        fprintf(stderr,
            "[grav debug] shortest dt_src rank=%d block=%d cell=(%d,%d,%d) "
            "x=(%.17e,%.17e,%.17e) rho=%.17e etot=%.17e "
            "v=(%.17e,%.17e,%.17e) accel=%.17e dt_src=%.17e\n",
            rank, block->id, i, j, k, x1, x2, x3, rho, u[PRJ_CONS_ETOT],
            u[PRJ_CONS_MOM1] * inv_rho, u[PRJ_CONS_MOM2] * inv_rho,
            u[PRJ_CONS_MOM3] * inv_rho, accel, dt_src_local);
#endif
        *dt_src = dt_src_local;
    }
}

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
    prj_riemann_set_mesh(mesh);
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (prj_timeint_local_block(block)) {
            prj_flux_update(eos, rad, block, block->W, block->eosvar, block->flux, 0);
        }
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (prj_timeint_local_block(block)) {
            prj_riemann_flux_send(block);
        }
    }
    PRJ_TIMER_STOP(timer, "flux_stage1");
#if PRJ_MHD
    PRJ_TIMER_START(timer, "flux_stage1");
    prj_timeint_mhd_update_mesh_emf(mesh, prj_timeint_stage1_array);
    PRJ_TIMER_STOP(timer, "flux_stage1");
#endif
    PRJ_TIMER_START(timer, "flux_stage1");
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (prj_timeint_local_block(block)) {
            int i;
            int j;
            int k;

            prj_flux_div(block->flux, block->area, block->vol, block->dUdt);
            prj_src_update(eos, block, block->W, block->dUdt);
#if PRJ_MHD
            prj_timeint_mhd_update_bf(block, 1, dt);
#endif
            for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
                for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                    for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                        double w[PRJ_NVAR_PRIM];
                        double u[PRJ_NVAR_CONS];
                        double u1[PRJ_NVAR_CONS];
                        int v;

                        prj_timeint_cell_prim(block->W, i, j, k, w);
                        prj_eos_prim2cons(eos, w, u);
                        prj_timeint_update_dt_src(block, u, i, j, k, dt_src);
                        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
                            u1[v] = u[v] + dt * block->dUdt[VIDX(v, i, j, k)];
                        }
#if PRJ_MHD
                        prj_timeint_mhd_set_cons_b_from_bf(block, block->Bf1, i, j, k, u1);
#endif
#if PRJ_NRAD > 0
                        {
                            double T_cell;
                            double lapse_cell = prj_timeint_cell_lapse(block, i, j, k);
                            /* Energy-space-flux part of SR redshift (Eqs. 21a/21b),
                             * applied to u1 before the stiff matter-coupling step.
                             * Stage1: closure built from W (cell-centred state at
                             * the start of the step), full dt weight. */
                            prj_rad_freq_flux_apply(rad, block, block->W, u1, i, j, k, lapse_cell, dt);
                            prj_rad_energy_update(rad, eos, u1, dt, lapse_cell, &T_cell);
                            prj_rad_momentum_update(rad, eos, u1, dt, lapse_cell, T_cell);
                        }
#endif
                        prj_eos_cons2prim(eos, u1, w);
                        prj_timeint_apply_eint_floor(mesh, u1, w);
                        prj_timeint_cell_cons_store(block->U, i, j, k, u1);
                        prj_timeint_cell_prim_store(block->W1, i, j, k, w);
                    }
                }
            }
        }
    }
    PRJ_TIMER_STOP(timer, "flux_stage1");
    PRJ_TIMER_START(timer, "ghost_fill_stage1");
    prj_eos_fill_active_cells(mesh, eos, 2);
    prj_boundary_fill_ghosts(mesh, bc, 2);
    prj_eos_fill_mesh(mesh, eos, 2);
#if PRJ_MHD
    prj_boundary_fill_bf(mesh, bc, 1);
#endif
    PRJ_TIMER_STOP(timer, "ghost_fill_stage1");
#if PRJ_USE_GRAVITY
    prj_gravity_monopole_reduce(mesh, 2);
    prj_gravity_monopole_integrate(mesh);
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
    prj_riemann_set_mesh(mesh);
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (prj_timeint_local_block(block)) {
            prj_flux_update(eos, rad, block, block->W1, block->eosvar, block->flux, 1);
        }
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (prj_timeint_local_block(block)) {
            prj_riemann_flux_send(block);
        }
    }
    PRJ_TIMER_STOP(timer, "flux_stage2");
#if PRJ_MHD
    PRJ_TIMER_START(timer, "flux_stage2");
    prj_timeint_mhd_update_mesh_emf(mesh, prj_timeint_stage2_array);
    PRJ_TIMER_STOP(timer, "flux_stage2");
#endif
    PRJ_TIMER_START(timer, "flux_stage2");
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (prj_timeint_local_block(block)) {
            int i;
            int j;
            int k;

            prj_flux_div(block->flux, block->area, block->vol, block->dUdt);
            prj_src_update(eos, block, block->W1, block->dUdt);
#if PRJ_MHD
            prj_timeint_mhd_update_bf(block, 2, dt);
#endif
            for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
                for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                    for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                        double w[PRJ_NVAR_PRIM];
                        double u[PRJ_NVAR_CONS];
                        double u1[PRJ_NVAR_CONS];
                        int v;

                        prj_timeint_cell_prim(block->W, i, j, k, w);
                        prj_eos_prim2cons(eos, w, u);
                        prj_timeint_cell_prim(block->W1, i, j, k, w);
                        prj_eos_prim2cons(eos, w, u1);
                        prj_timeint_update_dt_src(block, u, i, j, k, dt_src);
                        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
                            u[v] = 0.5 * u[v] + 0.5 * (u1[v] + dt * block->dUdt[VIDX(v, i, j, k)]);
                        }
#if PRJ_MHD
                        prj_timeint_mhd_set_cons_b_from_bf(block, block->Bf, i, j, k, u);
#endif
#if PRJ_NRAD > 0
                        {
                            double T_cell;
                            double lapse_cell = prj_timeint_cell_lapse(block, i, j, k);
                            /* Energy-space-flux part of SR redshift (Eqs. 21a/21b),
                             * applied to the post-average u before the stiff step.
                             * Stage2: closure from W1 (post-stage1 state); the 0.5·dt
                             * weight matches the RK2-Heun mixing of dUdt above. */
                            prj_rad_freq_flux_apply(rad, block, block->W1, u, i, j, k, lapse_cell, 0.5 * dt);
                            prj_rad_energy_update(rad, eos, u, dt, lapse_cell, &T_cell);
                            prj_rad_momentum_update(rad, eos, u, dt, lapse_cell, T_cell);
                        }
#endif
                        prj_eos_cons2prim(eos, u, w);
                        prj_timeint_apply_eint_floor(mesh, u, w);
                        prj_timeint_cell_cons_store(block->U, i, j, k, u);
                        prj_timeint_cell_prim_store(block->W, i, j, k, w);
                    }
                }
            }
        }
    }
    PRJ_TIMER_STOP(timer, "flux_stage2");
    PRJ_TIMER_START(timer, "ghost_fill_stage2");
    prj_eos_fill_active_cells(mesh, eos, 1);
    prj_boundary_fill_ghosts(mesh, bc, 1);
    prj_eos_fill_mesh(mesh, eos, 1);
#if PRJ_MHD
    prj_boundary_fill_bf(mesh, bc, 0);
#endif
    PRJ_TIMER_STOP(timer, "ghost_fill_stage2");
#if PRJ_USE_GRAVITY
    prj_gravity_monopole_reduce(mesh, 1);
    prj_gravity_monopole_integrate(mesh);
#endif
    PRJ_TIMER_STOP(timer, "stage2");
}

void prj_timeint_step(prj_mesh *mesh, const prj_coord *coord, const prj_bc *bc, prj_eos *eos,
    prj_rad *rad, double dt, double *dt_src, prj_timer *timer)
{
    prj_timeint_stage1(mesh, coord, bc, eos, rad, dt, dt_src, timer);
    prj_timeint_stage2(mesh, coord, bc, eos, rad, dt, dt_src, timer);
}
