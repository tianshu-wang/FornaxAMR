#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "prj.h"

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

#if PRJ_MHD
static void prj_timeint_mhd_fail(const char *message)
{
    fprintf(stderr, "%s\n", message);
    abort();
}

static double prj_timeint_cell_emf_prim(const double *w, int dir)
{
    if (w == 0) {
        prj_timeint_mhd_fail("prj_timeint_cell_emf_prim: primitive state is null");
    }
    if (dir == X1DIR) {
        return w[PRJ_PRIM_V2] * w[PRJ_PRIM_B3] - w[PRJ_PRIM_V3] * w[PRJ_PRIM_B2];
    }
    if (dir == X2DIR) {
        return w[PRJ_PRIM_V3] * w[PRJ_PRIM_B1] - w[PRJ_PRIM_V1] * w[PRJ_PRIM_B3];
    }
    if (dir == X3DIR) {
        return w[PRJ_PRIM_V1] * w[PRJ_PRIM_B2] - w[PRJ_PRIM_V2] * w[PRJ_PRIM_B1];
    }
    prj_timeint_mhd_fail("prj_timeint_cell_emf_prim: invalid direction");
    return 0.0;
}

static double prj_timeint_cell_emf(const double *src, int dir, int i, int j, int k)
{
    double w[PRJ_NVAR_PRIM];

    prj_timeint_cell_prim(src, i, j, k, w);
    return prj_timeint_cell_emf_prim(w, dir);
}

static double prj_timeint_face_emf(const prj_block *block, int face_dir, int emf_dir, int i, int j, int k)
{
    if (block == 0) {
        prj_timeint_mhd_fail("prj_timeint_face_emf: block is null");
    }
    if (face_dir == X1DIR) {
        if (emf_dir == X2DIR && block->vB1[X1DIR] != 0) {
            return block->vB1[X1DIR][IDX(i, j, k)];
        }
        if (emf_dir == X3DIR && block->vB2[X1DIR] != 0) {
            return block->vB2[X1DIR][IDX(i, j, k)];
        }
    } else if (face_dir == X2DIR) {
        if (emf_dir == X3DIR && block->vB1[X2DIR] != 0) {
            return block->vB1[X2DIR][IDX(i, j, k)];
        }
        if (emf_dir == X1DIR && block->vB2[X2DIR] != 0) {
            return block->vB2[X2DIR][IDX(i, j, k)];
        }
    } else if (face_dir == X3DIR) {
        if (emf_dir == X1DIR && block->vB1[X3DIR] != 0) {
            return block->vB1[X3DIR][IDX(i, j, k)];
        }
        if (emf_dir == X2DIR && block->vB2[X3DIR] != 0) {
            return block->vB2[X3DIR][IDX(i, j, k)];
        }
    }
    prj_timeint_mhd_fail("prj_timeint_face_emf: invalid face/emf component request");
    return 0.0;
}

static double prj_timeint_face_normal_velocity(const prj_block *block, int face_dir, int i, int j, int k)
{
    if (block == 0 || block->v_riemann[face_dir] == 0) {
        prj_timeint_mhd_fail("prj_timeint_face_normal_velocity: face velocity is not allocated");
    }
    return block->v_riemann[face_dir][face_dir * PRJ_BLOCK_NCELLS + IDX(i, j, k)];
}

static void prj_timeint_build_block_emf(prj_block *block, const double *Wstage)
{
    int dir;
    int i;
    int j;
    int k;

    if (block == 0 || Wstage == 0) {
        prj_timeint_mhd_fail("prj_timeint_build_block_emf: invalid input");
    }
    for (dir = 0; dir < 3; ++dir) {
        if (block->emf[dir] == 0) {
            prj_timeint_mhd_fail("prj_timeint_build_block_emf: emf storage is not allocated");
        }
        for (i = 0; i < PRJ_BLOCK_NCELLS; ++i) {
            block->emf[dir][i] = 0.0;
        }
    }

    for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
        for (j = 0; j <= PRJ_BLOCK_SIZE; ++j) {
            for (k = 0; k <= PRJ_BLOCK_SIZE; ++k) {
                double emf_face[4];
                double emf_cell[4];
                double face_velocity[4];

                emf_face[0] = prj_timeint_face_emf(block, X3DIR, X1DIR, i, j - 1, k);
                emf_face[1] = prj_timeint_face_emf(block, X2DIR, X1DIR, i, j, k - 1);
                emf_face[2] = prj_timeint_face_emf(block, X3DIR, X1DIR, i, j, k);
                emf_face[3] = prj_timeint_face_emf(block, X2DIR, X1DIR, i, j, k);
                face_velocity[0] = prj_timeint_face_normal_velocity(block, X3DIR, i, j - 1, k);
                face_velocity[1] = prj_timeint_face_normal_velocity(block, X2DIR, i, j, k - 1);
                face_velocity[2] = prj_timeint_face_normal_velocity(block, X3DIR, i, j, k);
                face_velocity[3] = prj_timeint_face_normal_velocity(block, X2DIR, i, j, k);
                emf_cell[0] = prj_timeint_cell_emf(Wstage, X1DIR, i, j - 1, k - 1);
                emf_cell[1] = prj_timeint_cell_emf(Wstage, X1DIR, i, j, k - 1);
                emf_cell[2] = prj_timeint_cell_emf(Wstage, X1DIR, i, j, k);
                emf_cell[3] = prj_timeint_cell_emf(Wstage, X1DIR, i, j - 1, k);
                block->emf[X1DIR][IDX(i, j, k)] = prj_mhd_emf_upwind(emf_face, emf_cell, face_velocity);
            }
        }
    }

    for (i = 0; i <= PRJ_BLOCK_SIZE; ++i) {
        for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
            for (k = 0; k <= PRJ_BLOCK_SIZE; ++k) {
                double emf_face[4];
                double emf_cell[4];
                double face_velocity[4];

                emf_face[0] = prj_timeint_face_emf(block, X1DIR, X2DIR, i, j, k - 1);
                emf_face[1] = prj_timeint_face_emf(block, X3DIR, X2DIR, i - 1, j, k);
                emf_face[2] = prj_timeint_face_emf(block, X1DIR, X2DIR, i, j, k);
                emf_face[3] = prj_timeint_face_emf(block, X3DIR, X2DIR, i, j, k);
                face_velocity[0] = prj_timeint_face_normal_velocity(block, X1DIR, i, j, k - 1);
                face_velocity[1] = prj_timeint_face_normal_velocity(block, X3DIR, i - 1, j, k);
                face_velocity[2] = prj_timeint_face_normal_velocity(block, X1DIR, i, j, k);
                face_velocity[3] = prj_timeint_face_normal_velocity(block, X3DIR, i, j, k);
                emf_cell[0] = prj_timeint_cell_emf(Wstage, X2DIR, i - 1, j, k - 1);
                emf_cell[1] = prj_timeint_cell_emf(Wstage, X2DIR, i - 1, j, k);
                emf_cell[2] = prj_timeint_cell_emf(Wstage, X2DIR, i, j, k);
                emf_cell[3] = prj_timeint_cell_emf(Wstage, X2DIR, i, j, k - 1);
                block->emf[X2DIR][IDX(i, j, k)] = prj_mhd_emf_upwind(emf_face, emf_cell, face_velocity);
            }
        }
    }

    for (i = 0; i <= PRJ_BLOCK_SIZE; ++i) {
        for (j = 0; j <= PRJ_BLOCK_SIZE; ++j) {
            for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                double emf_face[4];
                double emf_cell[4];
                double face_velocity[4];

                emf_face[0] = prj_timeint_face_emf(block, X2DIR, X3DIR, i - 1, j, k);
                emf_face[1] = prj_timeint_face_emf(block, X1DIR, X3DIR, i, j - 1, k);
                emf_face[2] = prj_timeint_face_emf(block, X2DIR, X3DIR, i, j, k);
                emf_face[3] = prj_timeint_face_emf(block, X1DIR, X3DIR, i, j, k);
                face_velocity[0] = prj_timeint_face_normal_velocity(block, X2DIR, i - 1, j, k);
                face_velocity[1] = prj_timeint_face_normal_velocity(block, X1DIR, i, j - 1, k);
                face_velocity[2] = prj_timeint_face_normal_velocity(block, X2DIR, i, j, k);
                face_velocity[3] = prj_timeint_face_normal_velocity(block, X1DIR, i, j, k);
                emf_cell[0] = prj_timeint_cell_emf(Wstage, X3DIR, i - 1, j - 1, k);
                emf_cell[1] = prj_timeint_cell_emf(Wstage, X3DIR, i, j - 1, k);
                emf_cell[2] = prj_timeint_cell_emf(Wstage, X3DIR, i, j, k);
                emf_cell[3] = prj_timeint_cell_emf(Wstage, X3DIR, i - 1, j, k);
                block->emf[X3DIR][IDX(i, j, k)] = prj_mhd_emf_upwind(emf_face, emf_cell, face_velocity);
#if PRJ_MHD_DEBUG
                if ((block->id == 1 && i == 15 && j == PRJ_BLOCK_SIZE) ||
                    (block->id == 5 && i == 15 && j == 0)) {
                    fprintf(stderr,
                        "[mhd-edge-debug] block=%d edge=(%d,%d,%d) "
                        "face_emf=(%.17g,%.17g,%.17g,%.17g) "
                        "cell_emf=(%.17g,%.17g,%.17g,%.17g) "
                        "face_v=(%.17g,%.17g,%.17g,%.17g) emf=%.17g\n",
                        block->id, i, j, k,
                        emf_face[0], emf_face[1], emf_face[2], emf_face[3],
                        emf_cell[0], emf_cell[1], emf_cell[2], emf_cell[3],
                        face_velocity[0], face_velocity[1], face_velocity[2], face_velocity[3],
                        block->emf[X3DIR][IDX(i, j, k)]);
                }
#endif
            }
        }
    }
}

static void prj_timeint_copy_bface_stage(prj_block *block, int dst_stage, int src_stage)
{
    int dir;
    int idx;

    if (block == 0) {
        prj_timeint_mhd_fail("prj_timeint_copy_bface_stage: block is null");
    }
    for (dir = 0; dir < 3; ++dir) {
        double *dst = prj_mhd_bface_stage(block, dst_stage, dir);
        const double *src = prj_mhd_bface_stage_const(block, src_stage, dir);

        if (dst == 0 || src == 0) {
            prj_timeint_mhd_fail("prj_timeint_copy_bface_stage: face storage is not allocated");
        }
        for (idx = 0; idx < PRJ_BLOCK_NCELLS; ++idx) {
            dst[idx] = src[idx];
        }
    }
}

static void prj_timeint_update_block_bface(prj_block *block, int dst_stage, int src_stage, double dt)
{
    double *dstface[3];
    const double *srcface[3];
    int i;
    int j;
    int k;

    if (block == 0) {
        prj_timeint_mhd_fail("prj_timeint_update_block_bface: block is null");
    }
    if (block->emf[0] == 0 || block->emf[1] == 0 || block->emf[2] == 0) {
        prj_timeint_mhd_fail("prj_timeint_update_block_bface: emf storage is not allocated");
    }
    prj_timeint_copy_bface_stage(block, dst_stage, src_stage);
    dstface[0] = prj_mhd_bface_stage(block, dst_stage, X1DIR);
    dstface[1] = prj_mhd_bface_stage(block, dst_stage, X2DIR);
    dstface[2] = prj_mhd_bface_stage(block, dst_stage, X3DIR);
    srcface[0] = prj_mhd_bface_stage_const(block, src_stage, X1DIR);
    srcface[1] = prj_mhd_bface_stage_const(block, src_stage, X2DIR);
    srcface[2] = prj_mhd_bface_stage_const(block, src_stage, X3DIR);

    for (i = 0; i <= PRJ_BLOCK_SIZE; ++i) {
        for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
            for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                dstface[X1DIR][IDX(i, j, k)] =
                    srcface[X1DIR][IDX(i, j, k)] +
                    dt * ((block->emf[X3DIR][IDX(i, j + 1, k)] - block->emf[X3DIR][IDX(i, j, k)]) / block->dx[1] -
                        (block->emf[X2DIR][IDX(i, j, k + 1)] - block->emf[X2DIR][IDX(i, j, k)]) / block->dx[2]);
            }
        }
    }
    for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
        for (j = 0; j <= PRJ_BLOCK_SIZE; ++j) {
            for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                dstface[X2DIR][IDX(i, j, k)] =
                    srcface[X2DIR][IDX(i, j, k)] +
                    dt * ((block->emf[X1DIR][IDX(i, j, k + 1)] - block->emf[X1DIR][IDX(i, j, k)]) / block->dx[2] -
                        (block->emf[X3DIR][IDX(i + 1, j, k)] - block->emf[X3DIR][IDX(i, j, k)]) / block->dx[0]);
            }
        }
    }
    for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
        for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
            for (k = 0; k <= PRJ_BLOCK_SIZE; ++k) {
                dstface[X3DIR][IDX(i, j, k)] =
                    srcface[X3DIR][IDX(i, j, k)] +
                    dt * ((block->emf[X2DIR][IDX(i + 1, j, k)] - block->emf[X2DIR][IDX(i, j, k)]) / block->dx[0] -
                        (block->emf[X1DIR][IDX(i, j + 1, k)] - block->emf[X1DIR][IDX(i, j, k)]) / block->dx[1]);
            }
        }
    }
}
#endif

static void prj_timeint_apply_eint_floor(const prj_mesh *mesh, double *u, double *w)
{
    double rho;
    double kinetic;
    double magnetic;

    if (mesh == 0 || u == 0 || w == 0 || mesh->E_floor <= 0.0) {
        return;
    }

    rho = w[PRJ_PRIM_RHO];
    if (rho <= 0.0 || w[PRJ_PRIM_EINT] >= mesh->E_floor) {
        return;
    }

    kinetic = prj_eos_kinetic_energy_density_prim(w);
    magnetic = prj_eos_magnetic_energy_density_prim(w);
    w[PRJ_PRIM_EINT] = mesh->E_floor;
    u[PRJ_CONS_ETOT] = rho * mesh->E_floor + kinetic + magnetic;
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
#if PRJ_MHD
                    {
                        double rho = w[PRJ_PRIM_RHO];
                        double a2 = q[PRJ_EOS_GAMMA] * q[PRJ_EOS_PRESSURE] / rho;
                        double b2 = (w[PRJ_PRIM_B1] * w[PRJ_PRIM_B1] +
                                     w[PRJ_PRIM_B2] * w[PRJ_PRIM_B2] +
                                     w[PRJ_PRIM_B3] * w[PRJ_PRIM_B3]) / rho;
                        double ab2 = a2 + b2;
                        double bn2_1 = w[PRJ_PRIM_B1] * w[PRJ_PRIM_B1] / rho;
                        double bn2_2 = w[PRJ_PRIM_B2] * w[PRJ_PRIM_B2] / rho;
                        double bn2_3 = w[PRJ_PRIM_B3] * w[PRJ_PRIM_B3] / rho;
                        double disc1 = ab2 * ab2 - 4.0 * a2 * bn2_1;
                        double disc2 = ab2 * ab2 - 4.0 * a2 * bn2_2;
                        double disc3 = ab2 * ab2 - 4.0 * a2 * bn2_3;
                        double cf1 = sqrt(0.5 * (ab2 + sqrt(disc1 > 0.0 ? disc1 : 0.0)));
                        double cf2 = sqrt(0.5 * (ab2 + sqrt(disc2 > 0.0 ? disc2 : 0.0)));
                        double cf3 = sqrt(0.5 * (ab2 + sqrt(disc3 > 0.0 ? disc3 : 0.0)));
                        denom =
                            (fabs(w[PRJ_PRIM_V1]) + cf1) / block->dx[0] +
                            (fabs(w[PRJ_PRIM_V2]) + cf2) / block->dx[1] +
                            (fabs(w[PRJ_PRIM_V3]) + cf3) / block->dx[2];
                    }
#else
                    denom =
                        (fabs(w[PRJ_PRIM_V1]) + cs) / block->dx[0] +
                        (fabs(w[PRJ_PRIM_V2]) + cs) / block->dx[1] +
                        (fabs(w[PRJ_PRIM_V3]) + cs) / block->dx[2];
#endif
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

void prj_timeint_stage1(prj_mesh *mesh, const prj_coord *coord, const prj_bc *bc, prj_eos *eos, prj_rad *rad, double dt)
{
    int bidx;
#if PRJ_MHD
    prj_mpi *mpi = prj_mpi_current();
#endif

    (void)coord;
    prj_eos_fill_active_cells(mesh, eos, 1);
    prj_boundary_fill_ghosts(mesh, bc, 1);
#if PRJ_MHD
    prj_boundary_fill_ghosts_bf(mesh, bc, 1);
#endif
    prj_eos_fill_mesh(mesh, eos, 1);
    prj_gravity_monopole_reduce(mesh, 1);
    prj_gravity_monopole_integrate(mesh);
    prj_riemann_set_mesh(mesh);
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (prj_timeint_local_block(block)) {
            prj_flux_update(eos, rad, block, block->W, block->eosvar, block->flux);
        }
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (prj_timeint_local_block(block)) {
            prj_riemann_flux_send(block);
        }
    }
#if PRJ_MHD
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (prj_timeint_local_block(block)) {
            prj_timeint_build_block_emf(block, block->W);
        }
    }
    prj_mhd_emf_send(mesh);
    if (mpi != 0) {
        prj_mpi_exchange_emf(mesh, mpi);
    }
#if PRJ_MHD_DEBUG
    prj_mhd_debug_check_emf(mesh);
#endif
#endif
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (prj_timeint_local_block(block)) {
            int i;
            int j;
            int k;

            prj_flux_div(block->flux, block->area, block->vol, block->dUdt);
            prj_src_update(eos, block, block->W, block->dUdt);
            for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
                for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                    for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                        double w[PRJ_NVAR_PRIM];
                        double u1[PRJ_NVAR_CONS];
                        int v;

                        prj_timeint_cell_prim(block->W, i, j, k, w);
                        prj_eos_prim2cons(eos, w, u1);
                        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
                            u1[v] += dt * block->dUdt[VIDX(v, i, j, k)];
                            block->U[VIDX(v, i, j, k)] = u1[v];
                        }
                    }
                }
            }
#if PRJ_MHD
            prj_timeint_update_block_bface(block, 2, 1, dt);
            prj_mhd_bf2bc(block, 2);
#endif
            for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
                for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                    for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                        double w[PRJ_NVAR_PRIM];
                        double u1[PRJ_NVAR_CONS];
                        int v;

                        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
                            u1[v] = block->U[VIDX(v, i, j, k)];
                        }
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
}

void prj_timeint_stage2(prj_mesh *mesh, const prj_coord *coord, const prj_bc *bc, prj_eos *eos, prj_rad *rad, double dt)
{
    int bidx;
#if PRJ_MHD
    prj_mpi *mpi = prj_mpi_current();
#endif

    (void)coord;
    prj_eos_fill_active_cells(mesh, eos, 2);
    prj_boundary_fill_ghosts(mesh, bc, 2);
#if PRJ_MHD
    prj_boundary_fill_ghosts_bf(mesh, bc, 2);
#endif
    prj_eos_fill_mesh(mesh, eos, 2);
    prj_gravity_monopole_reduce(mesh, 2);
    prj_gravity_monopole_integrate(mesh);
    prj_riemann_set_mesh(mesh);
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (prj_timeint_local_block(block)) {
            prj_flux_update(eos, rad, block, block->W1, block->eosvar, block->flux);
        }
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (prj_timeint_local_block(block)) {
            prj_riemann_flux_send(block);
        }
    }
#if PRJ_MHD
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (prj_timeint_local_block(block)) {
            prj_timeint_build_block_emf(block, block->W1);
        }
    }
    prj_mhd_emf_send(mesh);
    if (mpi != 0) {
        prj_mpi_exchange_emf(mesh, mpi);
    }
#if PRJ_MHD_DEBUG
    prj_mhd_debug_check_emf(mesh);
#endif
#endif
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (prj_timeint_local_block(block)) {
            int i;
            int j;
            int k;

            prj_flux_div(block->flux, block->area, block->vol, block->dUdt);
            prj_src_update(eos, block, block->W1, block->dUdt);
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
                        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
                            u[v] = 0.5 * u[v] + 0.5 * (u1[v] + dt * block->dUdt[VIDX(v, i, j, k)]);
                            block->U[VIDX(v, i, j, k)] = u[v];
                        }
                    }
                }
            }
#if PRJ_MHD
            prj_timeint_update_block_bface(block, 1, 2, dt);
            prj_mhd_bf2bc(block, 1);
#endif
            for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
                for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                    for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                        double w[PRJ_NVAR_PRIM];
                        double u[PRJ_NVAR_CONS];
                        int v;

                        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
                            u[v] = block->U[VIDX(v, i, j, k)];
                        }
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
    prj_eos_fill_active_cells(mesh, eos, 1);
}

void prj_timeint_step(prj_mesh *mesh, const prj_coord *coord, const prj_bc *bc, prj_eos *eos, prj_rad *rad, double dt)
{
    prj_timeint_stage1(mesh, coord, bc, eos, rad, dt);
    prj_timeint_stage2(mesh, coord, bc, eos, rad, dt);
#if PRJ_MHD_DEBUG
    prj_mhd_debug_check_divergence(mesh, 1);
#endif
}
