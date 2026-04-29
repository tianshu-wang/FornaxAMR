#include <math.h>

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

static void prj_timeint_apply_eint_floor(const prj_mesh *mesh, double *u, double *w)
{
    double rho;
    double kinetic;

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
    w[PRJ_PRIM_EINT] = mesh->E_floor;
    u[PRJ_CONS_ETOT] = rho * (mesh->E_floor + kinetic);
}

static void prj_timeint_update_dt_src(const prj_block *block, const double *u, int i, int j, int k, double *dt_src)
{
    double dt_src_local;

    if (dt_src == 0) {
        return;
    }

    dt_src_local = 0.02 * u[PRJ_CONS_ETOT] /
        (fabs(block->dUdt[VIDX(PRJ_CONS_ETOT, i, j, k)]) + 1.0e-50);
    if (dt_src_local < *dt_src) {
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
    prj_rad *rad, double dt, double *dt_src)
{
    int bidx;

    (void)coord;
    prj_eos_fill_active_cells(mesh, eos, 1);
    prj_boundary_fill_ghosts(mesh, bc, 1);
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
                        double u[PRJ_NVAR_CONS];
                        double u1[PRJ_NVAR_CONS];
                        int v;

                        prj_timeint_cell_prim(block->W, i, j, k, w);
                        prj_eos_prim2cons(eos, w, u);
                        prj_timeint_update_dt_src(block, u, i, j, k, dt_src);
                        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
                            u1[v] = u[v] + dt * block->dUdt[VIDX(v, i, j, k)];
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

void prj_timeint_stage2(prj_mesh *mesh, const prj_coord *coord, const prj_bc *bc, prj_eos *eos,
    prj_rad *rad, double dt, double *dt_src)
{
    int bidx;

    (void)coord;
    prj_eos_fill_active_cells(mesh, eos, 2);
    prj_boundary_fill_ghosts(mesh, bc, 2);
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
                        prj_timeint_update_dt_src(block, u, i, j, k, dt_src);
                        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
                            u[v] = 0.5 * u[v] + 0.5 * (u1[v] + dt * block->dUdt[VIDX(v, i, j, k)]);
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

void prj_timeint_step(prj_mesh *mesh, const prj_coord *coord, const prj_bc *bc, prj_eos *eos,
    prj_rad *rad, double dt, double *dt_src)
{
    prj_timeint_stage1(mesh, coord, bc, eos, rad, dt, dt_src);
    prj_timeint_stage2(mesh, coord, bc, eos, rad, dt, dt_src);
}
