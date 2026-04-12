#include <math.h>

#include "prj.h"

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

double prj_timeint_calc_dt(const prj_mesh *mesh, prj_eos *eos, double cfl)
{
    double dt_min = 1.0e99;
    int bidx;

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
                    double dt_cell;

                    prj_timeint_cell_prim(block->W, i, j, k, w);
                    prj_eos_rey(eos, w[PRJ_PRIM_RHO], w[PRJ_PRIM_EINT], w[PRJ_PRIM_YE], q);
                    denom =
                        (fabs(w[PRJ_PRIM_V1]) + sqrt(q[PRJ_EOS_GAMMA] * q[PRJ_EOS_PRESSURE] / w[PRJ_PRIM_RHO])) / block->dx[0] +
                        (fabs(w[PRJ_PRIM_V2]) + sqrt(q[PRJ_EOS_GAMMA] * q[PRJ_EOS_PRESSURE] / w[PRJ_PRIM_RHO])) / block->dx[1] +
                        (fabs(w[PRJ_PRIM_V3]) + sqrt(q[PRJ_EOS_GAMMA] * q[PRJ_EOS_PRESSURE] / w[PRJ_PRIM_RHO])) / block->dx[2];
                    dt_cell = cfl / denom;
                    if (PRJ_NRAD_VAR > 0) {
                        double dx_min = block->dx[0];

                        if (block->dx[1] < dx_min) {
                            dx_min = block->dx[1];
                        }
                        if (block->dx[2] < dx_min) {
                            dx_min = block->dx[2];
                        }
                        if (cfl * dx_min / PRJ_CLIGHT < dt_cell) {
                            dt_cell = cfl * dx_min / PRJ_CLIGHT;
                        }
                    }
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

    (void)coord;
    prj_gravity_monopole_reduce(mesh, eos);
    prj_gravity_monopole_integrate(mesh);
    prj_boundary_fill_ghosts(mesh, bc, 1);
    prj_riemann_set_mesh(mesh);
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (prj_timeint_local_block(block)) {
            prj_flux_update(eos, rad, block->W, block->flux);
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
            prj_src_update(mesh, eos, block->W, block->U, block->dUdt);
            for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
                for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                    for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                        double w[PRJ_NVAR_PRIM];
                        double u[PRJ_NVAR_CONS];
                        double u1[PRJ_NVAR_CONS];
                        int v;

                        prj_timeint_cell_prim(block->W, i, j, k, w);
                        prj_eos_prim2cons(eos, w, u);
                        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
                            u1[v] = u[v] + dt * block->dUdt[VIDX(v, i, j, k)];
                        }
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

    (void)coord;
    prj_gravity_monopole_reduce(mesh, eos);
    prj_gravity_monopole_integrate(mesh);
    prj_boundary_fill_ghosts(mesh, bc, 2);
    prj_riemann_set_mesh(mesh);
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (prj_timeint_local_block(block)) {
            prj_flux_update(eos, rad, block->W1, block->flux);
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
            prj_src_update(mesh, eos, block->W1, block->U, block->dUdt);
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
                        }
                        prj_eos_cons2prim(eos, u, w);
                        prj_timeint_apply_eint_floor(mesh, u, w);
                        prj_timeint_cell_cons_store(block->U, i, j, k, u);
                        prj_timeint_cell_prim_store(block->W, i, j, k, w);
                    }
                }
            }
        }
    }
}

void prj_timeint_step(prj_mesh *mesh, const prj_coord *coord, const prj_bc *bc, prj_eos *eos, prj_rad *rad, double dt)
{
    prj_timeint_stage1(mesh, coord, bc, eos, rad, dt);
    prj_timeint_stage2(mesh, coord, bc, eos, rad, dt);
}
