#include <math.h>

#include "prj.h"

#if PRJ_NRAD > 0
static double prj_timeint_cell_lapse(const prj_block *block, int i, int j, int k)
{
    const prj_grav_mono *gm = prj_gravity_active_monopole();
    double x1 = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
    double x2 = block->xmin[1] + ((double)j + 0.5) * block->dx[1];
    double x3 = block->xmin[2] + ((double)k + 0.5) * block->dx[2];
    double r = sqrt(x1 * x1 + x2 * x2 + x3 * x3);
    int idx;

    if (gm == 0 || gm->nbins <= 0 || gm->lapse == 0) {
        return 1.0;
    }
    if (r <= gm->rf[0]) {
        return gm->lapse[0];
    }
    for (idx = 0; idx < gm->nbins; ++idx) {
        double r0 = gm->rf[idx];
        double r1 = gm->rf[idx + 1];

        if (r <= r1) {
            double w = (r - r0) / (r1 - r0);

            return (1.0 - w) * gm->lapse[idx] + w * gm->lapse[idx + 1];
        }
    }
    return gm->lapse[gm->nbins];
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
    prj_boundary_fill_ghosts(mesh, bc, 1);
    prj_eos_fill_mesh(mesh, eos, 1);
    prj_gravity_monopole_reduce(mesh, 1);
    prj_gravity_monopole_integrate(mesh);
    prj_riemann_set_mesh(mesh);
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (prj_timeint_local_block(block)) {
            prj_flux_update(eos, rad, block->W, block->eosvar, block->flux);
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
#if PRJ_NRAD > 0
                        prj_rad_implicit_update(rad, eos, u1, dt,
                            prj_timeint_cell_lapse(block, i, j, k));
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

    (void)coord;
    prj_boundary_fill_ghosts(mesh, bc, 2);
    prj_eos_fill_mesh(mesh, eos, 2);
    prj_gravity_monopole_reduce(mesh, 2);
    prj_gravity_monopole_integrate(mesh);
    prj_riemann_set_mesh(mesh);
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (prj_timeint_local_block(block)) {
            prj_flux_update(eos, rad, block->W1, block->eosvar, block->flux);
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
#if PRJ_NRAD > 0
                        prj_rad_implicit_update(rad, eos, u, dt,
                            prj_timeint_cell_lapse(block, i, j, k));
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
}

void prj_timeint_step(prj_mesh *mesh, const prj_coord *coord, const prj_bc *bc, prj_eos *eos, prj_rad *rad, double dt)
{
    prj_timeint_stage1(mesh, coord, bc, eos, rad, dt);
    prj_timeint_stage2(mesh, coord, bc, eos, rad, dt);
}
