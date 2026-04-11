#include <math.h>
#include <stdio.h>

#include "prj.h"

static double prj_total_cons(const prj_mesh *mesh, int var)
{
    double total = 0.0;
    int bidx;

    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        const prj_block *block = &mesh->blocks[bidx];
        int i;
        int j;
        int k;

        if (block->id < 0 || block->active != 1) {
            continue;
        }
        for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
            for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    total += block->U[VIDX(var, i, j, k)] * block->vol;
                }
            }
        }
    }
    return total;
}

static void prj_store_prim(double *dst, int i, int j, int k, const double *w)
{
    int v;

    for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
        dst[VIDX(v, i, j, k)] = w[v];
    }
}

static void prj_store_cons(double *dst, int i, int j, int k, const double *u)
{
    int v;

    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        dst[VIDX(v, i, j, k)] = u[v];
    }
}

static double prj_total_boundary_flux(const prj_mesh *mesh, int var)
{
    double total = 0.0;
    int bidx;

    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        const prj_block *block = &mesh->blocks[bidx];
        int j;
        int k;
        int i;

        if (block->id < 0 || block->active != 1) {
            continue;
        }
        if (fabs(block->xmin[0] - mesh->coord.x1min) < 1.0e-12) {
            for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    total -= block->area[X1DIR] * block->flux[X1DIR][VIDX(var, 0, j, k)];
                }
            }
        }
        if (fabs(block->xmax[0] - mesh->coord.x1max) < 1.0e-12) {
            for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    total += block->area[X1DIR] * block->flux[X1DIR][VIDX(var, PRJ_BLOCK_SIZE, j, k)];
                }
            }
        }
        if (fabs(block->xmin[1] - mesh->coord.x2min) < 1.0e-12) {
            for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    total -= block->area[X2DIR] * block->flux[X2DIR][VIDX(var, i, 0, k)];
                }
            }
        }
        if (fabs(block->xmax[1] - mesh->coord.x2max) < 1.0e-12) {
            for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    total += block->area[X2DIR] * block->flux[X2DIR][VIDX(var, i, PRJ_BLOCK_SIZE, k)];
                }
            }
        }
        if (fabs(block->xmin[2] - mesh->coord.x3min) < 1.0e-12) {
            for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
                for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                    total -= block->area[X3DIR] * block->flux[X3DIR][VIDX(var, i, j, 0)];
                }
            }
        }
        if (fabs(block->xmax[2] - mesh->coord.x3max) < 1.0e-12) {
            for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
                for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                    total += block->area[X3DIR] * block->flux[X3DIR][VIDX(var, i, j, PRJ_BLOCK_SIZE)];
                }
            }
        }
    }
    return total;
}

static void prj_prepare_boundary_fluxes(prj_mesh *mesh, prj_eos *eos, prj_rad *rad, const prj_bc *bc)
{
    int bidx;

    prj_boundary_fill_ghosts(mesh, bc, 1);
    prj_riemann_set_mesh(mesh);
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (block->id >= 0 && block->active == 1) {
            prj_flux_update(eos, rad, block->W, block->flux);
        }
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (block->id >= 0 && block->active == 1) {
            prj_riemann_flux_send(block);
        }
    }
}

static void prj_fill_smooth_pulse(prj_mesh *mesh)
{
    int bidx;

    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        int i;
        int j;
        int k;

        if (block->id < 0 || block->active != 1) {
            continue;
        }
        for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
            for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
                for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                    double x1 = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
                    double x2 = block->xmin[1] + ((double)j + 0.5) * block->dx[1];
                    double x3 = block->xmin[2] + ((double)k + 0.5) * block->dx[2];
                    double pulse = exp(-8.0 * (x1 * x1 + x2 * x2 + x3 * x3));
                    double W[PRJ_NVAR_PRIM];
                    double U[PRJ_NVAR_CONS];

                    W[PRJ_PRIM_RHO] = 1.0 + 0.1 * pulse;
                    W[PRJ_PRIM_V1] = 1.0e-3;
                    W[PRJ_PRIM_V2] = -8.0e-4;
                    W[PRJ_PRIM_V3] = 6.0e-4;
                    W[PRJ_PRIM_EINT] = 1.0 / W[PRJ_PRIM_RHO];
                    W[PRJ_PRIM_YE] = 0.1;
                    prj_store_prim(block->W, i, j, k, W);
                    prj_eos_prim2cons((prj_eos *)0, W, U);
                    prj_store_cons(block->U, i, j, k, U);
                }
            }
        }
    }
}

static void prj_fill_step_function(prj_mesh *mesh)
{
    int bidx;

    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        int i;
        int j;
        int k;

        if (block->id < 0 || block->active != 1) {
            continue;
        }
        for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
            for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
                for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                    double x1 = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
                    double W[PRJ_NVAR_PRIM];
                    double U[PRJ_NVAR_CONS];

                    W[PRJ_PRIM_RHO] = (x1 > -0.25 && x1 < 0.25) ? 2.0 : 1.0;
                    W[PRJ_PRIM_V1] = 1.0e-3;
                    W[PRJ_PRIM_V2] = -8.0e-4;
                    W[PRJ_PRIM_V3] = 6.0e-4;
                    W[PRJ_PRIM_EINT] = 0.75 / W[PRJ_PRIM_RHO];
                    W[PRJ_PRIM_YE] = 0.1;
                    prj_store_prim(block->W, i, j, k, W);
                    prj_eos_prim2cons((prj_eos *)0, W, U);
                    prj_store_cons(block->U, i, j, k, U);
                }
            }
        }
    }
}

int main(void)
{
    prj_coord coord = {-1.0, 1.0, -1.0, 1.0, -1.0, 1.0};
    prj_bc bc = {
        PRJ_BC_OUTFLOW, PRJ_BC_OUTFLOW,
        PRJ_BC_OUTFLOW, PRJ_BC_OUTFLOW,
        PRJ_BC_OUTFLOW, PRJ_BC_OUTFLOW
    };
    prj_eos eos;
    prj_rad rad;
    prj_mesh mesh;
    prj_mesh amr_mesh;
    double mass0;
    double energy0;
    double mass1;
    double energy1;
    double dt;
    double boundary_mass;
    double boundary_energy;
    double expected_mass;
    double expected_energy;
    int step;

    eos.filename[0] = '\0';
    prj_rad_init(&rad);

    if (prj_mesh_init(&mesh, 1, 1, 1, 0, &coord) != 0) {
        return 1;
    }
    prj_fill_smooth_pulse(&mesh);
    mass0 = prj_total_cons(&mesh, PRJ_CONS_RHO);
    energy0 = prj_total_cons(&mesh, PRJ_CONS_ETOT);
    expected_mass = mass0;
    expected_energy = energy0;
    for (step = 0; step < 10; ++step) {
        dt = prj_timeint_calc_dt(&mesh, &eos, 1.0e-6);
        prj_prepare_boundary_fluxes(&mesh, &eos, &rad, &bc);
        boundary_mass = prj_total_boundary_flux(&mesh, PRJ_CONS_RHO);
        boundary_energy = prj_total_boundary_flux(&mesh, PRJ_CONS_ETOT);
        printf("uniform step %d boundary fluxes %.12e %.12e\n", step, boundary_mass, boundary_energy);
        expected_mass += dt * boundary_mass;
        expected_energy += dt * boundary_energy;
        prj_timeint_step(&mesh, &coord, &bc, &eos, &rad, dt);
    }
    mass1 = prj_total_cons(&mesh, PRJ_CONS_RHO);
    energy1 = prj_total_cons(&mesh, PRJ_CONS_ETOT);
    if (fabs((mass1 - expected_mass) / mass0) > 1.0e-12) {
        prj_mesh_destroy(&mesh);
        return 2;
    }
    if (fabs((energy1 - expected_energy) / energy0) > 1.0e-12) {
        prj_mesh_destroy(&mesh);
        return 3;
    }
    prj_mesh_destroy(&mesh);

    if (prj_mesh_init(&amr_mesh, 2, 2, 2, 1, &coord) != 0) {
        return 4;
    }
    prj_fill_step_function(&amr_mesh);
    prj_amr_refine_block(&amr_mesh, 7);
    prj_fill_step_function(&amr_mesh);
    mass0 = prj_total_cons(&amr_mesh, PRJ_CONS_RHO);
    energy0 = prj_total_cons(&amr_mesh, PRJ_CONS_ETOT);
    expected_mass = mass0;
    expected_energy = energy0;
    for (step = 0; step < 10; ++step) {
        dt = prj_timeint_calc_dt(&amr_mesh, &eos, 1.0e-6);
        prj_prepare_boundary_fluxes(&amr_mesh, &eos, &rad, &bc);
        boundary_mass = prj_total_boundary_flux(&amr_mesh, PRJ_CONS_RHO);
        boundary_energy = prj_total_boundary_flux(&amr_mesh, PRJ_CONS_ETOT);
        printf("amr step %d boundary fluxes %.12e %.12e\n", step, boundary_mass, boundary_energy);
        expected_mass += dt * boundary_mass;
        expected_energy += dt * boundary_energy;
        prj_timeint_step(&amr_mesh, &coord, &bc, &eos, &rad, dt);
    }
    mass1 = prj_total_cons(&amr_mesh, PRJ_CONS_RHO);
    energy1 = prj_total_cons(&amr_mesh, PRJ_CONS_ETOT);
    if (fabs((mass1 - expected_mass) / mass0) > 1.0e-10) {
        prj_mesh_destroy(&amr_mesh);
        return 5;
    }
    if (fabs((energy1 - expected_energy) / energy0) > 1.0e-10) {
        prj_mesh_destroy(&amr_mesh);
        return 6;
    }
    prj_mesh_destroy(&amr_mesh);
    return 0;
}
