#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>

#include "prj.h"

static int prj_local_active_block(const prj_block *block)
{
    prj_mpi *mpi = prj_mpi_current();

    return block != 0 && block->id >= 0 && block->active == 1 &&
        (mpi == 0 || block->rank == mpi->rank);
}

static int prj_rank_is_root(void)
{
    prj_mpi *mpi = prj_mpi_current();

    return mpi == 0 || mpi->rank == 0;
}

static double prj_total_group_energy(const prj_mesh *mesh, int field, int group)
{
    double total = 0.0;
    int var = PRJ_CONS_RAD_E(field, group);
    int bidx;

    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        const prj_block *block = &mesh->blocks[bidx];
        int i;
        int j;
        int k;

        if (!prj_local_active_block(block)) {
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
    return prj_mpi_global_sum(total);
}

static double prj_total_boundary_flux(const prj_mesh *mesh, int var)
{
    double total = 0.0;
    int bidx;

    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        const prj_block *block = &mesh->blocks[bidx];
        int i;
        int j;
        int k;

        if (!prj_local_active_block(block)) {
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
    return prj_mpi_global_sum(total);
}

static void prj_prepare_boundary_fluxes(prj_mesh *mesh, prj_eos *eos, prj_rad *rad, const prj_bc *bc)
{
    int bidx;

    prj_boundary_fill_ghosts(mesh, bc, 1);
    prj_riemann_set_mesh(mesh);
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (prj_local_active_block(block)) {
            prj_flux_update(eos, rad, block->W, block->flux);
        }
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (prj_local_active_block(block)) {
            prj_riemann_flux_send(block);
        }
    }
}

static void prj_store_state(prj_block *block, int i, int j, int k, const double *W, const double *U)
{
    int v;

    for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
        block->W[VIDX(v, i, j, k)] = W[v];
        block->W1[VIDX(v, i, j, k)] = W[v];
    }
    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        block->U[VIDX(v, i, j, k)] = U[v];
    }
}

static void prj_fill_radiation_problem(prj_sim *sim)
{
    int bidx;

    for (bidx = 0; bidx < sim->mesh.nblocks; ++bidx) {
        prj_block *block = &sim->mesh.blocks[bidx];
        int i;
        int j;
        int k;

        if (!prj_local_active_block(block)) {
            continue;
        }
        for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
            for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
                for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                    double x = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
                    double y = block->xmin[1] + ((double)j + 0.5) * block->dx[1];
                    double z = block->xmin[2] + ((double)k + 0.5) * block->dx[2];
                    double r = sqrt(x * x + y * y + z * z);
                    double W[PRJ_NVAR_PRIM];
                    double U[PRJ_NVAR_CONS];
                    int field;
                    int group;

                    for (field = 0; field < PRJ_NRAD; ++field) {
                        for (group = 0; group < PRJ_NEGROUP; ++group) {
                            W[PRJ_PRIM_RAD_E(field, group)] = r < 0.1 * (double)(group + 1) ? 1.0 : 0.0;
                            W[PRJ_PRIM_RAD_F1(field, group)] = 0.0;
                            W[PRJ_PRIM_RAD_F2(field, group)] = 0.0;
                            W[PRJ_PRIM_RAD_F3(field, group)] = 0.0;
                        }
                    }
                    W[PRJ_PRIM_RHO] = 1.0;
                    W[PRJ_PRIM_V1] = 0.0;
                    W[PRJ_PRIM_V2] = 0.0;
                    W[PRJ_PRIM_V3] = 0.0;
                    W[PRJ_PRIM_EINT] = 1.0e-3 / (((5.0 / 3.0) - 1.0) * W[PRJ_PRIM_RHO]);
                    W[PRJ_PRIM_YE] = 0.1;
                    prj_eos_prim2cons(&sim->eos, W, U);
                    prj_store_state(block, i, j, k, W, U);
                }
            }
        }
    }
}

static double prj_boundary_radiation_energy(const prj_mesh *mesh)
{
    double total = 0.0;
    int bidx;

    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        const prj_block *block = &mesh->blocks[bidx];
        int field;
        int group;
        int i;
        int j;
        int k;

        if (!prj_local_active_block(block)) {
            continue;
        }
        for (field = 0; field < PRJ_NRAD; ++field) {
            for (group = 0; group < PRJ_NEGROUP; ++group) {
                int var = PRJ_CONS_RAD_E(field, group);

                if (fabs(block->xmin[0] - mesh->coord.x1min) < 1.0e-12) {
                    for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                        for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                            total += block->U[VIDX(var, 0, j, k)];
                        }
                    }
                }
                if (fabs(block->xmax[0] - mesh->coord.x1max) < 1.0e-12) {
                    for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                        for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                            total += block->U[VIDX(var, PRJ_BLOCK_SIZE - 1, j, k)];
                        }
                    }
                }
                if (fabs(block->xmin[1] - mesh->coord.x2min) < 1.0e-12) {
                    for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
                        for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                            total += block->U[VIDX(var, i, 0, k)];
                        }
                    }
                }
                if (fabs(block->xmax[1] - mesh->coord.x2max) < 1.0e-12) {
                    for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
                        for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                            total += block->U[VIDX(var, i, PRJ_BLOCK_SIZE - 1, k)];
                        }
                    }
                }
                if (fabs(block->xmin[2] - mesh->coord.x3min) < 1.0e-12) {
                    for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
                        for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                            total += block->U[VIDX(var, i, j, 0)];
                        }
                    }
                }
                if (fabs(block->xmax[2] - mesh->coord.x3max) < 1.0e-12) {
                    for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
                        for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                            total += block->U[VIDX(var, i, j, PRJ_BLOCK_SIZE - 1)];
                        }
                    }
                }
            }
        }
    }
    return total;
}

int main(int argc, char **argv)
{
    prj_sim sim;
    prj_mpi mpi;
    double boundary_energy;
    double energy0[PRJ_NEGROUP];
    double expected[PRJ_NEGROUP];
    double energy1[PRJ_NEGROUP];
    int group;

    memset(&sim, 0, sizeof(sim));
    sim.coord.x1min = -2.0;
    sim.coord.x1max = 2.0;
    sim.coord.x2min = -2.0;
    sim.coord.x2max = 2.0;
    sim.coord.x3min = -2.0;
    sim.coord.x3max = 2.0;
    sim.bc.bc_x1_inner = PRJ_BC_OUTFLOW;
    sim.bc.bc_x1_outer = PRJ_BC_OUTFLOW;
    sim.bc.bc_x2_inner = PRJ_BC_OUTFLOW;
    sim.bc.bc_x2_outer = PRJ_BC_OUTFLOW;
    sim.bc.bc_x3_inner = PRJ_BC_OUTFLOW;
    sim.bc.bc_x3_outer = PRJ_BC_OUTFLOW;
    sim.cfl = 0.4;
    sim.t_end = 1.0e99;
    sim.max_steps = 10;
    sim.output_interval = 10;
    sim.restart_interval = -1;
    sim.amr_interval = -1;
    strcpy(sim.output_dir, "output/dump");

    prj_mpi_init(&argc, &argv, &mpi);
    if (prj_mesh_init(&sim.mesh, 4, 4, 4, 0, &sim.coord) != 0) {
        prj_mpi_finalize();
        return 1;
    }
    prj_amr_init_neighbors(&sim.mesh);
    prj_mpi_decompose(&sim.mesh);
    prj_mpi_prepare(&sim.mesh, &mpi);
    prj_rad_init(&sim.rad);
    prj_fill_radiation_problem(&sim);
    if (prj_rank_is_root()) {
        mkdir("output", 0777);
    }
    prj_boundary_fill_ghosts(&sim.mesh, &sim.bc, 1);
    prj_io_write_dump(&sim.mesh, sim.output_dir, sim.step);
    for (group = 0; group < PRJ_NEGROUP; ++group) {
        energy0[group] = prj_total_group_energy(&sim.mesh, 0, group);
        expected[group] = energy0[group];
    }

    boundary_energy = prj_boundary_radiation_energy(&sim.mesh);
    while (sim.step < sim.max_steps) {
        sim.dt = prj_timeint_calc_dt(&sim.mesh, &sim.eos, sim.cfl);
        if (sim.time + sim.dt > sim.t_end) {
            sim.dt = sim.t_end - sim.time;
        }
        if (sim.dt <= 0.0) {
            prj_mesh_destroy(&sim.mesh);
            prj_mpi_finalize();
            return 2;
        }
        prj_prepare_boundary_fluxes(&sim.mesh, &sim.eos, &sim.rad, &sim.bc);
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            expected[group] += sim.dt * prj_total_boundary_flux(&sim.mesh, PRJ_CONS_RAD_E(0, group));
        }
        prj_timeint_step(&sim.mesh, &sim.coord, &sim.bc, &sim.eos, &sim.rad, sim.dt);
        sim.time += sim.dt;
        sim.step += 1;
        prj_boundary_fill_ghosts(&sim.mesh, &sim.bc, 1);
        if (sim.output_interval > 0 && sim.step % sim.output_interval == 0) {
            prj_io_write_dump(&sim.mesh, sim.output_dir, sim.step);
        }
        boundary_energy = prj_boundary_radiation_energy(&sim.mesh);
    }

    if (sim.step % 10 != 0) {
        prj_io_write_dump(&sim.mesh, sim.output_dir, sim.step);
    }
    for (group = 0; group < PRJ_NEGROUP; ++group) {
        energy1[group] = prj_total_group_energy(&sim.mesh, 0, group);
    }
    if (prj_rank_is_root()) {
        printf("radiation boundary energy after %d steps: %.12e\n", sim.step, boundary_energy);
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            double denom = fabs(energy0[group]) > 1.0e-12 ? fabs(energy0[group]) : 1.0;
            double rel_err = fabs(energy1[group] - expected[group]) / denom;

            printf("group=%02d energy0=%.12e expected=%.12e energy1=%.12e rel_err=%.12e\n",
                group, energy0[group], expected[group], energy1[group], rel_err);
        }
    }
    for (group = 0; group < PRJ_NEGROUP; ++group) {
        double denom = fabs(energy0[group]) > 1.0e-12 ? fabs(energy0[group]) : 1.0;
        double rel_err = fabs(energy1[group] - expected[group]) / denom;

        if (rel_err > 1.0e-8) {
            prj_mesh_destroy(&sim.mesh);
            prj_mpi_finalize();
            return 3;
        }
    }
    prj_mesh_destroy(&sim.mesh);
    prj_mpi_finalize();
    return 0;
}
