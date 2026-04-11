#include <math.h>
#include <stdio.h>
#include <stdlib.h>
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

static double prj_total_cons(const prj_mesh *mesh, int var)
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

static int prj_has_bad_numbers(const prj_mesh *mesh)
{
    double local_bad = 0.0;
    int bidx;

    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        const prj_block *block = &mesh->blocks[bidx];
        int i;
        int j;
        int k;
        int v;

        if (!prj_local_active_block(block)) {
            continue;
        }
        for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
            for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
                for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                    for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                        if (!isfinite(block->W[VIDX(v, i, j, k)])) {
                            local_bad = 1.0;
                            break;
                        }
                    }
                    if (local_bad > 0.0) {
                        break;
                    }
                }
                if (local_bad > 0.0) {
                    break;
                }
            }
            if (local_bad > 0.0) {
                break;
            }
        }
        if (local_bad > 0.0) {
            break;
        }
    }
    return prj_mpi_global_sum(local_bad) > 0.0 ? 1 : 0;
}

static int prj_run_case(int mode, int argc, char **argv)
{
    prj_sim sim;
    prj_mpi mpi;
    double mass0;
    double energy0;
    double mass;
    double energy;

    memset(&sim, 0, sizeof(sim));
    sim.mesh.root_nx[0] = 2;
    sim.mesh.root_nx[1] = 2;
    sim.mesh.root_nx[2] = 2;
    sim.mesh.max_level = mode == 2 ? 1 : 0;
    prj_problem_shock1d(&sim);
    prj_mpi_init(&argc, &argv, &mpi);
    prj_mpi_decompose(&sim.mesh);
    prj_mpi_prepare(&sim.mesh, &mpi);
    sim.amr_interval = mode == 2 ? 1 : -1;
    if (prj_rank_is_root()) {
        mkdir("output", 0777);
    }

    prj_boundary_fill_ghosts(&sim.mesh, &sim.bc, 1);
    mass0 = prj_total_cons(&sim.mesh, PRJ_CONS_RHO);
    energy0 = prj_total_cons(&sim.mesh, PRJ_CONS_ETOT);

    while (sim.time < sim.t_end && sim.step < sim.max_steps) {
        sim.dt = prj_timeint_calc_dt(&sim.mesh, &sim.eos, sim.cfl);
        if (sim.time + sim.dt > sim.t_end) {
            sim.dt = sim.t_end - sim.time;
        }
        if (sim.dt <= 0.0) {
            prj_mesh_destroy(&sim.mesh);
            return 2;
        }
        prj_timeint_step(&sim.mesh, &sim.coord, &sim.bc, &sim.eos, &sim.rad, sim.dt);
        sim.time += sim.dt;
        sim.step += 1;

        if (sim.amr_interval > 0 && sim.step % sim.amr_interval == 0) {
            prj_boundary_fill_ghosts(&sim.mesh, &sim.bc, 1);
            prj_amr_adapt(&sim.mesh, &sim.eos);
            prj_mpi_rebalance(&sim.mesh);
            prj_boundary_fill_ghosts(&sim.mesh, &sim.bc, 1);
        }

        mass = prj_total_cons(&sim.mesh, PRJ_CONS_RHO);
        energy = prj_total_cons(&sim.mesh, PRJ_CONS_ETOT);
        if (prj_rank_is_root()) {
            printf("step=%d time=%.6e dt=%.6e blocks=%d mass_err=%.6e energy_err=%.6e\n",
                sim.step, sim.time, sim.dt, prj_mesh_count_active(&sim.mesh),
                fabs((mass - mass0) / mass0), fabs((energy - energy0) / energy0));
        }

        if (prj_has_bad_numbers(&sim.mesh)) {
            prj_mesh_destroy(&sim.mesh);
            prj_mpi_finalize();
            return 3;
        }
        if (sim.output_interval > 0 && sim.step % sim.output_interval == 0) {
            prj_boundary_fill_ghosts(&sim.mesh, &sim.bc, 1);
            prj_io_write_dump(&sim.mesh, sim.output_dir, sim.step);
        }
        if (sim.restart_interval > 0 && sim.step % sim.restart_interval == 0) {
            prj_io_write_restart(&sim.mesh, sim.time, sim.step);
        }
    }

    prj_mesh_destroy(&sim.mesh);
    prj_mpi_finalize();
    return 0;
}

int main(int argc, char **argv)
{
    int mode = 1;

    if (argc > 1) {
        mode = atoi(argv[1]);
    }
    if (mode < 1 || mode > 2) {
        return 1;
    }
    return prj_run_case(mode, argc, argv);
}
