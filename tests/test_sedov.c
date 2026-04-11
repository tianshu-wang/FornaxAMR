#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

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

static int prj_world_is_mpi(void)
{
    prj_mpi *mpi = prj_mpi_current();

    return mpi != 0 && mpi->totrank > 1;
}

static double prj_global_max_double(double local_value)
{
#if defined(PRJ_ENABLE_MPI)
    double neg_local = -local_value;
    double neg_global = prj_mpi_min_dt(neg_local);

    return -neg_global;
#else
    return local_value;
#endif
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

static double prj_blast_radius(const prj_mesh *mesh, double cx, double cy, double cz)
{
    const double ambient_pressure = 1.0e-3;
    double max_radius = 0.0;
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
                    double rho = block->W[VIDX(PRJ_PRIM_RHO, i, j, k)];
                    double eint = block->W[VIDX(PRJ_PRIM_EINT, i, j, k)];
                    double pressure = ((5.0 / 3.0) - 1.0) * rho * eint;

                    if (pressure > 2.0 * ambient_pressure) {
                        double x = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
                        double y = block->xmin[1] + ((double)j + 0.5) * block->dx[1];
                        double z = block->xmin[2] + ((double)k + 0.5) * block->dx[2];
                        double dx = x - cx;
                        double dy = y - cy;
                        double dz = z - cz;
                        double radius = sqrt(dx * dx + dy * dy + dz * dz);

                        if (radius > max_radius) {
                            max_radius = radius;
                        }
                    }
                }
            }
        }
    }
    return prj_global_max_double(max_radius);
}

static int prj_run_case(int mode, int argc, char **argv)
{
    prj_sim sim;
    prj_mpi mpi;
    double mass0;
    double energy0;
    double mass;
    double energy;
    int do_visualize;
    int max_steps;
    double blast_radius;
    double center[3];
    int wrote_output;
    int wrote_restart;

    memset(&sim, 0, sizeof(sim));
    sim.mesh.root_nx[0] = mode == 5 ? 8 : 4;
    sim.mesh.root_nx[1] = mode == 5 ? 8 : 4;
    sim.mesh.root_nx[2] = mode == 5 ? 8 : 4;
    sim.mesh.max_level = (mode == 2 || mode == 4 || mode == 5) ? 1 : 0;

    if (mode == 5) {
        prj_problem_sedov_offcenter(&sim);
        center[0] = 0.2;
        center[1] = 0.2;
        center[2] = 0.2;
    } else {
        prj_problem_sedov(&sim);
        center[0] = 0.0;
        center[1] = 0.0;
        center[2] = 0.0;
    }
    sim.mesh.amr_pressure_reference = 1.0;
    prj_mpi_init(&argc, &argv, &mpi);
    prj_mpi_decompose(&sim.mesh);
    prj_mpi_prepare(&sim.mesh, &mpi);

    if (prj_rank_is_root()) {
        mkdir("output", 0777);
    }
    if (mode == 2 || mode == 4 || mode == 5) {
        sim.amr_interval = 1;
    } else {
        sim.amr_interval = -1;
    }
    if (mode == 4) {
        sim.output_interval = 10;
        sim.restart_interval = 100;
        do_visualize = 0;
        max_steps = 10000;
    } else if (mode >= 3) {
        sim.output_interval = 10;
        sim.restart_interval = 100;
        do_visualize = 1;
        max_steps = sim.max_steps;
    } else {
        do_visualize = 0;
        max_steps = mode == 2 ? 5 : 10;
    }
    if (prj_world_is_mpi()) {
        do_visualize = 0;
    }
    wrote_output = 0;
    wrote_restart = 0;

    prj_boundary_fill_ghosts(&sim.mesh, &sim.bc, 1);
    mass0 = prj_total_cons(&sim.mesh, PRJ_CONS_RHO);
    energy0 = prj_total_cons(&sim.mesh, PRJ_CONS_ETOT);
    blast_radius = prj_blast_radius(&sim.mesh, center[0], center[1], center[2]);
    while (sim.step < max_steps && (mode < 3 || blast_radius < 1.8)) {
        sim.dt = prj_timeint_calc_dt(&sim.mesh, &sim.eos, sim.cfl);
        if (sim.time + sim.dt > sim.t_end) {
            sim.dt = sim.t_end - sim.time;
        }
        if (sim.dt <= 0.0) {
            prj_mesh_destroy(&sim.mesh);
            prj_mpi_finalize();
            return 5;
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
        blast_radius = prj_blast_radius(&sim.mesh, center[0], center[1], center[2]);
        if (prj_rank_is_root()) {
            printf("step=%d time=%.6e dt=%.6e blocks=%d mass_err=%.6e energy_err=%.6e energy=%.12e\n",
                sim.step, sim.time, sim.dt,
                prj_mesh_count_active(&sim.mesh),
                fabs((mass - mass0) / mass0), fabs((energy - energy0) / energy0), energy);
        }
        if (prj_has_bad_numbers(&sim.mesh)) {
            prj_mesh_destroy(&sim.mesh);
            prj_mpi_finalize();
            return 1;
        }
        if (mode <= 3 && fabs((energy - energy0) / energy0) > 1.0e-5) {
            prj_mesh_destroy(&sim.mesh);
            prj_mpi_finalize();
            return 2;
        }
        if (mode <= 2 && fabs((mass - mass0) / mass0) > 1.0e-5) {
            prj_mesh_destroy(&sim.mesh);
            prj_mpi_finalize();
            return 3;
        }
        if (sim.output_interval > 0 && sim.step % sim.output_interval == 0) {
            prj_io_write_dump(&sim.mesh, sim.output_dir, sim.step);
            wrote_output = 1;
        }
        if (sim.restart_interval > 0 && sim.step % sim.restart_interval == 0) {
            prj_io_write_restart(&sim.mesh, sim.time, sim.step);
            wrote_restart = 1;
        }
    }
    if (mode >= 3 && blast_radius < 1.8) {
        prj_mesh_destroy(&sim.mesh);
        prj_mpi_finalize();
        return 6;
    }
    if (do_visualize) {
        if (system("OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 KMP_INIT_AT_FORK=FALSE KMP_AFFINITY=disabled KMP_USE_SHM=0 KMP_DUPLICATE_LIB_OK=TRUE MKL_THREADING_LAYER=SEQUENTIAL MPLCONFIGDIR=/tmp/prj-mpl /Users/tianshuw/opt/anaconda3/bin/python analysis/visualize.py pressure") != 0) {
            prj_mesh_destroy(&sim.mesh);
            prj_mpi_finalize();
            return 4;
        }
        if (!wrote_output || access("output/plots/pressure/xy-00010.png", F_OK) != 0) {
            prj_mesh_destroy(&sim.mesh);
            prj_mpi_finalize();
            return 7;
        }
    }
    if (sim.restart_interval > 0 && sim.step >= sim.restart_interval && !wrote_restart) {
        prj_mesh_destroy(&sim.mesh);
        prj_mpi_finalize();
        return 8;
    }
    prj_mesh_destroy(&sim.mesh);
    prj_mpi_finalize();
    return 0;
}

int main(int argc, char *argv[])
{
    int mode = 1;

    if (argc > 1) {
        mode = atoi(argv[1]);
    }
    if (mode < 1 || mode > 5) {
        return 10;
    }
    return prj_run_case(mode, argc, argv);
}
