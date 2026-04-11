#include <math.h>
#include <stdio.h>
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

static void prj_report_dt_diagnostics(const prj_mesh *mesh, prj_eos *eos, int rank)
{
    double dt_min = 1.0e99;
    double worst_rho = 0.0;
    double worst_eint = 0.0;
    double worst_ye = 0.0;
    double worst_p = 0.0;
    double worst_gamma = 0.0;
    double worst_cs = 0.0;
    double worst_v1 = 0.0;
    double worst_v2 = 0.0;
    double worst_v3 = 0.0;
    double worst_x = 0.0;
    double worst_y = 0.0;
    double worst_z = 0.0;
    int worst_block = -1;
    int worst_i = 0;
    int worst_j = 0;
    int worst_k = 0;
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
                    double q[PRJ_EOS_NQUANT];
                    double rho = block->W[VIDX(PRJ_PRIM_RHO, i, j, k)];
                    double eint = block->W[VIDX(PRJ_PRIM_EINT, i, j, k)];
                    double ye = block->W[VIDX(PRJ_PRIM_YE, i, j, k)];
                    double v1 = block->W[VIDX(PRJ_PRIM_V1, i, j, k)];
                    double v2 = block->W[VIDX(PRJ_PRIM_V2, i, j, k)];
                    double v3 = block->W[VIDX(PRJ_PRIM_V3, i, j, k)];
                    double denom;
                    double dt_cell;

                    prj_eos_rey(eos, rho, eint, ye, q);
                    denom =
                        (fabs(v1) + sqrt(q[PRJ_EOS_GAMMA] * q[PRJ_EOS_PRESSURE] / rho)) / block->dx[0] +
                        (fabs(v2) + sqrt(q[PRJ_EOS_GAMMA] * q[PRJ_EOS_PRESSURE] / rho)) / block->dx[1] +
                        (fabs(v3) + sqrt(q[PRJ_EOS_GAMMA] * q[PRJ_EOS_PRESSURE] / rho)) / block->dx[2];
                    dt_cell = 0.4 / denom;
                    if (dt_cell < dt_min || !isfinite(dt_cell)) {
                        dt_min = dt_cell;
                        worst_rho = rho;
                        worst_eint = eint;
                        worst_ye = ye;
                        worst_p = q[PRJ_EOS_PRESSURE];
                        worst_gamma = q[PRJ_EOS_GAMMA];
                        worst_cs = sqrt(q[PRJ_EOS_GAMMA] * q[PRJ_EOS_PRESSURE] / rho);
                        worst_v1 = v1;
                        worst_v2 = v2;
                        worst_v3 = v3;
                        worst_x = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
                        worst_y = block->xmin[1] + ((double)j + 0.5) * block->dx[1];
                        worst_z = block->xmin[2] + ((double)k + 0.5) * block->dx[2];
                        worst_block = block->id;
                        worst_i = i;
                        worst_j = j;
                        worst_k = k;
                    }
                }
            }
        }
    }

    fprintf(stderr,
        "rank=%d dt_diag dt_min=%g block=%d cell=(%d,%d,%d) pos=(%g,%g,%g) rho=%g eint=%g ye=%g p=%g gamma=%g cs=%g v=(%g,%g,%g)\n",
        rank, dt_min, worst_block, worst_i, worst_j, worst_k, worst_x, worst_y, worst_z,
        worst_rho, worst_eint, worst_ye, worst_p, worst_gamma, worst_cs, worst_v1, worst_v2, worst_v3);
}

static int prj_has_bad_numbers(const prj_mesh *mesh)
{
    double local_bad = 0.0;
    int bidx;

    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        const prj_block *block = &mesh->blocks[bidx];
        int v;
        int i;
        int j;
        int k;

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

int main(int argc, char **argv)
{
    prj_sim sim;
    prj_mpi mpi;
    double mass0;
    double mass;
    int local_blocks;

    memset(&sim, 0, sizeof(sim));
    sim.mesh.root_nx[0] = 4;
    sim.mesh.root_nx[1] = 4;
    sim.mesh.root_nx[2] = 4;
    sim.mesh.max_level = 2;
    prj_problem_cc(&sim);
    prj_mpi_init(&argc, &argv, &mpi);
    prj_mpi_decompose(&sim.mesh);
    prj_mpi_prepare(&sim.mesh, &mpi);
    if (prj_rank_is_root()) {
        mkdir("output", 0777);
    }

    if (sim.mesh.nblocks <= 0 || sim.mesh.max_level != 2) {
        prj_mesh_destroy(&sim.mesh);
        prj_mpi_finalize();
        return 1;
    }
    if (prj_mesh_count_active(&sim.mesh) <= 64) {
        prj_mesh_destroy(&sim.mesh);
        prj_mpi_finalize();
        return 2;
    }

    prj_boundary_fill_ghosts(&sim.mesh, &sim.bc, 1);
    prj_io_write_dump(&sim.mesh, sim.output_dir, sim.step);
    mass0 = prj_total_cons(&sim.mesh, PRJ_CONS_RHO);
    while (sim.step < 100) {
        sim.dt = prj_timeint_calc_dt(&sim.mesh, &sim.eos, sim.cfl);
        if (sim.dt <= 0.0 || !isfinite(sim.dt)) {
            local_blocks = 0;
            {
                int bidx;

                for (bidx = 0; bidx < sim.mesh.nblocks; ++bidx) {
                    if (prj_local_active_block(&sim.mesh.blocks[bidx])) {
                        local_blocks += 1;
                    }
                }
            }
            fprintf(stderr, "rank=%d invalid dt=%g step=%d local_blocks=%d\n",
                mpi.rank, sim.dt, sim.step, local_blocks);
            prj_report_dt_diagnostics(&sim.mesh, &sim.eos, mpi.rank);
            prj_mesh_destroy(&sim.mesh);
            prj_mpi_finalize();
            return 3;
        }
        prj_timeint_step(&sim.mesh, &sim.coord, &sim.bc, &sim.eos, &sim.rad, sim.dt);
        sim.time += sim.dt;
        sim.step += 1;
        if (sim.amr_interval > 0 && sim.step % sim.amr_interval == 0) {
            prj_boundary_fill_ghosts(&sim.mesh, &sim.bc, 1);
            prj_amr_adapt(&sim.mesh, &sim.eos);
            prj_boundary_fill_ghosts(&sim.mesh, &sim.bc, 1);
        }
        if (prj_has_bad_numbers(&sim.mesh)) {
            prj_mesh_destroy(&sim.mesh);
            prj_mpi_finalize();
            return 4;
        }
        if (sim.output_interval > 0 && sim.step % sim.output_interval == 0) {
            prj_boundary_fill_ghosts(&sim.mesh, &sim.bc, 1);
            prj_io_write_dump(&sim.mesh, sim.output_dir, sim.step);
        }
    }

    mass = prj_total_cons(&sim.mesh, PRJ_CONS_RHO);
    if (!(mass > 0.0) || !(mass0 > 0.0)) {
        prj_mesh_destroy(&sim.mesh);
        prj_mpi_finalize();
        return 5;
    }
    if (prj_rank_is_root() && (access("output/dump_00010.h5", F_OK) != 0 ||
        access("output/dump_00100.h5", F_OK) != 0)) {
        prj_mesh_destroy(&sim.mesh);
        prj_mpi_finalize();
        return 6;
    }

    prj_mesh_destroy(&sim.mesh);
    prj_mpi_finalize();
    return 0;
}
