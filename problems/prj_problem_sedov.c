#include <math.h>
#include <string.h>

#include "prj.h"

static int prj_problem_local_block(const prj_block *block)
{
    prj_mpi *mpi = prj_mpi_current();

    return block != 0 && block->id >= 0 && block->active == 1 &&
        (mpi == 0 || block->rank == mpi->rank);
}

static void prj_problem_store_cell(prj_block *block, int i, int j, int k, const double *W, const double *U)
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

static double prj_problem_ball_overlap_fraction(const prj_block *block, int i, int j, int k,
    double cx, double cy, double cz, double radius)
{
    const int nsample = 8;
    int si;
    int sj;
    int sk;
    int inside = 0;
    int total = nsample * nsample * nsample;

    for (si = 0; si < nsample; ++si) {
        for (sj = 0; sj < nsample; ++sj) {
            for (sk = 0; sk < nsample; ++sk) {
                double x = block->xmin[0] + ((double)i + ((double)si + 0.5) / (double)nsample) * block->dx[0];
                double y = block->xmin[1] + ((double)j + ((double)sj + 0.5) / (double)nsample) * block->dx[1];
                double z = block->xmin[2] + ((double)k + ((double)sk + 0.5) / (double)nsample) * block->dx[2];
                double dx = x - cx;
                double dy = y - cy;
                double dz = z - cz;

                if (dx * dx + dy * dy + dz * dz < radius * radius) {
                    inside += 1;
                }
            }
        }
    }

    return (double)inside / (double)total;
}

static void prj_problem_fill_ambient(prj_sim *sim, double rho, double pressure)
{
    double eint = pressure / ((5.0 / 3.0) - 1.0) / rho;
    int bidx;

    for (bidx = 0; bidx < sim->mesh.nblocks; ++bidx) {
        prj_block *block = &sim->mesh.blocks[bidx];
        int i;
        int j;
        int k;

        if (!prj_problem_local_block(block)) {
            continue;
        }
        for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
            for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
                for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                    double W[PRJ_NVAR_PRIM] = {0.0};
                    double U[PRJ_NVAR_CONS] = {0.0};

                    W[PRJ_PRIM_RHO] = rho;
                    W[PRJ_PRIM_V1] = 0.0;
                    W[PRJ_PRIM_V2] = 0.0;
                    W[PRJ_PRIM_V3] = 0.0;
                    W[PRJ_PRIM_EINT] = eint;
                    W[PRJ_PRIM_YE] = 0.1;
                    prj_eos_prim2cons(&sim->eos, W, U);
                    prj_problem_store_cell(block, i, j, k, W, U);
                }
            }
        }
    }
}

static void prj_problem_inject_energy(prj_sim *sim, double cx, double cy, double cz)
{
    const double injection_radius = 0.05;
    int bidx;
    int selected = 0;
    double best_dist[8];
    prj_block *best_block[8];
    int best_i[8];
    int best_j[8];
    int best_k[8];
    int n;

    for (n = 0; n < 8; ++n) {
        best_dist[n] = 1.0e99;
        best_block[n] = 0;
        best_i[n] = 0;
        best_j[n] = 0;
        best_k[n] = 0;
    }

    for (bidx = 0; bidx < sim->mesh.nblocks; ++bidx) {
        prj_block *block = &sim->mesh.blocks[bidx];
        int i;
        int j;
        int k;

        if (!prj_problem_local_block(block)) {
            continue;
        }
        for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
            for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    double x = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
                    double y = block->xmin[1] + ((double)j + 0.5) * block->dx[1];
                    double z = block->xmin[2] + ((double)k + 0.5) * block->dx[2];
                    double r = sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy) + (z - cz) * (z - cz));

                    if (r < injection_radius) {
                        selected += 1;
                    }
                    for (n = 0; n < 8; ++n) {
                        if (r < best_dist[n]) {
                            int m;

                            for (m = 7; m > n; --m) {
                                best_dist[m] = best_dist[m - 1];
                                best_block[m] = best_block[m - 1];
                                best_i[m] = best_i[m - 1];
                                best_j[m] = best_j[m - 1];
                                best_k[m] = best_k[m - 1];
                            }
                            best_dist[n] = r;
                            best_block[n] = block;
                            best_i[n] = i;
                            best_j[n] = j;
                            best_k[n] = k;
                            break;
                        }
                    }
                }
            }
        }
    }

    if (selected == 0) {
        double weights[8];
        double weight_sum = 0.0;

        selected = 8;
        for (n = 0; n < 8; ++n) {
            if (best_block[n] != 0) {
                weights[n] = prj_problem_ball_overlap_fraction(
                    best_block[n], best_i[n], best_j[n], best_k[n], cx, cy, cz, injection_radius);
                weight_sum += weights[n];
            } else {
                weights[n] = 0.0;
            }
        }
        if (weight_sum <= 0.0) {
            for (n = 0; n < 8; ++n) {
                weights[n] = best_block[n] != 0 ? 1.0 : 0.0;
                weight_sum += weights[n];
            }
        }
        for (n = 0; n < 8; ++n) {
            if (best_block[n] != 0) {
                double cell_energy_density = (weights[n] / weight_sum) / best_block[n]->vol;

                best_block[n]->U[VIDX(PRJ_CONS_ETOT, best_i[n], best_j[n], best_k[n])] += cell_energy_density;
            }
        }
    } else {
        for (bidx = 0; bidx < sim->mesh.nblocks; ++bidx) {
            prj_block *block = &sim->mesh.blocks[bidx];
            int i;
            int j;
            int k;

            if (!prj_problem_local_block(block)) {
                continue;
            }
            for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
                for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                    for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                        double x = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
                        double y = block->xmin[1] + ((double)j + 0.5) * block->dx[1];
                        double z = block->xmin[2] + ((double)k + 0.5) * block->dx[2];
                        double r = sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy) + (z - cz) * (z - cz));

                        if (r < injection_radius) {
                            double cell_energy_density = (1.0 / (double)selected) / block->vol;

                            block->U[VIDX(PRJ_CONS_ETOT, i, j, k)] += cell_energy_density;
                        }
                    }
                }
            }
        }
    }

    for (bidx = 0; bidx < sim->mesh.nblocks; ++bidx) {
        prj_block *block = &sim->mesh.blocks[bidx];
        int i;
        int j;
        int k;

        if (!prj_problem_local_block(block)) {
            continue;
        }
        for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
            for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
                for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                    double U[PRJ_NVAR_CONS] = {0.0};
                    double W[PRJ_NVAR_PRIM] = {0.0};
                    int v;

                    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
                        U[v] = block->U[VIDX(v, i, j, k)];
                    }
                    prj_eos_cons2prim(&sim->eos, U, W);
                    for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
                        block->W[VIDX(v, i, j, k)] = W[v];
                        block->W1[VIDX(v, i, j, k)] = W[v];
                    }
                }
            }
        }
    }
}

void prj_problem_sedov(prj_sim *sim)
{
    if (prj_mesh_init(&sim->mesh, sim->mesh.root_nx[0], sim->mesh.root_nx[1], sim->mesh.root_nx[2],
        sim->mesh.max_level, &sim->coord) != 0) {
        return;
    }
    prj_problem_fill_ambient(sim, 1.0, 1.0e-3);
    prj_problem_inject_energy(sim, 0.0, 0.0, 0.0);
    prj_mhd_init(sim);
}
