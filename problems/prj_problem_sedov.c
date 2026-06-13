#include <math.h>
#include <string.h>

#include "prj.h"

#define PRJ_PROBLEM_SEDOV_INJECTION_RADIUS 0.05

static int prj_problem_local_block(const prj_block *block)
{
    return block != 0 && block->id >= 0 && block->active == 1 &&
        block->W != 0 && block->W1 != 0 && block->U != 0;
}

static int prj_problem_block_overlaps_ball(const prj_block *block,
    double cx, double cy, double cz, double radius)
{
    double dx = 0.0;
    double dy = 0.0;
    double dz = 0.0;

    if (block == 0 || block->id < 0 || block->active != 1) {
        return 0;
    }
    if (block->xmin[0] > cx) {
        dx = block->xmin[0] - cx;
    } else if (block->xmax[0] < cx) {
        dx = cx - block->xmax[0];
    }
    if (block->xmin[1] > cy) {
        dy = block->xmin[1] - cy;
    } else if (block->xmax[1] < cy) {
        dy = cy - block->xmax[1];
    }
    if (block->xmin[2] > cz) {
        dz = block->xmin[2] - cz;
    } else if (block->xmax[2] < cz) {
        dz = cz - block->xmax[2];
    }
    return dx * dx + dy * dy + dz * dz < radius * radius;
}

static void prj_problem_ensure_data_allocated(prj_mesh *mesh)
{
    int bidx;

    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (block->id >= 0 && block->active == 1 && block->W == 0) {
            prj_block_alloc_data(block);
            prj_block_setup_geometry(block, &mesh->coord);
        }
    }
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

static void prj_problem_refine_injection_region(prj_sim *sim, const prj_mpi *mpi,
    double cx, double cy, double cz, double rho, double pressure)
{
    int level;
    int target_level = sim->mesh.max_level;

    if (target_level <= 0) {
        return;
    }

    for (level = 0; level < target_level; ++level) {
        int bidx;
        int any_marked = 0;
        int refined;

        for (bidx = 0; bidx < sim->mesh.nblocks; ++bidx) {
            prj_block *block = &sim->mesh.blocks[bidx];

            if (block->id < 0 || block->active != 1 || block->level != level) {
                continue;
            }
            block->refine_flag = 0;
            if (prj_problem_block_overlaps_ball(block, cx, cy, cz,
                    PRJ_PROBLEM_SEDOV_INJECTION_RADIUS)) {
                block->refine_flag = 1;
                any_marked = 1;
            }
        }

        if (any_marked == 0) {
            break;
        }

        prj_amr_enforce_two_to_one(&sim->mesh, mpi);
        refined = prj_amr_refine_marked_blocks(&sim->mesh, mpi);
        if (refined == 0) {
            break;
        }
        prj_amr_init_neighbors(&sim->mesh);
        prj_problem_ensure_data_allocated(&sim->mesh);
        prj_problem_fill_ambient(sim, rho, pressure);
    }
    prj_mesh_update_max_active_level(&sim->mesh);
}

static void prj_problem_inject_energy(prj_sim *sim, double cx, double cy, double cz)
{
    const double injection_radius = PRJ_PROBLEM_SEDOV_INJECTION_RADIUS;
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

void prj_problem_sedov(prj_sim *sim, prj_mpi *mpi)
{
    const double rho = 1.0;
    const double pressure = 1.0e-3;

    if (prj_mesh_init(&sim->mesh, sim->mesh.root_nx[0], sim->mesh.root_nx[1], sim->mesh.root_nx[2],
        sim->mesh.max_level, &sim->coord) != 0) {
        return;
    }
    prj_problem_fill_ambient(sim, rho, pressure);
    prj_problem_refine_injection_region(sim, mpi, 0.0, 0.0, 0.0, rho, pressure);
    prj_problem_inject_energy(sim, 0.0, 0.0, 0.0);
    prj_mpi_decompose(&(sim->mesh), mpi);
    prj_mpi_prepare(&(sim->mesh), mpi);
    prj_mhd_init(sim, mpi);
}
