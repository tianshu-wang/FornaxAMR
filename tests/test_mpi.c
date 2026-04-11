#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "prj.h"

static double test_value(double x, double y, double z)
{
    return 1.0 + 0.5 * x - 0.25 * y + 0.75 * z;
}

static void fill_gradient(prj_mesh *mesh)
{
    prj_mpi *mpi = prj_mpi_current();
    int bidx;

    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        int i;
        int j;
        int k;

        if (block->id < 0 || block->active != 1 ||
            (mpi != 0 && block->rank != mpi->rank) || block->W == 0 || block->W1 == 0) {
            continue;
        }
        for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
            for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
                for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                    double x = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
                    double y = block->xmin[1] + ((double)j + 0.5) * block->dx[1];
                    double z = block->xmin[2] + ((double)k + 0.5) * block->dx[2];
                    double rho = test_value(x, y, z);

                    block->W[VIDX(PRJ_PRIM_RHO, i, j, k)] = rho;
                    block->W[VIDX(PRJ_PRIM_V1, i, j, k)] = 0.0;
                    block->W[VIDX(PRJ_PRIM_V2, i, j, k)] = 0.0;
                    block->W[VIDX(PRJ_PRIM_V3, i, j, k)] = 0.0;
                    block->W[VIDX(PRJ_PRIM_EINT, i, j, k)] = 1.0;
                    block->W[VIDX(PRJ_PRIM_YE, i, j, k)] = 0.1;
                    block->W1[VIDX(PRJ_PRIM_RHO, i, j, k)] = rho;
                    block->W1[VIDX(PRJ_PRIM_V1, i, j, k)] = 0.0;
                    block->W1[VIDX(PRJ_PRIM_V2, i, j, k)] = 0.0;
                    block->W1[VIDX(PRJ_PRIM_V3, i, j, k)] = 0.0;
                    block->W1[VIDX(PRJ_PRIM_EINT, i, j, k)] = 1.0;
                    block->W1[VIDX(PRJ_PRIM_YE, i, j, k)] = 0.1;
                }
            }
        }
    }
}

static int expected_remote_ghost(const prj_mesh *mesh, const prj_block *target, int i, int j, int k, double *expected)
{
    int n;
    double x = target->xmin[0] + ((double)i + 0.5) * target->dx[0];
    double y = target->xmin[1] + ((double)j + 0.5) * target->dx[1];
    double z = target->xmin[2] + ((double)k + 0.5) * target->dx[2];

    for (n = 0; n < 56; ++n) {
        int nid = target->slot[n].id;
        const prj_block *neighbor;

        if (nid < 0 || nid >= mesh->nblocks) {
            continue;
        }
        neighbor = &mesh->blocks[nid];
        if (neighbor->id < 0 || neighbor->active != 1 || neighbor->rank == target->rank) {
            continue;
        }
        if (x > neighbor->xmin[0] - 1.0e-12 && x < neighbor->xmax[0] + 1.0e-12 &&
            y > neighbor->xmin[1] - 1.0e-12 && y < neighbor->xmax[1] + 1.0e-12 &&
            z > neighbor->xmin[2] - 1.0e-12 && z < neighbor->xmax[2] + 1.0e-12) {
            *expected = test_value(x, y, z);
            return 1;
        }
    }
    return 0;
}

int main(int argc, char **argv)
{
    prj_coord coord;
    prj_mesh mesh;
    prj_bc bc;
    prj_mpi mpi;
    int local_blocks = 0;
    int total_blocks = 0;
    int bidx;
    double local_volume = 0.0;
    double global_volume;
    int checked_remote_ghost = 0;

    memset(&mesh, 0, sizeof(mesh));
    memset(&bc, 0, sizeof(bc));
    coord.x1min = -1.0;
    coord.x1max = 1.0;
    coord.x2min = -1.0;
    coord.x2max = 1.0;
    coord.x3min = -1.0;
    coord.x3max = 1.0;

    prj_mpi_init(&argc, &argv, &mpi);
    if (mpi.totrank != 3) {
        if (mpi.rank == 0) {
            fprintf(stderr, "test_mpi expects 3 ranks\n");
        }
        prj_mpi_finalize();
        return 1;
    }
    if (prj_mesh_init(&mesh, 4, 4, 4, 1, &coord) != 0) {
        prj_mpi_finalize();
        return 2;
    }

    prj_amr_init_neighbors(&mesh);
    prj_mpi_decompose(&mesh);
    prj_mpi_prepare(&mesh, &mpi);
    for (bidx = 0; bidx < mesh.nblocks; ++bidx) {
        if (mesh.blocks[bidx].id >= 0 && mesh.blocks[bidx].active == 1) {
            total_blocks += 1;
            if (mesh.blocks[bidx].rank == mpi.rank) {
                local_blocks += 1;
                local_volume += (mesh.blocks[bidx].xmax[0] - mesh.blocks[bidx].xmin[0]) *
                    (mesh.blocks[bidx].xmax[1] - mesh.blocks[bidx].xmin[1]) *
                    (mesh.blocks[bidx].xmax[2] - mesh.blocks[bidx].xmin[2]);
            }
            if (mesh.blocks[bidx].rank < 0 || mesh.blocks[bidx].rank >= mpi.totrank) {
                fprintf(stderr, "invalid owner rank\n");
                prj_mesh_destroy(&mesh);
                prj_mpi_finalize();
                return 3;
            }
        }
    }
    if ((int)prj_mpi_global_sum((double)local_blocks) != total_blocks) {
        fprintf(stderr, "ownership accounting failed\n");
        prj_mesh_destroy(&mesh);
        prj_mpi_finalize();
        return 4;
    }

    fill_gradient(&mesh);
    prj_boundary_fill_ghosts(&mesh, &bc, 1);
    for (bidx = 0; bidx < mesh.nblocks; ++bidx) {
        prj_block *block = &mesh.blocks[bidx];
        int i;
        int j;
        int k;

        if (block->id < 0 || block->active != 1 || block->rank != mpi.rank) {
            continue;
        }
        for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
            for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
                for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                    double expected;
                    double actual;

                    if (i >= 0 && i < PRJ_BLOCK_SIZE &&
                        j >= 0 && j < PRJ_BLOCK_SIZE &&
                        k >= 0 && k < PRJ_BLOCK_SIZE) {
                        continue;
                    }
                    if (!expected_remote_ghost(&mesh, block, i, j, k, &expected)) {
                        continue;
                    }
                    actual = block->W[VIDX(PRJ_PRIM_RHO, i, j, k)];
                    if (fabs(actual - expected) > 1.0e-10) {
                        fprintf(stderr, "ghost exchange mismatch on rank %d\n", mpi.rank);
                        prj_mesh_destroy(&mesh);
                        prj_mpi_finalize();
                        return 5;
                    }
                    checked_remote_ghost = 1;
                }
            }
        }
    }
    if ((int)prj_mpi_global_sum((double)checked_remote_ghost) == 0) {
        fprintf(stderr, "no remote ghost cells were validated\n");
        prj_mesh_destroy(&mesh);
        prj_mpi_finalize();
        return 6;
    }

    global_volume = prj_mpi_global_sum(local_volume);
    if (fabs(global_volume - 8.0) > 1.0e-10) {
        fprintf(stderr, "global sum mismatch: %g\n", global_volume);
        prj_mesh_destroy(&mesh);
        prj_mpi_finalize();
        return 7;
    }

    prj_amr_refine_block(&mesh, 0);
    prj_amr_init_neighbors(&mesh);
    prj_mpi_rebalance(&mesh);
    prj_mpi_prepare(&mesh, &mpi);

    local_blocks = 0;
    for (bidx = 0; bidx < mesh.nblocks; ++bidx) {
        if (mesh.blocks[bidx].id >= 0 && mesh.blocks[bidx].active == 1 && mesh.blocks[bidx].rank == mpi.rank) {
            local_blocks += 1;
        }
    }
    total_blocks = (int)prj_mpi_global_sum((double)local_blocks);
    if (total_blocks <= 64) {
        fprintf(stderr, "rebalance test did not add refined blocks\n");
        prj_mesh_destroy(&mesh);
        prj_mpi_finalize();
        return 8;
    }
    prj_amr_coarsen_block(&mesh, 0);
    prj_amr_init_neighbors(&mesh);

    local_blocks = 0;
    local_volume = 0.0;
    for (bidx = 0; bidx < mesh.nblocks; ++bidx) {
        if (mesh.blocks[bidx].id >= 0 && mesh.blocks[bidx].active == 1) {
            if (mesh.blocks[bidx].rank == mpi.rank) {
                local_blocks += 1;
                local_volume += (mesh.blocks[bidx].xmax[0] - mesh.blocks[bidx].xmin[0]) *
                    (mesh.blocks[bidx].xmax[1] - mesh.blocks[bidx].xmin[1]) *
                    (mesh.blocks[bidx].xmax[2] - mesh.blocks[bidx].xmin[2]);
            }
        }
    }
    total_blocks = (int)prj_mpi_global_sum((double)local_blocks);
    if (total_blocks != 64) {
        fprintf(stderr, "coarsen test did not restore coarse active blocks\n");
        prj_mesh_destroy(&mesh);
        prj_mpi_finalize();
        return 9;
    }

    global_volume = prj_mpi_global_sum(local_volume);
    if (fabs(global_volume - 8.0) > 1.0e-10) {
        fprintf(stderr, "global volume mismatch after coarsen: %g\n", global_volume);
        prj_mesh_destroy(&mesh);
        prj_mpi_finalize();
        return 10;
    }

    prj_mesh_destroy(&mesh);
    prj_mpi_finalize();
    return 0;
}
