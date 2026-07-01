#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(PRJ_ENABLE_MPI)
#include <mpi.h>
#endif

#include "prj.h"

static void die(const char *msg)
{
    fprintf(stderr, "test_mesh_r_com: %s\n", msg);
    exit(1);
}

static double expected_r_com(const prj_block *block, const double x_com[3], int i, int j, int k)
{
    double x1 = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
    double x2 = block->xmin[1] + ((double)j + 0.5) * block->dx[1];
    double x3 = block->xmin[2] + ((double)k + 0.5) * block->dx[2];
    double dx1 = x1 - x_com[0];
    double dx2 = x2 - x_com[1];
    double dx3 = x3 - x_com[2];

    return sqrt(dx1 * dx1 + dx2 * dx2 + dx3 * dx3);
}

static void assert_r_com(const prj_block *block, const double x_com[3], int i, int j, int k)
{
    int cache_idx = prj_block_cache_index(i, j, k);
    double expected = expected_r_com(block, x_com, i, j, k);
    double actual = block->r_com[cache_idx];
    double tol = 1.0e-12 * fmax(1.0, expected);

    if (!isfinite(actual) || fabs(actual - expected) > tol) {
        fprintf(stderr,
            "test_mesh_r_com: mismatch at (%d,%d,%d): got %.17e expected %.17e tol %.3e\n",
            i, j, k, actual, expected, tol);
        exit(1);
    }
}

int main(int argc, char **argv)
{
    prj_coord coord;
    prj_mesh mesh;
    prj_block *block;
    double x_com[3] = {0.25, -0.5, 0.75};
    double x_com_shifted[3] = {-0.75, 0.5, 2.0};
    double before;
    int cache_idx;
    int d;

#if defined(PRJ_ENABLE_MPI)
    MPI_Init(&argc, &argv);
#else
    (void)argc;
    (void)argv;
#endif

    memset(&mesh, 0, sizeof(mesh));
    coord.x1min = -1.0;
    coord.x1max = 3.0;
    coord.x2min = -2.0;
    coord.x2max = 2.0;
    coord.x3min = 0.0;
    coord.x3max = 8.0;

    if (prj_mesh_init(&mesh, 1, 1, 1, 0, &coord, 0) != 0) {
        die("mesh init failed");
    }
    if (mesh.nblocks != 1) {
        die("unexpected root block count");
    }
    block = &mesh.blocks[0];
    if (block->r_com == 0) {
        die("block r_com was not allocated");
    }

    for (d = 0; d < 3; ++d) {
        mesh.x_com[d] = x_com[d];
    }
    prj_mesh_update_r_com(&mesh);

    assert_r_com(block, x_com, 0, 0, 0);
    assert_r_com(block, x_com, PRJ_BLOCK_SIZE - 1, PRJ_BLOCK_SIZE - 1, PRJ_BLOCK_SIZE - 1);
#if PRJ_NGHOST > 0
    assert_r_com(block, x_com, -PRJ_NGHOST, 0, 0);
    assert_r_com(block, x_com,
        PRJ_BLOCK_SIZE + PRJ_NGHOST - 1,
        PRJ_BLOCK_SIZE + PRJ_NGHOST - 1,
        PRJ_BLOCK_SIZE + PRJ_NGHOST - 1);
#endif

    cache_idx = prj_block_cache_index(0, 0, 0);
    before = block->r_com[cache_idx];
    for (d = 0; d < 3; ++d) {
        mesh.x_com[d] = x_com_shifted[d];
    }
    prj_mesh_update_r_com(&mesh);
    assert_r_com(block, x_com_shifted, 0, 0, 0);
    if (fabs(block->r_com[cache_idx] - before) <= 1.0e-12) {
        die("r_com did not change after x_com changed");
    }

    prj_mesh_destroy(&mesh);
    printf("test_mesh_r_com: ok\n");

#if defined(PRJ_ENABLE_MPI)
    MPI_Finalize();
#endif
    return 0;
}
