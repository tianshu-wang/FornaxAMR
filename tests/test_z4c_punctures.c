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
    fprintf(stderr, "test_z4c_punctures: %s\n", msg);
    exit(1);
}

static void assert_close(const char *name, double got, double expected, double tol)
{
    if (!isfinite(got) || fabs(got - expected) > tol) {
        fprintf(stderr, "test_z4c_punctures: %s got %.17e expected %.17e tol %.3e\n",
            name, got, expected, tol);
        exit(1);
    }
}

#if PRJ_DYNAMIC_GR
static void setup_mesh(prj_mesh *mesh)
{
    prj_coord coord;

    memset(mesh, 0, sizeof(*mesh));
    prj_z4c_init_params(&mesh->z4c_params);
    mesh->use_dynamic_gr = 1;
    mesh->max_blocks = 16;
    coord.x1min = -4.0;
    coord.x1max = 4.0;
    coord.x2min = -4.0;
    coord.x2max = 4.0;
    coord.x3min = -4.0;
    coord.x3max = 4.0;
    if (prj_mesh_init(mesh, 1, 1, 1, 0, &coord, 0) != 0) {
        die("mesh init failed");
    }
}

static void check_single_puncture(void)
{
    prj_mesh mesh;
    prj_block *block;
    double centers[1][3] = {{0.0, 0.0, 0.0}};
    double masses[1] = {1.0};
    double momenta[1][3] = {{0.0, 0.0, 0.0}};
    double x[3];
    double r;
    double psi;
    double *z;

    setup_mesh(&mesh);
    prj_z4c_init_punctures(&mesh, 0, 1, centers, masses, momenta, 0.0);
    if (mesh.z4c_initialized == 0) {
        die("puncture initializer did not mark mesh initialized");
    }

    block = &mesh.blocks[0];
    z = prj_block_z4c_stage(block, 0);
    x[0] = block->xmin[0] + 0.5 * block->dx[0];
    x[1] = block->xmin[1] + 0.5 * block->dx[1];
    x[2] = block->xmin[2] + 0.5 * block->dx[2];
    r = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
    psi = 1.0 + 0.5 / r;

    assert_close("single chi", z[Z4CIDX(PRJ_Z4C_CHI, 0, 0, 0)], pow(psi, -4.0), 1.0e-14);
    assert_close("single alpha", z[Z4CIDX(PRJ_Z4C_ALPHA, 0, 0, 0)], pow(psi, -2.0), 1.0e-14);
    assert_close("single gxx", z[Z4CIDX(PRJ_Z4C_GXX, 0, 0, 0)], 1.0, 0.0);
    assert_close("single Axx", z[Z4CIDX(PRJ_Z4C_AXX, 0, 0, 0)], 0.0, 0.0);
    prj_mesh_destroy(&mesh);
}

static void check_binary_puncture_momentum(void)
{
    prj_mesh mesh;
    prj_block *block;
    double centers[2][3] = {{1.0, 0.0, 0.0}, {-1.0, 0.0, 0.0}};
    double masses[2] = {0.5, 0.5};
    double momenta[2][3] = {{0.0, -0.1, 0.0}, {0.0, 0.1, 0.0}};
    double *z;
    double max_a = 0.0;
    int i, j, k;

    setup_mesh(&mesh);
    prj_z4c_init_punctures(&mesh, 0, 2, centers, masses, momenta, 0.0);
    block = &mesh.blocks[0];
    z = prj_block_z4c_stage(block, 0);

    for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
        for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
            for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                double axx = fabs(z[Z4CIDX(PRJ_Z4C_AXX, i, j, k)]);
                double axy = fabs(z[Z4CIDX(PRJ_Z4C_AXY, i, j, k)]);
                double ayy = fabs(z[Z4CIDX(PRJ_Z4C_AYY, i, j, k)]);

                max_a = PRJ_MAX(max_a, PRJ_MAX(axx, PRJ_MAX(axy, ayy)));
            }
        }
    }
    if (!(max_a > 0.0)) {
        die("binary Bowen-York momentum did not seed Aij");
    }
    prj_mesh_destroy(&mesh);
}
#endif

int main(int argc, char **argv)
{
#if defined(PRJ_ENABLE_MPI)
    MPI_Init(&argc, &argv);
#else
    (void)argc;
    (void)argv;
#endif
#if PRJ_DYNAMIC_GR
    check_single_puncture();
    check_binary_puncture_momentum();
#endif
    printf("test_z4c_punctures: ok\n");
#if defined(PRJ_ENABLE_MPI)
    MPI_Finalize();
#endif
    return 0;
}
