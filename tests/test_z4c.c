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
    fprintf(stderr, "test_z4c: %s\n", msg);
    exit(1);
}

static void assert_close(const char *name, double got, double expected, double tol)
{
    if (!isfinite(got) || fabs(got - expected) > tol) {
        fprintf(stderr, "test_z4c: %s got %.17e expected %.17e tol %.3e\n",
            name, got, expected, tol);
        exit(1);
    }
}

static void check_enum_names(void)
{
    static const char *const expected[PRJ_NZ4C] = {
        "z4c_chi",
        "z4c_gxx", "z4c_gxy", "z4c_gxz", "z4c_gyy", "z4c_gyz", "z4c_gzz",
        "z4c_Khat",
        "z4c_Axx", "z4c_Axy", "z4c_Axz", "z4c_Ayy", "z4c_Ayz", "z4c_Azz",
        "z4c_Gamx", "z4c_Gamy", "z4c_Gamz",
        "z4c_Theta",
        "z4c_alpha",
        "z4c_betax", "z4c_betay", "z4c_betaz"
    };
    int v;

    if (PRJ_NZ4C != 22) {
        die("unexpected PRJ_NZ4C");
    }
    for (v = 0; v < PRJ_NZ4C; ++v) {
        if (strcmp(prj_z4c_var_name(v), expected[v]) != 0) {
            die("unexpected Z4c variable name/order");
        }
    }
}

#if PRJ_DYNAMIC_GR
static void check_flat_state(void)
{
    prj_mesh mesh;
    prj_coord coord;
    prj_block *block;
    double *z;
    double *rhs;
    double dt;
    int i, j, k, v;

    memset(&mesh, 0, sizeof(mesh));
    prj_z4c_init_params(&mesh.z4c_params);
    mesh.use_dynamic_gr = 1;
    coord.x1min = 0.0;
    coord.x1max = (double)PRJ_BLOCK_SIZE * 3.0e10;
    coord.x2min = 0.0;
    coord.x2max = (double)PRJ_BLOCK_SIZE * 3.0e10;
    coord.x3min = 0.0;
    coord.x3max = (double)PRJ_BLOCK_SIZE * 3.0e10;

    if (prj_mesh_init(&mesh, 1, 1, 1, 0, &coord, 0) != 0) {
        die("mesh init failed");
    }
    block = &mesh.blocks[0];
    if (block->z4c == 0 || block->z4c_rhs == 0) {
        die("Z4c block storage missing");
    }

    dt = prj_z4c_calc_dt_seconds(&mesh, 0, 0.5);
    assert_close("dt_z4c", dt, 0.5, 1.0e-3);

    prj_z4c_init_mesh_flat(&mesh, 0);
    z = prj_block_z4c_stage(block, 0);
    for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
        for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
            for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                assert_close("chi", z[Z4CIDX(PRJ_Z4C_CHI, i, j, k)], 1.0, 0.0);
                assert_close("gxx", z[Z4CIDX(PRJ_Z4C_GXX, i, j, k)], 1.0, 0.0);
                assert_close("gyy", z[Z4CIDX(PRJ_Z4C_GYY, i, j, k)], 1.0, 0.0);
                assert_close("gzz", z[Z4CIDX(PRJ_Z4C_GZZ, i, j, k)], 1.0, 0.0);
                assert_close("alpha", z[Z4CIDX(PRJ_Z4C_ALPHA, i, j, k)], 1.0, 0.0);
                for (v = 0; v < PRJ_NZ4C; ++v) {
                    if (v != PRJ_Z4C_CHI && v != PRJ_Z4C_GXX &&
                        v != PRJ_Z4C_GYY && v != PRJ_Z4C_GZZ &&
                        v != PRJ_Z4C_ALPHA) {
                        assert_close("flat zero", z[Z4CIDX(v, i, j, k)], 0.0, 0.0);
                    }
                }
            }
        }
    }

    prj_z4c_compute_rhs(&mesh, 0, 0, 0, 0.0);
    rhs = prj_block_z4c_rhs_stage(block, 0);
    for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
        for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
            for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                for (v = 0; v < PRJ_NZ4C; ++v) {
                    assert_close("flat rhs", rhs[Z4CIDX(v, i, j, k)], 0.0, 1.0e-50);
                }
            }
        }
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
    check_enum_names();
#if PRJ_DYNAMIC_GR
    check_flat_state();
#endif
    printf("test_z4c: ok\n");
#if defined(PRJ_ENABLE_MPI)
    MPI_Finalize();
#endif
    return 0;
}
