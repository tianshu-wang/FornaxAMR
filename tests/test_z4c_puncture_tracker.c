#include <math.h>
#include <stdio.h>
#include <string.h>

#if defined(PRJ_ENABLE_MPI)
#include <mpi.h>
#endif

#include "prj.h"

#if PRJ_DYNAMIC_GR
static int close_enough(double a, double b)
{
    return fabs(a - b) <= 2.0e-13 * fmax(1.0, fmax(fabs(a), fabs(b)));
}

int main(int argc, char **argv)
{
    prj_mesh mesh;
    prj_mpi mpi;
    prj_coord coord = {-2.0, 2.0, -2.0, 2.0, -2.0, 2.0};
    const double initial[3][3] = {
        {0.2, -0.3, 0.4},
        {-0.7, 0.6, -0.5},
        {1.1, -1.0, 0.9}
    };
    prj_block *block;
    int stage;
    int i, j, k, d, p;
    int status = 0;

#if defined(PRJ_ENABLE_MPI)
    MPI_Init(&argc, &argv);
#else
    (void)argc;
    (void)argv;
#endif
    memset(&mesh, 0, sizeof(mesh));
    memset(&mpi, 0, sizeof(mpi));
    mesh.max_blocks = 8;
    mesh.min_dx = 0.0;
    prj_z4c_init_params(&mesh.z4c_params);
#if defined(PRJ_ENABLE_MPI)
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi.rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi.totrank);
#else
    mpi.rank = 0;
    mpi.totrank = 1;
#endif
    if (prj_mesh_init(&mesh, 1, 1, 1, 0, &coord, 0) != 0) {
        fprintf(stderr, "tracker test mesh initialization failed\n");
        status = 1;
        goto done;
    }
    block = &mesh.blocks[0];
    block->rank = 0;
    for (stage = 0; stage < PRJ_BLOCK_NSTAGES; ++stage) {
        double *z = prj_block_z4c_stage(block, stage);

        for (i = -PRJ_NGHOST_Z4C; i < PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C; ++i) {
            for (j = -PRJ_NGHOST_Z4C; j < PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C; ++j) {
                for (k = -PRJ_NGHOST_Z4C; k < PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C; ++k) {
                    double x[3] = {
                        block->xmin[0] + ((double)i + 0.5) * block->dx[0],
                        block->xmin[1] + ((double)j + 0.5) * block->dx[1],
                        block->xmin[2] + ((double)k + 0.5) * block->dx[2]
                    };

                    for (d = 0; d < 3; ++d) {
                        z[Z4CIDX(PRJ_Z4C_BETAX + d, i, j, k)] = x[d];
                    }
                }
            }
        }
    }
    prj_z4c_puncture_tracker_init(&mesh, 3, initial);
    if (prj_z4c_puncture_count(&mesh) != 3) {
        fprintf(stderr, "tracker count mismatch\n");
        status = 2;
        goto done;
    }
    prj_z4c_puncture_update(&mesh, &mpi, 1, 0, 1.0, 0, 0.0, 0, 0.1);
    for (p = 0; p < 3; ++p) {
        const double *stage1 = mesh.z4c_puncture_positions_cm +
            3U * (size_t)mesh.z4c_puncture_count + 3U * (size_t)p;

        for (d = 0; d < 3; ++d) {
            if (!close_enough(stage1[d], 0.9 * initial[p][d])) {
                fprintf(stderr, "tracker interpolation/update mismatch\n");
                status = 3;
                goto done;
            }
        }
    }
    prj_z4c_puncture_blend(&mesh, 1, 0.25);
    for (p = 0; p < 3; ++p) {
        double position[3];

        if (!prj_z4c_puncture_position(&mesh, p, position)) {
            status = 4;
            goto done;
        }
        for (d = 0; d < 3; ++d) {
            if (!close_enough(position[d], 0.975 * initial[p][d])) {
                fprintf(stderr, "tracker stage blend mismatch\n");
                status = 5;
                goto done;
            }
        }
    }
done:
    if (mesh.blocks != 0) {
        prj_mesh_destroy(&mesh);
    }
#if defined(PRJ_ENABLE_MPI)
    MPI_Finalize();
#endif
    return status;
}
#else
int main(void) { return 0; }
#endif
