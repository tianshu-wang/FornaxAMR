#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include <hdf5.h>

#include "prj.h"

static void prj_fill_block_data(prj_mesh *mesh)
{
    int bidx;

    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        int i;
        int j;
        int k;

        if (block->id < 0 || block->active != 1) {
            continue;
        }
        for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
            for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
                for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                    block->W[VIDX(PRJ_PRIM_RHO, i, j, k)] = (double)(block->id + 1) + 0.001 * (double)(i + 2 * j + 3 * k);
                    block->W[VIDX(PRJ_PRIM_V1, i, j, k)] = 0.1 * (double)(block->id + 1);
                    block->W[VIDX(PRJ_PRIM_V2, i, j, k)] = 0.2 * (double)(block->id + 1);
                    block->W[VIDX(PRJ_PRIM_V3, i, j, k)] = 0.3 * (double)(block->id + 1);
                    block->W[VIDX(PRJ_PRIM_EINT, i, j, k)] = 1.0 + 0.05 * (double)(block->id + 1);
                    block->W[VIDX(PRJ_PRIM_YE, i, j, k)] = 0.01 * (double)(block->id + 1);
                }
            }
        }
    }
}

static int prj_compare_mesh_data(const prj_mesh *a, const prj_mesh *b)
{
    int bidx;

    if (a->nblocks != b->nblocks) {
        return 1;
    }
    for (bidx = 0; bidx < a->nblocks; ++bidx) {
        const prj_block *ba = &a->blocks[bidx];
        const prj_block *bb = &b->blocks[bidx];
        int v;
        int i;
        int j;
        int k;
        int n;

        if (ba->id != bb->id || ba->active != bb->active || ba->level != bb->level) {
            return 2;
        }
        for (n = 0; n < 56; ++n) {
            if (ba->slot[n].id != bb->slot[n].id ||
                ba->slot[n].rank != bb->slot[n].rank ||
                memcmp(ba->slot[n].xmin, bb->slot[n].xmin, sizeof(ba->slot[n].xmin)) != 0 ||
                memcmp(ba->slot[n].xmax, bb->slot[n].xmax, sizeof(ba->slot[n].xmax)) != 0 ||
                memcmp(ba->slot[n].dx, bb->slot[n].dx, sizeof(ba->slot[n].dx)) != 0) {
                return 3;
            }
        }
        if (ba->id < 0 || ba->active != 1) {
            continue;
        }
        for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
            for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
                for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                    for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                        if (memcmp(&ba->W[VIDX(v, i, j, k)], &bb->W[VIDX(v, i, j, k)], sizeof(double)) != 0) {
                            return 4;
                        }
                    }
                }
            }
        }
    }
    return 0;
}

static int prj_check_dump_file(const char *filename, int expected_blocks)
{
    hid_t file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t dset = H5Dopen2(file, "density", H5P_DEFAULT);
    hid_t space = H5Dget_space(dset);
    hsize_t dims[4];

    H5Sget_simple_extent_dims(space, dims, 0);
    H5Sclose(space);
    H5Dclose(dset);
    H5Fclose(file);
    if ((int)dims[0] != expected_blocks) {
        return 1;
    }
    if ((int)dims[1] != PRJ_BLOCK_SIZE || (int)dims[2] != PRJ_BLOCK_SIZE || (int)dims[3] != PRJ_BLOCK_SIZE) {
        return 2;
    }
    return 0;
}

int main(void)
{
    prj_coord coord = {-1.0, 1.0, -1.0, 1.0, -1.0, 1.0};
    prj_eos eos;
    prj_mesh mesh;
    prj_mesh mesh_in;
    double time;
    int step;
    int status;

    eos.filename[0] = '\0';
    mkdir("output", 0777);
    mkdir("output/plots", 0777);

    if (prj_mesh_init(&mesh, 2, 2, 2, 1, &coord) != 0) {
        return 1;
    }
    prj_fill_block_data(&mesh);
    prj_io_write_restart(&mesh, 1.25, 3);
    prj_io_read_restart(&mesh_in, &eos, "output/restart_00000003.h5", &time, &step);
    status = prj_compare_mesh_data(&mesh, &mesh_in);
    if (status != 0 || fabs(time - 1.25) > 1.0e-12 || step != 3) {
        prj_mesh_destroy(&mesh);
        prj_mesh_destroy(&mesh_in);
        return 2;
    }
    prj_mesh_destroy(&mesh_in);

    prj_io_write_dump(&mesh, "output/dump", 2);
    if (prj_check_dump_file("output/dump_00002.h5", 8) != 0) {
        prj_mesh_destroy(&mesh);
        return 3;
    }

    prj_amr_refine_block(&mesh, 7);
    prj_fill_block_data(&mesh);
    prj_io_write_restart(&mesh, 2.5, 4);
    prj_io_read_restart(&mesh_in, &eos, "output/restart_00000004.h5", &time, &step);
    status = prj_compare_mesh_data(&mesh, &mesh_in);
    if (status != 0 || fabs(time - 2.5) > 1.0e-12 || step != 4) {
        prj_mesh_destroy(&mesh);
        prj_mesh_destroy(&mesh_in);
        return 4;
    }
    prj_mesh_destroy(&mesh_in);

    prj_io_write_dump(&mesh, "output/dump", 4);
    if (prj_check_dump_file("output/dump_00004.h5", prj_mesh_count_active(&mesh)) != 0) {
        prj_mesh_destroy(&mesh);
        return 5;
    }

    if (system("OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 KMP_INIT_AT_FORK=FALSE KMP_AFFINITY=disabled KMP_USE_SHM=0 KMP_DUPLICATE_LIB_OK=TRUE MKL_THREADING_LAYER=SEQUENTIAL MPLCONFIGDIR=/tmp/prj-mpl /Users/tianshuw/opt/anaconda3/bin/python analysis/visualize.py pressure") != 0) {
        prj_mesh_destroy(&mesh);
        return 6;
    }
    if (access("output/plots/pressure/xy-00004.png", F_OK) != 0 ||
        access("output/plots/pressure/yz-00004.png", F_OK) != 0 ||
        access("output/plots/pressure/xz-00004.png", F_OK) != 0) {
        prj_mesh_destroy(&mesh);
        return 7;
    }

    prj_mesh_destroy(&mesh);
    return 0;
}
