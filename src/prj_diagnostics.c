#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(PRJ_ENABLE_MPI)
#include <mpi.h>
#endif

#include "prj.h"

#define PRJ_DQDT_NCOMP 6

static int prj_diagnostics_is_root_rank(const prj_mpi *mpi)
{
    return mpi == 0 || mpi->rank == 0;
}

static int prj_diagnostics_block_is_local_active(const prj_mpi *mpi, const prj_block *block)
{
    return block != 0 && block->id >= 0 && block->active == 1 &&
        block->W != 0 && (mpi == 0 || block->rank == mpi->rank);
}

void prj_diagnostics_truncate_dqdt(const prj_mpi *mpi)
{
    FILE *fp;

    if (!prj_diagnostics_is_root_rank(mpi)) {
        return;
    }
    fp = fopen("output/dqdt.txt", "w");
    if (fp == 0) {
        fprintf(stderr, "failed to truncate output/dqdt.txt: %s\n", strerror(errno));
        return;
    }
    fclose(fp);
}

static void prj_diagnostics_calc_dqdt_local(const prj_mesh *mesh, const prj_mpi *mpi,
    double dqdt[PRJ_DQDT_NCOMP])
{
    int bidx;
    int n;

    for (n = 0; n < PRJ_DQDT_NCOMP; ++n) {
        dqdt[n] = 0.0;
    }
    if (mesh == 0) {
        return;
    }

    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        const prj_block *block = &mesh->blocks[bidx];
        int i;
        int j;
        int k;

        if (!prj_diagnostics_block_is_local_active(mpi, block)) {
            continue;
        }
        for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
            for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    const double rho = block->W[WIDX(PRJ_PRIM_RHO, i, j, k)];
                    const double v0 = block->W[WIDX(PRJ_PRIM_V1, i, j, k)];
                    const double v1 = block->W[WIDX(PRJ_PRIM_V2, i, j, k)];
                    const double v2 = block->W[WIDX(PRJ_PRIM_V3, i, j, k)];
                    const double x0 = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
                    const double x1 = block->xmin[1] + ((double)j + 0.5) * block->dx[1];
                    const double x2 = block->xmin[2] + ((double)k + 0.5) * block->dx[2];
                    const double dm = rho * block->vol;
                    const double vdotx = v0 * x0 + v1 * x1 + v2 * x2;
                    const double trace_part = vdotx / 3.0;

                    dqdt[0] += 2.0 * dm * (v0 * x0 - trace_part);
                    dqdt[1] += dm * (v0 * x1 + v1 * x0);
                    dqdt[2] += dm * (v0 * x2 + v2 * x0);
                    dqdt[3] += 2.0 * dm * (v1 * x1 - trace_part);
                    dqdt[4] += dm * (v1 * x2 + v2 * x1);
                    dqdt[5] += 2.0 * dm * (v2 * x2 - trace_part);
                }
            }
        }
    }
}

static void prj_diagnostics_write_dqdt_root(double time, const double dqdt[PRJ_DQDT_NCOMP])
{
    FILE *fp = fopen("output/dqdt.txt", "a");

    if (fp == 0) {
        fprintf(stderr, "failed to append output/dqdt.txt: %s\n", strerror(errno));
        return;
    }
    fprintf(fp, "%.17e %.17e %.17e %.17e %.17e %.17e %.17e\n",
        time, dqdt[0], dqdt[1], dqdt[2], dqdt[3], dqdt[4], dqdt[5]);
    fclose(fp);
}

void prj_diagnostics_write_dqdt(const prj_mesh *mesh, const prj_mpi *mpi, double time)
{
    double local[PRJ_DQDT_NCOMP];
    double global[PRJ_DQDT_NCOMP];
    int n;

    prj_diagnostics_calc_dqdt_local(mesh, mpi, local);

#if defined(PRJ_ENABLE_MPI)
    if (mpi != 0 && mpi->totrank > 1) {
        double *gathered = 0;

        if (mpi->rank == 0) {
            gathered = (double *)calloc((size_t)mpi->totrank * PRJ_DQDT_NCOMP, sizeof(*gathered));
            if (gathered == 0) {
                fprintf(stderr, "failed to allocate dqdt gather buffer\n");
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            }
        }
        MPI_Gather(local, PRJ_DQDT_NCOMP, MPI_DOUBLE,
            gathered, PRJ_DQDT_NCOMP, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (mpi->rank == 0) {
            int rank;

            for (n = 0; n < PRJ_DQDT_NCOMP; ++n) {
                global[n] = 0.0;
            }
            for (rank = 0; rank < mpi->totrank; ++rank) {
                for (n = 0; n < PRJ_DQDT_NCOMP; ++n) {
                    global[n] += gathered[(size_t)rank * PRJ_DQDT_NCOMP + (size_t)n];
                }
            }
            free(gathered);
            prj_diagnostics_write_dqdt_root(time, global);
        }
        return;
    }
#else
    (void)mpi;
#endif

    for (n = 0; n < PRJ_DQDT_NCOMP; ++n) {
        global[n] = local[n];
    }
    if (prj_diagnostics_is_root_rank(mpi)) {
        prj_diagnostics_write_dqdt_root(time, global);
    }
}
