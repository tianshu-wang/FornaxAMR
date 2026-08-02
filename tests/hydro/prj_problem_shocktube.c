/* 1D Sod shock tube (run as a thin 3D pencil: root_nx2 = root_nx3 = 1).
 *
 * Initial discontinuity at x1 = 0.5:
 *   left  (x1 < 0.5): rho = 1.0,   P = 1.0
 *   right (x1 > 0.5): rho = 0.125, P = 0.1
 * velocity zero everywhere, ideal gas gamma = 5/3.  All boundaries outflow.
 */
#include <string.h>

#include "prj.h"

#define PRJ_ST_GAMMA (5.0 / 3.0)

/* Fill every resident active block, independent of the pre-decomposition rank
 * assignment.  main.c runs prj_mpi_decompose AFTER this init, and decompose
 * migrates each block's data from its OLD owner; if only the matching rank
 * filled a block the new owner would receive zeros.  Filling on all ranks (as
 * prj_problem_sedov does) keeps the data replicated so migration is safe. */
static int prj_problem_local_block(const prj_mpi *mpi, const prj_block *block)
{
    (void)mpi;
    return block != 0 && block->id >= 0 && block->active == 1;
}

static void prj_problem_st_ic(double x1, double x2, double x3, double *data)
{
    double rho;
    double pressure;

    (void)x2;
    (void)x3;
    if (x1 < 0.5) {
        rho = 1.0;
        pressure = 1.0;
    } else {
        rho = 0.125;
        pressure = 0.1;
    }
    data[PRJ_PRIM_RHO] = rho;
    data[PRJ_PRIM_V1] = 0.0;
    data[PRJ_PRIM_V2] = 0.0;
    data[PRJ_PRIM_V3] = 0.0;
    data[PRJ_PRIM_EINT] = pressure / ((PRJ_ST_GAMMA - 1.0) * rho);
    data[PRJ_PRIM_YE] = 0.1;
#if PRJ_MHD
    data[PRJ_PRIM_B1] = 0.0;
    data[PRJ_PRIM_B2] = 0.0;
    data[PRJ_PRIM_B3] = 0.0;
#endif
}

static void prj_problem_store_cell(prj_block *block, int i, int j, int k, const double *W, const double *U)
{
    prj_block_store_prim_cell(block, 0, i, j, k, W);
    prj_block_store_prim_cell(block, 1, i, j, k, W);
    prj_block_store_cons_cell(block, i, j, k, U);
}

void prj_problem_shocktube(prj_sim *sim, prj_mpi *mpi)
{
    int bidx;

    if (prj_mesh_init(&sim->mesh, sim->mesh.root_nx[0], sim->mesh.root_nx[1], sim->mesh.root_nx[2],
        sim->mesh.max_level, &sim->coord, 0) != 0) {
        return;
    }

    for (bidx = 0; bidx < sim->mesh.nblocks; ++bidx) {
        prj_block *block = &sim->mesh.blocks[bidx];
        int i;
        int j;
        int k;

        if (!prj_problem_local_block(mpi, block)) {
            continue;
        }
        for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
            for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
                for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                    double x1 = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
                    double x2 = block->xmin[1] + ((double)j + 0.5) * block->dx[1];
                    double x3 = block->xmin[2] + ((double)k + 0.5) * block->dx[2];
                    double W[PRJ_NVAR_PRIM] = {0.0};
                    double U[PRJ_NVAR_CONS] = {0.0};

                    prj_problem_st_ic(x1, x2, x3, W);
                    prj_eos_prim2cons(&sim->eos, W, U);
                    prj_problem_store_cell(block, i, j, k, W, U);
                }
            }
        }
    }
    prj_mhd_init(sim, mpi);
}
