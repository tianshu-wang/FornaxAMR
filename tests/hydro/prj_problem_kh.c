/* 2D Kelvin-Helmholtz instability (run as a thin 3D slab: root_nx3 = 1).
 *
 * Two shear interfaces at y = 0.25 and y = 0.75 separate a dense inner stream
 * (rho = 2, v1 = +0.5) from a light outer stream (rho = 1, v1 = -0.5) in
 * uniform pressure.  A small localized v2 perturbation seeds the rollup.  The
 * domain is doubly periodic in x and y (bc_x1/bc_x2 = periodic); the ignorable
 * z-axis stays outflow.  With max_level > 0 the grid pre-refines to the shear
 * layers via the velocity AMR estimator (prj_problem_fill_until_amr_converged).
 */
#include <math.h>
#include <string.h>

#include "prj.h"

#define PRJ_KH_GAMMA (5.0 / 3.0)
#define PRJ_KH_PI 3.14159265358979323846

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

static unsigned long long prj_problem_mesh_signature(const prj_mesh *mesh)
{
    unsigned long long sig = 1469598103934665603ULL;
    int bidx;
    int oct;

    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        const prj_block *block = &mesh->blocks[bidx];

        sig ^= (unsigned long long)(unsigned int)(block->id + 3);
        sig *= 1099511628211ULL;
        sig ^= (unsigned long long)(unsigned int)(block->level + 7);
        sig *= 1099511628211ULL;
        sig ^= (unsigned long long)(unsigned int)(block->active + 11);
        sig *= 1099511628211ULL;
        sig ^= (unsigned long long)(unsigned int)(block->parent + 13);
        sig *= 1099511628211ULL;
        for (oct = 0; oct < 8; ++oct) {
            sig ^= (unsigned long long)(unsigned int)(block->children[oct] + 17 + oct);
            sig *= 1099511628211ULL;
        }
    }
    return sig;
}

static void prj_problem_kh_ic(double x1, double x2, double x3, double *data)
{
    const double amp = 0.1;
    const double sigma = 0.05;
    const double pressure = 2.5;
    double y = x2;
    int inner = (y > 0.25 && y < 0.75);
    double rho = inner ? 2.0 : 1.0;
    double v1 = inner ? 0.5 : -0.5;
    double d0 = (y - 0.25) / sigma;
    double d1 = (y - 0.75) / sigma;
    double v2 = amp * sin(4.0 * PRJ_KH_PI * x1) *
        (exp(-0.5 * d0 * d0) + exp(-0.5 * d1 * d1));

    (void)x3;
    data[PRJ_PRIM_RHO] = rho;
    data[PRJ_PRIM_V1] = v1;
    data[PRJ_PRIM_V2] = v2;
    data[PRJ_PRIM_V3] = 0.0;
    data[PRJ_PRIM_EINT] = pressure / ((PRJ_KH_GAMMA - 1.0) * rho);
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

static void prj_problem_fill_mesh(prj_sim *sim, prj_mpi *mpi)
{
    int bidx;

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

                    prj_problem_kh_ic(x1, x2, x3, W);
                    prj_eos_prim2cons(&sim->eos, W, U);
                    prj_problem_store_cell(block, i, j, k, W, U);
                }
            }
        }
    }
}

static void prj_problem_fill_until_amr_converged(prj_sim *sim, prj_mpi *mpi)
{
    unsigned long long prev_sig;
    unsigned long long next_sig;

    prj_problem_fill_mesh(sim, mpi);
    if (sim->mesh.max_level == 0) {
        return;
    }

    do {
        prev_sig = prj_problem_mesh_signature(&sim->mesh);
        prj_eos_fill_active_cells(&sim->mesh, &sim->eos, mpi, 1, PRJ_EOS_CTX_MAIN);
        prj_boundary_fill_ghosts(&sim->mesh, mpi, &sim->bc, 1);
        prj_eos_fill_mesh(&sim->mesh, &sim->eos, mpi, 1, PRJ_EOS_CTX_MAIN);
        prj_amr_adapt(&sim->mesh, &sim->eos, mpi);
        prj_problem_fill_mesh(sim, mpi);
        next_sig = prj_problem_mesh_signature(&sim->mesh);
    } while (next_sig != prev_sig);
}

void prj_problem_kh(prj_sim *sim, prj_mpi *mpi)
{
    if (prj_mesh_init(&sim->mesh, sim->mesh.root_nx[0], sim->mesh.root_nx[1], sim->mesh.root_nx[2],
        sim->mesh.max_level, &sim->coord, 0) != 0) {
        return;
    }
    prj_problem_fill_until_amr_converged(sim, mpi);
    prj_mhd_init(sim, mpi);
}
