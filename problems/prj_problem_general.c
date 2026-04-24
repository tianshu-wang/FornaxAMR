#include <string.h>

#include "prj.h"

static int prj_problem_local_block(const prj_block *block)
{
    prj_mpi *mpi = prj_mpi_current();

    return block != 0 && block->id >= 0 && block->active == 1 &&
        (mpi == 0 || block->rank == mpi->rank);
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

void prj_problem_initial_condition(double x1, double x2, double x3, double *data)
{
    (void)x1;
    (void)x2;
    (void)x3;
    prj_fill(data, PRJ_NVAR_PRIM, 0.0);
    data[PRJ_PRIM_RHO] = 1.0;
    data[PRJ_PRIM_V1] = 0.0;
    data[PRJ_PRIM_V2] = 0.0;
    data[PRJ_PRIM_V3] = 0.0;
    data[PRJ_PRIM_EINT] = 1.0;
    data[PRJ_PRIM_YE] = 0.1;
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

static void prj_problem_fill_mesh(prj_sim *sim, void (*init_fn)(double, double, double, double *))
{
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
                    double x1 = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
                    double x2 = block->xmin[1] + ((double)j + 0.5) * block->dx[1];
                    double x3 = block->xmin[2] + ((double)k + 0.5) * block->dx[2];
                    double W[PRJ_NVAR_PRIM];
                    double U[PRJ_NVAR_CONS];

                    prj_fill(W, PRJ_NVAR_PRIM, 0.0);
                    init_fn(x1, x2, x3, W);
                    prj_eos_prim2cons(&sim->eos, W, U);
                    prj_problem_store_cell(block, i, j, k, W, U);
                }
            }
        }
    }
}

static void prj_problem_fill_until_amr_converged(prj_sim *sim)
{
    unsigned long long prev_sig;
    unsigned long long next_sig;

    prj_problem_fill_mesh(sim, prj_problem_initial_condition);
    if (sim->mesh.max_level == 0) {
        return;
    }

    do {
        prev_sig = prj_problem_mesh_signature(&sim->mesh);
        prj_eos_fill_active_cells(&sim->mesh, &sim->eos, 1);
        prj_boundary_fill_ghosts(&sim->mesh, &sim->bc, 1);
        prj_eos_fill_mesh(&sim->mesh, &sim->eos, 1);
        prj_amr_adapt(&sim->mesh, &sim->eos);
        prj_problem_fill_mesh(sim, prj_problem_initial_condition);
        next_sig = prj_problem_mesh_signature(&sim->mesh);
    } while (next_sig != prev_sig);
}

void prj_problem_general(prj_sim *sim)
{
    if (prj_mesh_init(&sim->mesh, sim->mesh.root_nx[0], sim->mesh.root_nx[1], sim->mesh.root_nx[2],
        sim->mesh.max_level, &sim->coord) != 0) {
        return;
    }
    prj_problem_fill_until_amr_converged(sim);
    prj_mhd_init(sim);
}
