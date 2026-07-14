#include <math.h>

#include "prj.h"

static unsigned long long prj_perturb_splitmix64(unsigned long long *state)
{
    unsigned long long z;

    *state += 0x9E3779B97F4A7C15ULL;
    z = *state;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

static double prj_perturb_uniform01(unsigned long long *state)
{
    unsigned long long r = prj_perturb_splitmix64(state) >> 11;
    return (double)r * (1.0 / 9007199254740992.0);
}

static double prj_perturb_standard_normal(unsigned long long *state)
{
    double u;
    double v;

    do {
        u = prj_perturb_uniform01(state);
    } while (u <= 1.0e-300);
    v = prj_perturb_uniform01(state);
    return sqrt(-2.0 * log(u)) * cos(2.0 * M_PI * v);
}

/* Build a cell-specific RNG state by hashing the user seed together with
 * the block id and the cell coordinates. The resulting perturbation is
 * tied to the cell's identity in the mesh rather than to the MPI rank,
 * so the same configuration produces the same perturbation regardless of
 * how blocks are distributed across ranks. */
static unsigned long long prj_perturb_cell_seed(unsigned long long seed,
    int block_id, int i, int j, int k)
{
    unsigned long long state = seed;

    (void)prj_perturb_splitmix64(&state);
    state ^= (unsigned long long)(unsigned int)block_id;
    (void)prj_perturb_splitmix64(&state);
    state ^= (unsigned long long)(unsigned int)i;
    (void)prj_perturb_splitmix64(&state);
    state ^= (unsigned long long)(unsigned int)j;
    (void)prj_perturb_splitmix64(&state);
    state ^= (unsigned long long)(unsigned int)k;
    return state;
}

void prj_set_perturbation(prj_mesh *mesh, prj_eos *eos, const prj_mpi *mpi,
    double gaussian_norm, unsigned long long seed)
{
    int bidx;

    if (mesh == 0 || gaussian_norm == 0.0) {
        return;
    }

    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        int i;
        int j;
        int k;

        if (block->id < 0 || block->active != 1) {
            continue;
        }
        if (mpi != 0 && block->rank != mpi->rank) {
            continue;
        }
        for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
            for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    unsigned long long state;
                    double gauss;
                    double factor;
                    double W[PRJ_NVAR_PRIM];
                    double U[PRJ_NVAR_CONS];

                    state = prj_perturb_cell_seed(seed, block->id, i, j, k);
                    gauss = prj_perturb_standard_normal(&state);
                    factor = 1.0 + gaussian_norm * gauss;

                    prj_block_load_prim_cell_const(block, 0, i, j, k, W);
                    W[PRJ_PRIM_RHO] *= factor;
                    prj_block_set_prim_value(block, 0, PRJ_PRIM_RHO, i, j, k, W[PRJ_PRIM_RHO]);
                    prj_block_set_prim_value(block, 1, PRJ_PRIM_RHO, i, j, k, W[PRJ_PRIM_RHO]);

                    prj_block_load_cons_cell_const(block, i, j, k, U);
                    prj_eos_cell_prim2cons(eos, mesh, block, 0, i, j, k, W, U,
                        PRJ_EOS_CTX_MAIN);
                    prj_block_store_cons_cell(block, i, j, k, U);
                }
            }
        }
    }
}
