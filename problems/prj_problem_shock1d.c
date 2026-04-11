#include <string.h>

#include "prj.h"

static int prj_problem_local_block(const prj_block *block)
{
    prj_mpi *mpi = prj_mpi_current();

    return block != 0 && block->id >= 0 && block->active == 1 &&
        (mpi == 0 || block->rank == mpi->rank);
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

void prj_problem_shock1d(prj_sim *sim)
{
    int root_nx1 = sim->mesh.root_nx[0];
    int root_nx2 = sim->mesh.root_nx[1];
    int root_nx3 = sim->mesh.root_nx[2];
    int max_level = sim->mesh.max_level;
    prj_coord coord = sim->coord;
    prj_bc bc = sim->bc;
    double cfl = sim->cfl;
    double t_end = sim->t_end;
    int max_steps = sim->max_steps;
    int output_interval = sim->output_interval;
    int restart_interval = sim->restart_interval;
    int amr_interval = sim->amr_interval;
    char output_dir[sizeof(sim->output_dir)];
    double amr_refine_thresh = sim->mesh.amr_refine_thresh;
    double amr_derefine_thresh = sim->mesh.amr_derefine_thresh;
    double amr_pressure_reference = sim->mesh.amr_pressure_reference;
    int bidx;

    strncpy(output_dir, sim->output_dir, sizeof(output_dir) - 1);
    output_dir[sizeof(output_dir) - 1] = '\0';
    memset(sim, 0, sizeof(*sim));
    sim->coord = coord;
    sim->bc = bc;
    sim->cfl = cfl;
    sim->t_end = t_end;
    sim->max_steps = max_steps;
    sim->output_interval = output_interval;
    sim->restart_interval = restart_interval;
    sim->amr_interval = amr_interval;
    strncpy(sim->output_dir, output_dir, sizeof(sim->output_dir) - 1);
    sim->output_dir[sizeof(sim->output_dir) - 1] = '\0';
    if (prj_mesh_init(&sim->mesh, root_nx1, root_nx2, root_nx3, max_level, &sim->coord) != 0) {
        return;
    }
    sim->mesh.amr_refine_thresh = amr_refine_thresh;
    sim->mesh.amr_derefine_thresh = amr_derefine_thresh;
    sim->mesh.amr_pressure_reference = amr_pressure_reference;

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
                    double W[PRJ_NVAR_PRIM];
                    double U[PRJ_NVAR_CONS];

                    W[PRJ_PRIM_RHO] = 1.0;
                    W[PRJ_PRIM_V1] = 0.0;
                    W[PRJ_PRIM_V2] = 0.0;
                    W[PRJ_PRIM_V3] = 0.0;
                    W[PRJ_PRIM_EINT] = 1.0 / (((5.0 / 3.0) - 1.0) * W[PRJ_PRIM_RHO]);
                    W[PRJ_PRIM_YE] = 0.1;
                    prj_eos_prim2cons(&sim->eos, W, U);
                    prj_problem_store_cell(block, i, j, k, W, U);
                }
            }
        }
    }
}
