#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(PRJ_ENABLE_MPI)
#include <mpi.h>
#endif

#include "prj.h"

static void prj_problem_z4c_fail(const char *message)
{
    fprintf(stderr, "%s\n", message);
#if defined(PRJ_ENABLE_MPI)
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#endif
    exit(EXIT_FAILURE);
}

#if PRJ_DYNAMIC_GR
static int prj_problem_z4c_local_block(const prj_mpi *mpi, const prj_block *block)
{
    return block != 0 && block->id >= 0 && block->active == 1 &&
        (mpi == 0 || block->rank == mpi->rank);
}

static void prj_problem_z4c_require_enabled(const prj_sim *sim, const char *name)
{
    (void)name;
    if (sim == 0 || sim->mesh.use_dynamic_gr == 0) {
        prj_problem_z4c_fail("Z4c puncture problems require use_dynamic_gr=1 in the parameter file");
    }
}

static void prj_problem_z4c_fill_ambient(prj_sim *sim, const prj_mpi *mpi)
{
    int bidx;

    for (bidx = 0; bidx < sim->mesh.nblocks; ++bidx) {
        prj_block *block = &sim->mesh.blocks[bidx];
        int i, j, k;

        if (!prj_problem_z4c_local_block(mpi, block)) {
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
                    int stage;
                    int v;

                    prj_problem_initial_condition(x1, x2, x3, W);
                    prj_eos_prim2cons(&sim->eos, W, U);
                    for (stage = 0; stage < PRJ_BLOCK_NSTAGES; ++stage) {
                        prj_block_store_prim_cell(block, stage, i, j, k, W);
                    }
                    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
                        block->U[VIDX(v, i, j, k)] = U[v];
                    }
                }
            }
        }
    }
}

static void prj_problem_z4c_prepare_mesh(prj_sim *sim, prj_mpi *mpi)
{
    if (prj_mesh_init(&sim->mesh, sim->mesh.root_nx[0], sim->mesh.root_nx[1],
        sim->mesh.root_nx[2], sim->mesh.max_level, &sim->coord, 1) != 0) {
        prj_problem_z4c_fail("Z4c puncture mesh initialization failed");
    }
    prj_mpi_decompose(&sim->mesh, mpi);
    prj_mpi_prepare(&sim->mesh, mpi);
    prj_problem_z4c_fill_ambient(sim, mpi);
}

void prj_problem_z4c_one_puncture(prj_sim *sim, prj_mpi *mpi)
{
    double centers[1][3];
    double masses[1];
    double momenta[1][3] = {{0.0, 0.0, 0.0}};

    prj_problem_z4c_require_enabled(sim, "z4c_one_puncture");
    prj_problem_z4c_prepare_mesh(sim, mpi);

    centers[0][0] = sim->mesh.z4c_params.puncture_center_cm[0];
    centers[0][1] = sim->mesh.z4c_params.puncture_center_cm[1];
    centers[0][2] = sim->mesh.z4c_params.puncture_center_cm[2];
    masses[0] = sim->mesh.z4c_params.puncture_mass_cm;

    prj_z4c_init_punctures(&sim->mesh, mpi, 1, centers, masses, momenta,
        sim->mesh.z4c_params.puncture_floor_radius_cm);
    prj_mhd_init(sim, mpi);
}

void prj_problem_z4c_two_puncture(prj_sim *sim, prj_mpi *mpi)
{
    const double b = sim->mesh.z4c_params.puncture_half_separation_cm;
    double centers[2][3] = {
        { b, 0.0, 0.0 },
        {-b, 0.0, 0.0 }
    };
    double masses[2] = {
        sim->mesh.z4c_params.puncture_mass_plus_cm,
        sim->mesh.z4c_params.puncture_mass_minus_cm
    };
    double momenta[2][3];

    prj_problem_z4c_require_enabled(sim, "z4c_two_puncture");
    prj_problem_z4c_prepare_mesh(sim, mpi);

    memcpy(momenta[0], sim->mesh.z4c_params.puncture_momentum_plus_cm,
        sizeof(momenta[0]));
    memcpy(momenta[1], sim->mesh.z4c_params.puncture_momentum_minus_cm,
        sizeof(momenta[1]));

    prj_z4c_init_punctures(&sim->mesh, mpi, 2, centers, masses, momenta,
        sim->mesh.z4c_params.puncture_floor_radius_cm);
    prj_mhd_init(sim, mpi);
}

#else

void prj_problem_z4c_one_puncture(prj_sim *sim, prj_mpi *mpi)
{
    (void)sim;
    (void)mpi;
    prj_problem_z4c_fail("z4c_one_puncture requires rebuilding with DYNAMIC_GR=1");
}

void prj_problem_z4c_two_puncture(prj_sim *sim, prj_mpi *mpi)
{
    (void)sim;
    (void)mpi;
    prj_problem_z4c_fail("z4c_two_puncture requires rebuilding with DYNAMIC_GR=1");
}

#endif
