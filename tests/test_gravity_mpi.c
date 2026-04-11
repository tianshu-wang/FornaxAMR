#include <math.h>
#include <string.h>

#include "prj.h"

#define PRJ_TEST_PI 3.14159265358979323846

static double prj_abs_double(double x)
{
    return x < 0.0 ? -x : x;
}

static int prj_local_active_block(const prj_block *block)
{
    prj_mpi *mpi = prj_mpi_current();

    return block != 0 && block->id >= 0 && block->active == 1 &&
        (mpi == 0 || block->rank == mpi->rank);
}

static void prj_fill_uniform_sphere(prj_sim *sim, double radius, double total_mass)
{
    double rho0 = total_mass / ((4.0 / 3.0) * PRJ_TEST_PI * radius * radius * radius);
    int bidx;

    for (bidx = 0; bidx < sim->mesh.nblocks; ++bidx) {
        prj_block *block = &sim->mesh.blocks[bidx];
        int i;
        int j;
        int k;

        if (!prj_local_active_block(block)) {
            continue;
        }
        for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
            for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
                for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                    double x1 = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
                    double x2 = block->xmin[1] + ((double)j + 0.5) * block->dx[1];
                    double x3 = block->xmin[2] + ((double)k + 0.5) * block->dx[2];
                    double r = sqrt(x1 * x1 + x2 * x2 + x3 * x3);
                    double W[PRJ_NVAR_PRIM];
                    double U[PRJ_NVAR_CONS];
                    int v;

                    for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
                        W[v] = 0.0;
                    }
                    W[PRJ_PRIM_RHO] = r <= radius ? rho0 : 0.0;
                    W[PRJ_PRIM_EINT] = 1.0;
                    W[PRJ_PRIM_YE] = 0.1;
                    prj_eos_prim2cons(&sim->eos, W, U);
                    for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
                        block->W[VIDX(v, i, j, k)] = W[v];
                    }
                    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
                        block->U[VIDX(v, i, j, k)] = U[v];
                    }
                }
            }
        }
    }
}

static double prj_interp_accel(const prj_grav_mono *grav_mono, double r)
{
    int idx;

    if (r <= 0.5 * (grav_mono->rf[0] + grav_mono->rf[1])) {
        return grav_mono->accel[0];
    }
    for (idx = 0; idx < grav_mono->nbins - 1; ++idx) {
        double r0 = 0.5 * (grav_mono->rf[idx] + grav_mono->rf[idx + 1]);
        double r1 = 0.5 * (grav_mono->rf[idx + 1] + grav_mono->rf[idx + 2]);

        if (r <= r1) {
            double w = (r - r0) / (r1 - r0);

            return (1.0 - w) * grav_mono->accel[idx] + w * grav_mono->accel[idx + 1];
        }
    }
    return grav_mono->accel[grav_mono->nbins - 1];
}

int main(int argc, char **argv)
{
    const double radius = 1.3;
    const double total_mass = 1.0;
    const double rho0 = total_mass / ((4.0 / 3.0) * PRJ_TEST_PI * radius * radius * radius);
    const double inner_radii[3] = {0.25, 0.55, 0.95};
    const double outer_radii[3] = {1.45, 1.70, 1.90};
    prj_sim sim;
    prj_mpi mpi;
    const prj_grav_mono *grav_mono;
    int i;

    memset(&sim, 0, sizeof(sim));
    sim.coord.x1min = -2.0;
    sim.coord.x1max = 2.0;
    sim.coord.x2min = -2.0;
    sim.coord.x2max = 2.0;
    sim.coord.x3min = -2.0;
    sim.coord.x3max = 2.0;

    prj_mpi_init(&argc, &argv, &mpi);
    if (prj_mesh_init(&sim.mesh, 8, 8, 8, 0, &sim.coord) != 0) {
        prj_mpi_finalize();
        return 1;
    }
    prj_amr_init_neighbors(&sim.mesh);
    prj_mpi_decompose(&sim.mesh);
    prj_mpi_prepare(&sim.mesh, &mpi);
    prj_gravity_init(&sim);
    prj_fill_uniform_sphere(&sim, radius, total_mass);
    prj_gravity_monopole_reduce(&sim.mesh);
    prj_gravity_monopole_integrate(&sim.mesh);
    grav_mono = prj_gravity_active_monopole();
    if (grav_mono == 0) {
        prj_gravity_free(&sim.monopole_grav);
        prj_mesh_destroy(&sim.mesh);
        prj_mpi_finalize();
        return 2;
    }

    if (grav_mono->ms[grav_mono->nbins - 1] <= 0.0 ||
        prj_abs_double((grav_mono->ms[grav_mono->nbins - 1] - total_mass) / total_mass) > 0.03) {
        prj_gravity_free(&sim.monopole_grav);
        prj_mesh_destroy(&sim.mesh);
        prj_mpi_finalize();
        return 3;
    }

    for (i = 0; i < 3; ++i) {
        double r = inner_radii[i];
        double expected = -4.0 * PRJ_TEST_PI * PRJ_GNEWT * rho0 * r / 3.0;
        double actual = prj_interp_accel(grav_mono, r);

        if (expected != 0.0 && prj_abs_double((actual - expected) / expected) > 0.15) {
            prj_gravity_free(&sim.monopole_grav);
            prj_mesh_destroy(&sim.mesh);
            prj_mpi_finalize();
            return 4;
        }
    }

    for (i = 0; i < 3; ++i) {
        double r = outer_radii[i];
        double expected = -PRJ_GNEWT * total_mass / (r * r);
        double actual = prj_interp_accel(grav_mono, r);

        if (expected != 0.0 && prj_abs_double((actual - expected) / expected) > 0.20) {
            prj_gravity_free(&sim.monopole_grav);
            prj_mesh_destroy(&sim.mesh);
            prj_mpi_finalize();
            return 5;
        }
    }

    prj_gravity_free(&sim.monopole_grav);
    prj_mesh_destroy(&sim.mesh);
    prj_mpi_finalize();
    return 0;
}
