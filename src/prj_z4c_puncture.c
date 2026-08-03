#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(PRJ_ENABLE_MPI)
#include <mpi.h>
#endif

#include "prj.h"

#if PRJ_DYNAMIC_GR

static void prj_z4c_puncture_fail(const char *message)
{
    fprintf(stderr, "%s\n", message);
#if defined(PRJ_ENABLE_MPI)
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#endif
    exit(EXIT_FAILURE);
}

static int prj_z4c_puncture_is_root(const prj_mpi *mpi)
{
    return mpi == 0 || mpi->rank == 0;
}

static size_t prj_z4c_puncture_stage_offset(const prj_mesh *mesh, int stage)
{
    if (mesh == 0 || stage < 0 || stage >= PRJ_BLOCK_NSTAGES) {
        prj_z4c_puncture_fail("prj_z4c_puncture: invalid tracker stage");
    }
    return (size_t)stage * (size_t)mesh->z4c_puncture_count * 3U;
}

static double *prj_z4c_puncture_stage(prj_mesh *mesh, int stage)
{
    if (mesh == 0 || mesh->z4c_puncture_count <= 0 ||
        mesh->z4c_puncture_positions_cm == 0) {
        return 0;
    }
    return mesh->z4c_puncture_positions_cm +
        prj_z4c_puncture_stage_offset(mesh, stage);
}

static const double *prj_z4c_puncture_stage_const(const prj_mesh *mesh, int stage)
{
    if (mesh == 0 || mesh->z4c_puncture_count <= 0 ||
        mesh->z4c_puncture_positions_cm == 0) {
        return 0;
    }
    return mesh->z4c_puncture_positions_cm +
        prj_z4c_puncture_stage_offset(mesh, stage);
}

int prj_z4c_puncture_count(const prj_mesh *mesh)
{
    return mesh != 0 ? mesh->z4c_puncture_count : 0;
}

int prj_z4c_puncture_position(const prj_mesh *mesh, int puncture,
    double position_cm[3])
{
    const double *positions = prj_z4c_puncture_stage_const(mesh, 0);

    if (positions == 0 || position_cm == 0 || puncture < 0 ||
        puncture >= mesh->z4c_puncture_count) {
        return 0;
    }
    memcpy(position_cm, positions + 3U * (size_t)puncture,
        3U * sizeof(*position_cm));
    return 1;
}

const double *prj_z4c_puncture_positions(const prj_mesh *mesh)
{
    return prj_z4c_puncture_stage_const(mesh, 0);
}

void prj_z4c_puncture_tracker_free(prj_mesh *mesh)
{
    if (mesh == 0) {
        return;
    }
    free(mesh->z4c_puncture_positions_cm);
    mesh->z4c_puncture_positions_cm = 0;
    mesh->z4c_puncture_count = 0;
}

void prj_z4c_puncture_tracker_init(prj_mesh *mesh, int npunctures,
    const double positions_cm[][3])
{
    size_t values_per_stage;
    size_t total_values;
    int stage;
    int p;

    if (mesh == 0 || npunctures < 0 || (npunctures > 0 && positions_cm == 0)) {
        prj_z4c_puncture_fail("prj_z4c_puncture_tracker_init: invalid puncture data");
    }
#if TIME_INTEGRATION == PRJ_TIMEINT_IMEX
    if (npunctures > 0) {
        prj_z4c_puncture_fail(
            "Z4c puncture tracking requires RK2 or an eSSPRK time integrator");
    }
#endif
    prj_z4c_puncture_tracker_free(mesh);
    if (npunctures == 0) {
        return;
    }
    if ((size_t)npunctures > SIZE_MAX / 3U) {
        prj_z4c_puncture_fail("prj_z4c_puncture_tracker_init: size overflow");
    }
    values_per_stage = 3U * (size_t)npunctures;
    if (values_per_stage > SIZE_MAX / (size_t)PRJ_BLOCK_NSTAGES) {
        prj_z4c_puncture_fail("prj_z4c_puncture_tracker_init: stage size overflow");
    }
    total_values = values_per_stage * (size_t)PRJ_BLOCK_NSTAGES;
    if (total_values > SIZE_MAX / sizeof(double)) {
        prj_z4c_puncture_fail("prj_z4c_puncture_tracker_init: allocation size overflow");
    }
    mesh->z4c_puncture_positions_cm =
        (double *)prj_calloc(total_values, sizeof(double));
    if (mesh->z4c_puncture_positions_cm == 0) {
        prj_z4c_puncture_fail("prj_z4c_puncture_tracker_init: allocation failed");
    }
    mesh->z4c_puncture_count = npunctures;
    for (p = 0; p < npunctures; ++p) {
        int d;

        for (d = 0; d < 3; ++d) {
            if (!isfinite(positions_cm[p][d])) {
                prj_z4c_puncture_fail(
                    "prj_z4c_puncture_tracker_init: non-finite initial position");
            }
        }
    }
    for (stage = 0; stage < PRJ_BLOCK_NSTAGES; ++stage) {
        memcpy(prj_z4c_puncture_stage(mesh, stage), positions_cm,
            values_per_stage * sizeof(double));
    }
}

void prj_z4c_puncture_save_stage(prj_mesh *mesh, int dst_stage, int src_stage)
{
    double *dst;
    const double *src;

    if (prj_z4c_puncture_count(mesh) == 0) {
        return;
    }
    dst = prj_z4c_puncture_stage(mesh, dst_stage);
    src = prj_z4c_puncture_stage_const(mesh, src_stage);
    memcpy(dst, src, 3U * (size_t)mesh->z4c_puncture_count * sizeof(*dst));
}

static int prj_z4c_puncture_block_contains(const prj_mesh *mesh,
    const prj_block *block, const double position[3])
{
    const double domain_max[3] = {
        mesh->coord.x1max, mesh->coord.x2max, mesh->coord.x3max
    };
    int d;

    for (d = 0; d < 3; ++d) {
        if (position[d] < block->xmin[d] || position[d] > block->xmax[d]) {
            return 0;
        }
        if (position[d] == block->xmax[d] && block->xmax[d] != domain_max[d]) {
            return 0;
        }
    }
    return 1;
}

static void prj_z4c_puncture_interp_weights(const prj_block *block, int dir,
    double coordinate, int *base, double weights[PRJ_RECON_Z4C_ORDER])
{
    const int n = PRJ_RECON_Z4C_ORDER;
    double cell_coordinate =
        (coordinate - (block->xmin[dir] + 0.5 * block->dx[dir])) /
        block->dx[dir];
    int first = (int)floor(cell_coordinate) - (n / 2 - 1);
    int q;

    if (first < 0) {
        first = 0;
    }
    if (first > PRJ_BLOCK_SIZE - n) {
        first = PRJ_BLOCK_SIZE - n;
    }
    *base = first;
    for (q = 0; q < n; ++q) {
        double weight = 1.0;
        int r;

        for (r = 0; r < n; ++r) {
            if (r != q) {
                weight *= (cell_coordinate - (double)(first + r)) /
                    (double)(q - r);
            }
        }
        weights[q] = weight;
    }
}

static void prj_z4c_puncture_interp_shift(const prj_block *block, int stage,
    const double position[3], double beta[3])
{
    const int n = PRJ_RECON_Z4C_ORDER;
    const double *z = prj_block_z4c_stage_const(block, stage);
    double wx[PRJ_RECON_Z4C_ORDER];
    double wy[PRJ_RECON_Z4C_ORDER];
    double wz[PRJ_RECON_Z4C_ORDER];
    int ibase, jbase, kbase;
    int i, j, k, d;

    if (z == 0) {
        prj_z4c_puncture_fail("prj_z4c_puncture: missing Z4c stage data");
    }
    prj_z4c_puncture_interp_weights(block, 0, position[0], &ibase, wx);
    prj_z4c_puncture_interp_weights(block, 1, position[1], &jbase, wy);
    prj_z4c_puncture_interp_weights(block, 2, position[2], &kbase, wz);
    beta[0] = beta[1] = beta[2] = 0.0;
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            for (k = 0; k < n; ++k) {
                double weight = wx[i] * wy[j] * wz[k];

                for (d = 0; d < 3; ++d) {
                    beta[d] += weight * z[Z4CIDX(PRJ_Z4C_BETAX + d,
                        ibase + i, jbase + j, kbase + k)];
                }
            }
        }
    }
}

static void prj_z4c_puncture_shift(const prj_mesh *mesh, const prj_mpi *mpi,
    int state_stage, const double *positions, double *beta)
{
    int *local_owners;
    int *global_owners;
    double *local_beta;
    int p;
    int bidx;
    size_t nvalues = 3U * (size_t)mesh->z4c_puncture_count;

    local_owners = (int *)prj_calloc((size_t)mesh->z4c_puncture_count,
        sizeof(*local_owners));
    global_owners = (int *)prj_calloc((size_t)mesh->z4c_puncture_count,
        sizeof(*global_owners));
    local_beta = (double *)prj_calloc(nvalues, sizeof(*local_beta));
    if (local_owners == 0 || global_owners == 0 || local_beta == 0) {
        prj_z4c_puncture_fail("prj_z4c_puncture: reduction allocation failed");
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        const prj_block *block = &mesh->blocks[bidx];

        if (block->id < 0 || block->active != 1 || block->z4c == 0 ||
            (mpi != 0 && block->rank != mpi->rank)) {
            continue;
        }
        for (p = 0; p < mesh->z4c_puncture_count; ++p) {
            const double *position = positions + 3U * (size_t)p;

            if (prj_z4c_puncture_block_contains(mesh, block, position)) {
                double value[3];
                int d;

                prj_z4c_puncture_interp_shift(block, state_stage, position, value);
                local_owners[p] += 1;
                for (d = 0; d < 3; ++d) {
                    local_beta[3U * (size_t)p + (size_t)d] += value[d];
                }
            }
        }
    }
#if defined(PRJ_ENABLE_MPI)
    if (mpi != 0 && mpi->totrank > 1) {
        size_t offset = 0;

        MPI_Allreduce(local_owners, global_owners, mesh->z4c_puncture_count,
            MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        while (offset < nvalues) {
            int count = nvalues - offset > (size_t)INT_MAX ?
                INT_MAX : (int)(nvalues - offset);

            MPI_Allreduce(local_beta + offset, beta + offset, count,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            offset += (size_t)count;
        }
    } else
#endif
    {
        memcpy(global_owners, local_owners,
            (size_t)mesh->z4c_puncture_count * sizeof(*global_owners));
        memcpy(beta, local_beta, nvalues * sizeof(*beta));
    }
    for (p = 0; p < mesh->z4c_puncture_count; ++p) {
        int d;

        if (global_owners[p] != 1) {
            prj_z4c_puncture_fail(
                "prj_z4c_puncture: puncture does not have exactly one active-block owner");
        }
        for (d = 0; d < 3; ++d) {
            if (!isfinite(beta[3U * (size_t)p + (size_t)d])) {
                prj_z4c_puncture_fail("prj_z4c_puncture: non-finite interpolated shift");
            }
        }
    }
    free(local_owners);
    free(global_owners);
    free(local_beta);
}

void prj_z4c_puncture_update(prj_mesh *mesh, const prj_mpi *mpi,
    int dst_stage, int a_stage, double a_weight,
    int b_stage, double b_weight, int rhs_state_stage, double dtau_cm)
{
    const double *a;
    const double *b;
    const double *rhs_positions;
    double *dst;
    double *beta;
    size_t nvalues;
    size_t q;

    if (prj_z4c_puncture_count(mesh) == 0) {
        return;
    }
    if (PRJ_FIX_Z4C) {
        prj_z4c_puncture_save_stage(mesh, dst_stage, 0);
        return;
    }
    if (!isfinite(dtau_cm)) {
        prj_z4c_puncture_fail("prj_z4c_puncture_update: non-finite timestep");
    }
    a = prj_z4c_puncture_stage_const(mesh, a_stage);
    b = prj_z4c_puncture_stage_const(mesh, b_stage);
    rhs_positions = prj_z4c_puncture_stage_const(mesh, rhs_state_stage);
    dst = prj_z4c_puncture_stage(mesh, dst_stage);
    nvalues = 3U * (size_t)mesh->z4c_puncture_count;
    beta = (double *)prj_calloc(nvalues, sizeof(*beta));
    if (beta == 0) {
        prj_z4c_puncture_fail("prj_z4c_puncture_update: allocation failed");
    }
    prj_z4c_puncture_shift(mesh, mpi, rhs_state_stage, rhs_positions, beta);
    for (q = 0; q < nvalues; ++q) {
        dst[q] = a_weight * a[q] + b_weight * b[q] - dtau_cm * beta[q];
        if (!isfinite(dst[q])) {
            prj_z4c_puncture_fail("prj_z4c_puncture_update: non-finite position");
        }
    }
    free(beta);
}

void prj_z4c_puncture_blend(prj_mesh *mesh, int saved_stage,
    double saved_weight)
{
    double *current;
    const double *saved;
    size_t nvalues;
    size_t q;

    if (prj_z4c_puncture_count(mesh) == 0) {
        return;
    }
    if (PRJ_FIX_Z4C) {
        return;
    }
    current = prj_z4c_puncture_stage(mesh, 0);
    saved = prj_z4c_puncture_stage_const(mesh, saved_stage);
    nvalues = 3U * (size_t)mesh->z4c_puncture_count;
    for (q = 0; q < nvalues; ++q) {
        current[q] = saved_weight * saved[q] +
            (1.0 - saved_weight) * current[q];
    }
}

static void prj_z4c_puncture_write(const prj_mesh *mesh, const prj_mpi *mpi,
    double time_seconds, const char *mode, int write_header)
{
    FILE *file;
    const double *positions;
    int p;

    if (!prj_z4c_puncture_is_root(mpi) || prj_z4c_puncture_count(mesh) == 0) {
        return;
    }
    file = fopen("output/z4c_punctures.txt", mode);
    if (file == 0) {
        fprintf(stderr, "failed to open output/z4c_punctures.txt: %s\n",
            strerror(errno));
        prj_z4c_puncture_fail("prj_z4c_puncture: trajectory output failed");
    }
    if (write_header) {
        fprintf(file, "# time_s puncture_id x1_cm x2_cm x3_cm\n");
    }
    positions = prj_z4c_puncture_stage_const(mesh, 0);
    for (p = 0; p < mesh->z4c_puncture_count; ++p) {
        const double *x = positions + 3U * (size_t)p;

        fprintf(file, "%.17e %d %.17e %.17e %.17e\n",
            time_seconds, p, x[0], x[1], x[2]);
    }
    fclose(file);
}

void prj_z4c_puncture_output_init(const prj_mesh *mesh, const prj_mpi *mpi,
    double time_seconds, int restart)
{
    if (prj_z4c_puncture_count(mesh) == 0) {
        return;
    }
    if (restart) {
        return;
    }
    prj_z4c_puncture_write(mesh, mpi, time_seconds, "w", 1);
}

void prj_z4c_puncture_output_append(const prj_mesh *mesh, const prj_mpi *mpi,
    double time_seconds)
{
    prj_z4c_puncture_write(mesh, mpi, time_seconds, "a", 0);
}

#else

int prj_z4c_puncture_count(const prj_mesh *mesh) { (void)mesh; return 0; }
int prj_z4c_puncture_position(const prj_mesh *mesh, int puncture,
    double position_cm[3])
{ (void)mesh; (void)puncture; (void)position_cm; return 0; }
const double *prj_z4c_puncture_positions(const prj_mesh *mesh)
{ (void)mesh; return 0; }
void prj_z4c_puncture_tracker_init(prj_mesh *mesh, int npunctures,
    const double positions_cm[][3])
{ (void)mesh; (void)npunctures; (void)positions_cm; }
void prj_z4c_puncture_tracker_free(prj_mesh *mesh) { (void)mesh; }
void prj_z4c_puncture_save_stage(prj_mesh *mesh, int dst_stage, int src_stage)
{ (void)mesh; (void)dst_stage; (void)src_stage; }
void prj_z4c_puncture_update(prj_mesh *mesh, const prj_mpi *mpi,
    int dst_stage, int a_stage, double a_weight,
    int b_stage, double b_weight, int rhs_state_stage, double dtau_cm)
{ (void)mesh; (void)mpi; (void)dst_stage; (void)a_stage; (void)a_weight;
  (void)b_stage; (void)b_weight; (void)rhs_state_stage; (void)dtau_cm; }
void prj_z4c_puncture_blend(prj_mesh *mesh, int saved_stage,
    double saved_weight)
{ (void)mesh; (void)saved_stage; (void)saved_weight; }
void prj_z4c_puncture_output_init(const prj_mesh *mesh, const prj_mpi *mpi,
    double time_seconds, int restart)
{ (void)mesh; (void)mpi; (void)time_seconds; (void)restart; }
void prj_z4c_puncture_output_append(const prj_mesh *mesh, const prj_mpi *mpi,
    double time_seconds)
{ (void)mesh; (void)mpi; (void)time_seconds; }

#endif
