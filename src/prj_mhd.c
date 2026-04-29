#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "prj.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#if PRJ_MHD
static inline void prj_mhd_fail(const char *message)
{
    fprintf(stderr, "%s\n", message);
    exit(EXIT_FAILURE);
}

static inline int prj_mhd_local_block(const prj_block *block)
{
    prj_mpi *mpi = prj_mpi_current();

    return block != 0 && block->id >= 0 && block->active == 1 &&
        (mpi == 0 || block->rank == mpi->rank);
}

static inline void prj_mhd_check_block_storage(const prj_block *block)
{
    int d;

    if (block == 0) {
        prj_mhd_fail("prj_mhd: block is null");
    }
    if (block->W == 0 || block->W1 == 0 || block->U == 0) {
        prj_mhd_fail("prj_mhd: missing cell-centered block storage");
    }
    for (d = 0; d < 3; ++d) {
        if (block->Bf[d] == 0 || block->Bf1[d] == 0 || block->emf[d] == 0) {
            prj_mhd_fail("prj_mhd: missing staggered MHD block storage");
        }
    }
}

static inline double prj_mhd_edge_coord(const prj_block *block, int axis, int component, int index)
{
    double offset = axis == component ? 0.5 : 0.0;

    return block->xmin[axis] + ((double)index + offset) * block->dx[axis];
}

static inline void prj_mhd_vector_potential(int init_type, double B_norm, double B_scale,
    double x, double y, double z, double A[3])
{
    double factor;
    double taper = 1.0;

    if (!isfinite(B_norm)) {
        prj_mhd_fail("prj_mhd_init: B_norm is not finite");
    }
    if (init_type == PRJ_MHD_INIT_DIPOLE_CORE) {
        double r;

        if (!isfinite(B_scale) || B_scale <= 0.0) {
            prj_mhd_fail("prj_mhd_init: dipole B_scale must be finite and positive");
        }
        r = sqrt(x * x + y * y + z * z);
        taper = 1.0 / (1.0 + pow(r / B_scale, 3.0));
    } else if (init_type != PRJ_MHD_INIT_UNIFORM) {
        prj_mhd_fail("prj_mhd_init: unknown magnetic field initialization type");
    }

    factor = 0.5 * B_norm / sqrt(4.0 * M_PI) * taper;
    A[0] = -y * factor;
    A[1] = x * factor;
    A[2] = 0.0;
}

static inline void prj_mhd_store_vector_potential(prj_block *block, int init_type,
    double B_norm, double B_scale)
{
    int d;
    int i;
    int j;
    int k;

    for (d = 0; d < 3; ++d) {
        prj_fill(block->emf[d], (size_t)PRJ_BLOCK_NCELLS, 0.0);
    }

    for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
        for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
            for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                for (d = 0; d < 3; ++d) {
                    double x = prj_mhd_edge_coord(block, X1DIR, d, i);
                    double y = prj_mhd_edge_coord(block, X2DIR, d, j);
                    double z = prj_mhd_edge_coord(block, X3DIR, d, k);
                    double A[3];

                    prj_mhd_vector_potential(init_type, B_norm, B_scale, x, y, z, A);
                    block->emf[d][IDX(i, j, k)] = A[d];
                }
            }
        }
    }
}

static inline void prj_mhd_curl_a_to_bf(prj_block *block)
{
    int d;
    int i;
    int j;
    int k;
    int ilo = -PRJ_NGHOST;
    int ihi = PRJ_BLOCK_SIZE + PRJ_NGHOST - 1;

    for (d = 0; d < 3; ++d) {
        prj_fill(block->Bf[d], (size_t)PRJ_BLOCK_NCELLS, 0.0);
    }

    for (i = ilo; i < ihi; ++i) {
        for (j = ilo; j < ihi; ++j) {
            for (k = ilo; k < ihi; ++k) {
                double B1;
                double B2;
                double B3;

                B1 = (block->emf[X3DIR][IDX(i, j + 1, k)] - block->emf[X3DIR][IDX(i, j, k)]) / block->dx[1] -
                    (block->emf[X2DIR][IDX(i, j, k + 1)] - block->emf[X2DIR][IDX(i, j, k)]) / block->dx[2];
                B2 = (block->emf[X1DIR][IDX(i, j, k + 1)] - block->emf[X1DIR][IDX(i, j, k)]) / block->dx[2] -
                    (block->emf[X3DIR][IDX(i + 1, j, k)] - block->emf[X3DIR][IDX(i, j, k)]) / block->dx[0];
                B3 = (block->emf[X2DIR][IDX(i + 1, j, k)] - block->emf[X2DIR][IDX(i, j, k)]) / block->dx[0] -
                    (block->emf[X1DIR][IDX(i, j + 1, k)] - block->emf[X1DIR][IDX(i, j, k)]) / block->dx[1];

                if (!isfinite(B1) || !isfinite(B2) || !isfinite(B3)) {
                    prj_mhd_fail("prj_mhd_init: non-finite face-centered magnetic field");
                }
                block->Bf[X1DIR][IDX(i, j, k)] = B1;
                block->Bf[X2DIR][IDX(i, j, k)] = B2;
                block->Bf[X3DIR][IDX(i, j, k)] = B3;
            }
        }
    }
}

static inline void prj_mhd_copy_bf_to_bf1(prj_block *block)
{
    int d;
    int n;

    for (d = 0; d < 3; ++d) {
        for (n = 0; n < PRJ_BLOCK_NCELLS; ++n) {
            block->Bf1[d][n] = block->Bf[d][n];
        }
    }
}

void prj_mhd_bf2bc(prj_eos *eos, prj_block *block, int use_bf1)
{
    double *W;
    double *src[3];
    int i;
    int j;
    int k;
    int d;

    if (eos == 0) {
        prj_mhd_fail("prj_mhd_bf2bc: eos is null");
    }
    prj_mhd_check_block_storage(block);

    W = use_bf1 != 0 ? block->W1 : block->W;
    for (d = 0; d < 3; ++d) {
        src[d] = use_bf1 != 0 ? block->Bf1[d] : block->Bf[d];
    }

    for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
        for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
            for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                double wcell[PRJ_NVAR_PRIM];
                double ucell[PRJ_NVAR_CONS];
                double b1;
                double b2;
                double b3;
                int v;

                b1 = 0.5 * (src[X1DIR][IDX(i, j, k)] + src[X1DIR][IDX(i + 1, j, k)]);
                b2 = 0.5 * (src[X2DIR][IDX(i, j, k)] + src[X2DIR][IDX(i, j + 1, k)]);
                b3 = 0.5 * (src[X3DIR][IDX(i, j, k)] + src[X3DIR][IDX(i, j, k + 1)]);
                if (!isfinite(b1) || !isfinite(b2) || !isfinite(b3)) {
                    prj_mhd_fail("prj_mhd_bf2bc: non-finite cell-centered magnetic field");
                }

                for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
                    wcell[v] = W[VIDX(v, i, j, k)];
                }
                wcell[PRJ_PRIM_B1] = b1;
                wcell[PRJ_PRIM_B2] = b2;
                wcell[PRJ_PRIM_B3] = b3;
                prj_eos_prim2cons(eos, wcell, ucell);
                for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
                    W[VIDX(v, i, j, k)] = wcell[v];
                }
                for (v = 0; v < PRJ_NVAR_CONS; ++v) {
                    block->U[VIDX(v, i, j, k)] = ucell[v];
                }
            }
        }
    }
}

double prj_mhd_emf_upwind(prj_block *block, int dir, int i, int j, int k,
    const double emf_face[4], const double emf_cell[4], const double v_norm[4])
{
    double emf;
    int n;

    if (block == 0 || dir < 0 || dir >= 3 || block->emf[dir] == 0) {
        prj_mhd_fail("prj_mhd_emf_upwind: invalid emf destination");
    }
    if (i < -PRJ_NGHOST || i >= PRJ_BLOCK_SIZE + PRJ_NGHOST ||
        j < -PRJ_NGHOST || j >= PRJ_BLOCK_SIZE + PRJ_NGHOST ||
        k < -PRJ_NGHOST || k >= PRJ_BLOCK_SIZE + PRJ_NGHOST) {
        prj_mhd_fail("prj_mhd_emf_upwind: edge index is out of storage bounds");
    }
    if (emf_face == 0 || emf_cell == 0 || v_norm == 0) {
        prj_mhd_fail("prj_mhd_emf_upwind: null input array");
    }
    for (n = 0; n < 4; ++n) {
        if (!isfinite(emf_face[n]) || !isfinite(emf_cell[n]) || !isfinite(v_norm[n])) {
            prj_mhd_fail("prj_mhd_emf_upwind: non-finite input");
        }
    }

    emf = 0.25 * (emf_face[0] + emf_face[1] + emf_face[2] + emf_face[3]);
    if (v_norm[0] > 0.0) {
        emf += 0.25 * (emf_face[3] - emf_cell[3]);
    } else if (v_norm[0] < 0.0) {
        emf += 0.25 * (emf_face[1] - emf_cell[0]);
    } else {
        emf += 0.125 * (emf_face[3] - emf_cell[3] + emf_face[1] - emf_cell[0]);
    }
    if (v_norm[1] > 0.0) {
        emf += 0.25 * (emf_face[2] - emf_cell[1]);
    } else if (v_norm[1] < 0.0) {
        emf += 0.25 * (emf_face[0] - emf_cell[0]);
    } else {
        emf += 0.125 * (emf_face[2] - emf_cell[1] + emf_face[0] - emf_cell[0]);
    }
    if (v_norm[2] > 0.0) {
        emf += 0.25 * (emf_face[3] - emf_cell[2]);
    } else if (v_norm[2] < 0.0) {
        emf += 0.25 * (emf_face[1] - emf_cell[1]);
    } else {
        emf += 0.125 * (emf_face[3] - emf_cell[2] + emf_face[1] - emf_cell[1]);
    }
    if (v_norm[3] > 0.0) {
        emf += 0.25 * (emf_face[2] - emf_cell[2]);
    } else if (v_norm[3] < 0.0) {
        emf += 0.25 * (emf_face[0] - emf_cell[3]);
    } else {
        emf += 0.125 * (emf_face[2] - emf_cell[2] + emf_face[0] - emf_cell[3]);
    }

    if (!isfinite(emf)) {
        prj_mhd_fail("prj_mhd_emf_upwind: non-finite result");
    }
    block->emf[dir][IDX(i, j, k)] = emf;
    return emf;
}

void prj_mhd_init(prj_sim *sim)
{
    int bidx;

    if (sim == 0) {
        prj_mhd_fail("prj_mhd_init: sim is null");
    }
    if (sim->mhd_init_type != PRJ_MHD_INIT_UNIFORM &&
        sim->mhd_init_type != PRJ_MHD_INIT_DIPOLE_CORE) {
        prj_mhd_fail("prj_mhd_init: unknown magnetic field initialization type");
    }

    for (bidx = 0; bidx < sim->mesh.nblocks; ++bidx) {
        prj_block *block = &sim->mesh.blocks[bidx];
        int d;

        if (!prj_mhd_local_block(block)) {
            continue;
        }
        prj_mhd_check_block_storage(block);
        prj_mhd_store_vector_potential(block, sim->mhd_init_type, sim->mhd_B_norm, sim->mhd_B_scale);
        prj_mhd_curl_a_to_bf(block);
        prj_mhd_copy_bf_to_bf1(block);
        for (d = 0; d < 3; ++d) {
            prj_fill(block->emf[d], (size_t)PRJ_BLOCK_NCELLS, 0.0);
        }
        prj_mhd_bf2bc(&sim->eos, block, 0);
        prj_mhd_bf2bc(&sim->eos, block, 1);
    }
}
#endif
