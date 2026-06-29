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

static inline int prj_mhd_local_block(const prj_mpi *mpi, const prj_block *block)
{
    return block != 0 && block->id >= 0 && block->active == 1 &&
        (mpi == 0 || block->rank == mpi->rank);
}

static inline int prj_mhd_initialized_storage_block(const prj_block *block)
{
    return block != 0 && block->id >= 0 && block->active == 1 &&
        block->W != 0 && block->Bf[0] != 0 && block->emf[0] != 0;
}

static inline void prj_mhd_check_block_storage(const prj_block *block)
{
    int d;

    if (block == 0) {
        prj_mhd_fail("prj_mhd: block is null");
    }
    if (block->W == 0) {
        prj_mhd_fail("prj_mhd: missing cell-centered block storage");
    }
    for (d = 0; d < 3; ++d) {
        if (block->Bf[d] == 0 || block->emf[d] == 0) {
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
        prj_fill(block->emf[d], (size_t)PRJ_BLOCK_NEDGES, 0.0);
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
                    block->emf[d][EDGE_IDX(d, i, j, k)] = A[d];
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
        prj_fill(block->Bf[d], (size_t)PRJ_BLOCK_NFACES, 0.0);
    }

    for (i = ilo; i < ihi; ++i) {
        for (j = ilo; j < ihi; ++j) {
            for (k = ilo; k < ihi; ++k) {
                double B1;
                double B2;
                double B3;

                B1 = (block->emf[X3DIR][EDGE_IDX(X3DIR, i, j + 1, k)] - block->emf[X3DIR][EDGE_IDX(X3DIR, i, j, k)]) / block->dx[1] -
                     (block->emf[X2DIR][EDGE_IDX(X2DIR, i, j, k + 1)] - block->emf[X2DIR][EDGE_IDX(X2DIR, i, j, k)]) / block->dx[2];
                B2 = (block->emf[X1DIR][EDGE_IDX(X1DIR, i, j, k + 1)] - block->emf[X1DIR][EDGE_IDX(X1DIR, i, j, k)]) / block->dx[2] -
                     (block->emf[X3DIR][EDGE_IDX(X3DIR, i + 1, j, k)] - block->emf[X3DIR][EDGE_IDX(X3DIR, i, j, k)]) / block->dx[0];
                B3 = (block->emf[X2DIR][EDGE_IDX(X2DIR, i + 1, j, k)] - block->emf[X2DIR][EDGE_IDX(X2DIR, i, j, k)]) / block->dx[0] -
                     (block->emf[X1DIR][EDGE_IDX(X1DIR, i, j + 1, k)] - block->emf[X1DIR][EDGE_IDX(X1DIR, i, j, k)]) / block->dx[1];

                if (!isfinite(B1) || !isfinite(B2) || !isfinite(B3)) {
                    prj_mhd_fail("prj_mhd_init: non-finite face-centered magnetic field");
                }
                block->Bf[X1DIR][FACE_IDX(X1DIR, i, j, k)] = B1;
                block->Bf[X2DIR][FACE_IDX(X2DIR, i, j, k)] = B2;
                block->Bf[X3DIR][FACE_IDX(X3DIR, i, j, k)] = B3;
            }
        }
    }
}

static inline void prj_mhd_copy_bf_to_bf1(prj_block *block)
{
    int d;
    int n;

    for (d = 0; d < 3; ++d) {
        double *bf1 = prj_block_bf_stage(block, d, 1);

        for (n = 0; n < PRJ_BLOCK_NFACES; ++n) {
            bf1[n] = block->Bf[d][n];
        }
    }
}

static inline void prj_mhd_check_bf_storage(const prj_block *block)
{
    int d;

    if (block == 0) {
        prj_mhd_fail("prj_mhd_check_bf_storage: block is null");
    }
    for (d = 0; d < 3; ++d) {
        if (block->face_fidelity[d] == 0) {
            prj_mhd_fail("prj_mhd_check_bf_storage: missing face fidelity storage");
        }
        if (block->Bf[d] == 0) {
            prj_mhd_fail("prj_mhd_check_bf_storage: missing face-centered magnetic field storage");
        }
    }
}

static inline double *prj_mhd_bf_array(prj_block *block, int dir, int use_bf1)
{
    return prj_block_bf_stage(block, dir, use_bf1 != 0 ? 1 : 0);
}

static inline int prj_mhd_face_storage_index_ok(int dir, int i, int j, int k)
{
    int idx[3] = {i, j, k};
    int d;

    if (dir < 0 || dir >= 3) {
        return 0;
    }
    for (d = 0; d < 3; ++d) {
        int lo = -PRJ_NGHOST;
        int hi = dir == d ? PRJ_BLOCK_SIZE + PRJ_NGHOST : PRJ_BLOCK_SIZE + PRJ_NGHOST - 1;
        if (idx[d] < lo || idx[d] > hi) {
            return 0;
        }
    }
    return 1;
}

static inline void prj_mhd_write_bf_face(prj_block *block, int use_bf1,
    int dir, int i, int j, int k, double value, int fidelity)
{
    int idx;
    double *dst;

    prj_mhd_check_bf_storage(block);
    if (!prj_mhd_face_storage_index_ok(dir, i, j, k)) {
        fprintf(stderr, "prj_mhd_write_bf_face: invalid face index dir=%d i=%d j=%d k=%d\n",
            dir, i, j, k);
        exit(EXIT_FAILURE);
    }
    if (!isfinite(value) || fidelity < PRJ_MHD_FIDELITY_NONE ||
        fidelity > PRJ_MHD_FIDELITY_FINER) {
        prj_mhd_fail("prj_mhd_write_bf_face: invalid value or fidelity");
    }
    idx = FACE_IDX(dir, i, j, k);
    if (fidelity < block->face_fidelity[dir][idx]) {
        return;
    }
    dst = prj_mhd_bf_array(block, dir, use_bf1);
    dst[idx] = value;
    if (fidelity > block->face_fidelity[dir][idx]) {
        block->face_fidelity[dir][idx] = fidelity;
    }
}

static inline double prj_mhd_buf_read(const double *buf,
    const int buf_lo[3], const int buf_n[3], int i, int j, int k)
{
    return buf[((i - buf_lo[0]) * buf_n[1] + (j - buf_lo[1])) * buf_n[2] +
               (k - buf_lo[2])];
}

static inline double prj_mhd_buf_sign_half(int bit)
{
    return bit == 0 ? -1.0 : 1.0;
}

static inline double prj_mhd_interp_x1_buf(const double *buf,
    const int buf_lo[3], const int buf_n[3], const double dx[3],
    int i, int j, int k, int fine_j, int fine_k, double area, int use_BJ)
{
    double stencil[9];
    double target[2];
    double value;
    int dj;
    int dk;

    (void)dx;
    for (dj = -1; dj <= 1; ++dj) {
        for (dk = -1; dk <= 1; ++dk) {
            stencil[prj_reconstruct_stencil2_index(dj, dk)] =
                prj_mhd_buf_read(buf, buf_lo, buf_n, i, j + dj, k + dk);
        }
    }
    target[0] = 0.25 * prj_mhd_buf_sign_half(fine_j);
    target[1] = 0.25 * prj_mhd_buf_sign_half(fine_k);
    value = prj_reconstruct_face_for_prolongate(stencil, target, use_BJ);

    return value * area;
}

static inline double prj_mhd_interp_x2_buf(const double *buf,
    const int buf_lo[3], const int buf_n[3], const double dx[3],
    int i, int j, int k, int fine_i, int fine_k, double area, int use_BJ)
{
    double stencil[9];
    double target[2];
    double value;
    int di;
    int dk;

    (void)dx;
    for (di = -1; di <= 1; ++di) {
        for (dk = -1; dk <= 1; ++dk) {
            stencil[prj_reconstruct_stencil2_index(di, dk)] =
                prj_mhd_buf_read(buf, buf_lo, buf_n, i + di, j, k + dk);
        }
    }
    target[0] = 0.25 * prj_mhd_buf_sign_half(fine_i);
    target[1] = 0.25 * prj_mhd_buf_sign_half(fine_k);
    value = prj_reconstruct_face_for_prolongate(stencil, target, use_BJ);

    return value * area;
}

static inline double prj_mhd_interp_x3_buf(const double *buf,
    const int buf_lo[3], const int buf_n[3], const double dx[3],
    int i, int j, int k, int fine_i, int fine_j, double area, int use_BJ)
{
    double stencil[9];
    double target[2];
    double value;
    int di;
    int dj;

    (void)dx;
    for (di = -1; di <= 1; ++di) {
        for (dj = -1; dj <= 1; ++dj) {
            stencil[prj_reconstruct_stencil2_index(di, dj)] =
                prj_mhd_buf_read(buf, buf_lo, buf_n, i + di, j + dj, k);
        }
    }
    target[0] = 0.25 * prj_mhd_buf_sign_half(fine_i);
    target[1] = 0.25 * prj_mhd_buf_sign_half(fine_j);
    value = prj_reconstruct_face_for_prolongate(stencil, target, use_BJ);

    return value * area;
}

static inline void prj_mhd_compute_buffer_inner_fluxes(double u[3][2][2],
    double v[2][3][2], double w[2][2][3], double dx1, double dx2,
    double dx3)
{
    double Uxx = 0.0;
    double Vyy = 0.0;
    double Wzz = 0.0;
    double Uxyz = 0.0;
    double Vxyz = 0.0;
    double Wxyz = 0.0;
    double dx1s = dx1 * dx1;
    double dx2s = dx2 * dx2;
    double dx3s = dx3 * dx3;
    int i;
    int j;
    int k;

    for (i = 0; i < 2; ++i) {
        double si = prj_mhd_buf_sign_half(i);
        for (j = 0; j < 2; ++j) {
            double sj = prj_mhd_buf_sign_half(j);
            for (k = 0; k < 2; ++k) {
                double sk = prj_mhd_buf_sign_half(k);
                Uxx += si * sj * v[i][2 * j][k] + si * sk * w[i][j][2 * k];
                Vyy += sj * sk * w[i][j][2 * k] + sj * si * u[2 * i][j][k];
                Wzz += sk * si * u[2 * i][j][k] + sk * sj * v[i][2 * j][k];
                Uxyz += si * sj * sk * u[2 * i][j][k] / (dx2s + dx3s);
                Vxyz += si * sj * sk * v[i][2 * j][k] / (dx1s + dx3s);
                Wxyz += si * sj * sk * w[i][j][2 * k] / (dx1s + dx2s);
            }
        }
    }
    Uxx *= 0.125;
    Vyy *= 0.125;
    Wzz *= 0.125;
    Uxyz *= 0.125;
    Vxyz *= 0.125;
    Wxyz *= 0.125;

    for (j = 0; j < 2; ++j) {
        double sj = prj_mhd_buf_sign_half(j);
        for (k = 0; k < 2; ++k) {
            double sk = prj_mhd_buf_sign_half(k);
            u[1][j][k] = 0.5 * (u[2][j][k] + u[0][j][k]) +
                Uxx + sk * dx3s * Vxyz + sj * dx2s * Wxyz;
        }
    }
    for (i = 0; i < 2; ++i) {
        double si = prj_mhd_buf_sign_half(i);
        for (k = 0; k < 2; ++k) {
            double sk = prj_mhd_buf_sign_half(k);
            v[i][1][k] = 0.5 * (v[i][2][k] + v[i][0][k]) +
                Vyy + si * dx1s * Wxyz + sk * dx3s * Uxyz;
        }
    }
    for (i = 0; i < 2; ++i) {
        double si = prj_mhd_buf_sign_half(i);
        for (j = 0; j < 2; ++j) {
            double sj = prj_mhd_buf_sign_half(j);
            w[i][j][1] = 0.5 * (w[i][j][2] + w[i][j][0]) +
                Wzz + sj * dx2s * Uxyz + si * dx1s * Vxyz;
        }
    }
}

void prj_mhd_prolong_bf_from_buffer(const double *buf[3],
    int buf_lo[3][3], int buf_n[3][3], const double coarse_dx[3],
    prj_block *fine, int ci, int cj, int ck, int fi, int fj, int fk,
    int use_bf1, int use_BJ)
{
    double u[3][2][2];
    double v[2][3][2];
    double w[2][2][3];
    double *dst[3];
    double area_u;
    double area_v;
    double area_w;
    int i;
    int j;
    int k;
    int d;

    prj_mhd_check_bf_storage(fine);

    for (d = 0; d < 3; ++d) {
        dst[d] = prj_mhd_bf_array(fine, d, use_bf1);
    }

    area_u = fine->dx[1] * fine->dx[2];
    area_v = fine->dx[0] * fine->dx[2];
    area_w = fine->dx[0] * fine->dx[1];

    for (j = 0; j < 2; ++j) {
        for (k = 0; k < 2; ++k) {
            double flux;
            int idx;

            flux = prj_mhd_interp_x1_buf(buf[X1DIR], buf_lo[X1DIR],
                buf_n[X1DIR], coarse_dx, ci, cj, ck, j, k, area_u, use_BJ);
            idx = FACE_IDX(X1DIR, fi, fj + j, fk + k);
            if (fine->face_fidelity[X1DIR][idx] > PRJ_MHD_FIDELITY_COARSER) {
                u[0][j][k] = dst[X1DIR][idx] * area_u;
            } else {
                u[0][j][k] = flux;
                dst[X1DIR][idx] = flux / area_u;
                fine->face_fidelity[X1DIR][idx] = PRJ_MHD_FIDELITY_COARSER;
            }

            flux = prj_mhd_interp_x1_buf(buf[X1DIR], buf_lo[X1DIR],
                buf_n[X1DIR], coarse_dx, ci + 1, cj, ck, j, k, area_u, use_BJ);
            idx = FACE_IDX(X1DIR, fi + 2, fj + j, fk + k);
            if (fine->face_fidelity[X1DIR][idx] > PRJ_MHD_FIDELITY_COARSER) {
                u[2][j][k] = dst[X1DIR][idx] * area_u;
            } else {
                u[2][j][k] = flux;
                dst[X1DIR][idx] = flux / area_u;
                fine->face_fidelity[X1DIR][idx] = PRJ_MHD_FIDELITY_COARSER;
            }
        }
    }

    for (i = 0; i < 2; ++i) {
        for (k = 0; k < 2; ++k) {
            double flux;
            int idx;

            flux = prj_mhd_interp_x2_buf(buf[X2DIR], buf_lo[X2DIR],
                buf_n[X2DIR], coarse_dx, ci, cj, ck, i, k, area_v, use_BJ);
            idx = FACE_IDX(X2DIR, fi + i, fj, fk + k);
            if (fine->face_fidelity[X2DIR][idx] > PRJ_MHD_FIDELITY_COARSER) {
                v[i][0][k] = dst[X2DIR][idx] * area_v;
            } else {
                v[i][0][k] = flux;
                dst[X2DIR][idx] = flux / area_v;
                fine->face_fidelity[X2DIR][idx] = PRJ_MHD_FIDELITY_COARSER;
            }

            flux = prj_mhd_interp_x2_buf(buf[X2DIR], buf_lo[X2DIR],
                buf_n[X2DIR], coarse_dx, ci, cj + 1, ck, i, k, area_v, use_BJ);
            idx = FACE_IDX(X2DIR, fi + i, fj + 2, fk + k);
            if (fine->face_fidelity[X2DIR][idx] > PRJ_MHD_FIDELITY_COARSER) {
                v[i][2][k] = dst[X2DIR][idx] * area_v;
            } else {
                v[i][2][k] = flux;
                dst[X2DIR][idx] = flux / area_v;
                fine->face_fidelity[X2DIR][idx] = PRJ_MHD_FIDELITY_COARSER;
            }
        }
    }

    for (i = 0; i < 2; ++i) {
        for (j = 0; j < 2; ++j) {
            double flux;
            int idx;

            flux = prj_mhd_interp_x3_buf(buf[X3DIR], buf_lo[X3DIR],
                buf_n[X3DIR], coarse_dx, ci, cj, ck, i, j, area_w, use_BJ);
            idx = FACE_IDX(X3DIR, fi + i, fj + j, fk);
            if (fine->face_fidelity[X3DIR][idx] > PRJ_MHD_FIDELITY_COARSER) {
                w[i][j][0] = dst[X3DIR][idx] * area_w;
            } else {
                w[i][j][0] = flux;
                dst[X3DIR][idx] = flux / area_w;
                fine->face_fidelity[X3DIR][idx] = PRJ_MHD_FIDELITY_COARSER;
            }

            flux = prj_mhd_interp_x3_buf(buf[X3DIR], buf_lo[X3DIR],
                buf_n[X3DIR], coarse_dx, ci, cj, ck + 1, i, j, area_w, use_BJ);
            idx = FACE_IDX(X3DIR, fi + i, fj + j, fk + 2);
            if (fine->face_fidelity[X3DIR][idx] > PRJ_MHD_FIDELITY_COARSER) {
                w[i][j][2] = dst[X3DIR][idx] * area_w;
            } else {
                w[i][j][2] = flux;
                dst[X3DIR][idx] = flux / area_w;
                fine->face_fidelity[X3DIR][idx] = PRJ_MHD_FIDELITY_COARSER;
            }
        }
    }

    prj_mhd_compute_buffer_inner_fluxes(u, v, w,
        fine->dx[0], fine->dx[1], fine->dx[2]);

    for (j = 0; j < 2; ++j) {
        for (k = 0; k < 2; ++k) {
            prj_mhd_write_bf_face(fine, use_bf1, X1DIR,
                fi + 1, fj + j, fk + k,
                u[1][j][k] / area_u, PRJ_MHD_FIDELITY_COARSER);
        }
    }
    for (i = 0; i < 2; ++i) {
        for (k = 0; k < 2; ++k) {
            prj_mhd_write_bf_face(fine, use_bf1, X2DIR,
                fi + i, fj + 1, fk + k,
                v[i][1][k] / area_v, PRJ_MHD_FIDELITY_COARSER);
        }
    }
    for (i = 0; i < 2; ++i) {
        for (j = 0; j < 2; ++j) {
            prj_mhd_write_bf_face(fine, use_bf1, X3DIR,
                fi + i, fj + j, fk + 1,
                w[i][j][1] / area_w, PRJ_MHD_FIDELITY_COARSER);
        }
    }
}

static inline void prj_mhd_check_edge_storage_index(const char *label, int dir, int i, int j, int k)
{
    int idx[3] = {i, j, k};
    int d;

    for (d = 0; d < 3; ++d) {
        int lo = -PRJ_NGHOST;
        int hi = dir == d ? PRJ_BLOCK_SIZE + PRJ_NGHOST - 1 : PRJ_BLOCK_SIZE + PRJ_NGHOST;
        if (idx[d] < lo || idx[d] > hi) {
            fprintf(stderr,
                "%s: index out of edge storage bounds dir=%d i=%d j=%d k=%d (axis %d valid [%d, %d])\n",
                label, dir, i, j, k, d, lo, hi);
            exit(EXIT_FAILURE);
        }
    }
}

static void prj_mhd_bf2bc_impl(prj_eos *eos, prj_block *block, int use_bf1,
    int use_cell_derived_mask)
{
    double *W;
    double *src[3];
    int i;
    int j;
    int k;
    int d;

    (void)eos;
    prj_mhd_check_block_storage(block);

    W = prj_block_prim_stage(block, use_bf1 != 0 ? 1 : 0);
    for (d = 0; d < 3; ++d) {
        src[d] = prj_block_bf_stage(block, d, use_bf1 != 0 ? 1 : 0);
    }

    for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
        for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
            for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                double b1;
                double b2;
                double b3;

                if (use_cell_derived_mask != 0 &&
                    block->cell_derived_done != 0 &&
                    block->cell_derived_done[IDX(i, j, k)] != 0) {
                    continue;
                }
                b1 = 0.5 * (src[X1DIR][FACE_IDX(X1DIR, i, j, k)] + src[X1DIR][FACE_IDX(X1DIR, i + 1, j, k)]);
                b2 = 0.5 * (src[X2DIR][FACE_IDX(X2DIR, i, j, k)] + src[X2DIR][FACE_IDX(X2DIR, i, j + 1, k)]);
                b3 = 0.5 * (src[X3DIR][FACE_IDX(X3DIR, i, j, k)] + src[X3DIR][FACE_IDX(X3DIR, i, j, k + 1)]);
                if (!isfinite(b1) || !isfinite(b2) || !isfinite(b3)) {
                    prj_mhd_fail("prj_mhd_bf2bc: non-finite cell-centered magnetic field");
                }

                W[WIDX(PRJ_PRIM_B1, i, j, k)] = b1;
                W[WIDX(PRJ_PRIM_B2, i, j, k)] = b2;
                W[WIDX(PRJ_PRIM_B3, i, j, k)] = b3;
            }
        }
    }
}

void prj_mhd_bf2bc(prj_eos *eos, prj_block *block, int use_bf1)
{
    prj_mhd_bf2bc_impl(eos, block, use_bf1, 1);
}

void prj_mhd_bf2bc_all(prj_eos *eos, prj_block *block, int use_bf1)
{
    prj_mhd_bf2bc_impl(eos, block, use_bf1, 0);
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
    block->emf[dir][EDGE_IDX(dir, i, j, k)] = emf;
    return emf;
}

static inline int prj_mhd_edge_axis_active_max(int dir, int axis)
{
    return dir == axis ? PRJ_BLOCK_SIZE - 1 : PRJ_BLOCK_SIZE;
}

static inline int prj_mhd_edge_active(int dir, int i, int j, int k)
{
    int idx[3];
    int d;

    idx[0] = i;
    idx[1] = j;
    idx[2] = k;
    for (d = 0; d < 3; ++d) {
        if (idx[d] < 0 || idx[d] > prj_mhd_edge_axis_active_max(dir, d)) {
            return 0;
        }
    }
    return 1;
}

static inline void prj_mhd_check_emf_storage(const prj_block *block)
{
    int d;

    if (block == 0) {
        prj_mhd_fail("prj_mhd_emf_send: missing edge fidelity storage");
    }
    for (d = 0; d < 3; ++d) {
        if (block->edge_fidelity[d] == 0) {
            prj_mhd_fail("prj_mhd_emf_send: missing edge fidelity storage");
        }
        if (block->emf[d] == 0) {
            prj_mhd_fail("prj_mhd_emf_send: missing emf storage");
        }
    }
}

static inline int prj_mhd_edge_point_inside(const prj_block *block, const double x[3])
{
    const double tol = 1.0e-12;
    int d;

    if (block == 0) {
        return 0;
    }
    for (d = 0; d < 3; ++d) {
        if (x[d] < block->xmin[d] - tol || x[d] > block->xmax[d] + tol) {
            return 0;
        }
    }
    return 1;
}

static inline int prj_mhd_nearest_int(double x)
{
    return x >= 0.0 ? (int)(x + 0.5) : (int)(x - 0.5);
}

static inline double prj_mhd_restrict_emf_value(const prj_block *fine, int dir,
    const double x[3])
{
    const double tol = 1.0e-8;
    int idx[3] = {0, 0, 0};
    int edge_base = 0;
    double sum = 0.0;
    int n;
    int d;

    prj_mhd_check_emf_storage(fine);
    for (d = 0; d < 3; ++d) {
        double q;

        if (fine->dx[d] <= 0.0) {
            prj_mhd_fail("prj_mhd_emf_send: invalid fine cell size");
        }
        q = (x[d] - fine->xmin[d]) / fine->dx[d];
        if (d == dir) {
            double r = q - 0.5;

            edge_base = prj_mhd_nearest_int(r - 0.5);
            if (fabs(r - ((double)edge_base + 0.5)) > tol) {
                prj_mhd_fail("prj_mhd_emf_send: coarse edge is not centered on two fine edges");
            }
            idx[d] = edge_base;
        } else {
            int edge_idx = prj_mhd_nearest_int(q);

            if (fabs(q - (double)edge_idx) > tol) {
                prj_mhd_fail("prj_mhd_emf_send: coarse edge is not fine-aligned");
            }
            idx[d] = edge_idx;
        }
    }
    for (n = 0; n < 2; ++n) {
        int eidx[3];
        double value;

        eidx[0] = idx[0];
        eidx[1] = idx[1];
        eidx[2] = idx[2];
        eidx[dir] = edge_base + n;
        prj_mhd_check_edge_storage_index("prj_mhd_emf_send: fine edge", dir, eidx[0], eidx[1], eidx[2]);
        value = fine->emf[dir][EDGE_IDX(dir, eidx[0], eidx[1], eidx[2])];
        if (!isfinite(value)) {
            prj_mhd_fail("prj_mhd_emf_send: non-finite fine emf");
        }
        sum += value;
    }
    return 0.5 * sum;
}

static inline void prj_mhd_write_emf_edge(prj_block *block, int dir,
    int i, int j, int k, double value, int fidelity)
{
    int idx;

    prj_mhd_check_emf_storage(block);
    if (dir < 0 || dir >= 3 || !isfinite(value) ||
        fidelity < PRJ_MHD_FIDELITY_NONE || fidelity > PRJ_MHD_FIDELITY_FINER) {
        prj_mhd_fail("prj_mhd_emf_send: invalid emf write");
    }
    prj_mhd_check_edge_storage_index("prj_mhd_emf_send: destination edge", dir, i, j, k);
    idx = EDGE_IDX(dir, i, j, k);
    if (fidelity < block->edge_fidelity[dir][idx]) {
        return;
    }
    block->emf[dir][idx] = value;
    if (fidelity > block->edge_fidelity[dir][idx]) {
        block->edge_fidelity[dir][idx] = fidelity;
    }
}

static inline void prj_mhd_restrict_emf_to_coarse(const prj_block *fine,
    prj_block *coarse)
{
    int dir;

    prj_mhd_check_emf_storage(fine);
    prj_mhd_check_emf_storage(coarse);

    int axis[3] = {0,0,0};
    int touch = 0;
    int d;

    for (d = 0; d < 3; ++d) {
        if (fabs(fine->xmax[d] - coarse->xmin[d]) < 1.0e-12*fine->dx[d]) {
            axis[d] = 1;
            touch += 1;
        } else if (fabs(coarse->xmax[d] - fine->xmin[d]) < 1.0e-12*fine->dx[d]) {
            axis[d] = -1;
            touch += 1;
        } else {
            axis[d] = 0;
        }
    }
    if (touch==0||touch>2){return;}


    int it_recv[3] = {-100,-100,-100};
    int it_send[3] = {-100,-100,-100};

    for (d = 0; d < 3; ++d) {
        if (axis[d]==1) {
            it_send[d] = PRJ_BLOCK_SIZE;
            it_recv[d] = 0;
        }
        if (axis[d]==-1) {
            it_send[d] = 0;
            it_recv[d] = PRJ_BLOCK_SIZE;
        }
    }

    for (dir = 0; dir < 3; ++dir) {
        if (axis[dir]!=0) {continue;}
        int i_offset = ((fine->xmin[dir]+fine->xmax[dir]) < (coarse->xmin[dir]+coarse->xmax[dir])) ? 0 : PRJ_BLOCK_SIZE/2;

        int tan0 = (dir+1)%3;
        int tan1 = (dir+2)%3;

        
        if (axis[tan0]!=0&&axis[tan1]!=0) {
            int i;
            for (i = 0; i<PRJ_BLOCK_SIZE; i+=2) {
                double value = 0;
                it_send[dir] = i;
                value += fine->emf[dir][EDGE_IDX(dir, it_send[0], it_send[1], it_send[2])];   
                it_send[dir] = i+1;
                value += fine->emf[dir][EDGE_IDX(dir, it_send[0], it_send[1], it_send[2])];   

                value *= 0.5;

                it_recv[dir] = i/2 + i_offset;
                prj_mhd_write_emf_edge(coarse, dir, it_recv[0], it_recv[1], it_recv[2], value,
                    PRJ_MHD_FIDELITY_FINER);
            }
        } else if (axis[tan0]!=0) {
            int i;
            int j;
            int j_offset = ((fine->xmin[tan1]+fine->xmax[tan1]) < (coarse->xmin[tan1]+coarse->xmax[tan1])) ? 0 : PRJ_BLOCK_SIZE/2;
            for (i = 0; i<PRJ_BLOCK_SIZE; i+=2) {
                for (j = 0; j<=PRJ_BLOCK_SIZE; j+=2) {
                    it_send[tan1] = j;
                    it_recv[tan1] = j/2 + j_offset;

                    double value = 0;
                    it_send[dir] = i;
                    value += fine->emf[dir][EDGE_IDX(dir, it_send[0], it_send[1], it_send[2])];   
                    it_send[dir] = i+1;
                    value += fine->emf[dir][EDGE_IDX(dir, it_send[0], it_send[1], it_send[2])];   

                    value *= 0.5;

                    it_recv[dir] = i/2 + i_offset;
                    prj_mhd_write_emf_edge(coarse, dir, it_recv[0], it_recv[1], it_recv[2], value,
                        PRJ_MHD_FIDELITY_FINER);
                }
            }
        } else if (axis[tan1]!=0) {
            int i;
            int j;
            int j_offset = ((fine->xmin[tan0]+fine->xmax[tan0]) < (coarse->xmin[tan0]+coarse->xmax[tan0])) ? 0 : PRJ_BLOCK_SIZE/2;
            for (i = 0; i<PRJ_BLOCK_SIZE; i+=2) {
                for (j = 0; j<=PRJ_BLOCK_SIZE; j+=2) {
                    it_send[tan0] = j;
                    it_recv[tan0] = j/2 + j_offset;

                    double value = 0;
                    it_send[dir] = i;
                    value += fine->emf[dir][EDGE_IDX(dir, it_send[0], it_send[1], it_send[2])];   
                    it_send[dir] = i+1;
                    value += fine->emf[dir][EDGE_IDX(dir, it_send[0], it_send[1], it_send[2])];   

                    value *= 0.5;
                    
                    it_recv[dir] = i/2 + i_offset;
                    prj_mhd_write_emf_edge(coarse, dir, it_recv[0], it_recv[1], it_recv[2], value,
                        PRJ_MHD_FIDELITY_FINER);
                }
            }
        } else {
            fprintf(stderr,"Unknown emf edge type\n"); 
            exit(1);
        }
    }
}

static inline void prj_mhd_init_edge_fidelity(prj_mesh *mesh, const prj_mpi *mpi)
{
    int bidx;

    if (mesh == 0) {
        prj_mhd_fail("prj_mhd_emf_send: mesh is null");
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        int dir;

        if (!prj_mhd_local_block(mpi, block)) {
            continue;
        }
        prj_mhd_check_emf_storage(block);
        for (dir = 0; dir < 3; ++dir) {
            int i;
            int j;
            int k;
            int n;

            for (n = 0; n < PRJ_BLOCK_NEDGES; ++n) {
                block->edge_fidelity[dir][n] = PRJ_MHD_FIDELITY_NONE;
            }

            for (i = 0; i <= prj_mhd_edge_axis_active_max(dir, 0); ++i) {
                for (j = 0; j <= prj_mhd_edge_axis_active_max(dir, 1); ++j) {
                    for (k = 0; k <= prj_mhd_edge_axis_active_max(dir, 2); ++k) {
                        block->edge_fidelity[dir][EDGE_IDX(dir, i, j, k)] = PRJ_MHD_FIDELITY_SAME;
                    }
                }
            }
        }
    }
}

void prj_mhd_emf_send(prj_mesh *mesh, const prj_mpi *mpi)
{
    int bidx;

    prj_mhd_init_edge_fidelity(mesh, mpi);
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        int n;

        if (!prj_mhd_local_block(mpi, block)) {
            continue;
        }
        for (n = 0; n < 56; ++n) {
            int nid = block->slot[n].id;
            prj_block *neighbor;

            if (nid < 0 || nid >= mesh->nblocks || block->slot[n].rel_level >= 0) {
                continue;
            }
            neighbor = &mesh->blocks[nid];
            if (!prj_mhd_local_block(mpi, neighbor) || neighbor->rank != block->rank) {
                continue;
            }
            prj_mhd_restrict_emf_to_coarse(block, neighbor);
        }
    }
}

void prj_mhd_debug_check_divb(const prj_mesh *mesh, const prj_mpi *mpi, int use_bf1)
{
    int bidx;

    if (mesh == 0) {
        prj_mhd_fail("prj_mhd_debug_check_divb: mesh is null");
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        const prj_block *block = &mesh->blocks[bidx];
        const double *bf[3];
        int i;
        int j;
        int k;
        int d;

        if (!prj_mhd_local_block(mpi, block)) {
            continue;
        }
        prj_mhd_check_bf_storage(block);
        for (d = 0; d < 3; ++d) {
            bf[d] = prj_block_bf_stage_const(block, d, use_bf1 != 0 ? 1 : 0);
        }
        for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
            for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    double divb =
                        (bf[X1DIR][FACE_IDX(X1DIR, i + 1, j, k)] - bf[X1DIR][FACE_IDX(X1DIR, i, j, k)]) / block->dx[0] +
                        (bf[X2DIR][FACE_IDX(X2DIR, i, j + 1, k)] - bf[X2DIR][FACE_IDX(X2DIR, i, j, k)]) / block->dx[1] +
                        (bf[X3DIR][FACE_IDX(X3DIR, i, j, k + 1)] - bf[X3DIR][FACE_IDX(X3DIR, i, j, k)]) / block->dx[2];

                    if (!isfinite(divb) || fabs(divb) > 1.0e-10) {
                        fprintf(stderr,
                            "prj_mhd_debug_check_divb: block=%d cell=(%d,%d,%d) divB=%.17e\n",
                            block->id, i, j, k, divb);
                        exit(EXIT_FAILURE);
                    }
                }
            }
        }
    }
}

void prj_mhd_debug_check_emf(const prj_mesh *mesh, const prj_mpi *mpi)
{
    int bidx;

    if (mesh == 0) {
        prj_mhd_fail("prj_mhd_debug_check_emf: mesh is null");
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        const prj_block *block = &mesh->blocks[bidx];
        int n;

        if (!prj_mhd_local_block(mpi, block)) {
            continue;
        }
        prj_mhd_check_emf_storage(block);
        for (n = 0; n < 56; ++n) {
            int nid = block->slot[n].id;
            const prj_block *neighbor;
            int dir;

            if (nid < 0 || nid >= mesh->nblocks) {
                continue;
            }
            neighbor = &mesh->blocks[nid];
            if (!prj_mhd_local_block(mpi, neighbor)) {
                continue;
            }
            prj_mhd_check_emf_storage(neighbor);
            if (block->slot[n].rel_level == 0) {
                for (dir = 0; dir < 3; ++dir) {
                    int i;
                    int j;
                    int k;

                    for (i = 0; i <= prj_mhd_edge_axis_active_max(dir, 0); ++i) {
                        for (j = 0; j <= prj_mhd_edge_axis_active_max(dir, 1); ++j) {
                            for (k = 0; k <= prj_mhd_edge_axis_active_max(dir, 2); ++k) {
                                double x[3];
                                int idx[3];
                                double q;
                                int d;
                                double a;
                                double b;

                                x[0] = prj_mhd_edge_coord(block, 0, dir, i);
                                x[1] = prj_mhd_edge_coord(block, 1, dir, j);
                                x[2] = prj_mhd_edge_coord(block, 2, dir, k);
                                if (!prj_mhd_edge_point_inside(neighbor, x)) {
                                    continue;
                                }
                                for (d = 0; d < 3; ++d) {
                                    double offset = d == dir ? 0.5 : 0.0;

                                    q = (x[d] - neighbor->xmin[d]) / neighbor->dx[d] - offset;
                                    idx[d] = prj_mhd_nearest_int(q);
                                    if (fabs(q - (double)idx[d]) > 1.0e-8) {
                                        prj_mhd_fail("prj_mhd_debug_check_emf: non-aligned same-level edge");
                                    }
                                }
                                if (!prj_mhd_edge_active(dir, idx[0], idx[1], idx[2])) {
                                    continue;
                                }
                                a = block->emf[dir][EDGE_IDX(dir, i, j, k)];
                                b = neighbor->emf[dir][EDGE_IDX(dir, idx[0], idx[1], idx[2])];
                                if (!isfinite(a) || !isfinite(b) || fabs(a - b) > 1.0e-10) {
                                    fprintf(stderr,
                                        "prj_mhd_debug_check_emf: same-level mismatch blocks %d/%d dir=%d value=(%.17e, %.17e) index=(%d, %d, %d), (%d, %d, %d) fidelity=(%d, %d)\n",
                                        block->id, neighbor->id, dir, a, b, i,j,k,idx[0],idx[1],idx[2],block->edge_fidelity[dir][EDGE_IDX(dir, i, j, k)], neighbor->edge_fidelity[dir][EDGE_IDX(dir, idx[0], idx[1], idx[2])]);
                                    exit(EXIT_FAILURE);
                                }
                            }
                        }
                    }
                }
            } else if (block->slot[n].rel_level < 0) {
                for (dir = 0; dir < 3; ++dir) {
                    int i;
                    int j;
                    int k;

                    for (i = 0; i <= prj_mhd_edge_axis_active_max(dir, 0); ++i) {
                        for (j = 0; j <= prj_mhd_edge_axis_active_max(dir, 1); ++j) {
                            for (k = 0; k <= prj_mhd_edge_axis_active_max(dir, 2); ++k) {
                                double x[3];
                                double restricted;
                                double coarse_value;

                                x[0] = prj_mhd_edge_coord(neighbor, 0, dir, i);
                                x[1] = prj_mhd_edge_coord(neighbor, 1, dir, j);
                                x[2] = prj_mhd_edge_coord(neighbor, 2, dir, k);
                                if (!prj_mhd_edge_point_inside(block, x)) {
                                    continue;
                                }
                                restricted = prj_mhd_restrict_emf_value(block, dir, x);
                                coarse_value = neighbor->emf[dir][EDGE_IDX(dir, i, j, k)];
                                if (!isfinite(coarse_value) ||
                                    fabs(restricted - coarse_value) > 1.0e-10) {
                                    fprintf(stderr,
                                        "prj_mhd_debug_check_emf: coarse/fine mismatch fine=%d coarse=%d dir=%d value=(%.17e, %.17e) index=(%d, %d, %d)\n",
                                        block->id, neighbor->id, dir, restricted, coarse_value, i, j, k);
                                    exit(EXIT_FAILURE);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void prj_mhd_init(prj_sim *sim, prj_mpi *mpi)
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

        if (!prj_mhd_initialized_storage_block(block)) {
            continue;
        }
        prj_mhd_check_block_storage(block);
        prj_mhd_store_vector_potential(block, sim->mhd_init_type, sim->mhd_B_norm, sim->mhd_B_scale);
    }

    prj_mhd_emf_send(&sim->mesh, mpi);
    prj_mpi_exchange_emf(&sim->mesh, mpi);
#if PRJ_MHD_DEBUG
    prj_mhd_debug_check_emf(&sim->mesh, mpi);
#endif

    for (bidx = 0; bidx < sim->mesh.nblocks; ++bidx) {
        prj_block *block = &sim->mesh.blocks[bidx];
        int d;

        if (!prj_mhd_initialized_storage_block(block)) {
            continue;
        }
        prj_mhd_check_block_storage(block);
        prj_mhd_curl_a_to_bf(block);
        prj_mhd_copy_bf_to_bf1(block);
        for (d = 0; d < 3; ++d) {
            prj_fill(block->emf[d], (size_t)PRJ_BLOCK_NEDGES, 0.0);
        }
        prj_mhd_bf2bc_all(&sim->eos, block, 0);
        prj_mhd_bf2bc_all(&sim->eos, block, 1);
    }
#if PRJ_MHD_DEBUG
    prj_mhd_debug_check_divb(&sim->mesh, mpi, 0);
#endif
}
#endif
