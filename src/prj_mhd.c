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

static inline int prj_mhd_initialized_storage_block(const prj_block *block)
{
    return block != 0 && block->id >= 0 && block->active == 1 &&
        block->W != 0 && block->W1 != 0 && block->U != 0 &&
        block->Bf[0] != 0 && block->Bf1[0] != 0 && block->emf[0] != 0;
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

enum {
    PRJ_MHD_FACE_FROM_COARSER = PRJ_MHD_FIDELITY_COARSER
};

static inline void prj_mhd_check_bf_storage(const prj_block *block)
{
    int d;

    if (block == 0) {
        prj_mhd_fail("prj_mhd_bf_prolongate: block is null");
    }
    for (d = 0; d < 3; ++d) {
        if (block->face_fidelity[d] == 0) {
            prj_mhd_fail("prj_mhd_bf_prolongate: missing face fidelity storage");
        }
        if (block->Bf[d] == 0 || block->Bf1[d] == 0) {
            prj_mhd_fail("prj_mhd_bf_prolongate: missing face-centered magnetic field storage");
        }
    }
}

static inline void prj_mhd_check_storage_index(const char *label, int i, int j, int k)
{
    if (i < -PRJ_NGHOST || i >= PRJ_BLOCK_SIZE + PRJ_NGHOST ||
        j < -PRJ_NGHOST || j >= PRJ_BLOCK_SIZE + PRJ_NGHOST ||
        k < -PRJ_NGHOST || k >= PRJ_BLOCK_SIZE + PRJ_NGHOST) {
        fprintf(stderr,
            "%s: index out of storage bounds i=%d j=%d k=%d (valid [%d, %d])\n",
            label, i, j, k, -PRJ_NGHOST, PRJ_BLOCK_SIZE + PRJ_NGHOST - 1);
        exit(EXIT_FAILURE);
    }
}

static inline double prj_mhd_sign_half(int bit)
{
    if (bit == 0) {
        return -1.0;
    }
    if (bit == 1) {
        return 1.0;
    }
    prj_mhd_fail("prj_mhd_bf_prolongate: expected fine half index 0 or 1");
    return 0.0;
}

static inline double prj_mhd_minmod_slope(double left, double center, double right, double dx)
{
    double sl;
    double sr;

    if (!isfinite(left) || !isfinite(center) || !isfinite(right) ||
        !isfinite(dx) || dx <= 0.0) {
        prj_mhd_fail("prj_mhd_bf_prolongate: invalid minmod stencil");
    }
    sl = (center - left) / dx;
    sr = (right - center) / dx;
    if (!isfinite(sl) || !isfinite(sr)) {
        prj_mhd_fail("prj_mhd_bf_prolongate: non-finite minmod slope");
    }
    if (sl * sr <= 0.0) {
        return 0.0;
    }
    if (fabs(sl) < fabs(sr)) {
        return sl;
    }
    return sr;
}

static inline double prj_mhd_read_bf_face(const double *bf,
    const char *label, int i, int j, int k)
{
    double value;

    if (bf == 0) {
        prj_mhd_fail("prj_mhd_bf_prolongate: null face-centered source");
    }
    prj_mhd_check_storage_index(label, i, j, k);
    value = bf[IDX(i, j, k)];
    if (!isfinite(value)) {
        fprintf(stderr, "%s: non-finite face-centered magnetic field at i=%d j=%d k=%d\n",
            label, i, j, k);
        exit(EXIT_FAILURE);
    }
    return value;
}

static inline double prj_mhd_interp_x1_face_flux(const prj_block *coarse,
    const double *bf, int i, int j, int k, int fine_j, int fine_k, double area)
{
    double base;
    double sy;
    double sz;
    double value;

    base = prj_mhd_read_bf_face(bf, "prj_mhd_bf_prolongate: x1 face", i, j, k);
    sy = prj_mhd_minmod_slope(
        prj_mhd_read_bf_face(bf, "prj_mhd_bf_prolongate: x1 y-stencil", i, j - 1, k),
        base,
        prj_mhd_read_bf_face(bf, "prj_mhd_bf_prolongate: x1 y-stencil", i, j + 1, k),
        coarse->dx[1]);
    sz = prj_mhd_minmod_slope(
        prj_mhd_read_bf_face(bf, "prj_mhd_bf_prolongate: x1 z-stencil", i, j, k - 1),
        base,
        prj_mhd_read_bf_face(bf, "prj_mhd_bf_prolongate: x1 z-stencil", i, j, k + 1),
        coarse->dx[2]);
    value = base +
        sy * (0.25 * prj_mhd_sign_half(fine_j) * coarse->dx[1]) +
        sz * (0.25 * prj_mhd_sign_half(fine_k) * coarse->dx[2]);
    if (!isfinite(value) || !isfinite(area) || area <= 0.0) {
        prj_mhd_fail("prj_mhd_bf_prolongate: invalid x1 interpolated flux");
    }
    return value * area;
}

static inline double prj_mhd_interp_x2_face_flux(const prj_block *coarse,
    const double *bf, int i, int j, int k, int fine_i, int fine_k, double area)
{
    double base;
    double sx;
    double sz;
    double value;

    base = prj_mhd_read_bf_face(bf, "prj_mhd_bf_prolongate: x2 face", i, j, k);
    sx = prj_mhd_minmod_slope(
        prj_mhd_read_bf_face(bf, "prj_mhd_bf_prolongate: x2 x-stencil", i - 1, j, k),
        base,
        prj_mhd_read_bf_face(bf, "prj_mhd_bf_prolongate: x2 x-stencil", i + 1, j, k),
        coarse->dx[0]);
    sz = prj_mhd_minmod_slope(
        prj_mhd_read_bf_face(bf, "prj_mhd_bf_prolongate: x2 z-stencil", i, j, k - 1),
        base,
        prj_mhd_read_bf_face(bf, "prj_mhd_bf_prolongate: x2 z-stencil", i, j, k + 1),
        coarse->dx[2]);
    value = base +
        sx * (0.25 * prj_mhd_sign_half(fine_i) * coarse->dx[0]) +
        sz * (0.25 * prj_mhd_sign_half(fine_k) * coarse->dx[2]);
    if (!isfinite(value) || !isfinite(area) || area <= 0.0) {
        prj_mhd_fail("prj_mhd_bf_prolongate: invalid x2 interpolated flux");
    }
    return value * area;
}

static inline double prj_mhd_interp_x3_face_flux(const prj_block *coarse,
    const double *bf, int i, int j, int k, int fine_i, int fine_j, double area)
{
    double base;
    double sx;
    double sy;
    double value;

    base = prj_mhd_read_bf_face(bf, "prj_mhd_bf_prolongate: x3 face", i, j, k);
    sx = prj_mhd_minmod_slope(
        prj_mhd_read_bf_face(bf, "prj_mhd_bf_prolongate: x3 x-stencil", i - 1, j, k),
        base,
        prj_mhd_read_bf_face(bf, "prj_mhd_bf_prolongate: x3 x-stencil", i + 1, j, k),
        coarse->dx[0]);
    sy = prj_mhd_minmod_slope(
        prj_mhd_read_bf_face(bf, "prj_mhd_bf_prolongate: x3 y-stencil", i, j - 1, k),
        base,
        prj_mhd_read_bf_face(bf, "prj_mhd_bf_prolongate: x3 y-stencil", i, j + 1, k),
        coarse->dx[1]);
    value = base +
        sx * (0.25 * prj_mhd_sign_half(fine_i) * coarse->dx[0]) +
        sy * (0.25 * prj_mhd_sign_half(fine_j) * coarse->dx[1]);
    if (!isfinite(value) || !isfinite(area) || area <= 0.0) {
        prj_mhd_fail("prj_mhd_bf_prolongate: invalid x3 interpolated flux");
    }
    return value * area;
}

static inline void prj_mhd_outer_or_prolongated_flux(prj_block *fine,
    double *dst[3], int dir, int i, int j, int k, double area,
    double prolongated_flux, double *flux)
{
    int idx;
    double value;

    if (fine == 0 || dst == 0 || flux == 0 || dir < 0 || dir >= 3 || dst[dir] == 0) {
        prj_mhd_fail("prj_mhd_bf_prolongate: invalid fine face destination");
    }
    if (!isfinite(area) || area <= 0.0 || !isfinite(prolongated_flux)) {
        prj_mhd_fail("prj_mhd_bf_prolongate: invalid outer face flux");
    }
    prj_mhd_check_storage_index("prj_mhd_bf_prolongate: fine outer face", i, j, k);
    idx = IDX(i, j, k);
    if (fine->face_fidelity[dir][idx] > PRJ_MHD_FACE_FROM_COARSER) {
        value = dst[dir][idx];
        if (!isfinite(value)) {
            prj_mhd_fail("prj_mhd_bf_prolongate: non-finite existing fine face");
        }
        *flux = value * area;
        return;
    }
    *flux = prolongated_flux;
    dst[dir][idx] = prolongated_flux / area;
    fine->face_fidelity[dir][idx] = PRJ_MHD_FACE_FROM_COARSER;
}

static inline void prj_mhd_write_inner_flux(prj_block *fine, double *dst[3],
    int dir, int i, int j, int k, double area, double flux)
{
    int idx;

    if (fine == 0 || dst == 0 || dir < 0 || dir >= 3 || dst[dir] == 0) {
        prj_mhd_fail("prj_mhd_bf_prolongate: invalid fine inner face destination");
    }
    if (!isfinite(area) || area <= 0.0 || !isfinite(flux)) {
        prj_mhd_fail("prj_mhd_bf_prolongate: invalid inner face flux");
    }
    prj_mhd_check_storage_index("prj_mhd_bf_prolongate: fine inner face", i, j, k);
    idx = IDX(i, j, k);
    if (fine->face_fidelity[dir][idx] > PRJ_MHD_FACE_FROM_COARSER) {
        if (!isfinite(dst[dir][idx])) {
            prj_mhd_fail("prj_mhd_bf_prolongate: non-finite existing fine inner face");
        }
        return;
    }
    dst[dir][idx] = flux / area;
    fine->face_fidelity[dir][idx] = PRJ_MHD_FACE_FROM_COARSER;
}

static inline void prj_mhd_check_prolongation_geometry(const prj_block *coarse,
    const prj_block *fine, int ci, int cj, int ck, int fi, int fj, int fk)
{
    const double tol = 1.0e-10;
    int cidx[3];
    int fidx[3];
    int d;

    cidx[0] = ci;
    cidx[1] = cj;
    cidx[2] = ck;
    fidx[0] = fi;
    fidx[1] = fj;
    fidx[2] = fk;

    for (d = 0; d < 3; ++d) {
        double scale;
        double dx_err;
        double coarse_lo;
        double fine_lo;
        double align_err;

        if (!isfinite(coarse->dx[d]) || !isfinite(fine->dx[d]) ||
            coarse->dx[d] <= 0.0 || fine->dx[d] <= 0.0) {
            prj_mhd_fail("prj_mhd_bf_prolongate: invalid block cell size");
        }
        scale = fabs(coarse->dx[d]) + fabs(fine->dx[d]) + 1.0;
        dx_err = fabs(coarse->dx[d] - 2.0 * fine->dx[d]);
        if (dx_err > tol * scale) {
            prj_mhd_fail("prj_mhd_bf_prolongate: coarse/fine cell sizes are not 2:1");
        }
        coarse_lo = coarse->xmin[d] + (double)cidx[d] * coarse->dx[d];
        fine_lo = fine->xmin[d] + (double)fidx[d] * fine->dx[d];
        align_err = fabs(coarse_lo - fine_lo);
        scale = fabs(coarse_lo) + fabs(fine_lo) + fabs(coarse->dx[d]) + 1.0;
        if (align_err > tol * scale) {
            prj_mhd_fail("prj_mhd_bf_prolongate: coarse and fine patches are not aligned");
        }
    }
}

static inline void prj_mhd_compute_inner_fluxes(double u[3][2][2],
    double v[2][3][2], double w[2][2][3], double dx1, double dx2, double dx3)
{
    double Uxx = 0.0;
    double Vyy = 0.0;
    double Wzz = 0.0;
    double Uxyz = 0.0;
    double Vxyz = 0.0;
    double Wxyz = 0.0;
    double dx1s;
    double dx2s;
    double dx3s;
    int i;
    int j;
    int k;

    if (!isfinite(dx1) || !isfinite(dx2) || !isfinite(dx3) ||
        dx1 <= 0.0 || dx2 <= 0.0 || dx3 <= 0.0) {
        prj_mhd_fail("prj_mhd_bf_prolongate: invalid fine cell sizes");
    }
    dx1s = dx1 * dx1;
    dx2s = dx2 * dx2;
    dx3s = dx3 * dx3;

    for (i = 0; i < 2; ++i) {
        double si = prj_mhd_sign_half(i);

        for (j = 0; j < 2; ++j) {
            double sj = prj_mhd_sign_half(j);

            for (k = 0; k < 2; ++k) {
                double sk = prj_mhd_sign_half(k);

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
        double sj = prj_mhd_sign_half(j);

        for (k = 0; k < 2; ++k) {
            double sk = prj_mhd_sign_half(k);

            u[1][j][k] = 0.5 * (u[2][j][k] + u[0][j][k]) +
                Uxx + sk * dx3s * Vxyz + sj * dx2s * Wxyz;
        }
    }
    for (i = 0; i < 2; ++i) {
        double si = prj_mhd_sign_half(i);

        for (k = 0; k < 2; ++k) {
            double sk = prj_mhd_sign_half(k);

            v[i][1][k] = 0.5 * (v[i][2][k] + v[i][0][k]) +
                Vyy + si * dx1s * Wxyz + sk * dx3s * Uxyz;
        }
    }
    for (i = 0; i < 2; ++i) {
        double si = prj_mhd_sign_half(i);

        for (j = 0; j < 2; ++j) {
            double sj = prj_mhd_sign_half(j);

            w[i][j][1] = 0.5 * (w[i][j][2] + w[i][j][0]) +
                Wzz + sj * dx2s * Uxyz + si * dx1s * Vxyz;
        }
    }
}

void prj_mhd_bf_prolongate(const prj_block *coarse, prj_block *fine,
    int ci, int cj, int ck, int fi, int fj, int fk, int use_bf1)
{
    double u[3][2][2];
    double v[2][3][2];
    double w[2][2][3];
    const double *src[3];
    double *dst[3];
    double area_u;
    double area_v;
    double area_w;
    int i;
    int j;
    int k;
    int d;

    prj_mhd_check_bf_storage(coarse);
    prj_mhd_check_bf_storage(fine);
    prj_mhd_check_prolongation_geometry(coarse, fine, ci, cj, ck, fi, fj, fk);

    for (d = 0; d < 3; ++d) {
        src[d] = use_bf1 != 0 ? coarse->Bf1[d] : coarse->Bf[d];
        dst[d] = use_bf1 != 0 ? fine->Bf1[d] : fine->Bf[d];
    }

    area_u = fine->dx[1] * fine->dx[2];
    area_v = fine->dx[0] * fine->dx[2];
    area_w = fine->dx[0] * fine->dx[1];
    if (!isfinite(area_u) || !isfinite(area_v) || !isfinite(area_w) ||
        area_u <= 0.0 || area_v <= 0.0 || area_w <= 0.0) {
        prj_mhd_fail("prj_mhd_bf_prolongate: invalid fine face area");
    }

    for (j = 0; j < 2; ++j) {
        for (k = 0; k < 2; ++k) {
            prj_mhd_outer_or_prolongated_flux(fine, dst, X1DIR, fi, fj + j, fk + k, area_u,
                prj_mhd_interp_x1_face_flux(coarse, src[X1DIR], ci, cj, ck, j, k, area_u),
                &u[0][j][k]);
            prj_mhd_outer_or_prolongated_flux(fine, dst, X1DIR, fi + 2, fj + j, fk + k, area_u,
                prj_mhd_interp_x1_face_flux(coarse, src[X1DIR], ci + 1, cj, ck, j, k, area_u),
                &u[2][j][k]);
        }
    }
    for (i = 0; i < 2; ++i) {
        for (k = 0; k < 2; ++k) {
            prj_mhd_outer_or_prolongated_flux(fine, dst, X2DIR, fi + i, fj, fk + k, area_v,
                prj_mhd_interp_x2_face_flux(coarse, src[X2DIR], ci, cj, ck, i, k, area_v),
                &v[i][0][k]);
            prj_mhd_outer_or_prolongated_flux(fine, dst, X2DIR, fi + i, fj + 2, fk + k, area_v,
                prj_mhd_interp_x2_face_flux(coarse, src[X2DIR], ci, cj + 1, ck, i, k, area_v),
                &v[i][2][k]);
        }
    }
    for (i = 0; i < 2; ++i) {
        for (j = 0; j < 2; ++j) {
            prj_mhd_outer_or_prolongated_flux(fine, dst, X3DIR, fi + i, fj + j, fk, area_w,
                prj_mhd_interp_x3_face_flux(coarse, src[X3DIR], ci, cj, ck, i, j, area_w),
                &w[i][j][0]);
            prj_mhd_outer_or_prolongated_flux(fine, dst, X3DIR, fi + i, fj + j, fk + 2, area_w,
                prj_mhd_interp_x3_face_flux(coarse, src[X3DIR], ci, cj, ck + 1, i, j, area_w),
                &w[i][j][2]);
        }
    }

    prj_mhd_compute_inner_fluxes(u, v, w, fine->dx[0], fine->dx[1], fine->dx[2]);

    for (j = 0; j < 2; ++j) {
        for (k = 0; k < 2; ++k) {
            prj_mhd_write_inner_flux(fine, dst, X1DIR, fi + 1, fj + j, fk + k, area_u,
                u[1][j][k]);
        }
    }
    for (i = 0; i < 2; ++i) {
        for (k = 0; k < 2; ++k) {
            prj_mhd_write_inner_flux(fine, dst, X2DIR, fi + i, fj + 1, fk + k, area_v,
                v[i][1][k]);
        }
    }
    for (i = 0; i < 2; ++i) {
        for (j = 0; j < 2; ++j) {
            prj_mhd_write_inner_flux(fine, dst, X3DIR, fi + i, fj + j, fk + 1, area_w,
                w[i][j][1]);
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
        prj_mhd_check_storage_index("prj_mhd_emf_send: fine edge", eidx[0], eidx[1], eidx[2]);
        value = fine->emf[dir][IDX(eidx[0], eidx[1], eidx[2])];
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
    prj_mhd_check_storage_index("prj_mhd_emf_send: destination edge", i, j, k);
    idx = IDX(i, j, k);
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
    for (dir = 0; dir < 3; ++dir) {
        int i;
        int j;
        int k;

        for (i = 0; i <= prj_mhd_edge_axis_active_max(dir, 0); ++i) {
            for (j = 0; j <= prj_mhd_edge_axis_active_max(dir, 1); ++j) {
                for (k = 0; k <= prj_mhd_edge_axis_active_max(dir, 2); ++k) {
                    double x[3];
                    double value;

                    x[0] = prj_mhd_edge_coord(coarse, 0, dir, i);
                    x[1] = prj_mhd_edge_coord(coarse, 1, dir, j);
                    x[2] = prj_mhd_edge_coord(coarse, 2, dir, k);
                    if (!prj_mhd_edge_point_inside(fine, x)) {
                        continue;
                    }
                    value = prj_mhd_restrict_emf_value(fine, dir, x);
                    prj_mhd_write_emf_edge(coarse, dir, i, j, k, value,
                        PRJ_MHD_FIDELITY_FINER);
                }
            }
        }
    }
}

static inline void prj_mhd_init_edge_fidelity(prj_mesh *mesh)
{
    int bidx;

    if (mesh == 0) {
        prj_mhd_fail("prj_mhd_emf_send: mesh is null");
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        int dir;

        if (!prj_mhd_local_block(block)) {
            continue;
        }
        prj_mhd_check_emf_storage(block);
        for (dir = 0; dir < 3; ++dir) {
            int i;
            int j;
            int k;
            int n;

            for (n = 0; n < PRJ_BLOCK_NCELLS; ++n) {
                block->edge_fidelity[dir][n] = PRJ_MHD_FIDELITY_NONE;
            }

            for (i = 0; i <= prj_mhd_edge_axis_active_max(dir, 0); ++i) {
                for (j = 0; j <= prj_mhd_edge_axis_active_max(dir, 1); ++j) {
                    for (k = 0; k <= prj_mhd_edge_axis_active_max(dir, 2); ++k) {
                        block->edge_fidelity[dir][IDX(i, j, k)] = PRJ_MHD_FIDELITY_SAME;
                    }
                }
            }
        }
    }
}

void prj_mhd_emf_send(prj_mesh *mesh)
{
    int bidx;

    prj_mhd_init_edge_fidelity(mesh);
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        int n;

        if (!prj_mhd_local_block(block)) {
            continue;
        }
        for (n = 0; n < 56; ++n) {
            int nid = block->slot[n].id;
            prj_block *neighbor;

            if (nid < 0 || nid >= mesh->nblocks || block->slot[n].rel_level >= 0) {
                continue;
            }
            neighbor = &mesh->blocks[nid];
            if (!prj_mhd_local_block(neighbor) || neighbor->rank != block->rank) {
                continue;
            }
            prj_mhd_restrict_emf_to_coarse(block, neighbor);
        }
    }
}

void prj_mhd_debug_check_divb(const prj_mesh *mesh, int use_bf1)
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

        if (!prj_mhd_local_block(block)) {
            continue;
        }
        prj_mhd_check_bf_storage(block);
        for (d = 0; d < 3; ++d) {
            bf[d] = use_bf1 != 0 ? block->Bf1[d] : block->Bf[d];
        }
        for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
            for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    double divb =
                        (bf[X1DIR][IDX(i + 1, j, k)] - bf[X1DIR][IDX(i, j, k)]) / block->dx[0] +
                        (bf[X2DIR][IDX(i, j + 1, k)] - bf[X2DIR][IDX(i, j, k)]) / block->dx[1] +
                        (bf[X3DIR][IDX(i, j, k + 1)] - bf[X3DIR][IDX(i, j, k)]) / block->dx[2];

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

void prj_mhd_debug_check_emf(const prj_mesh *mesh)
{
    int bidx;

    if (mesh == 0) {
        prj_mhd_fail("prj_mhd_debug_check_emf: mesh is null");
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        const prj_block *block = &mesh->blocks[bidx];
        int n;

        if (!prj_mhd_local_block(block)) {
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
            if (!prj_mhd_local_block(neighbor)) {
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
                                a = block->emf[dir][IDX(i, j, k)];
                                b = neighbor->emf[dir][IDX(idx[0], idx[1], idx[2])];
                                if (!isfinite(a) || !isfinite(b) || fabs(a - b) > 1.0e-10) {
                                    fprintf(stderr,
                                        "prj_mhd_debug_check_emf: same-level mismatch blocks %d/%d dir=%d value=(%.17e, %.17e)\n",
                                        block->id, neighbor->id, dir, a, b);
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
                                coarse_value = neighbor->emf[dir][IDX(i, j, k)];
                                if (!isfinite(coarse_value) ||
                                    fabs(restricted - coarse_value) > 1.0e-10) {
                                    fprintf(stderr,
                                        "prj_mhd_debug_check_emf: coarse/fine mismatch fine=%d coarse=%d dir=%d value=(%.17e, %.17e)\n",
                                        block->id, neighbor->id, dir, restricted, coarse_value);
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

        if (!prj_mhd_initialized_storage_block(block)) {
            continue;
        }
        prj_mhd_check_block_storage(block);
        prj_mhd_store_vector_potential(block, sim->mhd_init_type, sim->mhd_B_norm, sim->mhd_B_scale);
    }

    prj_mhd_emf_send(&sim->mesh);
    prj_mpi_exchange_emf(&sim->mesh, prj_mpi_current());
#if PRJ_MHD_DEBUG
    prj_mhd_debug_check_emf(&sim->mesh);
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
            prj_fill(block->emf[d], (size_t)PRJ_BLOCK_NCELLS, 0.0);
        }
        prj_mhd_bf2bc(&sim->eos, block, 0);
        prj_mhd_bf2bc(&sim->eos, block, 1);
    }
#if PRJ_MHD_DEBUG
    prj_mhd_debug_check_divb(&sim->mesh, 0);
#endif
}
#endif
