#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "prj.h"

#if PRJ_USE_RADIATION_FSA && PRJ_USE_RADIAL_FRAME_FSA
static void die(const char *msg)
{
    fprintf(stderr, "test_fsa_rotation_geometry: %s\n", msg);
    exit(1);
}

static double dot3(const double a[3], const double b[3])
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static void cross3(const double a[3], const double b[3], double c[3])
{
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

static void normalize3(double a[3])
{
    double mag = sqrt(dot3(a, a));

    if (mag <= 0.0) {
        die("zero-length vector");
    }
    a[0] /= mag;
    a[1] /= mag;
    a[2] /= mag;
}

static void assert_close(const char *name, double got, double expect, double tol)
{
    if (fabs(got - expect) > tol) {
        fprintf(stderr,
            "test_fsa_rotation_geometry: %s got %.17e expect %.17e tol %.3e\n",
            name, got, expect, tol);
        exit(1);
    }
}

static void read_rotation(const prj_block *block, int i, int j, int k, double R[3][3])
{
    int row;
    int col;

    for (row = 0; row < 3; ++row) {
        for (col = 0; col < 3; ++col) {
            R[row][col] = block->rotation_matrix_fsa[
                PRJ_FSA_ROT_IDX(row, col, i, j, k)];
        }
    }
}

static void arc_center_dir(const prj_rad *rad, const prj_block *block,
    int arc, int i, int j, int k, double n[3])
{
    int c0 = rad->arc_neighbor[2 * arc];
    int c1 = rad->arc_neighbor[2 * arc + 1];
    double n0[3];
    double n1[3];
    int d;

    prj_rad_fsa_rotated_angle_dir(rad, block, c0, i, j, k, n0);
    prj_rad_fsa_rotated_angle_dir(rad, block, c1, i, j, k, n1);
    for (d = 0; d < 3; ++d) {
        n[d] = n0[d] + n1[d];
    }
    normalize3(n);
}

static void check_rotation_matrix(const prj_block *block, int i, int j, int k)
{
    double R[3][3];
    double pos[3];
    double rhat[3];
    double col[3][3];
    double c0xc1[3];
    double r;
    int a;
    int b;
    int d;

    read_rotation(block, i, j, k, R);
    pos[0] = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
    pos[1] = block->xmin[1] + ((double)j + 0.5) * block->dx[1];
    pos[2] = block->xmin[2] + ((double)k + 0.5) * block->dx[2];
    r = sqrt(dot3(pos, pos));
    for (d = 0; d < 3; ++d) {
        rhat[d] = pos[d] / r;
    }

    for (d = 0; d < 3; ++d) {
        assert_close("R ez maps to rhat", R[d][2], rhat[d], 1.0e-13);
    }

    for (a = 0; a < 3; ++a) {
        for (d = 0; d < 3; ++d) {
            col[a][d] = R[d][a];
        }
    }
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            assert_close("rotation column orthonormality", dot3(col[a], col[b]),
                a == b ? 1.0 : 0.0, 1.0e-13);
        }
    }
    cross3(col[0], col[1], c0xc1);
    assert_close("rotation determinant", dot3(c0xc1, col[2]), 1.0, 1.0e-13);
}

static void check_ang_geom_matches_finite_difference(const prj_rad *rad,
    const prj_block *block, int i, int j, int k)
{
    const int arc = 0;
    double n_center[3];
    double expected[3] = {0.0, 0.0, 0.0};
    double got[3];
    int axis;
    int d;

    arc_center_dir(rad, block, arc, i, j, k, n_center);
    for (axis = 0; axis < 3; ++axis) {
        int ip[3] = {i, j, k};
        int im[3] = {i, j, k};
        double n_plus[3];
        double n_minus[3];

        ip[axis] += 1;
        im[axis] -= 1;
        arc_center_dir(rad, block, arc, ip[0], ip[1], ip[2], n_plus);
        arc_center_dir(rad, block, arc, im[0], im[1], im[2], n_minus);
        for (d = 0; d < 3; ++d) {
            double deriv = (n_plus[d] - n_minus[d]) / (2.0 * block->dx[axis]);

            expected[d] += n_center[axis] * deriv;
        }
    }
    for (d = 0; d < 3; ++d) {
        got[d] = block->ang_geom_fsa[PRJ_FSA_ANG_GEOM_IDX(arc, d, i, j, k)] /
            PRJ_CLIGHT;
        assert_close("ang_geom finite-difference check", got[d], expected[d], 2.0e-6);
    }
}

int main(void)
{
    prj_rad rad;
    prj_block block;
    const int ic = 2;
    const int jc = 2;
    const int kc = 2;
    const double dx = 1.0e-4;
    const double target[3] = {0.7, -0.4, 1.1};
    int d;

    memset(&rad, 0, sizeof(rad));
    memset(&block, 0, sizeof(block));
    prj_rad_fsa_calculate_directions(&rad);

    block.id = 0;
    block.active = 1;
    block.rank = 0;
    for (d = 0; d < 3; ++d) {
        block.dx[d] = dx;
        block.xmin[d] = target[d] - ((double)ic + 0.5) * dx;
        block.xmax[d] = block.xmin[d] + (double)PRJ_BLOCK_SIZE * dx;
    }
    if (prj_block_alloc_data(&block) != 0) {
        die("block allocation failed");
    }
    prj_rad_fsa_refresh_block_geometry(&rad, &block);

    check_rotation_matrix(&block, ic, jc, kc);
    check_ang_geom_matches_finite_difference(&rad, &block, ic, jc, kc);

    prj_block_free_data(&block);
    prj_rad_fsa_free_geometry(&rad);
    printf("test_fsa_rotation_geometry: ok\n");
    return 0;
}
#else
int main(void)
{
    printf("test_fsa_rotation_geometry: skipped (radial-frame FSA=0)\n");
    return 0;
}
#endif
