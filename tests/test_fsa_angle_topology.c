#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "prj.h"

#if PRJ_USE_RADIATION_FSA
#define TEST_MAX_CELL_ARCS 6

static void die(const char *msg)
{
    fprintf(stderr, "test_fsa_angle_topology: %s\n", msg);
    exit(1);
}

static double dot3(const double a[3], const double b[3])
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
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

static int arc_has_cell(const prj_rad *rad, int arc, int cell)
{
    return rad->arc_neighbor[2 * arc] == cell ||
        rad->arc_neighbor[2 * arc + 1] == cell;
}

static int cell_has_arc(const prj_rad *rad, int cell, int arc)
{
    int s;

    for (s = 0; s < TEST_MAX_CELL_ARCS; ++s) {
        if (rad->cell_neighbor[cell * TEST_MAX_CELL_ARCS + s] == arc) {
            return 1;
        }
    }
    return 0;
}

static int cell_arc_count(const prj_rad *rad, int cell)
{
    int count = 0;
    int s;

    for (s = 0; s < TEST_MAX_CELL_ARCS; ++s) {
        int arc = rad->cell_neighbor[cell * TEST_MAX_CELL_ARCS + s];

        if (arc >= 0) {
            count += 1;
        }
    }
    return count;
}

static void check_first_cell_points_north(const prj_rad *rad)
{
    if (fabs(rad->n0[0][0]) > 1.0e-14 ||
        fabs(rad->n0[0][1]) > 1.0e-14 ||
        fabs(rad->n0[0][2] - 1.0) > 1.0e-14) {
        die("first angular cell is not aligned with +z");
    }
}

static void check_cells(const prj_rad *rad, int cell_degree[PRJ_NANGLE])
{
    const double pi = acos(-1.0);
    double omega_sum = 0.0;
    double omega_pent = 0.0;
    double omega_hex = 0.0;
    int npent = 0;
    int nhex = 0;
    int cell;

    for (cell = 0; cell < PRJ_NANGLE; ++cell) {
        int s;
        int degree = cell_arc_count(rad, cell);

        cell_degree[cell] = degree;
        if (!isfinite(rad->solid_angle[cell]) || rad->solid_angle[cell] <= 0.0) {
            die("invalid cell solid angle");
        }
        omega_sum += rad->solid_angle[cell];
        if (degree == 5) {
            npent += 1;
            omega_pent += rad->solid_angle[cell];
        } else if (degree == 6) {
            nhex += 1;
            omega_hex += rad->solid_angle[cell];
        } else {
            die("cell degree is neither five nor six");
        }
        for (s = 0; s < TEST_MAX_CELL_ARCS; ++s) {
            int arc = rad->cell_neighbor[cell * TEST_MAX_CELL_ARCS + s];

            if (arc < 0) {
                continue;
            }
            if (arc >= PRJ_NARC) {
                die("cell arc neighbor index out of range");
            }
            if (!arc_has_cell(rad, arc, cell)) {
                die("cell arc neighbor does not include this cell");
            }
        }
    }
    if (npent != 12) {
        die("expected exactly twelve pentagons");
    }
    if (nhex != PRJ_NANGLE - 12) {
        die("unexpected hexagon count");
    }
    if (fabs(omega_sum - 4.0 * pi) > 1.0e-12) {
        die("solid angles do not sum to 4*pi");
    }
    if (nhex > 0 && fabs(omega_pent / (double)npent - omega_hex / (double)nhex) < 1.0e-12) {
        die("pentagon and hexagon average solid angles are not distinct");
    }
}

static void check_arcs(const prj_rad *rad, const int cell_degree[PRJ_NANGLE])
{
    const double pi = acos(-1.0);
    int pent_pent = 0;
    int pent_hex = 0;
    int hex_hex = 0;
    int arc;

    for (arc = 0; arc < PRJ_NARC; ++arc) {
        int c0 = rad->arc_neighbor[2 * arc];
        int c1 = rad->arc_neighbor[2 * arc + 1];
        double vec[3];
        double nface[3];
        double mid[3];
        double delta[3];
        double delta_proj[3];
        double delta_dot_mid;
        double vec_norm;
        double nface_norm;
        double nface_err;
        double tangent_err;
        double alignment;
        int d;

        if (c0 < 0 || c0 >= PRJ_NANGLE || c1 < 0 || c1 >= PRJ_NANGLE || c0 == c1) {
            die("arc cell neighbor index out of range");
        }
        if (!cell_has_arc(rad, c0, arc) || !cell_has_arc(rad, c1, arc)) {
            die("arc cell neighbor does not list this arc");
        }
        if (!isfinite(rad->arc_angle[arc]) ||
            rad->arc_angle[arc] <= 0.0 || rad->arc_angle[arc] >= pi) {
            die("invalid arc angle");
        }

        if (cell_degree[c0] == 5 && cell_degree[c1] == 5) {
            pent_pent += 1;
        } else if ((cell_degree[c0] == 5 && cell_degree[c1] == 6) ||
                   (cell_degree[c0] == 6 && cell_degree[c1] == 5)) {
            pent_hex += 1;
        } else if (cell_degree[c0] == 6 && cell_degree[c1] == 6) {
            hex_hex += 1;
        } else {
            die("unexpected arc type");
        }

        for (d = 0; d < 3; ++d) {
            vec[d] = rad->arc_vec[3 * arc + d];
            nface[d] = rad->arc_nface[3 * arc + d];
            mid[d] = rad->n0[c0][d] + rad->n0[c1][d];
            delta[d] = rad->n0[c1][d] - rad->n0[c0][d];
        }
        normalize3(mid);
        vec_norm = sqrt(dot3(vec, vec));
        if (!isfinite(vec_norm) || fabs(vec_norm - 1.0) > 1.0e-12) {
            die("arc_vec is not a unit vector");
        }
        nface_norm = sqrt(dot3(nface, nface));
        if (!isfinite(nface_norm) || fabs(nface_norm - 1.0) > 1.0e-12) {
            die("arc_nface is not a unit vector");
        }
        nface_err = 0.0;
        for (d = 0; d < 3; ++d) {
            double err = fabs(nface[d] - mid[d]);

            if (err > nface_err) {
                nface_err = err;
            }
        }
        if (nface_err > 1.0e-12) {
            die("arc_nface is not the normalized angular-cell average");
        }
        tangent_err = fabs(dot3(vec, mid));
        if (tangent_err > 1.0e-12) {
            die("arc_vec is not tangent to the sphere");
        }
        delta_dot_mid = dot3(delta, mid);
        for (d = 0; d < 3; ++d) {
            delta_proj[d] = delta[d] - delta_dot_mid * mid[d];
        }
        normalize3(delta_proj);
        alignment = dot3(vec, delta_proj);
        if (fabs(alignment - 1.0) > 1.0e-12) {
            die("arc_vec is not aligned from first cell to second cell");
        }
    }

#if PRJ_N_ANGLE_LEV == 1
    if (pent_pent != PRJ_NARC || pent_hex != 0 || hex_hex != 0) {
        die("unexpected level-1 arc type counts");
    }
#else
    if (pent_pent != 0 || pent_hex != 60 || hex_hex != PRJ_NARC - 60) {
        die("unexpected general-level arc type counts");
    }
#endif
}

int main(void)
{
    prj_rad rad;
    int cell_degree[PRJ_NANGLE];

    memset(&rad, 0, sizeof(rad));
    prj_rad_fsa_calculate_directions(&rad);
    if (rad.arc_angle == 0 || rad.arc_vec == 0 ||
        rad.arc_nface == 0 || rad.arc_neighbor == 0 || rad.cell_neighbor == 0) {
        die("FSA angular topology arrays were not allocated");
    }

    check_first_cell_points_north(&rad);
    check_cells(&rad, cell_degree);
    check_arcs(&rad, cell_degree);
    prj_rad_fsa_free_geometry(&rad);

    printf("test_fsa_angle_topology: ok (N_ANGLE_LEV=%d, cells=%d, arcs=%d)\n",
        PRJ_N_ANGLE_LEV, PRJ_NANGLE, PRJ_NARC);
    return 0;
}
#else
int main(void)
{
    printf("test_fsa_angle_topology: skipped (RADIATION_FSA=0)\n");
    return 0;
}
#endif
