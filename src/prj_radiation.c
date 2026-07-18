#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "prj.h"

static double prj_rad_m1_chi_exact(double f)
{
    if (f <= 0.0) {
        return 1.0 / 3.0;
    }
    if (f >= 1.0) {
        return 1.0;
    }
    return (3.0 + 4.0 * f * f) / (5.0 + 2.0 * sqrt(4.0 - 3.0 * f * f));
}

#if PRJ_NRAD > 0
/* Levermore/Vaytet third-moment scalar q(f), evaluated through the equivalent
 * boost parameter beta = 3f / (2 + sqrt(4 - 3f^2)). */
static double prj_rad_levermore_q_factor_exact(double f)
{
    double a;
    double beta;
    double beta2;
    double beta4;
    double one_minus_beta2;

    if (f <= 0.0) {
        return 0.0;
    }
    if (f >= 1.0) {
        return 1.0;
    }

    a = sqrt(fmax(0.0, 4.0 - 3.0 * f * f));
    beta = 3.0 * f / (2.0 + a);
    if (beta < 5.0e-2) {
        beta2 = beta * beta;
        return 4.0 * beta * (1.0 / 5.0 + beta2 * (1.0 / 21.0 - beta2 / 315.0));
    }

    beta2 = beta * beta;
    beta4 = beta2 * beta2;
    one_minus_beta2 = 1.0 - beta2;
    return (9.0 * beta - 8.0 / beta + 3.0 / (beta2 * beta) -
        3.0 * one_minus_beta2 * one_minus_beta2 * one_minus_beta2 *
        atanh(beta) / beta4) / (3.0 + beta2);
}

static void prj_rad_init_closure(prj_rad *rad)
{
    int i;

    if (rad == 0) {
        return;
    }
    for (i = 0; i <= NCLOSURE; ++i) {
        double f = (double)i / (double)NCLOSURE;

        rad->chi[i] = prj_rad_m1_chi_exact(f);
        rad->q[i] = prj_rad_levermore_q_factor_exact(f);
    }
}

static int prj_rad_closure_ready(const prj_rad *rad)
{
    return rad != 0 && rad->chi[0] > 0.0 && rad->chi[NCLOSURE] > 0.0 &&
        rad->q[NCLOSURE] > 0.0;
}

static double prj_rad_closure_lookup(const double values[NCLOSURE + 1], double f)
{
    double scaled;
    double w;
    int idx;

    if (f <= 0.0) {
        return values[0];
    }
    if (f >= 1.0) {
        return values[NCLOSURE];
    }
    scaled = f * (double)NCLOSURE;
    idx = (int)scaled;
    w = scaled - (double)idx;
    return values[idx] + w * (values[idx + 1] - values[idx]);
}
#endif

static double prj_rad_m1_chi(const prj_rad *rad, double f)
{
#if PRJ_NRAD > 0
    return prj_rad_closure_lookup(rad->chi, f);
#else
    (void)rad;
    return 0;
#endif
}

#if PRJ_NRAD > 0
static double prj_rad_levermore_q_factor(const prj_rad *rad, double f)
{
    return prj_rad_closure_lookup(rad->q, f);
}
#endif

#if PRJ_USE_RADIATION_FSA
#define PRJ_RAD_FSA_ICOS_NVERT 12
#define PRJ_RAD_FSA_ICOS_NFACE 20
#define PRJ_RAD_FSA_MAX_CELL_VERTS 6
#define PRJ_RAD_FSA_NTRI (20 * (PRJ_N_ANGLE_LEV) * (PRJ_N_ANGLE_LEV))

typedef struct prj_rad_fsa_triangle {
    int v[3];
    double center[3];
} prj_rad_fsa_triangle;

typedef struct prj_rad_fsa_cell {
    int tri[PRJ_RAD_FSA_MAX_CELL_VERTS];
    int ntri;
} prj_rad_fsa_cell;

typedef struct prj_rad_fsa_arc {
    int cell[2];
    int tri[2];
    int ntri;
} prj_rad_fsa_arc;

static void prj_rad_fsa_fail(const char *msg)
{
    fprintf(stderr, "prj_rad_fsa_calculate_directions: %s\n", msg);
    exit(EXIT_FAILURE);
}

static double prj_rad_fsa_dot(const double a[3], const double b[3])
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static void prj_rad_fsa_cross(const double a[3], const double b[3], double c[3])
{
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

static double prj_rad_fsa_clamp_dot(double x)
{
    if (x < -1.0) {
        return -1.0;
    }
    if (x > 1.0) {
        return 1.0;
    }
    return x;
}

static void prj_rad_fsa_normalize(double a[3])
{
    double n = sqrt(prj_rad_fsa_dot(a, a));

    if (n <= 0.0) {
        prj_rad_fsa_fail("zero-length vector");
    }
    a[0] /= n;
    a[1] /= n;
    a[2] /= n;
}

static void prj_rad_fsa_rotate_vector_to_north(const double from[3], double v[3])
{
    double c = prj_rad_fsa_clamp_dot(from[2]);
    double axis[3] = {from[1], -from[0], 0.0};
    double s2 = prj_rad_fsa_dot(axis, axis);
    double s;
    double cross1[3];
    double axis_dot_v;
    double out[3];
    int d;

    if (c > 1.0 - 1.0e-14) {
        return;
    }
    if (c < -1.0 + 1.0e-14) {
        v[1] = -v[1];
        v[2] = -v[2];
        return;
    }

    s = sqrt(s2);
    for (d = 0; d < 3; ++d) {
        axis[d] /= s;
    }
    prj_rad_fsa_cross(axis, v, cross1);
    axis_dot_v = prj_rad_fsa_dot(axis, v);
    for (d = 0; d < 3; ++d) {
        out[d] = c * v[d] + s * cross1[d] + (1.0 - c) * axis_dot_v * axis[d];
    }
    for (d = 0; d < 3; ++d) {
        v[d] = out[d];
    }
}

static void prj_rad_fsa_rotate_icosahedron_to_north(
    double x[PRJ_RAD_FSA_ICOS_NVERT][3])
{
    double from[3];
    int i;
    int d;

    for (d = 0; d < 3; ++d) {
        from[d] = x[0][d];
    }
    prj_rad_fsa_normalize(from);

    for (i = 0; i < PRJ_RAD_FSA_ICOS_NVERT; ++i) {
        prj_rad_fsa_rotate_vector_to_north(from, x[i]);
        prj_rad_fsa_normalize(x[i]);
    }
    x[0][0] = 0.0;
    x[0][1] = 0.0;
    x[0][2] = 1.0;
}

static void prj_rad_fsa_init_icosahedron(double x[PRJ_RAD_FSA_ICOS_NVERT][3])
{
    const double phi = 0.5 * (1.0 + sqrt(5.0));
    double raw[PRJ_RAD_FSA_ICOS_NVERT][3] = {
        {-1.0,  phi,  0.0}, { 1.0,  phi,  0.0},
        {-1.0, -phi,  0.0}, { 1.0, -phi,  0.0},
        { 0.0, -1.0,  phi}, { 0.0,  1.0,  phi},
        { 0.0, -1.0, -phi}, { 0.0,  1.0, -phi},
        { phi,  0.0, -1.0}, { phi,  0.0,  1.0},
        {-phi,  0.0, -1.0}, {-phi,  0.0,  1.0}
    };
    int i;
    int d;

    for (i = 0; i < PRJ_RAD_FSA_ICOS_NVERT; ++i) {
        for (d = 0; d < 3; ++d) {
            x[i][d] = raw[i][d];
        }
        prj_rad_fsa_normalize(x[i]);
    }
    prj_rad_fsa_rotate_icosahedron_to_north(x);
}

static void prj_rad_fsa_flat_face_point(
    const double ico[PRJ_RAD_FSA_ICOS_NVERT][3], const int face[3],
    int i, int j, double p[3])
{
    const int nlev = PRJ_N_ANGLE_LEV;
    double w0 = (double)(nlev - i - j);
    double w1 = (double)i;
    double w2 = (double)j;
    int d;

    for (d = 0; d < 3; ++d) {
        p[d] = (w0 * ico[face[0]][d] + w1 * ico[face[1]][d] +
            w2 * ico[face[2]][d]) / (double)nlev;
    }
}

static int prj_rad_fsa_find_or_add_vertex(double vertices[PRJ_NANGLE][3],
    int *nvertices, const double p_in[3])
{
    const double tol2 = 1.0e-24;
    double p[3];
    int i;
    int d;

    for (d = 0; d < 3; ++d) {
        p[d] = p_in[d];
    }
    prj_rad_fsa_normalize(p);

    for (i = 0; i < *nvertices; ++i) {
        double dx = p[0] - vertices[i][0];
        double dy = p[1] - vertices[i][1];
        double dz = p[2] - vertices[i][2];

        if (dx * dx + dy * dy + dz * dz < tol2) {
            return i;
        }
    }

    if (*nvertices >= PRJ_NANGLE) {
        prj_rad_fsa_fail("too many angular cell centers");
    }
    i = *nvertices;
    for (d = 0; d < 3; ++d) {
        vertices[i][d] = p[d];
    }
    *nvertices = i + 1;
    return i;
}

static void prj_rad_fsa_add_cell_triangle(prj_rad_fsa_cell *cells,
    int vertex_id, int triangle_id)
{
    prj_rad_fsa_cell *cell = &cells[vertex_id];

    if (cell->ntri >= PRJ_RAD_FSA_MAX_CELL_VERTS) {
        prj_rad_fsa_fail("angular cell has more than six vertices");
    }
    cell->tri[cell->ntri] = triangle_id;
    cell->ntri += 1;
}

static void prj_rad_fsa_add_triangle(prj_rad_fsa_triangle *triangles,
    prj_rad_fsa_cell *cells, int *ntriangles,
    int v0, int v1, int v2, const double center_in[3])
{
    prj_rad_fsa_triangle *tri;
    int tri_id;
    int d;

    if (*ntriangles >= PRJ_RAD_FSA_NTRI) {
        prj_rad_fsa_fail("too many subdivided triangles");
    }

    tri_id = *ntriangles;
    tri = &triangles[tri_id];
    tri->v[0] = v0;
    tri->v[1] = v1;
    tri->v[2] = v2;
    for (d = 0; d < 3; ++d) {
        tri->center[d] = center_in[d];
    }
    prj_rad_fsa_normalize(tri->center);

    prj_rad_fsa_add_cell_triangle(cells, v0, tri_id);
    prj_rad_fsa_add_cell_triangle(cells, v1, tri_id);
    prj_rad_fsa_add_cell_triangle(cells, v2, tri_id);
    *ntriangles = tri_id + 1;
}

static void prj_rad_fsa_triangle_center(
    const double ico[PRJ_RAD_FSA_ICOS_NVERT][3], const int face[3],
    int i0, int j0, int i1, int j1, int i2, int j2, double center[3])
{
    double p0[3];
    double p1[3];
    double p2[3];
    int d;

    prj_rad_fsa_flat_face_point(ico, face, i0, j0, p0);
    prj_rad_fsa_flat_face_point(ico, face, i1, j1, p1);
    prj_rad_fsa_flat_face_point(ico, face, i2, j2, p2);
    prj_rad_fsa_normalize(p0);
    prj_rad_fsa_normalize(p1);
    prj_rad_fsa_normalize(p2);
    {
        double e01[3];
        double e02[3];
        double sum[3];

        for (d = 0; d < 3; ++d) {
            e01[d] = p0[d] - p1[d];
            e02[d] = p0[d] - p2[d];
            sum[d] = p0[d] + p1[d] + p2[d];
        }
        prj_rad_fsa_cross(e01, e02, center);
        if (prj_rad_fsa_dot(center, sum) < 0.0) {
            center[0] = -center[0];
            center[1] = -center[1];
            center[2] = -center[2];
        }
        prj_rad_fsa_normalize(center);
    }
}

static void prj_rad_fsa_find_or_add_arc(prj_rad_fsa_arc *arcs,
    int *narcs, int cell0, int cell1, int tri_id)
{
    int a;

    if (cell0 == cell1) {
        prj_rad_fsa_fail("degenerate angular arc");
    }
    if (cell1 < cell0) {
        int tmp = cell0;

        cell0 = cell1;
        cell1 = tmp;
    }

    for (a = 0; a < *narcs; ++a) {
        if (arcs[a].cell[0] == cell0 && arcs[a].cell[1] == cell1) {
            if (arcs[a].ntri >= 2) {
                prj_rad_fsa_fail("angular arc has more than two triangle neighbors");
            }
            arcs[a].tri[arcs[a].ntri] = tri_id;
            arcs[a].ntri += 1;
            return;
        }
    }

    if (*narcs >= PRJ_NARC) {
        prj_rad_fsa_fail("too many angular arcs");
    }
    a = *narcs;
    arcs[a].cell[0] = cell0;
    arcs[a].cell[1] = cell1;
    arcs[a].tri[0] = tri_id;
    arcs[a].tri[1] = -1;
    arcs[a].ntri = 1;
    *narcs = a + 1;
}

static void prj_rad_fsa_cell_add_arc(int *cell_neighbor, int cell, int arc)
{
    int s;

    for (s = 0; s < PRJ_RAD_FSA_MAX_CELL_VERTS; ++s) {
        int *slot = &cell_neighbor[cell * PRJ_RAD_FSA_MAX_CELL_VERTS + s];

        if (*slot == arc) {
            prj_rad_fsa_fail("duplicate angular cell arc");
        }
        if (*slot == -1) {
            *slot = arc;
            return;
        }
    }
    prj_rad_fsa_fail("angular cell has too many arc neighbors");
}

static void prj_rad_fsa_build_arcs(const double vertices[PRJ_NANGLE][3],
    const prj_rad_fsa_triangle *triangles, const prj_rad_fsa_cell *cells,
    int ntriangles, prj_rad *rad)
{
    static const int edge_pair[3][2] = {{0, 1}, {1, 2}, {2, 0}};
    prj_rad_fsa_arc *arcs;
    int narcs = 0;
    int idx;
    int tri_id;
    int a;
    int cell;

    arcs = (prj_rad_fsa_arc *)prj_calloc((size_t)PRJ_NARC, sizeof(*arcs));
    for (idx = 0; idx < PRJ_NARC; ++idx) {
        rad->arc_angle[idx] = 0.0;
        rad->arc_vec[3 * idx] = 0.0;
        rad->arc_vec[3 * idx + 1] = 0.0;
        rad->arc_vec[3 * idx + 2] = 0.0;
        rad->arc_nface[3 * idx] = 0.0;
        rad->arc_nface[3 * idx + 1] = 0.0;
        rad->arc_nface[3 * idx + 2] = 0.0;
        rad->arc_neighbor[2 * idx] = -1;
        rad->arc_neighbor[2 * idx + 1] = -1;
    }
    for (idx = 0; idx < PRJ_RAD_FSA_MAX_CELL_VERTS * PRJ_NANGLE; ++idx) {
        rad->cell_neighbor[idx] = -1;
    }

    for (tri_id = 0; tri_id < ntriangles; ++tri_id) {
        int e;

        for (e = 0; e < 3; ++e) {
            int c0 = triangles[tri_id].v[edge_pair[e][0]];
            int c1 = triangles[tri_id].v[edge_pair[e][1]];

            prj_rad_fsa_find_or_add_arc(arcs, &narcs, c0, c1, tri_id);
        }
    }
    if (narcs != PRJ_NARC) {
        prj_rad_fsa_fail("unexpected number of angular arcs");
    }

    for (a = 0; a < narcs; ++a) {
        const double *p0;
        const double *p1;
        int c0 = arcs[a].cell[0];
        int c1 = arcs[a].cell[1];
        double mid[3];
        double dn[3];
        double v[3];
        double v_dot_mid;
        double p_dot;
        int d;

        if (arcs[a].ntri != 2) {
            prj_rad_fsa_fail("angular arc does not have two triangle neighbors");
        }

        p0 = triangles[arcs[a].tri[0]].center;
        p1 = triangles[arcs[a].tri[1]].center;
        p_dot = prj_rad_fsa_clamp_dot(prj_rad_fsa_dot(p0, p1));
        rad->arc_angle[a] = acos(p_dot);
        if (rad->arc_angle[a] <= 0.0) {
            prj_rad_fsa_fail("non-positive angular arc length");
        }

        for (d = 0; d < 3; ++d) {
            mid[d] = p0[d] + p1[d];
            dn[d] = vertices[c1][d] - vertices[c0][d];
        }
        prj_rad_fsa_normalize(mid);
        v_dot_mid = prj_rad_fsa_dot(dn, mid);
        for (d = 0; d < 3; ++d) {
            v[d] = dn[d] - v_dot_mid * mid[d];
        }
        prj_rad_fsa_normalize(v);
        if (prj_rad_fsa_dot(v, dn) < 0.0) {
            v[0] = -v[0];
            v[1] = -v[1];
            v[2] = -v[2];
        }

        rad->arc_neighbor[2 * a] = c0;
        rad->arc_neighbor[2 * a + 1] = c1;
        for (d = 0; d < 3; ++d) {
            rad->arc_vec[3 * a + d] = v[d];
        }
        {
            /* The angular face normal nface = normalize(n0[c0] + n0[c1]) is
             * purely geometric (n0 == vertices), so precompute the unit vector
             * once instead of rebuilding it per spatial cell in the kernel. */
            double nf[3];

            for (d = 0; d < 3; ++d) {
                nf[d] = vertices[c0][d] + vertices[c1][d];
            }
            prj_rad_fsa_normalize(nf);
            for (d = 0; d < 3; ++d) {
                rad->arc_nface[3 * a + d] = nf[d];
            }
        }
        prj_rad_fsa_cell_add_arc(rad->cell_neighbor, c0, a);
        prj_rad_fsa_cell_add_arc(rad->cell_neighbor, c1, a);
    }

    for (cell = 0; cell < PRJ_NANGLE; ++cell) {
        int count = 0;
        int s;

        for (s = 0; s < PRJ_RAD_FSA_MAX_CELL_VERTS; ++s) {
            if (rad->cell_neighbor[cell * PRJ_RAD_FSA_MAX_CELL_VERTS + s] >= 0) {
                count += 1;
            }
        }
        if (count != cells[cell].ntri) {
            prj_rad_fsa_fail("angular cell arc count does not match valence");
        }
    }

    free(arcs);
}

void prj_rad_fsa_free_geometry(prj_rad *rad)
{
    if (rad == 0) {
        return;
    }
    free(rad->arc_angle);
    free(rad->arc_vec);
    free(rad->arc_nface);
    free(rad->arc_neighbor);
    free(rad->cell_neighbor);
    rad->arc_angle = 0;
    rad->arc_vec = 0;
    rad->arc_nface = 0;
    rad->arc_neighbor = 0;
    rad->cell_neighbor = 0;
}

static void prj_rad_fsa_build_grid(double vertices[PRJ_NANGLE][3],
    prj_rad_fsa_triangle *triangles, prj_rad_fsa_cell *cells,
    int *nvertices, int *ntriangles)
{
    static const int faces[PRJ_RAD_FSA_ICOS_NFACE][3] = {
        {0, 11, 5}, {0, 5, 1}, {0, 1, 7}, {0, 7, 10}, {0, 10, 11},
        {1, 5, 9}, {5, 11, 4}, {11, 10, 2}, {10, 7, 6}, {7, 1, 8},
        {3, 9, 4}, {3, 4, 2}, {3, 2, 6}, {3, 6, 8}, {3, 8, 9},
        {4, 9, 5}, {2, 4, 11}, {6, 2, 10}, {8, 6, 7}, {9, 8, 1}
    };
    const int nlev = PRJ_N_ANGLE_LEV;
    const int stride = nlev + 1;
    double ico[PRJ_RAD_FSA_ICOS_NVERT][3];
    int *local;
    int f;

    prj_rad_fsa_init_icosahedron(ico);
    local = (int *)prj_malloc((size_t)stride * (size_t)stride * sizeof(*local));
    *nvertices = 0;
    *ntriangles = 0;

    for (f = 0; f < PRJ_RAD_FSA_ICOS_NFACE; ++f) {
        int i;
        int j;

        for (i = 0; i <= nlev; ++i) {
            for (j = 0; j <= nlev - i; ++j) {
                double p[3];

                prj_rad_fsa_flat_face_point(ico, faces[f], i, j, p);
                local[i * stride + j] =
                    prj_rad_fsa_find_or_add_vertex(vertices, nvertices, p);
            }
        }

        for (i = 0; i < nlev; ++i) {
            for (j = 0; j < nlev - i; ++j) {
                double center[3];
                int v0 = local[i * stride + j];
                int v1 = local[(i + 1) * stride + j];
                int v2 = local[i * stride + j + 1];

                prj_rad_fsa_triangle_center(ico, faces[f],
                    i, j, i + 1, j, i, j + 1, center);
                prj_rad_fsa_add_triangle(triangles, cells, ntriangles,
                    v0, v1, v2, center);

                if (i + j < nlev - 1) {
                    int v3 = local[(i + 1) * stride + j + 1];

                    prj_rad_fsa_triangle_center(ico, faces[f],
                        i + 1, j, i + 1, j + 1, i, j + 1, center);
                    prj_rad_fsa_add_triangle(triangles, cells, ntriangles,
                        v1, v3, v2, center);
                }
            }
        }
    }

    free(local);
}

static double prj_rad_fsa_spherical_triangle_area(const double a[3],
    const double b[3], const double c[3])
{
    double bx_c[3];
    double det;
    double denom;

    prj_rad_fsa_cross(b, c, bx_c);
    det = prj_rad_fsa_dot(a, bx_c);
    denom = 1.0 + prj_rad_fsa_dot(a, b) + prj_rad_fsa_dot(b, c) +
        prj_rad_fsa_dot(c, a);
    return 2.0 * atan2(fabs(det), denom);
}

static double prj_rad_fsa_cell_solid_angle(const double center[3],
    const prj_rad_fsa_cell *cell, const prj_rad_fsa_triangle *triangles)
{
    double ref[3] = {0.0, 0.0, 1.0};
    double e1[3];
    double e2[3];
    double angles[PRJ_RAD_FSA_MAX_CELL_VERTS];
    int order[PRJ_RAD_FSA_MAX_CELL_VERTS];
    double area = 0.0;
    int m = cell->ntri;
    int i;

    if (m < 3 || m > PRJ_RAD_FSA_MAX_CELL_VERTS) {
        prj_rad_fsa_fail("invalid angular cell valence");
    }
    if (fabs(center[2]) > 0.9) {
        ref[0] = 1.0;
        ref[1] = 0.0;
        ref[2] = 0.0;
    }
    prj_rad_fsa_cross(ref, center, e1);
    prj_rad_fsa_normalize(e1);
    prj_rad_fsa_cross(center, e1, e2);

    for (i = 0; i < m; ++i) {
        const double *p = triangles[cell->tri[i]].center;

        order[i] = cell->tri[i];
        angles[i] = atan2(prj_rad_fsa_dot(p, e2), prj_rad_fsa_dot(p, e1));
    }
    for (i = 1; i < m; ++i) {
        double angle = angles[i];
        int tri = order[i];
        int j = i - 1;

        while (j >= 0 && angles[j] > angle) {
            angles[j + 1] = angles[j];
            order[j + 1] = order[j];
            --j;
        }
        angles[j + 1] = angle;
        order[j + 1] = tri;
    }

    for (i = 0; i < m; ++i) {
        const double *p0 = triangles[order[i]].center;
        const double *p1 = triangles[order[(i + 1) % m]].center;

        area += prj_rad_fsa_spherical_triangle_area(center, p0, p1);
    }
    return area;
}

void prj_rad_fsa_calculate_directions(prj_rad *rad)
{
    double (*vertices)[3];
    prj_rad_fsa_triangle *triangles;
    prj_rad_fsa_cell *cells;
    int nvertices;
    int ntriangles;
    int npent = 0;
    int nhex = 0;
    int n;
    int d;

    if (rad == 0) {
        return;
    }

    prj_rad_fsa_free_geometry(rad);
    rad->arc_angle = (double *)prj_malloc((size_t)PRJ_NARC * sizeof(*rad->arc_angle));
    rad->arc_vec = (double *)prj_malloc((size_t)3 * (size_t)PRJ_NARC *
        sizeof(*rad->arc_vec));
    rad->arc_nface = (double *)prj_malloc((size_t)3 * (size_t)PRJ_NARC *
        sizeof(*rad->arc_nface));
    rad->arc_neighbor = (int *)prj_malloc((size_t)2 * (size_t)PRJ_NARC *
        sizeof(*rad->arc_neighbor));
    rad->cell_neighbor = (int *)prj_malloc((size_t)PRJ_RAD_FSA_MAX_CELL_VERTS *
        (size_t)PRJ_NANGLE * sizeof(*rad->cell_neighbor));

    vertices = (double (*)[3])prj_malloc((size_t)PRJ_NANGLE * sizeof(*vertices));
    triangles = (prj_rad_fsa_triangle *)prj_malloc(
        (size_t)PRJ_RAD_FSA_NTRI * sizeof(*triangles));
    cells = (prj_rad_fsa_cell *)prj_calloc((size_t)PRJ_NANGLE, sizeof(*cells));
    prj_rad_fsa_build_grid(vertices, triangles, cells, &nvertices, &ntriangles);

    if (nvertices != PRJ_NANGLE) {
        prj_rad_fsa_fail("unexpected number of angular cells");
    }
    if (ntriangles != PRJ_RAD_FSA_NTRI) {
        prj_rad_fsa_fail("unexpected number of triangulated faces");
    }

    for (n = 0; n < PRJ_NANGLE; ++n) {
        if (cells[n].ntri == 5) {
            npent += 1;
        } else if (cells[n].ntri == 6) {
            nhex += 1;
        } else {
            prj_rad_fsa_fail("angular cell is not a pentagon or hexagon");
        }
        for (d = 0; d < 3; ++d) {
            rad->n0[n][d] = vertices[n][d];
        }
        rad->solid_angle[n] =
            prj_rad_fsa_cell_solid_angle(vertices[n], &cells[n], triangles);
    }
    if (npent != 12 || nhex != PRJ_NANGLE - 12) {
        prj_rad_fsa_fail("unexpected pentagon/hexagon count");
    }
    prj_rad_fsa_build_arcs(vertices, triangles, cells, ntriangles, rad);

    free(cells);
    free(triangles);
    free(vertices);
}

#if PRJ_USE_RADIAL_FRAME_FSA
static void prj_rad_fsa_set_rotation_axis_fallback(double qz, double R[9],
    double omega[3][3])
{
    int i;
    int j;

    for (i = 0; i < 9; ++i) {
        R[i] = 0.0;
    }
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            omega[i][j] = 0.0;
        }
    }

    if (qz >= 0.0) {
        R[0] = 1.0;
        R[4] = 1.0;
        R[8] = 1.0;
    } else {
        R[0] = 1.0;
        R[4] = -1.0;
        R[8] = -1.0;
    }
}

static void prj_rad_fsa_rotation_omega_at(double x1, double x2, double x3,
    double R[9], double omega[3][3])
{
    const double eps_s2 = 1.0e-28;
    double r2 = x1 * x1 + x2 * x2 + x3 * x3;
    double r;
    double inv_r;
    double qx;
    double qy;
    double qz;
    double s2;
    double s;
    double inv_s;
    double inv_s2;
    int i;
    int j;

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            omega[i][j] = 0.0;
        }
    }

    if (r2 <= 0.0) {
        prj_rad_fsa_set_rotation_axis_fallback(1.0, R, omega);
        return;
    }

    r = sqrt(r2);
    inv_r = 1.0 / r;
    qx = x1 * inv_r;
    qy = x2 * inv_r;
    qz = x3 * inv_r;
    s2 = qx * qx + qy * qy;
    if (s2 <= eps_s2) {
        prj_rad_fsa_set_rotation_axis_fallback(qz, R, omega);
        return;
    }
    s = sqrt(s2);
    inv_s = 1.0 / s;
    inv_s2 = 1.0 / s2;

    R[0] = qx * qz * inv_s;
    R[1] = -qy * inv_s;
    R[2] = qx;
    R[3] = qy * qz * inv_s;
    R[4] = qx * inv_s;
    R[5] = qy;
    R[6] = -s;
    R[7] = 0.0;
    R[8] = qz;

    omega[0][0] = -qx * qy * qz * inv_s2 * inv_r;
    omega[0][1] = qx * qx * qz * inv_s2 * inv_r;
    omega[0][2] = -qy * inv_s2 * inv_r;
    omega[1][0] = -qy * qy * qz * inv_s2 * inv_r;
    omega[1][1] = qx * qy * qz * inv_s2 * inv_r;
    omega[1][2] = qx * inv_s2 * inv_r;
    omega[2][0] = qy * inv_r;
    omega[2][1] = -qx * inv_r;
    omega[2][2] = 0.0;
}

static void prj_rad_fsa_mat_vec(const double R[9], const double a[3], double out[3])
{
    out[0] = R[0] * a[0] + R[1] * a[1] + R[2] * a[2];
    out[1] = R[3] * a[0] + R[4] * a[1] + R[5] * a[2];
    out[2] = R[6] * a[0] + R[7] * a[1] + R[8] * a[2];
}
#endif

void prj_rad_fsa_rotated_dir(const prj_block *block, int i, int j, int k,
    const double n0[3], double n[3])
{
#if PRJ_USE_RADIAL_FRAME_FSA
    int row;
#endif

    if (n == 0) {
        return;
    }
    if (n0 == 0) {
        n[0] = 0.0;
        n[1] = 0.0;
        n[2] = 0.0;
        return;
    }
#if PRJ_USE_RADIAL_FRAME_FSA
    if (block == 0 || block->rotation_matrix_fsa == 0) {
        n[0] = n0[0];
        n[1] = n0[1];
        n[2] = n0[2];
        return;
    }
    for (row = 0; row < 3; ++row) {
        n[row] =
            block->rotation_matrix_fsa[PRJ_FSA_ROT_IDX(row, 0, i, j, k)] * n0[0] +
            block->rotation_matrix_fsa[PRJ_FSA_ROT_IDX(row, 1, i, j, k)] * n0[1] +
            block->rotation_matrix_fsa[PRJ_FSA_ROT_IDX(row, 2, i, j, k)] * n0[2];
    }
#else
    (void)block;
    (void)i;
    (void)j;
    (void)k;
    n[0] = n0[0];
    n[1] = n0[1];
    n[2] = n0[2];
#endif
}

void prj_rad_fsa_rotated_angle_dir(const prj_rad *rad, const prj_block *block,
    int angle, int i, int j, int k, double n[3])
{
    if (rad == 0 || angle < 0 || angle >= PRJ_NANGLE) {
        if (n != 0) {
            n[0] = 0.0;
            n[1] = 0.0;
            n[2] = 0.0;
        }
        return;
    }
    prj_rad_fsa_rotated_dir(block, i, j, k, rad->n0[angle], n);
}

#if PRJ_USE_RADIAL_FRAME_FSA
static void prj_rad_fsa_store_rotation(prj_block *block, int i, int j, int k,
    const double R[9])
{
    int row;
    int col;

    if (block == 0 || block->rotation_matrix_fsa == 0) {
        return;
    }
    for (row = 0; row < 3; ++row) {
        for (col = 0; col < 3; ++col) {
            block->rotation_matrix_fsa[PRJ_FSA_ROT_IDX(row, col, i, j, k)] =
                R[3 * row + col];
        }
    }
}

static void prj_rad_fsa_store_ang_geom(prj_block *block, int arc, int i, int j, int k,
    const double geom[3])
{
    int d;

    if (block == 0 || block->ang_geom_fsa == 0) {
        return;
    }
    for (d = 0; d < 3; ++d) {
        block->ang_geom_fsa[PRJ_FSA_ANG_GEOM_IDX(arc, d, i, j, k)] = geom[d];
    }
}

void prj_rad_fsa_refresh_block_geometry(const prj_rad *rad, prj_block *block)
{
    int i;
    int j;
    int k;

    if (block == 0 || block->rotation_matrix_fsa == 0) {
        return;
    }

    for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
        for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
            for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                double x1 = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
                double x2 = block->xmin[1] + ((double)j + 0.5) * block->dx[1];
                double x3 = block->xmin[2] + ((double)k + 0.5) * block->dx[2];
                double R[9];
                double omega[3][3];
                int arc;

                prj_rad_fsa_rotation_omega_at(x1, x2, x3, R, omega);
                prj_rad_fsa_store_rotation(block, i, j, k, R);

                if (block->ang_geom_fsa == 0) {
                    continue;
                }
                for (arc = 0; arc < PRJ_NARC; ++arc) {
                    double geom[3] = {0.0, 0.0, 0.0};

                    if (rad != 0 && rad->arc_neighbor != 0) {
                        int c0 = rad->arc_neighbor[2 * arc];
                        int c1 = rad->arc_neighbor[2 * arc + 1];

                        if (c0 >= 0 && c0 < PRJ_NANGLE && c1 >= 0 && c1 < PRJ_NANGLE) {
                            double n0[3];
                            double n1[3];
                            double n_arc[3];
                            double mag;
                            int d;
                            int axis;

                            prj_rad_fsa_mat_vec(R, rad->n0[c0], n0);
                            prj_rad_fsa_mat_vec(R, rad->n0[c1], n1);
                            for (d = 0; d < 3; ++d) {
                                n_arc[d] = n0[d] + n1[d];
                            }
                            mag = sqrt(prj_rad_fsa_dot(n_arc, n_arc));
                            if (mag > 0.0) {
                                for (d = 0; d < 3; ++d) {
                                    n_arc[d] /= mag;
                                }
                                for (axis = 0; axis < 3; ++axis) {
                                    double cross[3];

                                    prj_rad_fsa_cross(omega[axis], n_arc, cross);
                                    for (d = 0; d < 3; ++d) {
                                        geom[d] += n_arc[axis] * cross[d];
                                    }
                                }
                                for (d = 0; d < 3; ++d) {
                                    geom[d] *= PRJ_CLIGHT;
                                }
                            }
                        }
                    }
                    prj_rad_fsa_store_ang_geom(block, arc, i, j, k, geom);
                }
            }
        }
    }
}

void prj_rad_fsa_refresh_mesh_geometry(const prj_rad *rad, prj_mesh *mesh,
    const prj_mpi *mpi)
{
    int bidx;

    if (mesh == 0) {
        return;
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (block->id < 0 || block->active != 1 || block->W_mhd == 0) {
            continue;
        }
        if (mpi != 0 && block->rank != mpi->rank) {
            continue;
        }
        prj_rad_fsa_refresh_block_geometry(rad, block);
    }
}
#endif
#endif

void prj_rad_init(prj_rad *rad)
{
#if PRJ_USE_RADIATION_FSA
    prj_rad_fsa_calculate_directions(rad);
#endif
#if PRJ_USE_RADIATION_M1
    prj_rad_init_closure(rad);
#endif
#if PRJ_USE_RADIATION_M1 || PRJ_USE_RADIATION_FSA
    prj_rad3_opac_init(rad);
#endif
#if PRJ_USE_RADIATION_M1 || PRJ_USE_RADIATION_FSA
    prj_rad_eleinel_init(rad);
#elif PRJ_NRAD == 0
    (void)rad;
#endif
}

void prj_rad_prim2cons(const double *W, double *U)
{
#if PRJ_USE_RADIATION_FSA
    int field;
    int group;
    int angle;

    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            for (angle = 0; angle < PRJ_NANGLE; ++angle) {
                U[PRJ_CONS_RAD_I(field, group, angle)] =
                    W[PRJ_PRIM_RAD_I(field, group, angle)];
            }
        }
    }
#elif PRJ_USE_RADIATION_M1
    int field;
    int group;

    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            U[PRJ_CONS_RAD_E(field, group)] = W[PRJ_PRIM_RAD_E(field, group)];
            U[PRJ_CONS_RAD_F1(field, group)] = W[PRJ_PRIM_RAD_F1(field, group)];
            U[PRJ_CONS_RAD_F2(field, group)] = W[PRJ_PRIM_RAD_F2(field, group)];
            U[PRJ_CONS_RAD_F3(field, group)] = W[PRJ_PRIM_RAD_F3(field, group)];
        }
    }
#else
    (void)W;
    (void)U;
#endif
}

void prj_rad_cons2prim(const double *U, double *W)
{
#if PRJ_USE_RADIATION_FSA
    int field;
    int group;
    int angle;

    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            for (angle = 0; angle < PRJ_NANGLE; ++angle) {
                W[PRJ_PRIM_RAD_I(field, group, angle)] =
                    U[PRJ_CONS_RAD_I(field, group, angle)];
            }
        }
    }
#elif PRJ_USE_RADIATION_M1
    int field;
    int group;

    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            W[PRJ_PRIM_RAD_E(field, group)] = U[PRJ_CONS_RAD_E(field, group)];
            W[PRJ_PRIM_RAD_F1(field, group)] = U[PRJ_CONS_RAD_F1(field, group)];
            W[PRJ_PRIM_RAD_F2(field, group)] = U[PRJ_CONS_RAD_F2(field, group)];
            W[PRJ_PRIM_RAD_F3(field, group)] = U[PRJ_CONS_RAD_F3(field, group)];
        }
    }
#else
    (void)U;
    (void)W;
#endif
}

/* Public M1 closure for the pressure tensor.  P^{ij} = E * D^{ij} with the
 * Levermore Eddington tensor D^{ij} = a δ^{ij} + b n^i n^j, n = F/|F|, and
 * χ(f) = (3 + 4f²)/(5 + 2√(4 - 3f²)), f = |F|/(c E).  Falls back to the
 * isotropic limit P^{ij} = (E/3) δ^{ij} when |F| or E vanishes. */
void prj_rad_m1_pressure(const prj_rad *rad, double E, double F1, double F2, double F3,
    double P[3][3])
{
    double Fmag;
    double cE;
    double f;
    double chi;
    double a_c;
    double b_c;
    double n[3];
    int a;
    int b;

    Fmag = sqrt(F1 * F1 + F2 * F2 + F3 * F3);
    cE = PRJ_CLIGHT * (E > 0.0 ? E : 0.0);

    if (cE <= 0.0 || Fmag <= 0.0) {
        double third = (E > 0.0 ? E : 0.0) / 3.0;

        for (a = 0; a < 3; ++a) {
            for (b = 0; b < 3; ++b) {
                P[a][b] = (a == b) ? third : 0.0;
            }
        }
        return;
    }

    f = Fmag / cE;
    if (f > 1.0) {
        f = 1.0;
    }
    chi = prj_rad_m1_chi(rad, f);
    a_c = 0.5 * (1.0 - chi);
    b_c = 0.5 * (3.0 * chi - 1.0);
    n[0] = F1 / Fmag;
    n[1] = F2 / Fmag;
    n[2] = F3 / Fmag;
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            P[a][b] = E * (a_c * (a == b ? 1.0 : 0.0) + b_c * n[a] * n[b]);
        }
    }
}

#if PRJ_DYNAMIC_GR && PRJ_USE_RADIATION_M1
static void prj_rad_gr_m1_raise_vec_ctx(const prj_rad_gr_m1_closure_ctx *ctx,
    const double vcov[3], double vcon[3])
{
    int a;
    int b;

    for (a = 0; a < 3; ++a) {
        vcon[a] = 0.0;
        for (b = 0; b < 3; ++b) {
            vcon[a] += ctx->gamma_inv[a][b] * vcov[b];
        }
    }
}

static void prj_rad_gr_m1_lower_vec_ctx(const prj_rad_gr_m1_closure_ctx *ctx,
    const double vcon[3], double vcov[3])
{
    int a;
    int b;

    for (a = 0; a < 3; ++a) {
        vcov[a] = 0.0;
        for (b = 0; b < 3; ++b) {
            vcov[a] += ctx->gamma[a][b] * vcon[b];
        }
    }
}

static double prj_rad_gr_m1_dot_con_con(const prj_rad_gr_m1_closure_ctx *ctx,
    const double v[3], const double w[3])
{
    double dot = 0.0;
    int a;
    int b;

    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            dot += ctx->gamma[a][b] * v[a] * w[b];
        }
    }
    return dot;
}

static void prj_rad_gr_m1_christoffel(const prj_rad_gr_m1_closure_ctx *ctx,
    double Gamma[3][3][3])
{
    int a;
    int b;
    int c;
    int d;

    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            for (c = 0; c < 3; ++c) {
                Gamma[a][b][c] = 0.0;
                for (d = 0; d < 3; ++d) {
                    Gamma[a][b][c] += 0.5 * ctx->gamma_inv[a][d] *
                        (ctx->dgamma[b][c][d] + ctx->dgamma[c][b][d] -
                            ctx->dgamma[d][b][c]);
                }
            }
        }
    }
}

#if PRJ_INCLUDE_RADIATION_VISCOSITY
/* Shear tensor + mean-free-path helpers: only used by the viscous closure path
 * (the -(4 lbar/15) sigma term). Not compiled when radiation viscosity is off. */
static void prj_rad_gr_m1_shear(const prj_rad_gr_m1_closure_ctx *ctx,
    const double vcon[3], const double vcov[3], double wlor,
    double sigma_con[3][3], double *divu_out, double *sigma2_out)
{
    double Gamma[3][3][3];
    double du_cov[3][3];
    double D_cov[3][3];
    double h_cov[3][3];
    double h_con[3][3];
    double bracket[3][3];
    double u_cov[3];
    double dbeta2[3];
    double dw[3];
    double divu = 0.0;
    double sigma_cov[3][3];
    double sigma2 = 0.0;
    int a;
    int b;
    int c;
    int d;

    memset(sigma_con, 0, 9 * sizeof(double));
    *divu_out = 0.0;
    *sigma2_out = 0.0;
    if (ctx == 0 || !ctx->have_shear) {
        return;
    }

    prj_rad_gr_m1_christoffel(ctx, Gamma);
    for (a = 0; a < 3; ++a) {
        u_cov[a] = wlor * vcov[a];
    }
    for (d = 0; d < 3; ++d) {
        dbeta2[d] = 0.0;
        for (a = 0; a < 3; ++a) {
            for (b = 0; b < 3; ++b) {
                dbeta2[d] += ctx->dgamma[d][a][b] * vcon[a] * vcon[b] +
                    2.0 * ctx->gamma[a][b] * vcon[a] * ctx->dvdx[d][b];
            }
        }
        dw[d] = 0.5 * wlor * wlor * wlor * dbeta2[d];
    }

    for (d = 0; d < 3; ++d) {
        for (a = 0; a < 3; ++a) {
            double dvcov = 0.0;

            for (b = 0; b < 3; ++b) {
                dvcov += ctx->dgamma[d][a][b] * vcon[b] +
                    ctx->gamma[a][b] * ctx->dvdx[d][b];
            }
            du_cov[d][a] = dw[d] * vcov[a] + wlor * dvcov;
        }
    }

    for (d = 0; d < 3; ++d) {
        for (a = 0; a < 3; ++a) {
            D_cov[d][a] = du_cov[d][a] - wlor * ctx->K_dd[d][a];
            for (b = 0; b < 3; ++b) {
                D_cov[d][a] -= Gamma[b][d][a] * u_cov[b];
            }
            divu += ctx->gamma_inv[d][a] * D_cov[d][a];
        }
    }

    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            h_cov[a][b] = ctx->gamma[a][b] + u_cov[a] * u_cov[b];
            h_con[a][b] = ctx->gamma_inv[a][b] +
                wlor * wlor * vcon[a] * vcon[b];
        }
    }

    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            bracket[a][b] = D_cov[a][b] + D_cov[b][a] -
                (2.0 / 3.0) * h_cov[a][b] * divu;
        }
    }
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            for (c = 0; c < 3; ++c) {
                for (d = 0; d < 3; ++d) {
                    sigma_con[a][b] += h_con[a][c] * h_con[b][d] *
                        bracket[c][d];
                }
            }
        }
    }

    memset(sigma_cov, 0, sizeof(sigma_cov));
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            for (c = 0; c < 3; ++c) {
                for (d = 0; d < 3; ++d) {
                    sigma_cov[a][b] += ctx->gamma[a][c] * ctx->gamma[b][d] *
                        sigma_con[c][d];
                }
            }
            sigma2 += sigma_con[a][b] * sigma_cov[a][b];
        }
    }
    if (!isfinite(divu)) {
        divu = 0.0;
    }
    if (!isfinite(sigma2) || sigma2 < 0.0) {
        sigma2 = 0.0;
    }
    *divu_out = divu;
    *sigma2_out = sigma2;
}

static double prj_rad_gr_m1_lbar_by_shear(double divu, double sigma2)
{
    if (sigma2 > 0.0 && divu > 0.0) {
        return PRJ_RAD_GR_M1_CSIGMA * sqrt(divu / sigma2);
    }
    return 1.0e300;
}

static double prj_rad_gr_m1_lbar_from_side(const prj_rad_gr_m1_closure_ctx *ctx,
    const prj_rad_gr_m1_side_data *side)
{
    double by_opacity = 1.0e300;
    double by_shear = side != 0 ? side->lbar_by_shear : 1.0e300;
    double opacity = ctx != 0 ? ctx->opacity : 0.0;

    if (isfinite(opacity) && opacity > 0.0) {
        by_opacity = 1.0 / opacity;
    }
    if (by_opacity < by_shear) {
        return by_opacity;
    }
    if (by_shear < 1.0e299) {
        return by_shear;
    }
    return by_opacity < 1.0e299 ? by_opacity : 0.0;
}
#endif /* PRJ_INCLUDE_RADIATION_VISCOSITY */

static const int prj_rad_gr_m1_sym_i[6] = {0, 1, 2, 0, 0, 1};
static const int prj_rad_gr_m1_sym_j[6] = {0, 1, 2, 1, 2, 2};

typedef struct prj_rad_gr_m1_pressure_data {
    const prj_rad *rad;
    const prj_rad_gr_m1_closure_ctx *ctx;
    double E;
    double Fcov[3];
    double Fcon[3];
    double Fhat_con[3];
    double Pthin[3][3];
    double A0[3][3];
    double vcon[3];
    double vcov[3];
    double u_cov[3];
    double h_mix[3][3];
    double H0[3];
    double Hcoef[3][6];
    double Jcoef[6];
    double wlor;
    double J0;
    double beta2;
    double f_euler;
    /* fbar-independent pieces of the 6x6 closure system, precomputed once per
     * solve so each root-finding evaluation (Newton or the Illinois fallback)
     * only rescales by the scalars thin_w/thick_w instead of rebuilding
     * A0/H0/vcon products every iteration. Evaluation order matches the inline
     * expressions it replaces so the result is bit-identical (no -ffast-math
     * reassociation). */
    double base_coef[6][6];
    double Pthin_sym[6];
    double Prad_sym[6];
    /* Max row-sum of |base_coef|. The per-fbar system is (I - thick_w*base_coef)
     * P = rhs; base_coef carries the fluid-velocity coupling and scales as
     * O(v^2/c^2), so for non-relativistic flows this matrix is I minus a tiny
     * perturbation. When |thick_w|*base_row_absmax is safely < 1 the inverse is
     * evaluated as a truncated Neumann series (a few branch-free 6x6 mat-vecs)
     * instead of a full Gaussian elimination; decided once per group here so the
     * per-iteration solve carries no norm/branch overhead. */
    double base_row_absmax;
} prj_rad_gr_m1_pressure_data;

static int prj_rad_gr_m1_solve6(double A[6][7], double x[6])
{
    int col;
    int row;

    /* Gaussian elimination with partial pivoting + back substitution. Only the
     * rows below each pivot are eliminated (forward sweep), then the upper-
     * triangular system is back-solved -- ~40% fewer flops than a Gauss-Jordan
     * sweep that also clears rows above every pivot. Solves the same system to
     * machine precision (result differs by rounding only). */
    for (col = 0; col < 6; ++col) {
        int pivot = col;
        double pivabs = fabs(A[col][col]);
        double inv_piv;

        for (row = col + 1; row < 6; ++row) {
            double v = fabs(A[row][col]);

            if (v > pivabs) {
                pivabs = v;
                pivot = row;
            }
        }
        if (!isfinite(pivabs) || pivabs < 1.0e-14) {
            return 0;
        }
        if (pivot != col) {
            int j;

            for (j = col; j < 7; ++j) {
                double tmp = A[col][j];

                A[col][j] = A[pivot][j];
                A[pivot][j] = tmp;
            }
        }
        inv_piv = 1.0 / A[col][col];
        for (row = col + 1; row < 6; ++row) {
            double factor = A[row][col] * inv_piv;
            int j;

            if (factor == 0.0) {
                continue;
            }
            for (j = col + 1; j < 7; ++j) {
                A[row][j] -= factor * A[col][j];
            }
        }
    }
    for (row = 5; row >= 0; --row) {
        double sum = A[row][6];
        int j;

        for (j = row + 1; j < 6; ++j) {
            sum -= A[row][j] * x[j];
        }
        x[row] = sum / A[row][row];
        if (!isfinite(x[row])) {
            return 0;
        }
    }
    return 1;
}

#if PRJ_INCLUDE_RADIATION_VISCOSITY
/* Per-side kinematics + shear for the GR M1 closure. Depends only on the fluid
 * velocity and the face geometry (NOT on E/F or opacity), so it is identical
 * for every energy group of a face side. Split out of prepare_pressure so the
 * interface flux can build it once per side instead of once per group; the
 * arithmetic (including operand order for the A0/h_mix/Jcoef/Hcoef terms) is
 * unchanged, so prepare_pressure below reproduces the old result bit-for-bit. */
void prj_rad_gr_m1_prepare_side(const prj_rad_gr_m1_closure_ctx *ctx,
    prj_rad_gr_m1_side_data *side)
{
    double beta2;
    int a;
    int b;
    int m;

    for (a = 0; a < 3; ++a) {
        side->vcon[a] = isfinite(ctx->vcon[a]) ? ctx->vcon[a] : 0.0;
    }
    beta2 = prj_rad_gr_m1_dot_con_con(ctx, side->vcon, side->vcon);
    if (!isfinite(beta2) || beta2 < 0.0) {
        beta2 = 0.0;
        side->vcon[0] = side->vcon[1] = side->vcon[2] = 0.0;
    }
    if (beta2 >= 1.0) {
        double scale = sqrt((1.0 - 1.0e-12) / beta2);

        for (a = 0; a < 3; ++a) {
            side->vcon[a] *= scale;
        }
        beta2 = 1.0 - 1.0e-12;
    }
    side->beta2 = beta2;
    prj_rad_gr_m1_lower_vec_ctx(ctx, side->vcon, side->vcov);
    side->wlor = 1.0 / sqrt(1.0 - beta2);
    for (a = 0; a < 3; ++a) {
        side->u_cov[a] = side->wlor * side->vcov[a];
    }

    prj_rad_gr_m1_shear(ctx, side->vcon, side->vcov, side->wlor,
        side->sigma_con, &side->divu, &side->sigma2);
    side->lbar_by_shear = prj_rad_gr_m1_lbar_by_shear(side->divu,
        side->sigma2);
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            side->A0_kin[a][b] =
                (ctx->gamma_inv[a][b] + 4.0 * side->vcon[a] * side->vcon[b]) /
                    3.0;
            side->h_mix[a][b] = (a == b ? 1.0 : 0.0) +
                side->wlor * side->wlor * side->vcon[a] * side->vcov[b];
        }
    }

    for (m = 0; m < 6; ++m) {
        int p = prj_rad_gr_m1_sym_i[m];
        int q = prj_rad_gr_m1_sym_j[m];

        side->Jcoef[m] = side->u_cov[p] * side->u_cov[q] *
            (p == q ? 1.0 : 2.0);
        for (a = 0; a < 3; ++a) {
            side->Hcoef[a][m] = -side->h_mix[a][p] * side->u_cov[q];
            if (p != q) {
                side->Hcoef[a][m] -= side->h_mix[a][q] * side->u_cov[p];
            }
        }
    }

    /* Velocity/metric-only parts of the closure coefficient matrix (see the
     * struct comment). base_coef_kin is the A0_kin + Hcoef*vcon contribution;
     * sigma_jcoef isolates the shear part that prepare_pressure scales by lbar.
     * Row sums (of |.|) bound ||base_coef||_inf for the Neumann guard. */
    for (m = 0; m < 6; ++m) {
        int i = prj_rad_gr_m1_sym_i[m];
        int j = prj_rad_gr_m1_sym_j[m];
        double kin_row = 0.0;
        double sig_row = 0.0;
        int n;

        for (n = 0; n < 6; ++n) {
            double kin = side->A0_kin[i][j] * side->Jcoef[n] +
                side->Hcoef[i][n] * side->vcon[j] +
                side->Hcoef[j][n] * side->vcon[i];
            double sig = side->sigma_con[i][j] * side->Jcoef[n];

            side->base_coef_kin[m][n] = kin;
            side->sigma_jcoef[m][n] = sig;
            kin_row += fabs(kin);
            sig_row += fabs(sig);
        }
        side->base_kin_rowsum[m] = kin_row;
        side->sigma_jcoef_rowsum[m] = sig_row;
    }
}
#else /* !PRJ_INCLUDE_RADIATION_VISCOSITY */
/* Inviscid per-side kinematics for the GR M1 closure: identical to the viscous
 * prepare_side above but with the radiation shear-viscosity machinery removed
 * (no shear tensor, no mean-free-path length, no sigma_jcoef). With the
 * -(4 lbar/15) sigma term gone, base_coef == base_coef_kin, so only that piece
 * is built here and its row sum is the exact ||base_coef||_inf. */
void prj_rad_gr_m1_prepare_side(const prj_rad_gr_m1_closure_ctx *ctx,
    prj_rad_gr_m1_side_data *side)
{
    double beta2;
    int a;
    int b;
    int m;

    for (a = 0; a < 3; ++a) {
        side->vcon[a] = isfinite(ctx->vcon[a]) ? ctx->vcon[a] : 0.0;
    }
    beta2 = prj_rad_gr_m1_dot_con_con(ctx, side->vcon, side->vcon);
    if (!isfinite(beta2) || beta2 < 0.0) {
        beta2 = 0.0;
        side->vcon[0] = side->vcon[1] = side->vcon[2] = 0.0;
    }
    if (beta2 >= 1.0) {
        double scale = sqrt((1.0 - 1.0e-12) / beta2);

        for (a = 0; a < 3; ++a) {
            side->vcon[a] *= scale;
        }
        beta2 = 1.0 - 1.0e-12;
    }
    side->beta2 = beta2;
    prj_rad_gr_m1_lower_vec_ctx(ctx, side->vcon, side->vcov);
    side->wlor = 1.0 / sqrt(1.0 - beta2);
    for (a = 0; a < 3; ++a) {
        side->u_cov[a] = side->wlor * side->vcov[a];
    }

    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            side->A0_kin[a][b] =
                (ctx->gamma_inv[a][b] + 4.0 * side->vcon[a] * side->vcon[b]) /
                    3.0;
            side->h_mix[a][b] = (a == b ? 1.0 : 0.0) +
                side->wlor * side->wlor * side->vcon[a] * side->vcov[b];
        }
    }

    for (m = 0; m < 6; ++m) {
        int p = prj_rad_gr_m1_sym_i[m];
        int q = prj_rad_gr_m1_sym_j[m];

        side->Jcoef[m] = side->u_cov[p] * side->u_cov[q] *
            (p == q ? 1.0 : 2.0);
        for (a = 0; a < 3; ++a) {
            side->Hcoef[a][m] = -side->h_mix[a][p] * side->u_cov[q];
            if (p != q) {
                side->Hcoef[a][m] -= side->h_mix[a][q] * side->u_cov[p];
            }
        }
    }

    for (m = 0; m < 6; ++m) {
        int i = prj_rad_gr_m1_sym_i[m];
        int j = prj_rad_gr_m1_sym_j[m];
        double kin_row = 0.0;
        int n;

        for (n = 0; n < 6; ++n) {
            double kin = side->A0_kin[i][j] * side->Jcoef[n] +
                side->Hcoef[i][n] * side->vcon[j] +
                side->Hcoef[j][n] * side->vcon[i];

            side->base_coef_kin[m][n] = kin;
            kin_row += fabs(kin);
        }
        side->base_kin_rowsum[m] = kin_row;
    }
}
#endif /* PRJ_INCLUDE_RADIATION_VISCOSITY */

#if PRJ_INCLUDE_RADIATION_VISCOSITY
static void prj_rad_gr_m1_prepare_pressure(const prj_rad *rad,
    const prj_rad_gr_m1_closure_ctx *ctx,
    const prj_rad_gr_m1_side_data *side, double E, const double Fcov_in[3],
    prj_rad_gr_m1_pressure_data *data)
{
    double F2;
    double Fmag;
    double cE;
    double FdotV;
    double Fdotu;
    double Q0;
    double lbar;
    int a;
    int b;
    int m;

    memset(data, 0, sizeof(*data));
    data->rad = rad;
    data->ctx = ctx;
    if (!isfinite(E) || E < 0.0) {
        E = 0.0;
    }
    data->E = E;
    for (a = 0; a < 3; ++a) {
        data->Fcov[a] = isfinite(Fcov_in[a]) ? Fcov_in[a] : 0.0;
    }
    prj_rad_gr_m1_raise_vec_ctx(ctx, data->Fcov, data->Fcon);
    F2 = 0.0;
    for (a = 0; a < 3; ++a) {
        F2 += data->Fcov[a] * data->Fcon[a];
    }
    if (!isfinite(F2) || F2 < 0.0) {
        F2 = 0.0;
    }
    Fmag = sqrt(F2);
    cE = PRJ_CLIGHT * E;
    if (Fmag > cE && Fmag > 0.0) {
        double scale = cE / Fmag;

        for (a = 0; a < 3; ++a) {
            data->Fcov[a] *= scale;
            data->Fcon[a] *= scale;
        }
        F2 *= scale * scale;
        Fmag = cE;
    }
    data->f_euler = cE > 0.0 ? Fmag / cE : 0.0;
    if (data->f_euler > 1.0) {
        data->f_euler = 1.0;
    }

    memset(data->Pthin, 0, sizeof(data->Pthin));
    if (E > 0.0 && F2 > 0.0) {
        for (a = 0; a < 3; ++a) {
            for (b = 0; b < 3; ++b) {
                data->Pthin[a][b] = E * data->Fcon[a] * data->Fcon[b] / F2;
            }
        }
    }

    /* Per-side kinematics were computed once in prj_rad_gr_m1_prepare_side. */
    for (a = 0; a < 3; ++a) {
        data->vcon[a] = side->vcon[a];
        data->vcov[a] = side->vcov[a];
        data->u_cov[a] = side->u_cov[a];
        data->Fhat_con[a] = data->Fcon[a] / PRJ_CLIGHT;
    }
    data->wlor = side->wlor;
    data->beta2 = side->beta2;
    memcpy(data->h_mix, side->h_mix, sizeof(data->h_mix));
    memcpy(data->Jcoef, side->Jcoef, sizeof(data->Jcoef));
    memcpy(data->Hcoef, side->Hcoef, sizeof(data->Hcoef));

    lbar = prj_rad_gr_m1_lbar_from_side(ctx, side);
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            data->A0[a][b] = side->A0_kin[a][b] -
                (4.0 * lbar / 15.0) * side->sigma_con[a][b];
        }
    }

    FdotV = 0.0;
    for (a = 0; a < 3; ++a) {
        FdotV += data->Fhat_con[a] * data->vcov[a];
    }
    Fdotu = data->wlor * FdotV;
    data->J0 = E * data->wlor * data->wlor - 2.0 * data->wlor * Fdotu;
    Q0 = E * data->wlor - Fdotu;
    for (a = 0; a < 3; ++a) {
        double h_n = -data->wlor * data->wlor * data->vcon[a];
        double h_F = data->Fhat_con[a] +
            data->wlor * data->wlor * data->vcon[a] * FdotV;

        data->H0[a] = Q0 * h_n + data->wlor * h_F;
    }

    /* Precompute the fbar-independent RHS terms (symmetric-storage Pthin/Prad).
     * The coefficient matrix itself is assembled just below from the per-side
     * pieces. */
    for (m = 0; m < 6; ++m) {
        int i = prj_rad_gr_m1_sym_i[m];
        int j = prj_rad_gr_m1_sym_j[m];

        data->Pthin_sym[m] = data->Pthin[i][j];
        data->Prad_sym[m] = data->A0[i][j] * data->J0 +
            data->H0[i] * data->vcon[j] +
            data->H0[j] * data->vcon[i];
    }

    /* base_coef = base_coef_kin - (4 lbar/15) * sigma_jcoef. Both matrices are
     * velocity/metric-only (E/F-independent) and were built once per side in
     * prepare_side; here only the scalar lbar (opacity-dependent, hence per
     * group) recombines them, replacing the full 6x6 triple-product build. The
     * distributed multiply rounds differently than the old fused form, so the
     * result matches to rounding (~1e-16), not bit-for-bit. */
    {
        double clbar = (4.0 * lbar / 15.0);

        data->base_row_absmax = 0.0;
        for (m = 0; m < 6; ++m) {
            int n;

            for (n = 0; n < 6; ++n) {
                data->base_coef[m][n] = side->base_coef_kin[m][n] -
                    clbar * side->sigma_jcoef[m][n];
            }
            /* Upper bound on the row sum of |base_coef| via the triangle
             * inequality; conservative for the Neumann guard (may fall back to
             * the full solve slightly more often, never less safely). */
            {
                double row = side->base_kin_rowsum[m] +
                    fabs(clbar) * side->sigma_jcoef_rowsum[m];

                if (row > data->base_row_absmax) {
                    data->base_row_absmax = row;
                }
            }
        }
    }
}
#else /* !PRJ_INCLUDE_RADIATION_VISCOSITY */
/* Inviscid closure assembly + E/F -> fluid-frame J/H conversion: identical to
 * the viscous prepare_pressure above but with the radiation-viscosity term
 * dropped. A0 == A0_kin (no lbar, no shear), so Prad and the coefficient matrix
 * carry no mean-free-path dependence, and base_coef == base_coef_kin. */
static void prj_rad_gr_m1_prepare_pressure(const prj_rad *rad,
    const prj_rad_gr_m1_closure_ctx *ctx,
    const prj_rad_gr_m1_side_data *side, double E, const double Fcov_in[3],
    prj_rad_gr_m1_pressure_data *data)
{
    double F2;
    double Fmag;
    double cE;
    double FdotV;
    double Fdotu;
    double Q0;
    int a;
    int b;
    int m;

    memset(data, 0, sizeof(*data));
    data->rad = rad;
    data->ctx = ctx;
    if (!isfinite(E) || E < 0.0) {
        E = 0.0;
    }
    data->E = E;
    for (a = 0; a < 3; ++a) {
        data->Fcov[a] = isfinite(Fcov_in[a]) ? Fcov_in[a] : 0.0;
    }
    prj_rad_gr_m1_raise_vec_ctx(ctx, data->Fcov, data->Fcon);
    F2 = 0.0;
    for (a = 0; a < 3; ++a) {
        F2 += data->Fcov[a] * data->Fcon[a];
    }
    if (!isfinite(F2) || F2 < 0.0) {
        F2 = 0.0;
    }
    Fmag = sqrt(F2);
    cE = PRJ_CLIGHT * E;
    if (Fmag > cE && Fmag > 0.0) {
        double scale = cE / Fmag;

        for (a = 0; a < 3; ++a) {
            data->Fcov[a] *= scale;
            data->Fcon[a] *= scale;
        }
        F2 *= scale * scale;
        Fmag = cE;
    }
    data->f_euler = cE > 0.0 ? Fmag / cE : 0.0;
    if (data->f_euler > 1.0) {
        data->f_euler = 1.0;
    }

    memset(data->Pthin, 0, sizeof(data->Pthin));
    if (E > 0.0 && F2 > 0.0) {
        for (a = 0; a < 3; ++a) {
            for (b = 0; b < 3; ++b) {
                data->Pthin[a][b] = E * data->Fcon[a] * data->Fcon[b] / F2;
            }
        }
    }

    /* Per-side kinematics were computed once in prj_rad_gr_m1_prepare_side. */
    for (a = 0; a < 3; ++a) {
        data->vcon[a] = side->vcon[a];
        data->vcov[a] = side->vcov[a];
        data->u_cov[a] = side->u_cov[a];
        data->Fhat_con[a] = data->Fcon[a] / PRJ_CLIGHT;
    }
    data->wlor = side->wlor;
    data->beta2 = side->beta2;
    memcpy(data->h_mix, side->h_mix, sizeof(data->h_mix));
    memcpy(data->Jcoef, side->Jcoef, sizeof(data->Jcoef));
    memcpy(data->Hcoef, side->Hcoef, sizeof(data->Hcoef));

    /* No radiation viscosity: A0 is the pure kinematic Eddington tensor. */
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            data->A0[a][b] = side->A0_kin[a][b];
        }
    }

    FdotV = 0.0;
    for (a = 0; a < 3; ++a) {
        FdotV += data->Fhat_con[a] * data->vcov[a];
    }
    Fdotu = data->wlor * FdotV;
    data->J0 = E * data->wlor * data->wlor - 2.0 * data->wlor * Fdotu;
    Q0 = E * data->wlor - Fdotu;
    for (a = 0; a < 3; ++a) {
        double h_n = -data->wlor * data->wlor * data->vcon[a];
        double h_F = data->Fhat_con[a] +
            data->wlor * data->wlor * data->vcon[a] * FdotV;

        data->H0[a] = Q0 * h_n + data->wlor * h_F;
    }

    for (m = 0; m < 6; ++m) {
        int i = prj_rad_gr_m1_sym_i[m];
        int j = prj_rad_gr_m1_sym_j[m];

        data->Pthin_sym[m] = data->Pthin[i][j];
        data->Prad_sym[m] = data->A0[i][j] * data->J0 +
            data->H0[i] * data->vcon[j] +
            data->H0[j] * data->vcon[i];
    }

    /* base_coef == base_coef_kin (no viscosity); base_kin_rowsum is then the
     * exact row-sum bound on ||base_coef||_inf for the Neumann guard. */
    data->base_row_absmax = 0.0;
    for (m = 0; m < 6; ++m) {
        int n;

        for (n = 0; n < 6; ++n) {
            data->base_coef[m][n] = side->base_coef_kin[m][n];
        }
        if (side->base_kin_rowsum[m] > data->base_row_absmax) {
            data->base_row_absmax = side->base_kin_rowsum[m];
        }
    }
}
#endif /* PRJ_INCLUDE_RADIATION_VISCOSITY */

/* Solve (I - thick_w*base_coef) sol = rhs by a truncated Neumann series when the
 * coupling is weak (caller guarantees ||thick_w*base_coef||_inf < 1). N carries
 * the O(v^2/c^2) velocity coupling, so sol = (I + N + N^2 + ...) rhs converges
 * geometrically; we accumulate terms until the next one is negligible relative
 * to the running sum. Replaces a 6x6 Gaussian elimination (with pivoting and 6
 * divisions) by a few branch-free mat-vecs -- the dominant win for the closure
 * in non-relativistic flows. */
static void prj_rad_gr_m1_neumann_solve(const double base_coef[6][6],
    double thick_w, const double rhs[6], double sol[6])
{
    double term[6];
    int k;
    int m;
    int n;

    for (m = 0; m < 6; ++m) {
        sol[m] = rhs[m];
        term[m] = rhs[m];
    }
    for (k = 0; k < 8; ++k) {
        double next[6];
        double tmax = 0.0;
        double smax = 0.0;

        for (m = 0; m < 6; ++m) {
            double s = 0.0;

            for (n = 0; n < 6; ++n) {
                s += thick_w * base_coef[m][n] * term[n];
            }
            next[m] = s;
        }
        for (m = 0; m < 6; ++m) {
            double a;

            sol[m] += next[m];
            term[m] = next[m];
            a = fabs(next[m]);
            if (a > tmax) {
                tmax = a;
            }
            a = fabs(sol[m]);
            if (a > smax) {
                smax = a;
            }
        }
        if (tmax <= 1.0e-15 * smax) {
            break;
        }
    }
}

static void prj_rad_gr_m1_pressure_for_fbar(
    const prj_rad_gr_m1_pressure_data *data, double fbar, double P[3][3])
{
    double amat[6][7];
    double rhs[6];
    double sol[6];
    double chi;
    double thin_w;
    double thick_w;
    int m;

    if (!isfinite(fbar) || fbar < 0.0) {
        fbar = 0.0;
    } else if (fbar > 1.0) {
        fbar = 1.0;
    }
    chi = prj_rad_m1_chi(data->rad, fbar);
    thin_w = 0.5 * (3.0 * chi - 1.0);
    thick_w = 1.5 * (1.0 - chi);
    for (m = 0; m < 6; ++m) {
        rhs[m] = thin_w * data->Pthin_sym[m] + thick_w * data->Prad_sym[m];
    }

    if (fabs(thick_w) * data->base_row_absmax < PRJ_RAD_GR_M1_NEUMANN_MAXNORM) {
        prj_rad_gr_m1_neumann_solve(data->base_coef, thick_w, rhs, sol);
    } else {
        /* Relativistic coupling: fall back to the full solve. Every column is
         * written below, so no memset of amat is needed. */
        for (m = 0; m < 6; ++m) {
            int n;

            for (n = 0; n < 6; ++n) {
                amat[m][n] =
                    (m == n ? 1.0 : 0.0) - thick_w * data->base_coef[m][n];
            }
            amat[m][6] = rhs[m];
        }
        if (!prj_rad_gr_m1_solve6(amat, sol)) {
            for (m = 0; m < 6; ++m) {
                sol[m] = rhs[m];
            }
        }
    }

    memset(P, 0, 9 * sizeof(double));
    for (m = 0; m < 6; ++m) {
        int i = prj_rad_gr_m1_sym_i[m];
        int j = prj_rad_gr_m1_sym_j[m];

        P[i][j] = sol[m];
        P[j][i] = sol[m];
    }
}

static double prj_rad_gr_m1_derived_fbar(
    const prj_rad_gr_m1_pressure_data *data, const double P[3][3])
{
    const prj_rad_gr_m1_closure_ctx *ctx = data->ctx;
    double R0;
    double Rcon[3];
    double J;
    double numerator;
    double denom;
    int a;
    int b;

    /* Eq. 3.29 is evaluated in the local Eulerian-normal frame:
     * M^{00}=E, M^{0i}=F^i/c, M^{ij}=P^{ij}.  Eq. 6.24 then uses
     * h_{alpha gamma} R^alpha R^gamma / J^2 with
     * R^alpha=M^{alpha beta}u_beta and J=M^{alpha beta}u_alpha u_beta. */
    R0 = -data->E * data->wlor;
    J = data->J0;
    for (a = 0; a < 3; ++a) {
        R0 += data->Fhat_con[a] * data->u_cov[a];
        Rcon[a] = -data->wlor * data->Fhat_con[a];
        for (b = 0; b < 3; ++b) {
            Rcon[a] += P[a][b] * data->u_cov[b];
            J += P[a][b] * data->u_cov[a] * data->u_cov[b];
        }
    }

    numerator = (data->wlor * data->wlor - 1.0) * R0 * R0;
    for (a = 0; a < 3; ++a) {
        numerator += -2.0 * data->wlor * data->u_cov[a] * R0 * Rcon[a];
        for (b = 0; b < 3; ++b) {
            numerator += (ctx->gamma[a][b] +
                data->u_cov[a] * data->u_cov[b]) * Rcon[a] * Rcon[b];
        }
    }
    if (!isfinite(numerator) || numerator < 0.0) {
        numerator = 0.0;
    }
    denom = fabs(J);
    if (!isfinite(denom) || denom <= 0.0) {
        return 0.0;
    }
    numerator = sqrt(numerator) / denom;
    if (!isfinite(numerator) || numerator < 0.0) {
        return 0.0;
    }
    if (numerator > 1.0) {
        return 1.0;
    }
    return numerator;
}

/* Newton-Raphson attempt on g(fbar) = derived_fbar(P(fbar)) - fbar. There is no
 * analytic Jacobian through the chi/q closure interpolation, so g' is estimated
 * by a one-sided finite difference (stepped inward at the fbar=1 edge), costing
 * two residual evaluations per iteration. Seeded from the Eulerian closure it
 * typically converges in ~3-4 iterations. Returns 1 and fills P / *fbar_out on
 * success; returns 0 (leaving P untouched) on any failure -- non-finite iterate,
 * a step leaving [0, 1], a degenerate derivative, stagnation short of tolerance,
 * or exhausting the iteration cap -- so the caller can fall back to the
 * guaranteed-convergent Illinois bracketing. */
static int prj_rad_gr_m1_pressure_newton(
    const prj_rad_gr_m1_pressure_data *data, double f0, double P[3][3],
    double *fbar_out)
{
    const double tol = PRJ_RAD_GR_M1_FBAR_TOL;
    double f = f0;
    int iter;

    if (!isfinite(f) || f < 0.0) {
        f = 0.0;
    } else if (f > 1.0) {
        f = 1.0;
    }

    for (iter = 0; iter < 15; ++iter) {
        double Pf[3][3];
        double Ph[3][3];
        double fderiv;
        double g;
        double h;
        double fh;
        double gh;
        double deriv;
        double fnew;

        prj_rad_gr_m1_pressure_for_fbar(data, f, Pf);
        fderiv = prj_rad_gr_m1_derived_fbar(data, Pf);
        g = fderiv - f;
        if (fabs(g) < tol) {
            memcpy(P, Pf, 9 * sizeof(double));
            if (fbar_out != 0) {
                *fbar_out = fderiv;
            }
            return 1;
        }

        h = 1.0e-7 * (1.0 + fabs(f));
        fh = f + h;
        if (fh > 1.0) {
            fh = f - h;
        }
        prj_rad_gr_m1_pressure_for_fbar(data, fh, Ph);
        gh = prj_rad_gr_m1_derived_fbar(data, Ph) - fh;

        deriv = (gh - g) / (fh - f);
        if (!isfinite(deriv) || fabs(deriv) < 1.0e-12) {
            return 0;
        }
        fnew = f - g / deriv;
        if (!isfinite(fnew) || fnew < 0.0 || fnew > 1.0) {
            return 0;
        }
        if (fabs(fnew - f) < tol) {
            /* Step has stalled without meeting the residual tolerance: hand off
             * to Illinois rather than accept a possibly spurious root. */
            return 0;
        }
        f = fnew;
    }
    return 0;
}

static void prj_rad_gr_m1_pressure_implicit(
    const prj_rad_gr_m1_pressure_data *data, double P[3][3],
    double *fbar_out)
{
    double Plo[3][3];
    double Phi[3][3];
    double Pmid[3][3];
    double flo;
    double fhi;
    double fmid;
    double glo;
    double ghi;
    double gmid;
    double lo = 0.0;
    double hi = 1.0;
    int iter;

    if (data->beta2 < PRJ_RAD_GR_M1_EULERIAN_FBAR_BETA *
            PRJ_RAD_GR_M1_EULERIAN_FBAR_BETA) {
        prj_rad_gr_m1_pressure_for_fbar(data, data->f_euler, P);
        if (fbar_out != 0) {
            *fbar_out = data->f_euler;
        }
        return;
    }

    /* Fast path: Newton-Raphson seeded from the Eulerian closure. On success it
     * fills P / fbar_out and returns; on any failure it leaves P untouched and
     * we drop through to the Illinois bracketing, which is guaranteed to bracket
     * and converge on the same root. */
    if (prj_rad_gr_m1_pressure_newton(data, data->f_euler, P, fbar_out)) {
        return;
    }

    /* Illinois fallback (reached only when Newton above bailed). Seed with the
     * Eulerian closure f_euler (leading-order root). If it is already within
     * tolerance we are done in one evaluation; otherwise it brackets the root
     * against a single endpoint on a much tighter interval than [0, 1], so
     * Illinois needs fewer iterations. Only if f_euler fails to bracket (rare,
     * non-monotone g) do we fall through to the full [0, 1] bracketing below.
     * Physically equivalent to that bracket, not bit-identical. */
    {
        double fe = data->f_euler;
        double Pe[3][3];
        double fe_derived;
        double ge;

        prj_rad_gr_m1_pressure_for_fbar(data, fe, Pe);
        fe_derived = prj_rad_gr_m1_derived_fbar(data, Pe);
        ge = fe_derived - fe;
        if (fabs(ge) < PRJ_RAD_GR_M1_FBAR_TOL) {
            memcpy(P, Pe, 9 * sizeof(double));
            if (fbar_out != 0) {
                *fbar_out = fe_derived;
            }
            return;
        }
        if (ge > 0.0) {
            /* root in (fe, hi]; confirm the hi=1 endpoint has g < 0 */
            prj_rad_gr_m1_pressure_for_fbar(data, hi, Phi);
            fhi = prj_rad_gr_m1_derived_fbar(data, Phi);
            ghi = fhi - hi;
            if (fabs(ghi) < PRJ_RAD_GR_M1_FBAR_TOL) {
                memcpy(P, Phi, 9 * sizeof(double));
                if (fbar_out != 0) {
                    *fbar_out = fhi;
                }
                return;
            }
            if (ghi < 0.0) {
                lo = fe;
                glo = ge;
                memcpy(Plo, Pe, 9 * sizeof(double));
                goto illinois;
            }
        } else {
            /* root in [lo, fe); confirm the lo=0 endpoint has g > 0 */
            prj_rad_gr_m1_pressure_for_fbar(data, lo, Plo);
            flo = prj_rad_gr_m1_derived_fbar(data, Plo);
            glo = flo - lo;
            if (fabs(glo) < PRJ_RAD_GR_M1_FBAR_TOL) {
                memcpy(P, Plo, 9 * sizeof(double));
                if (fbar_out != 0) {
                    *fbar_out = flo;
                }
                return;
            }
            if (glo > 0.0) {
                hi = fe;
                ghi = ge;
                memcpy(Phi, Pe, 9 * sizeof(double));
                goto illinois;
            }
        }
        /* f_euler did not bracket the root: reset and use full [0, 1]. */
        lo = 0.0;
        hi = 1.0;
    }

    prj_rad_gr_m1_pressure_for_fbar(data, lo, Plo);
    flo = prj_rad_gr_m1_derived_fbar(data, Plo);
    glo = flo - lo;
    if (fabs(glo) < PRJ_RAD_GR_M1_FBAR_TOL) {
        memcpy(P, Plo, 9 * sizeof(double));
        if (fbar_out != 0) {
            *fbar_out = flo;
        }
        return;
    }

    prj_rad_gr_m1_pressure_for_fbar(data, hi, Phi);
    fhi = prj_rad_gr_m1_derived_fbar(data, Phi);
    ghi = fhi - hi;
    if (fabs(ghi) < PRJ_RAD_GR_M1_FBAR_TOL) {
        memcpy(P, Phi, 9 * sizeof(double));
        if (fbar_out != 0) {
            *fbar_out = fhi;
        }
        return;
    }

    if (glo < 0.0 || ghi > 0.0) {
        if (fabs(glo) <= fabs(ghi)) {
            memcpy(P, Plo, 9 * sizeof(double));
            if (fbar_out != 0) {
                *fbar_out = flo;
            }
        } else {
            memcpy(P, Phi, 9 * sizeof(double));
            if (fbar_out != 0) {
                *fbar_out = fhi;
            }
        }
        return;
    }

    /* Illinois-modified regula falsi. The bracket [lo, hi] has glo > 0 and
     * ghi < 0 (guaranteed above), and every step keeps opposite signs at the
     * endpoints, so convergence is guaranteed like bisection -- but the secant
     * step gives superlinear convergence, typically ~6-10 iterations instead of
     * the ~40 bisection needs to reach hi-lo < 1e-13. The Illinois halving of a
     * repeatedly retained endpoint's g-value prevents the classic regula-falsi
     * stagnation where one endpoint never moves. Same root of
     * g(fbar) = derived_fbar(P(fbar)) - fbar to the same 1e-13 tolerance:
     * physically equivalent to the old bisection, not bit-identical. */
illinois:
    {
        int side = 0;

        for (iter = 0; iter < 80; ++iter) {
            double denom = ghi - glo;
            double mid;


            /* False-position estimate; fall back to bisection if the secant is
             * degenerate or would step outside the open bracket. */
            if (!isfinite(denom) || fabs(denom) < 1.0e-300) {
                mid = 0.5 * (lo + hi);
            } else {
                mid = (lo * ghi - hi * glo) / denom;
            }
            if (!isfinite(mid) || mid <= lo || mid >= hi) {
                mid = 0.5 * (lo + hi);
            }

            prj_rad_gr_m1_pressure_for_fbar(data, mid, Pmid);
            fmid = prj_rad_gr_m1_derived_fbar(data, Pmid);
            gmid = fmid - mid;
            if (fabs(gmid) < PRJ_RAD_GR_M1_FBAR_TOL || hi - lo < PRJ_RAD_GR_M1_FBAR_TOL) {
                memcpy(P, Pmid, 9 * sizeof(double));
                if (fbar_out != 0) {
                    *fbar_out = fmid;
                }
                return;
            }
            if (gmid > 0.0) {
                /* Same sign as glo: root lies in [mid, hi]; advance lo. */
                lo = mid;
                glo = gmid;
                memcpy(Plo, Pmid, 9 * sizeof(double));
                if (side == 1) {
                    ghi *= 0.5;   /* hi retained again -> down-weight it */
                }
                side = 1;
            } else {
                /* Same sign as ghi: root lies in [lo, mid]; advance hi. */
                hi = mid;
                ghi = gmid;
                memcpy(Phi, Pmid, 9 * sizeof(double));
                if (side == -1) {
                    glo *= 0.5;   /* lo retained again -> down-weight it */
                }
                side = -1;
            }
        }
    }

    (void)glo;
    (void)ghi;
    {
        double fbar = 0.5 * (lo + hi);

        prj_rad_gr_m1_pressure_for_fbar(data, fbar, P);
        if (fbar_out != 0) {
            *fbar_out = fbar;
        }
    }
}

static void prj_rad_gr_m1_decompose_m2(
    const prj_rad_gr_m1_pressure_data *data, const double P[3][3],
    double *J_out, double H[4], double L[4][4], double ucon[4],
    double ucov[4], double hcon[4][4], double hcov[4][4])
{
    const prj_rad_gr_m1_closure_ctx *ctx = data->ctx;
    double M2[4][4];
    double J = data->J0;
    int a;
    int b;
    int m;

    memset(M2, 0, sizeof(M2));
    M2[0][0] = data->E;
    for (a = 0; a < 3; ++a) {
        M2[0][a + 1] = data->Fhat_con[a];
        M2[a + 1][0] = data->Fhat_con[a];
        for (b = 0; b < 3; ++b) {
            M2[a + 1][b + 1] = P[a][b];
        }
    }

    for (m = 0; m < 6; ++m) {
        int i = prj_rad_gr_m1_sym_i[m];
        int j = prj_rad_gr_m1_sym_j[m];

        J += data->Jcoef[m] * P[i][j];
    }
    *J_out = J;

    memset(H, 0, 4 * sizeof(double));
    for (a = 0; a < 3; ++a) {
        H[a + 1] = data->H0[a];
        for (m = 0; m < 6; ++m) {
            int i = prj_rad_gr_m1_sym_i[m];
            int j = prj_rad_gr_m1_sym_j[m];

            H[a + 1] += data->Hcoef[a][m] * P[i][j];
        }
        H[0] += data->u_cov[a] * H[a + 1];
    }
    if (data->wlor > 0.0) {
        H[0] /= data->wlor;
    } else {
        H[0] = 0.0;
    }

    ucon[0] = data->wlor;
    ucov[0] = -data->wlor;
    for (a = 0; a < 3; ++a) {
        ucon[a + 1] = data->wlor * data->vcon[a];
        ucov[a + 1] = data->u_cov[a];
    }

    memset(hcon, 0, 16 * sizeof(double));
    memset(hcov, 0, 16 * sizeof(double));
    hcon[0][0] = -1.0 + ucon[0] * ucon[0];
    hcov[0][0] = -1.0 + ucov[0] * ucov[0];
    for (a = 0; a < 3; ++a) {
        hcon[0][a + 1] = ucon[0] * ucon[a + 1];
        hcon[a + 1][0] = hcon[0][a + 1];
        hcov[0][a + 1] = ucov[0] * ucov[a + 1];
        hcov[a + 1][0] = hcov[0][a + 1];
        for (b = 0; b < 3; ++b) {
            hcon[a + 1][b + 1] =
                ctx->gamma_inv[a][b] + ucon[a + 1] * ucon[b + 1];
            hcov[a + 1][b + 1] =
                ctx->gamma[a][b] + ucov[a + 1] * ucov[b + 1];
        }
    }

    for (a = 0; a < 4; ++a) {
        for (b = 0; b < 4; ++b) {
            L[a][b] = M2[a][b] - J * ucon[a] * ucon[b] -
                H[a] * ucon[b] - H[b] * ucon[a];
        }
    }
}

static void prj_rad_gr_m1_grad_cu_cov(
    const prj_rad_gr_m1_pressure_data *data, double grad_cu_cov[3][4])
{
    const prj_rad_gr_m1_closure_ctx *ctx = data->ctx;
    double Gamma[3][3][3];
    double du_cov[3][3];
    double D_cov[3][3];
    double dbeta2[3];
    double dw[3];
    int a;
    int b;
    int d;

    for (d = 0; d < 3; ++d) {
        for (a = 0; a < 4; ++a) {
            grad_cu_cov[d][a] = 0.0;
        }
    }
    if (ctx == 0) {
        return;
    }

    prj_rad_gr_m1_christoffel(ctx, Gamma);
    for (d = 0; d < 3; ++d) {
        dbeta2[d] = 0.0;
        for (a = 0; a < 3; ++a) {
            for (b = 0; b < 3; ++b) {
                dbeta2[d] += ctx->dgamma[d][a][b] * data->vcon[a] *
                    data->vcon[b] +
                    2.0 * ctx->gamma[a][b] * data->vcon[a] *
                    ctx->dvdx[d][b];
            }
        }
        dw[d] = 0.5 * data->wlor * data->wlor * data->wlor * dbeta2[d];
    }

    for (d = 0; d < 3; ++d) {
        for (a = 0; a < 3; ++a) {
            double dvcov = 0.0;

            for (b = 0; b < 3; ++b) {
                dvcov += ctx->dgamma[d][a][b] * data->vcon[b] +
                    ctx->gamma[a][b] * ctx->dvdx[d][b];
            }
            du_cov[d][a] = dw[d] * data->vcov[a] + data->wlor * dvcov;
        }
    }

    for (d = 0; d < 3; ++d) {
        for (a = 0; a < 3; ++a) {
            D_cov[d][a] = du_cov[d][a] - data->wlor * ctx->K_dd[d][a];
            for (b = 0; b < 3; ++b) {
                D_cov[d][a] -= Gamma[b][d][a] * data->u_cov[b];
            }
            grad_cu_cov[d][a + 1] = PRJ_CLIGHT * D_cov[d][a];
            grad_cu_cov[d][0] -= data->vcon[a] * grad_cu_cov[d][a + 1];
        }
    }
    /* TODO(GR-M1 frequency): include explicit temporal derivative support for
     * nabla_0(c u_beta) once stage-time fluid derivatives are available.  The
     * current GR transport equations only contract the spatial k index here. */
}

static void prj_rad_gr_m1_frequency_drifts(
    const prj_rad_gr_m1_pressure_data *data, const double P[3][3],
    double fbar, double *energy_drift, double momentum_drift_cov[3])
{
    const prj_rad_gr_m1_closure_ctx *ctx = data->ctx;
    double J;
    double H[4];
    double L[4][4];
    double ucon[4];
    double ucov[4];
    double hcon[4][4];
    double hcov[4][4];
    double Hnorm2 = 0.0;
    double inv_Hnorm3 = 0.0;
    double chi;
    double thin_w;
    double thick_w;
    double grad_cu_cov[3][4];
    double contract_con[4] = {0.0, 0.0, 0.0, 0.0};
    int alpha;
    int beta;
    int gamma;
    int b;
    int d;
    int i;

    prj_rad_gr_m1_decompose_m2(data, P, &J, H, L, ucon, ucov, hcon, hcov);
    prj_rad_gr_m1_grad_cu_cov(data, grad_cu_cov);
    for (alpha = 0; alpha < 4; ++alpha) {
        for (beta = 0; beta < 4; ++beta) {
            Hnorm2 += hcov[alpha][beta] * H[alpha] * H[beta];
        }
    }
    if (isfinite(Hnorm2) && Hnorm2 > 0.0) {
        inv_Hnorm3 = 1.0 / (Hnorm2 * sqrt(Hnorm2));
    }

    if (!isfinite(fbar) || fbar < 0.0) {
        fbar = 0.0;
    } else if (fbar > 1.0) {
        fbar = 1.0;
    }
    chi = prj_rad_m1_chi(data->rad, fbar);
    thin_w = 0.5 * (3.0 * chi - 1.0);
    thick_w = 1.5 * (1.0 - chi);

    for (alpha = 0; alpha < 4; ++alpha) {
        for (b = 0; b < 4; ++b) {
            beta = b;
            for (d = 0; d < 3; ++d) {
                double Nthin = 0.0;
                double Nthick;
                double N;
                double M3;

                gamma = d + 1;
                if (inv_Hnorm3 > 0.0) {
                    Nthin = J * H[alpha] * H[beta] * H[gamma] *
                        inv_Hnorm3;
                }
                Nthick = 0.2 * (H[alpha] * hcon[beta][gamma] +
                    H[beta] * hcon[alpha][gamma] +
                    H[gamma] * hcon[alpha][beta]);
                N = thin_w * Nthin + thick_w * Nthick;
                M3 = J * ucon[alpha] * ucon[beta] * ucon[gamma] +
                    H[alpha] * ucon[beta] * ucon[gamma] +
                    H[beta] * ucon[alpha] * ucon[gamma] +
                    H[gamma] * ucon[alpha] * ucon[beta] +
                    L[alpha][beta] * ucon[gamma] +
                    L[alpha][gamma] * ucon[beta] +
                    L[beta][gamma] * ucon[alpha] + N;

                contract_con[alpha] += M3 * grad_cu_cov[d][beta];
            }
        }
    }

    if (energy_drift != 0) {
        *energy_drift = contract_con[0];
    }
    for (i = 0; i < 3; ++i) {
        momentum_drift_cov[i] = 0.0;
        for (b = 0; b < 3; ++b) {
            /* M3 is assembled as a stress-tensor moment with M^{0i}=F^i/c.
             * Convert the spatial projection back to the evolved physical
             * flux F_i before the frequency-bin update. */
            momentum_drift_cov[i] += PRJ_CLIGHT * ctx->gamma[i][b] *
                contract_con[b + 1];
        }
    }
}

void prj_rad_gr_m1_pressure(const prj_rad *rad,
    const prj_rad_gr_m1_closure_ctx *ctx, double E, const double Fcov_in[3],
    double P[3][3])
{
    prj_rad_gr_m1_side_data side;

    if (P == 0) {
        return;
    }
    if (ctx == 0 || Fcov_in == 0) {
        prj_rad_m1_pressure(rad, E,
            Fcov_in != 0 ? Fcov_in[0] : 0.0,
            Fcov_in != 0 ? Fcov_in[1] : 0.0,
            Fcov_in != 0 ? Fcov_in[2] : 0.0, P);
        return;
    }
    prj_rad_gr_m1_prepare_side(ctx, &side);
    prj_rad_gr_m1_pressure_cached(rad, ctx, &side, E, Fcov_in, P);
}

/* Fast path for callers (interface flux) that solve the closure for many energy
 * groups at one point: the per-side kinematics `side` are built once via
 * prj_rad_gr_m1_prepare_side and reused, so only the E/F/opacity-dependent work
 * runs per group. */
void prj_rad_gr_m1_pressure_cached(const prj_rad *rad,
    const prj_rad_gr_m1_closure_ctx *ctx,
    const prj_rad_gr_m1_side_data *side, double E, const double Fcov_in[3],
    double P[3][3])
{
    prj_rad_gr_m1_pressure_data data;

    if (P == 0 || ctx == 0 || side == 0 || Fcov_in == 0) {
        return;
    }
    prj_rad_gr_m1_prepare_pressure(rad, ctx, side, E, Fcov_in, &data);
    prj_rad_gr_m1_pressure_implicit(&data, P, 0);
}

void prj_rad_gr_m1_pressure_fbar(const prj_rad *rad,
    const prj_rad_gr_m1_closure_ctx *ctx, double E, const double Fcov_in[3],
    double P[3][3], double *fbar_out)
{
    prj_rad_gr_m1_side_data side;

    if (fbar_out != 0) {
        *fbar_out = 0.0;
    }
    if (P == 0) {
        return;
    }
    if (ctx == 0 || Fcov_in == 0) {
        prj_rad_m1_pressure(rad, E,
            Fcov_in != 0 ? Fcov_in[0] : 0.0,
            Fcov_in != 0 ? Fcov_in[1] : 0.0,
            Fcov_in != 0 ? Fcov_in[2] : 0.0, P);
        return;
    }
    prj_rad_gr_m1_prepare_side(ctx, &side);
    prj_rad_gr_m1_pressure_fbar_cached(rad, ctx, &side, E, Fcov_in, P,
        fbar_out);
}

void prj_rad_gr_m1_pressure_fbar_cached(const prj_rad *rad,
    const prj_rad_gr_m1_closure_ctx *ctx,
    const prj_rad_gr_m1_side_data *side, double E, const double Fcov_in[3],
    double P[3][3], double *fbar_out)
{
    prj_rad_gr_m1_pressure_data data;

    if (fbar_out != 0) {
        *fbar_out = 0.0;
    }
    if (P == 0 || ctx == 0 || side == 0 || Fcov_in == 0) {
        return;
    }
    prj_rad_gr_m1_prepare_pressure(rad, ctx, side, E, Fcov_in, &data);
    prj_rad_gr_m1_pressure_implicit(&data, P, fbar_out);
}

/* --- Cell-centered closure (P, fbar) cache, shared source -> frequency flux --- */
#define PRJ_CLCACHE_NCELL (PRJ_BLOCK_SIZE * PRJ_BLOCK_SIZE * PRJ_BLOCK_SIZE)
#define PRJ_CLCACHE_NGRP (PRJ_NRAD * PRJ_NEGROUP)
static double g_clcache_P[PRJ_CLCACHE_NCELL][PRJ_CLCACHE_NGRP][3][3];
static double g_clcache_fbar[PRJ_CLCACHE_NCELL][PRJ_CLCACHE_NGRP];
static const void *g_clcache_block;
static const void *g_clcache_wrad;
static int g_clcache_stage = -1;
static int g_clcache_valid;

void prj_rad_gr_m1_closure_cache_begin(const void *block, const void *wrad,
    int stage)
{
    g_clcache_block = block;
    g_clcache_wrad = wrad;
    g_clcache_stage = stage;
    g_clcache_valid = 1;
}

void prj_rad_gr_m1_closure_cache_put(int active_cell, int grp,
    const double P[3][3], double fbar)
{
    if (!g_clcache_valid || active_cell < 0 ||
        active_cell >= PRJ_CLCACHE_NCELL || grp < 0 || grp >= PRJ_CLCACHE_NGRP) {
        return;
    }
    memcpy(g_clcache_P[active_cell][grp], P, 9 * sizeof(double));
    g_clcache_fbar[active_cell][grp] = fbar;
}

int prj_rad_gr_m1_closure_cache_get(const void *block, const void *wrad,
    int stage, int active_cell, int grp, double P[3][3], double *fbar_out)
{
    if (!g_clcache_valid || block != g_clcache_block ||
        wrad != g_clcache_wrad || stage != g_clcache_stage) {
        return 0;
    }
    if (active_cell < 0 || active_cell >= PRJ_CLCACHE_NCELL || grp < 0 ||
        grp >= PRJ_CLCACHE_NGRP) {
        return 0;
    }
    memcpy(P, g_clcache_P[active_cell][grp], 9 * sizeof(double));
    if (fbar_out != 0) {
        *fbar_out = g_clcache_fbar[active_cell][grp];
    }
    return 1;
}
#endif

#if PRJ_NRAD > 0
/* Contraction M[b] = sum_{a,c} Q^{abc} dvdx[a][c] of the Levermore third moment
 * with the velocity gradient, computed analytically so the 27 components of Q
 * never have to be materialised. With
 *   Q^{abc} = coef_nnn n^a n^b n^c + coef_mix (n^a d_bc + n^b d_ca + n^c d_ab),
 * the contraction collapses to
 *   M[b] = coef_nnn n[b] S + coef_mix (T1[b] + n[b] divv + T3[b]),
 * with S = n.dvdx.n, T1 = n^T dvdx, T3 = dvdx n, divv = tr(dvdx). Mathematically
 * identical to building Q and summing (validated to machine epsilon), with the
 * isotropic E<=0 / Fmag<=0 limit returning zero exactly as m1_third_moment does. */
static void prj_rad_m1_third_moment_contract(const prj_rad *rad, double E,
    double F1, double F2, double F3, const double dvdx[3][3], double M[3])
{
    double E_pos;
    double Fmag;
    double cE;
    double f;
    double q_fac;
    double n[3];
    double coef_nnn;
    double coef_mix;
    double divv;
    double T1[3];
    double T3[3];
    double S;
    int b;

    M[0] = 0.0;
    M[1] = 0.0;
    M[2] = 0.0;

    E_pos = E > 0.0 ? E : 0.0;
    Fmag = sqrt(F1 * F1 + F2 * F2 + F3 * F3);
    cE = PRJ_CLIGHT * E_pos;
    if (cE <= 0.0 || Fmag <= 0.0) {
        return;
    }

    f = Fmag / cE;
    if (f > 1.0) {
        f = 1.0;
    }
    q_fac = prj_rad_levermore_q_factor(rad, f);
    n[0] = F1 / Fmag;
    n[1] = F2 / Fmag;
    n[2] = F3 / Fmag;
    coef_nnn = 0.5 * cE * (5.0 * q_fac - 3.0 * f);
    coef_mix = 0.5 * cE * (f - q_fac);

    divv = dvdx[0][0] + dvdx[1][1] + dvdx[2][2];
    for (b = 0; b < 3; ++b) {
        T1[b] = n[0] * dvdx[0][b] + n[1] * dvdx[1][b] + n[2] * dvdx[2][b];
        T3[b] = dvdx[b][0] * n[0] + dvdx[b][1] * n[1] + dvdx[b][2] * n[2];
    }
    S = n[0] * T3[0] + n[1] * T3[1] + n[2] * T3[2];

    for (b = 0; b < 3; ++b) {
        M[b] = coef_nnn * n[b] * S + coef_mix * (T1[b] + n[b] * divv + T3[b]);
    }
}

static void prj_rad_m1_phys_flux_with_fluxmag(const prj_rad *rad, double E, double F1,
    double F2, double F3, double Fmag, double inv_Fmag, double f,
    double *fE, double *fF1, double *fF2, double *fF3)
{
    double chi;
    double n1;
    double n2;
    double n3;
    double D11;
    double D12;
    double D13;
    double c2;

    c2 = PRJ_CLIGHT * PRJ_CLIGHT;

    if (E <= 0.0 || Fmag <= 0.0) {
        /* isotropic: P = (E/3) I */
        *fE = F1;
        *fF1 = c2 * E / 3.0;
        *fF2 = 0.0;
        *fF3 = 0.0;
        return;
    }

    chi = prj_rad_m1_chi(rad, f);

    n1 = F1 * inv_Fmag;
    n2 = F2 * inv_Fmag;
    n3 = F3 * inv_Fmag;

    {
        double a = 0.5 * (1.0 - chi);
        double b = 0.5 * (3.0 * chi - 1.0);
        D11 = a + b * n1 * n1;
        D12 = b * n1 * n2;
        D13 = b * n1 * n3;
    }

    *fE = F1;
    *fF1 = c2 * E * D11;
    *fF2 = c2 * E * D12;
    *fF3 = c2 * E * D13;
}

static void prj_rad_enforce_flux_limit(double *E, double *F1, double *F2, double *F3,
    double *Fmag_out, double *inv_Fmag_out, double *f_out)
{
    double Fmag;
    double cE;
    double f;
    double scale;

    if (*E < 0.0) {
        *E = 0.0;
    }
    Fmag = sqrt((*F1) * (*F1) + (*F2) * (*F2) + (*F3) * (*F3));
    cE = PRJ_CLIGHT * (*E);
    if (Fmag > cE && Fmag > 0.0) {
        scale = cE / Fmag;
        *F1 *= scale;
        *F2 *= scale;
        *F3 *= scale;
        Fmag = cE;
    }
    if (cE > 0.0) {
        f = Fmag / cE;
        if (f > 1.0) {
            f = 1.0;
        }
    } else {
        f = 0.0;
    }
    *Fmag_out = Fmag;
    *inv_Fmag_out = (Fmag > 0.0) ? (1.0 / Fmag) : 0.0;
    *f_out = f;
}

void prj_rad_m1_wavespeeds_with_fluxmag(double E, double F1, double Fmag, double inv_Fmag,
    double f, double *lam_min, double *lam_max)
{
    double mu;
    double fsq;
    double ffac;
    double inv_ffac;
    double lterm;

    if (E <= 0.0 || Fmag <= 0.0) {
        *lam_min = -1.0 / sqrt(3.0);
        *lam_max = 1.0 / sqrt(3.0);
        return;
    }

    mu = F1 * inv_Fmag;

    fsq = f * f;
    ffac = sqrt(4.0 - 3.0 * fsq);
    inv_ffac = 1.0 / ffac;
    lterm = sqrt(fabs((2.0 / 3.0) * (4.0 - 3.0 * fsq - ffac) + 2.0 * mu * mu * (2.0 - fsq - ffac)));
    *lam_min = (mu * f - lterm) * inv_ffac;
    *lam_max = (mu * f + lterm) * inv_ffac;
    if (*lam_min < -1.0) {
        *lam_min = -1.0;
    }
    if (*lam_max > 1.0) {
        *lam_max = 1.0;
    }
}

void prj_rad_m1_wavespeeds(double E, double F1, double F2, double F3,
    double *lam_min, double *lam_max)
{
    double Fmag;
    double inv_Fmag;
    double cE;
    double f;

    cE = PRJ_CLIGHT * (E > 0.0 ? E : 0.0);
    Fmag = sqrt(F1 * F1 + F2 * F2 + F3 * F3);
    inv_Fmag = (Fmag > 0.0) ? (1.0 / Fmag) : 0.0;
    if (cE > 0.0) {
        f = Fmag / cE;
        if (f > 1.0) {
            f = 1.0;
        }
    } else {
        f = 0.0;
    }
    prj_rad_m1_wavespeeds_with_fluxmag(E, F1, Fmag, inv_Fmag, f, lam_min, lam_max);
}
#endif

void prj_rad_flux(const prj_rad *rad, const double *WL, const double *WR,
    double lapse, const double *chi_face,
    double dx_dir, double v_face, double *flux)
{
#if PRJ_USE_RADIATION_M1
    {
        int field;
        int group;

        for (field = 0; field < PRJ_NRAD; ++field) {
            for (group = 0; group < PRJ_NEGROUP; ++group) {
                int idx = field * PRJ_NEGROUP + group;
                double EL = WL[PRJ_PRIM_RAD_E(field, group)];
                double ER = WR[PRJ_PRIM_RAD_E(field, group)];
                double F1L = WL[PRJ_PRIM_RAD_F1(field, group)];
                double F2L = WL[PRJ_PRIM_RAD_F2(field, group)];
                double F3L = WL[PRJ_PRIM_RAD_F3(field, group)];
                double F1R = WR[PRJ_PRIM_RAD_F1(field, group)];
                double F2R = WR[PRJ_PRIM_RAD_F2(field, group)];
                double F3R = WR[PRJ_PRIM_RAD_F3(field, group)];
                double Fmag_L;
                double inv_Fmag_L;
                double f_L;
                double Fmag_R;
                double inv_Fmag_R;
                double f_R;
                double fLE;
                double fLF1;
                double fLF2;
                double fLF3;
                double fRE;
                double fRF1;
                double fRF2;
                double fRF3;
                double lamL_min;
                double lamL_max;
                double lamR_min;
                double lamR_max;
                double sL;
                double sR;
                double denom;
                double inv_denom;
                double chi_ext;
                double tau;
                double eps;
                double eps2;

                prj_rad_enforce_flux_limit(&EL, &F1L, &F2L, &F3L, &Fmag_L, &inv_Fmag_L, &f_L);
                prj_rad_enforce_flux_limit(&ER, &F1R, &F2R, &F3R, &Fmag_R, &inv_Fmag_R, &f_R);

                prj_rad_m1_phys_flux_with_fluxmag(rad, EL, F1L, F2L, F3L, Fmag_L, inv_Fmag_L, f_L,
                    &fLE, &fLF1, &fLF2, &fLF3);
                prj_rad_m1_phys_flux_with_fluxmag(rad, ER, F1R, F2R, F3R, Fmag_R, inv_Fmag_R, f_R,
                    &fRE, &fRF1, &fRF2, &fRF3);

                prj_rad_m1_wavespeeds_with_fluxmag(EL, F1L, Fmag_L, inv_Fmag_L, f_L,
                    &lamL_min, &lamL_max);
                prj_rad_m1_wavespeeds_with_fluxmag(ER, F1R, Fmag_R, inv_Fmag_R, f_R,
                    &lamR_min, &lamR_max);
                sL = PRJ_CLIGHT * (lamL_min < lamR_min ? lamL_min : lamR_min);
                sR = PRJ_CLIGHT * (lamL_max > lamR_max ? lamL_max : lamR_max);
                if (sL > 0.0) {
                    sL = 0.0;
                }
                if (sR < 0.0) {
                    sR = 0.0;
                }
                if (sR - sL < 1.0e-30) {
                    sL = -PRJ_CLIGHT;
                    sR = PRJ_CLIGHT;
                }
                denom = sR - sL;
                inv_denom = 1.0 / denom;

                chi_ext = chi_face[idx];
                tau = chi_ext * dx_dir;
                eps = 3.0 / (5.0 * tau + 1.0e-10);
                if (eps > 1.0) {
                    eps = 1.0;
                }
                eps2 = eps*eps;

                /* Equation 49 and 50 of Audit et al. 2002 */
                flux[PRJ_CONS_RAD_E(field, group)] = lapse *
                    (sR * fLE - sL * fRE + eps * sL * sR * (ER - EL)) * inv_denom;
                flux[PRJ_CONS_RAD_F1(field, group)] = lapse *
                    ((eps2*(sR * fLF1 - sL * fRF1) + eps * sL * sR * (F1R - F1L)) * inv_denom
                    +(1-eps2)*(fLF1+fRF1)*0.5);
                flux[PRJ_CONS_RAD_F2(field, group)] = lapse *
                    ((eps2*(sR * fLF2 - sL * fRF2) + eps * sL * sR * (F2R - F2L)) * inv_denom
                    +(1-eps2)*(fLF2+fRF2)*0.5);
                flux[PRJ_CONS_RAD_F3(field, group)] = lapse *
                    ((eps2*(sR * fLF3 - sL * fRF3) + eps * sL * sR * (F3R - F3L)) * inv_denom
                    +(1-eps2)*(fLF3+fRF3)*0.5);

                /* O(v/c) fluid advection: upwinded v_face * {E, F_i} term. */
                {
                    double E_up = v_face >= 0.0 ? EL : ER;
                    double F1_up = v_face >= 0.0 ? F1L : F1R;
                    double F2_up = v_face >= 0.0 ? F2L : F2R;
                    double F3_up = v_face >= 0.0 ? F3L : F3R;

                    flux[PRJ_CONS_RAD_E(field, group)] += lapse * v_face * E_up;
                    flux[PRJ_CONS_RAD_F1(field, group)] += lapse * v_face * F1_up;
                    flux[PRJ_CONS_RAD_F2(field, group)] += lapse * v_face * F2_up;
                    flux[PRJ_CONS_RAD_F3(field, group)] += lapse * v_face * F3_up;
                }
            }
        }
    }
#else
    (void)rad;
    (void)WL;
    (void)WR;
    (void)chi_face;
    (void)dx_dir;
    (void)v_face;
    (void)lapse;
    (void)flux;
#endif
}

#if PRJ_USE_RADIATION_FSA
/* Clamp negative angular intensities to zero.  A negative J is an unphysical
 * numerical undershoot (the FSA spatial advection is not strictly positivity-
 * preserving at a marginal multidimensional CFL), not real radiation.  We
 * discard it WITHOUT any matter back-reaction: the energy/Ye/momentum the gas
 * would appear to exchange for this J change is fictitious (it is not emission
 * or absorption), so only the radiation slots are touched here and the hydro
 * slots are left untouched.  Applied to the conserved state before the
 * radiation update so the transport and the matter coupling never see J < 0. */
void prj_rad_fsa_clamp_intensities(double *u)
{
    int field;
    int group;
    int angle;

    if (u == 0) {
        return;
    }
    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            for (angle = 0; angle < PRJ_NANGLE; ++angle) {
                int iv = PRJ_CONS_RAD_I(field, group, angle);

                if (u[iv] < 0.0) {
                    u[iv] = 0.0;
                }
            }
        }
    }
}

void prj_rad_flux_fsa(const prj_rad *rad, const prj_block *block,
    const double *WL, const double *WR, double lapse, int dir, double v_face,
    int il, int jl, int kl, int ir, int jr, int kr, double *flux)
{
    int field;
    int group;
    int angle;

    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            for (angle = 0; angle < PRJ_NANGLE; ++angle) {
                int v = PRJ_CONS_RAD_I(field, group, angle);
                double nL[3];
                double nR[3];
                double n_face_dir;
                double speed;
                double J_face;

                prj_rad_fsa_rotated_angle_dir(rad, block, angle, il, jl, kl, nL);
                prj_rad_fsa_rotated_angle_dir(rad, block, angle, ir, jr, kr, nR);
                n_face_dir = 0.5 * (nL[dir] + nR[dir]);
                speed = v_face + lapse * PRJ_CLIGHT * n_face_dir;
                J_face = speed >= 0.0 ? WL[PRJ_PRIM_RAD_I(field, group, angle)] :
                    WR[PRJ_PRIM_RAD_I(field, group, angle)];
                flux[v] = speed * J_face;
            }
        }
    }
}
#endif

#if PRJ_NRAD > 0
/* Building-block derivatives of the implicit residual w.r.t. (lnT, Ye), filled
 * on request by prj_rad_implicit_residuals(). */
typedef struct prj_rad_resid_deriv {
    double dlnkappa_dlnT[PRJ_NRAD * PRJ_NEGROUP];
    double dlnkappa_dYe[PRJ_NRAD * PRJ_NEGROUP];
    double dlneta_dlnT[PRJ_NRAD * PRJ_NEGROUP];
    double dlneta_dYe[PRJ_NRAD * PRJ_NEGROUP];
    double deint_dlnT;
    double deint_dYe;
} prj_rad_resid_deriv;

static void prj_rad_energy_failure_diagnostics(const char *reason,
    const double *u_input, double dt, double lapse, int iter, int maxiter,
    double res, double tol2, double T, double Ye, double F1, double F2,
    double rho, double Uint_old, double Ye_old)
{
    int v;

    fprintf(stderr, "prj_rad_energy_update: %s\n", reason);
    fprintf(stderr,
        "  solver: iter=%d maxiter=%d res=%.17e tol2=%.17e "
        "T=%.17e Ye=%.17e F1=%.17e F2=%.17e\n",
        iter, maxiter, res, tol2, T, Ye, F1, F2);
    fprintf(stderr,
        "  derived input: rho=%.17e Uint_old=%.17e Ye_old=%.17e\n",
        rho, Uint_old, Ye_old);
    fprintf(stderr,
        "  raw input: dt=%.17e lapse=%.17e PRJ_NVAR_CONS=%d "
        "PRJ_NRAD=%d PRJ_NEGROUP=%d PRJ_MHD=%d\n",
        dt, lapse, PRJ_NVAR_CONS, PRJ_NRAD, PRJ_NEGROUP, PRJ_MHD);
    fprintf(stderr, "  double u[PRJ_NVAR_CONS] = {\n");
    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        fprintf(stderr, "      [%d] = %.17e,\n", v, u_input[v]);
    }
    fprintf(stderr, "  };\n");
    fflush(stderr);
}

static void prj_rad_implicit_residuals(prj_rad *rad, prj_eos *eos, double *u,
    double dt, double lapse, double rho, double Uint_old, double Ye_old,
    const double *E_nu_old, double T, double Ye, double *F1, double *F2,
    double *E_nu_new_out, double *kappa_out, double *eta_out, prj_rad_resid_deriv *deriv)
{
    double kappa_local[PRJ_NRAD * PRJ_NEGROUP];
    double eta_local[PRJ_NRAD * PRJ_NEGROUP];
    double dlnkappa_dlnT_local[PRJ_NRAD * PRJ_NEGROUP];
    double dlnkappa_dYe_local[PRJ_NRAD * PRJ_NEGROUP];
    double dlneta_dlnT_local[PRJ_NRAD * PRJ_NEGROUP];
    double dlneta_dYe_local[PRJ_NRAD * PRJ_NEGROUP];
    double *kappa = kappa_out != 0 ? kappa_out : kappa_local;
    double *eta = eta_out != 0 ? eta_out : eta_local;
    double *dlnkappa_dlnT = deriv != 0 ? deriv->dlnkappa_dlnT : dlnkappa_dlnT_local;
    double *dlnkappa_dYe = deriv != 0 ? deriv->dlnkappa_dYe : dlnkappa_dYe_local;
    double *dlneta_dlnT = deriv != 0 ? deriv->dlneta_dlnT : dlneta_dlnT_local;
    double *dlneta_dYe = deriv != 0 ? deriv->dlneta_dYe : dlneta_dYe_local;
    double eint_new;
    double deint_dlnT;
    double deint_dYe;
    double Uint_new;
    double sum_dE = 0.0;
    double sum_dE_xe = 0.0;
    int nu;
    int g;

    (void)u;
    prj_rad3_opac_lookup_ke(rad, rho, T, Ye, kappa, eta,
        dlnkappa_dlnT, dlnkappa_dYe, dlneta_dlnT, dlneta_dYe);
    eint_new = prj_eos_rty_eint(eos, rho, T, Ye, &deint_dlnT, &deint_dYe, PRJ_EOS_CTX_MAIN);
    Uint_new = rho * eint_new;

    for (nu = 0; nu < PRJ_NRAD; ++nu) {
        for (g = 0; g < PRJ_NEGROUP; ++g) {
            int idx = nu * PRJ_NEGROUP + g;
            double num = E_nu_old[idx] + dt * lapse * eta[idx];
            double den = 1.0 + dt * lapse * PRJ_CLIGHT * kappa[idx];
            double Enew = num / den;
            double dE = Enew - E_nu_old[idx];

            sum_dE += dE;
            sum_dE_xe += dE * rad->x_e[nu][g];
            if (E_nu_new_out != 0) {
                E_nu_new_out[idx] = Enew;
            }
        }
    }

    /* sum_dE is a radiation-energy change in RAD_SCALE*erg units; multiply back
       to erg to balance the gas internal energy.  sum_dE_xe already carries
       RAD_SCALE through x_e, so the lepton residual needs no extra factor. */
    *F1 = Uint_new - Uint_old + sum_dE * RAD_SCALE;
    *F2 = rho * Ye - rho * Ye_old + sum_dE_xe;

    if (deriv != 0) {
        deriv->deint_dlnT = deint_dlnT;
        deriv->deint_dYe = deint_dYe;
    }
}

static void prj_rad_implicit_jacobian_from_deriv(const prj_rad *rad,
    const double *E_nu_old, const double *E_nu_new, const double *kappa,
    const prj_rad_resid_deriv *deriv, double dt, double lapse, double rho,
    double T, double *dFdT_1, double *dFdT_2, double *dFdY_1, double *dFdY_2)
{
    double dF1_dlnT;
    double dF2_dlnT;
    double dF1_dYe;
    double dF2_dYe;
    double dt_lapse;
    double inv_T;
    int nu;
    int g;

    dF1_dlnT = rho * deriv->deint_dlnT;
    dF2_dlnT = 0.0;
    dF1_dYe = rho * deriv->deint_dYe;
    dF2_dYe = rho;
    dt_lapse = dt * lapse;
    inv_T = T > 0.0 ? 1.0 / T : 0.0;

    for (nu = 0; nu < PRJ_NRAD; ++nu) {
        for (g = 0; g < PRJ_NEGROUP; ++g) {
            int idx = nu * PRJ_NEGROUP + g;
            double den = 1.0 + dt_lapse * PRJ_CLIGHT * kappa[idx];
            double eta = 0.0;
            double dE_dlnT = 0.0;
            double dE_dYe = 0.0;

            if (dt_lapse != 0.0) {
                eta = (E_nu_new[idx] * den - E_nu_old[idx]) / dt_lapse;
                dE_dlnT = dt_lapse *
                    (eta * deriv->dlneta_dlnT[idx] -
                        E_nu_new[idx] * PRJ_CLIGHT * kappa[idx] *
                            deriv->dlnkappa_dlnT[idx]) / den;
                dE_dYe = dt_lapse *
                    (eta * deriv->dlneta_dYe[idx] -
                        E_nu_new[idx] * PRJ_CLIGHT * kappa[idx] *
                            deriv->dlnkappa_dYe[idx]) / den;
            }

            dF1_dlnT += RAD_SCALE * dE_dlnT;
            dF2_dlnT += rad->x_e[nu][g] * dE_dlnT;
            dF1_dYe += RAD_SCALE * dE_dYe;
            dF2_dYe += rad->x_e[nu][g] * dE_dYe;
        }
    }

    *dFdT_1 = dF1_dlnT * inv_T;
    *dFdT_2 = dF2_dlnT * inv_T;
    *dFdY_1 = dF1_dYe;
    *dFdY_2 = dF2_dYe;
}

static void prj_rad_energy_update_impl(prj_rad *rad, prj_eos *eos, double *u, double dt, double lapse, double *final_temperature, double *final_ye, double *kappa_out, double *eta_out)
{
    double u_input[PRJ_NVAR_CONS];
    double E_nu_old[PRJ_NRAD * PRJ_NEGROUP];
    double E_nu_new[PRJ_NRAD * PRJ_NEGROUP];
    double last_kappa[PRJ_NRAD * PRJ_NEGROUP];
    double last_eta[PRJ_NRAD * PRJ_NEGROUP];
    /* eta is threaded out only when the caller asks for it (FSA); NULL keeps the
     * M1 path free of the extra per-iteration eta copy. */
    double *eta_capture = eta_out != 0 ? last_eta : 0;
    double rho;
    double KE;
    double Emag = 0.0;
    double Uint_old;
    double Ye_old;
    double eint_old;
    double eos_q[PRJ_EOS_NQUANT];
    double T;
    double Ye;
    double err_scale_1;
    double err_scale_2;
    double res_cur;
    double cached_F1 = 0.0;
    double cached_F2 = 0.0;
    double cached_res = 0.0;
    prj_rad_resid_deriv cached_deriv = {0};
    int have_cached_residual = 0;
    int have_final_residual = 0;
    int iter;
    int nu;
    int g;
    int v;
    const double alpha_ls = 1.0e-4;
    const double tol2 = rad->implicit_err_tol * rad->implicit_err_tol;

    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        u_input[v] = u[v];
    }

    rho = u[PRJ_CONS_RHO];
    KE = 0.5 * (u[PRJ_CONS_MOM1] * u[PRJ_CONS_MOM1] +
        u[PRJ_CONS_MOM2] * u[PRJ_CONS_MOM2] +
        u[PRJ_CONS_MOM3] * u[PRJ_CONS_MOM3]) / rho;
#if PRJ_MHD
    Emag = 0.5 * (u[PRJ_CONS_B1] * u[PRJ_CONS_B1] +
        u[PRJ_CONS_B2] * u[PRJ_CONS_B2] +
        u[PRJ_CONS_B3] * u[PRJ_CONS_B3]);
#endif
    Uint_old = u[PRJ_CONS_ETOT] - KE - Emag;
    Ye_old = u[PRJ_CONS_YE] / rho;
    eint_old = Uint_old / rho;

    for (nu = 0; nu < PRJ_NRAD; ++nu) {
        for (g = 0; g < PRJ_NEGROUP; ++g) {
            E_nu_old[nu * PRJ_NEGROUP + g] = u[PRJ_CONS_RAD_E(nu, g)];
        }
    }

    prj_eos_rey(eos, rho, eint_old, Ye_old, eos_q, PRJ_EOS_CTX_MAIN);
    T = eos_q[PRJ_EOS_TEMPERATURE];
    Ye = Ye_old;
    err_scale_1 = fabs(Uint_old) > 0.0 ? Uint_old : 1.0;
    err_scale_2 = fabs(rho * Ye_old) > 0.0 ? rho * Ye_old : 1.0;
    res_cur = 1.0e30;

    for (iter = 0; iter < rad->maxiter; ++iter) {
        double F1;
        double F2;
        double f1;
        double f2;
        prj_rad_resid_deriv deriv;
        prj_rad_resid_deriv trial_deriv;
        double dFdT_1;
        double dFdT_2;
        double dFdY_1;
        double dFdY_2;
        double J00;
        double J01;
        double J10;
        double J11;
        double r0;
        double r1;
        double col_scale0;
        double col_scale1;
        double det;
        double s0;
        double s1;
        double dT;
        double dY;
        double step_scale;
        double gradf0;
        double gradf1;
        double gradfdx;
        double lam;
        double lamold;
        double resold;
        double Ttrial;
        double Ytrial;
        double res_trial;
        double F1_trial;
        double F2_trial;
        int inner_iter;
        int accepted_trial;

        if (have_cached_residual) {
            F1 = cached_F1;
            F2 = cached_F2;
            res_cur = cached_res;
            deriv = cached_deriv;
            have_cached_residual = 0;
            have_final_residual = 1;
        } else {
            prj_rad_implicit_residuals(rad, eos, u, dt, lapse, rho, Uint_old, Ye_old,
                E_nu_old, T, Ye, &F1, &F2, E_nu_new, last_kappa, eta_capture, &deriv);
            f1 = F1 / err_scale_1;
            f2 = F2 / err_scale_2;
            res_cur = 0.5 * (f1 * f1 + f2 * f2);
            have_final_residual = 1;
        }

        if (res_cur < tol2) {
            break;
        }

        prj_rad_implicit_jacobian_from_deriv(rad, E_nu_old, E_nu_new, last_kappa,
            &deriv, dt, lapse, rho, T, &dFdT_1, &dFdT_2, &dFdY_1, &dFdY_2);

        /* grad(½||F||²) = J^T F  (unscaled, for line search). */
        gradf0 = dFdT_1 * F1 / (err_scale_1 * err_scale_1) +
            dFdT_2 * F2 / (err_scale_2 * err_scale_2);
        gradf1 = dFdY_1 * F1 / (err_scale_1 * err_scale_1) +
            dFdY_2 * F2 / (err_scale_2 * err_scale_2);

        /* Row-max equilibration on the augmented matrix [J | -F]. */
        J00 = dFdT_1; J01 = dFdY_1; r0 = F1;
        J10 = dFdT_2; J11 = dFdY_2; r1 = F2;
        {
            double row_max;

            row_max = fabs(J00) > fabs(J01) ? fabs(J00) : fabs(J01);
            if (row_max == 0.0) row_max = 1.0e-10;
            J00 /= row_max; J01 /= row_max; r0 /= row_max;

            row_max = fabs(J10) > fabs(J11) ? fabs(J10) : fabs(J11);
            if (row_max == 0.0) row_max = 1.0e-10;
            J10 /= row_max; J11 /= row_max; r1 /= row_max;
        }
        /* Column-max equilibration. */
        col_scale0 = fabs(J00) > fabs(J10) ? fabs(J00) : fabs(J10);
        col_scale1 = fabs(J01) > fabs(J11) ? fabs(J01) : fabs(J11);
        if (col_scale0 == 0.0) col_scale0 = 1.0e-16;
        if (col_scale1 == 0.0) col_scale1 = 1.0e-16;
        J00 /= col_scale0; J10 /= col_scale0;
        J01 /= col_scale1; J11 /= col_scale1;

        /* Negate RHS for Newton step. */
        r0 = -r0;
        r1 = -r1;

        /* 2x2 Cramer solve. */
        det = J00 * J11 - J01 * J10;
        if (fabs(det) < 1.0e-30) {
            prj_rad_energy_failure_diagnostics("singular Jacobian", u_input,
                dt, lapse, iter, rad->maxiter, res_cur, tol2, T, Ye, F1, F2,
                rho, Uint_old, Ye_old);
            fprintf(stderr, "  scaled Jacobian determinant: %.17e\n", det);
            fflush(stderr);
            exit(1);
        }
        s0 = (J11 * r0 - J01 * r1) / det;
        s1 = (-J10 * r0 + J00 * r1) / det;

        /* Undo column scaling. */
        dT = s0 / col_scale0;
        dY = s1 / col_scale1;

        /* Step limiter: cap at 3% of current values. */
        step_scale = 1.0;
        if (fabs(dT) > 0.03 * T) {
            step_scale = 0.03 * T / fabs(dT);
        }
        if (fabs(dY) > 0.03 * Ye && 0.03 * Ye / fabs(dY) < step_scale) {
            step_scale = 0.03 * Ye / fabs(dY);
        }
        dT *= step_scale;
        dY *= step_scale;

        /* Directional derivative for Armijo check. */
        gradfdx = gradf0 * dT + gradf1 * dY;

        /* Backtracking line search with cubic/quadratic interpolation. */
        lam = 1.0;
        lamold = 0.0;
        resold = 0.0;
        accepted_trial = 0;
        F1_trial = F1;
        F2_trial = F2;
        res_trial = res_cur;
        for (inner_iter = 0; inner_iter < 6; ++inner_iter) {
            double ft1;
            double ft2;

            Ttrial = T + lam * dT;
            Ytrial = Ye + lam * dY;
            prj_rad_implicit_residuals(rad, eos, u, dt, lapse, rho, Uint_old, Ye_old,
                E_nu_old, Ttrial, Ytrial, &F1_trial, &F2_trial, E_nu_new,
                last_kappa, eta_capture, &trial_deriv);
            ft1 = F1_trial / err_scale_1;
            ft2 = F2_trial / err_scale_2;
            res_trial = 0.5 * (ft1 * ft1 + ft2 * ft2);

            if (res_trial < tol2 || res_trial < res_cur + alpha_ls * lam * gradfdx) {
                accepted_trial = 1;
                break;
            }

            {
                double templam;

                if (inner_iter == 0) {
                    templam = -gradfdx / (2.0 * (res_trial - res_cur - gradfdx));
                } else {
                    double rhs1 = res_trial - res_cur - lam * gradfdx;
                    double rhs2 = resold - res_cur - lamold * gradfdx;
                    double a_c = (rhs1 / (lam * lam) - rhs2 / (lamold * lamold)) / (lam - lamold);
                    double b_c = (-lamold * rhs1 / (lam * lam) + lam * rhs2 / (lamold * lamold)) / (lam - lamold);

                    if (a_c == 0.0) {
                        templam = -gradfdx / (2.0 * b_c);
                    } else {
                        double disc = b_c * b_c - 3.0 * a_c * gradfdx;

                        templam = disc >= 0.0 ? (-b_c + sqrt(disc)) / (3.0 * a_c) : 0.5 * lam;
                    }
                    if (templam > 0.5 * lam) {
                        templam = 0.5 * lam;
                    }
                }
                lamold = lam;
                resold = res_trial;
                lam = templam > 0.1 * lam ? templam : 0.1 * lam;
            }
        }

        T = T + lam * dT;
        Ye = Ye + lam * dY;
        if (T <= 0.0) {
            T = 0.5 * (T - lam * dT);
            accepted_trial = 0;
            have_final_residual = 0;
        }
        if (accepted_trial) {
            cached_F1 = F1_trial;
            cached_F2 = F2_trial;
            cached_res = res_trial;
            cached_deriv = trial_deriv;
            have_cached_residual = 1;
            have_final_residual = 1;
        } else {
            have_final_residual = 0;
        }

        res_cur = res_trial;
        if (have_final_residual && res_cur < tol2) {
            break;
        }
    }
    if (iter == rad->maxiter) {
        double F1_final;
        double F2_final;
        double f1_final;
        double f2_final;

        prj_rad_implicit_residuals(rad, eos, u, dt, lapse, rho, Uint_old, Ye_old,
            E_nu_old, T, Ye, &F1_final, &F2_final, E_nu_new, last_kappa, eta_capture, 0);
        f1_final = F1_final / err_scale_1;
        f2_final = F2_final / err_scale_2;
        res_cur = 0.5 * (f1_final * f1_final + f2_final * f2_final);
        prj_rad_energy_failure_diagnostics("failed to converge", u_input,
            dt, lapse, iter, rad->maxiter, res_cur, tol2, T, Ye,
            F1_final, F2_final, rho, Uint_old, Ye_old);
        exit(1);
    }

    /* Final pass at converged (T, Ye) to populate E_nu_new if the accepted
     * line-search residual was not already evaluated at the final state. */
    {
        double eint_new;

        if (!have_final_residual) {
            double F1_final;
            double F2_final;

            prj_rad_implicit_residuals(rad, eos, u, dt, lapse, rho, Uint_old, Ye_old,
                E_nu_old, T, Ye, &F1_final, &F2_final, E_nu_new, last_kappa, eta_capture, 0);
        }
        prj_eos_rty(eos, rho, T, Ye, eos_q, PRJ_EOS_CTX_MAIN);
        eint_new = eos_q[PRJ_EOS_EINT];
        u[PRJ_CONS_ETOT] = rho * eint_new + KE + Emag;
        u[PRJ_CONS_YE] = rho * Ye;
        for (nu = 0; nu < PRJ_NRAD; ++nu) {
            for (g = 0; g < PRJ_NEGROUP; ++g) {
                u[PRJ_CONS_RAD_E(nu, g)] = E_nu_new[nu * PRJ_NEGROUP + g];
            }
        }
    }

    if (final_temperature != 0) {
        *final_temperature = T;
    }
    if (final_ye != 0) {
        *final_ye = Ye;
    }
    if (kappa_out != 0) {
        int i;
        for (i = 0; i < PRJ_NRAD * PRJ_NEGROUP; ++i) {
            kappa_out[i] = last_kappa[i];
        }
    }
    if (eta_out != 0) {
        int i;
        for (i = 0; i < PRJ_NRAD * PRJ_NEGROUP; ++i) {
            eta_out[i] = last_eta[i];
        }
    }
}

void prj_rad_energy_update(prj_rad *rad, prj_eos *eos, double *u, double dt, double lapse, double *final_temperature, double *kappa_out)
{
    prj_rad_energy_update_impl(rad, eos, u, dt, lapse, final_temperature, 0, kappa_out, 0);
}

#if PRJ_USE_RADIATION_FSA
/* Same converged implicit energy/lepton solve as prj_rad_energy_update, but also
 * returns the converged Ye and the emissivity eta at the converged (T, Ye).  The
 * FSA energy-momentum update reuses kappa and eta from here (and looks up only
 * sigma/delta at the same converged Ye) so the opacity is self-consistent and no
 * redundant kappa/eta interpolation is done. */
void prj_rad_energy_update_fsa(prj_rad *rad, prj_eos *eos, double *u, double dt, double lapse, double *final_temperature, double *final_ye, double *kappa_out, double *eta_out)
{
    prj_rad_energy_update_impl(rad, eos, u, dt, lapse, final_temperature, final_ye, kappa_out, eta_out);
}
#endif

void prj_rad_momentum_update(prj_rad *rad, prj_eos *eos, double *u, double dt, double lapse, double temperature, const double *kappa_in)
{
    double sigma[PRJ_NRAD * PRJ_NEGROUP];
    double delta[PRJ_NRAD * PRJ_NEGROUP];
    double rho;
    double Ye;
    double inv_c2;
    double dmom[3];
    double e_unchanged;
    int nu;
    int g;
    int d;

    (void)eos;

    rho = u[PRJ_CONS_RHO];
    Ye = u[PRJ_CONS_YE] / rho;
    inv_c2 = 1.0 / (PRJ_CLIGHT * PRJ_CLIGHT);

    prj_rad3_opac_lookup(rad, rho, temperature, Ye, 0, sigma, delta, 0);

    dmom[0] = 0.0;
    dmom[1] = 0.0;
    dmom[2] = 0.0;
    e_unchanged = u[PRJ_CONS_ETOT] - 0.5 * (u[PRJ_CONS_MOM1] * u[PRJ_CONS_MOM1] +
                                            u[PRJ_CONS_MOM2] * u[PRJ_CONS_MOM2] +
                                            u[PRJ_CONS_MOM3] * u[PRJ_CONS_MOM3]) / rho;

    for (nu = 0; nu < PRJ_NRAD; ++nu) {
        for (g = 0; g < PRJ_NEGROUP; ++g) {
            int idx = nu * PRJ_NEGROUP + g;
            double chi = kappa_in[idx] + sigma[idx] * (1.0 - delta[idx] / 3.0);
            double factor = 1.0 / (1.0 + dt * lapse * PRJ_CLIGHT * chi) - 1.0;
            int fi[3];

            fi[0] = PRJ_CONS_RAD_F1(nu, g);
            fi[1] = PRJ_CONS_RAD_F2(nu, g);
            fi[2] = PRJ_CONS_RAD_F3(nu, g);

            double F_old[PRJ_NDIM];
            for (d = 0; d < 3; ++d) {
                F_old[d] = u[fi[d]];
                double dF = F_old[d] * factor;
                u[fi[d]] = F_old[d] + dF;
                dmom[d] += dF * inv_c2;
            }

            double E = u[PRJ_CONS_RAD_E(nu, g)];
            double F1 = u[fi[0]];
            double F2 = u[fi[1]];
            double F3 = u[fi[2]];
            double Fmag = sqrt(F1 * F1 + F2 * F2 + F3 * F3);
            double cE = PRJ_CLIGHT * E;

            if (E > 0.0 && Fmag > cE) {
                double scale = cE / Fmag;
                u[fi[0]] = F1 * scale;
                u[fi[1]] = F2 * scale;
                u[fi[2]] = F3 * scale;
            }
        }
    }

    /* dmom/detot accumulate radiation-flux changes in RAD_SCALE*erg units;
       multiply back to erg for the gas momentum/energy back-reaction. */
    u[PRJ_CONS_MOM1] -= dmom[0] * RAD_SCALE;
    u[PRJ_CONS_MOM2] -= dmom[1] * RAD_SCALE;
    u[PRJ_CONS_MOM3] -= dmom[2] * RAD_SCALE;
    u[PRJ_CONS_ETOT] = e_unchanged + 0.5 * (u[PRJ_CONS_MOM1] * u[PRJ_CONS_MOM1] +
                                            u[PRJ_CONS_MOM2] * u[PRJ_CONS_MOM2] +
                                            u[PRJ_CONS_MOM3] * u[PRJ_CONS_MOM3]) / rho;
}

#if PRJ_USE_RADIATION_FSA
void prj_rad_energy_momentum_update_fsa(prj_rad *rad, const prj_block *block,
    int ic, int jc, int kc, prj_eos *eos, double *u, double dt, double lapse)
{
    double u_tmp[PRJ_NVAR_CONS];
    double kappa[PRJ_NRAD * PRJ_NEGROUP];
    double sigma[PRJ_NRAD * PRJ_NEGROUP];
    double delta[PRJ_NRAD * PRJ_NEGROUP];
    double emis[PRJ_NRAD * PRJ_NEGROUP];
    double final_temperature = 0.0;
    double final_ye = 0.0;
    double rho;
    double dt_lapse;
    double e_unchanged;
    double dmom[3] = {0.0, 0.0, 0.0};
    double dE_matter = 0.0;
    double dYe_matter = 0.0;
    int field;
    int group;
    int angle;
    int d;
    int v;
    const double four_pi = 4.0 * M_PI;

    if (rad == 0 || eos == 0 || u == 0) {
        return;
    }
    dt_lapse = dt * lapse;
    if (dt_lapse == 0.0) {
        return;
    }

    /* prj_rad_energy_update reads only the hydro slots and the per-group E slots
     * (prj_rad_implicit_residuals ignores u entirely), and the reconstruction
     * below writes every E/F moment slot.  The angular intensity slots of u_tmp
     * are never read, so copying only the hydro block is bit-identical and skips
     * ~NANGLE-4 dead doubles per group. */
    for (v = 0; v < PRJ_NHYDRO; ++v) {
        u_tmp[v] = u[v];
    }

    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            double E = 0.0;
            double first_moment[3] = {0.0, 0.0, 0.0};

            for (angle = 0; angle < PRJ_NANGLE; ++angle) {
                int iv = PRJ_CONS_RAD_I(field, group, angle);
                double J = u[iv];
                double n[3];

                E += J;
                prj_rad_fsa_rotated_angle_dir(rad, block, angle, ic, jc, kc, n);
                for (d = 0; d < 3; ++d) {
                    first_moment[d] += J * n[d];
                }
            }

            u_tmp[PRJ_CONS_RAD_E(field, group)] = E;
            u_tmp[PRJ_CONS_RAD_F1(field, group)] = PRJ_CLIGHT * first_moment[0];
            u_tmp[PRJ_CONS_RAD_F2(field, group)] = PRJ_CLIGHT * first_moment[1];
            u_tmp[PRJ_CONS_RAD_F3(field, group)] = PRJ_CLIGHT * first_moment[2];
        }
    }

    /* kappa and emis(eta) come back from the converged solve, evaluated at the
     * exact converged (T, Ye); reuse them and look up only the scattering pair at
     * that same converged Ye (not u_tmp[YE]/rho, which round-trips through rho and
     * would perturb kappa/emis relative to sigma/delta). */
    prj_rad_energy_update_fsa(rad, eos, u_tmp, dt, lapse, &final_temperature, &final_ye, kappa, emis);

    rho = u_tmp[PRJ_CONS_RHO];
    prj_rad3_opac_lookup(rad, rho, final_temperature, final_ye,
        0, sigma, delta, 0);

    rho = u[PRJ_CONS_RHO];
    e_unchanged = u[PRJ_CONS_ETOT] - 0.5 *
        (u[PRJ_CONS_MOM1] * u[PRJ_CONS_MOM1] +
         u[PRJ_CONS_MOM2] * u[PRJ_CONS_MOM2] +
         u[PRJ_CONS_MOM3] * u[PRJ_CONS_MOM3]) / rho;

    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            int idx = field * PRJ_NEGROUP + group;
            double J_old[PRJ_NANGLE];
            double J_abs[PRJ_NANGLE];
            double E_abs = 0.0;
            double E_matter_group = 0.0;
            double scatter_rate = sigma[idx] * (1.0 - delta[idx] / 3.0);
            double scatter_den = 1.0 + dt_lapse * PRJ_CLIGHT * scatter_rate;
            double scatter_fac = 1.0 / scatter_den;
            /* den depends only on the group (idx), not the angle: hoist the
             * NANGLE-fold redundant recomputation out of the angle loop. */
            double den = 1.0 + dt_lapse * PRJ_CLIGHT * kappa[idx];

            for (angle = 0; angle < PRJ_NANGLE; ++angle) {
                int iv = PRJ_CONS_RAD_I(field, group, angle);
                double domega = rad->solid_angle[angle];

                J_old[angle] = u[iv];
                J_abs[angle] = (J_old[angle] +
                    dt_lapse * emis[idx] * domega / four_pi) / den;
                E_abs += J_abs[angle];
            }

            for (angle = 0; angle < PRJ_NANGLE; ++angle) {
                int iv = PRJ_CONS_RAD_I(field, group, angle);
                double J_iso = rad->solid_angle[angle] * E_abs / four_pi;
                double J_new = J_iso + (J_abs[angle] - J_iso) * scatter_fac;
                double dE = J_old[angle] - J_new;
                double n[3];

                u[iv] = J_new;
                E_matter_group += dE;
                prj_rad_fsa_rotated_angle_dir(rad, block, angle, ic, jc, kc, n);
                for (d = 0; d < 3; ++d) {
                    dmom[d] += dE * n[d] / PRJ_CLIGHT;
                }
            }

            dE_matter += E_matter_group;
            dYe_matter += rad->x_e[field][group] * E_matter_group;
        }
    }

    e_unchanged += dE_matter * RAD_SCALE;
    u[PRJ_CONS_YE] += dYe_matter;
    u[PRJ_CONS_MOM1] += dmom[0] * RAD_SCALE;
    u[PRJ_CONS_MOM2] += dmom[1] * RAD_SCALE;
    u[PRJ_CONS_MOM3] += dmom[2] * RAD_SCALE;
    u[PRJ_CONS_ETOT] = e_unchanged + 0.5 *
        (u[PRJ_CONS_MOM1] * u[PRJ_CONS_MOM1] +
         u[PRJ_CONS_MOM2] * u[PRJ_CONS_MOM2] +
         u[PRJ_CONS_MOM3] * u[PRJ_CONS_MOM3]) / rho;
}
#endif

#if PRJ_USE_RADIATION_FSA && DO_FFC
/* Fermi constant used by the fast-flavor-conversion rate [MeV cm^3]. */
#define PRJ_FFC_GF 8.958e-44

/* Fast flavor conversion of neutrinos.  Species 0 = nu_e (J1), 1 = nubar_e (J2),
 * 2 = nu_x (J3, one of the four heavy-lepton flavors, hence the /4 factors).
 * All species must share the same energy grid (checked at init).  J in the
 * physics formulae is the physical intensity, i.e. the stored intensity times
 * RAD_SCALE; the mixing and BGK relaxation are linear so they run directly in
 * stored units, while Ip/Im (the growth rate) carry the RAD_SCALE factor. */
void prj_rad_ffc_fsa(prj_rad *rad, double *u, double dt)
{
    /* Ip/Im coefficient: sqrt(2) G_F / (hbar c), with G_F converted to erg cm^3
     * and hbar c = HPLANCK*CLIGHT/(2*pi).  The stored->physical RAD_SCALE for J
     * is folded in here so Ip/Im come out as physical inverse lengths. */
    const double ffc_coeff = sqrt(2.0) * (PRJ_FFC_GF * PRJ_MEV_TO_ERG)
        / (PRJ_HPLANCK * PRJ_CLIGHT / (2.0 * M_PI)) * RAD_SCALE;
    double Ip = 0.0;
    double Im = 0.0;
    double rate;
    double decay;
    int g;
    int angle;

    if (rad == 0 || u == 0 || dt <= 0.0) {
        return;
    }

    /* Ip = sum over (group, angle) of max(J1-J2,0)/erg * coeff, Im likewise with
     * max(J2-J1,0).  erg is the group-center energy (shared by all species). */
    for (g = 0; g < PRJ_NEGROUP; ++g) {
        double inv_erg = 1.0 / rad->egroup_erg[0][g];

        for (angle = 0; angle < PRJ_NANGLE; ++angle) {
            double d = u[PRJ_CONS_RAD_I(0, g, angle)] - u[PRJ_CONS_RAD_I(1, g, angle)];

            if (d > 0.0) {
                Ip += d * inv_erg;
            } else if (d < 0.0) {
                Im += (-d) * inv_erg;
            }
        }
    }
    Ip *= ffc_coeff;
    Im *= ffc_coeff;

    rate = sqrt(Ip * Im);
    decay = exp(-rate * PRJ_CLIGHT * dt);

    for (angle = 0; angle < PRJ_NANGLE; ++angle) {
        double eln = 0.0;
        double P;

        /* Per-angle electron lepton number sets which side of the crossing this
         * angle is on; the survival probability is shared by all energy groups. */
        for (g = 0; g < PRJ_NEGROUP; ++g) {
            eln += (u[PRJ_CONS_RAD_I(0, g, angle)] - u[PRJ_CONS_RAD_I(1, g, angle)])
                / rad->egroup_erg[0][g];
        }

        if (Ip > Im) {
            /* J1>J2 side depleted to 1 - 2/3 Im/Ip; J2>J1 side fully mixed. */
            P = (eln > 0.0) ? (1.0 - (2.0 / 3.0) * Im / Ip) : (1.0 / 3.0);
        } else {
            double ratio = (Im > 0.0) ? (Ip / Im) : 0.0;

            P = (eln > 0.0) ? (1.0 / 3.0) : (1.0 - (2.0 / 3.0) * ratio);
        }

        for (g = 0; g < PRJ_NEGROUP; ++g) {
            int i1 = PRJ_CONS_RAD_I(0, g, angle);
            int i2 = PRJ_CONS_RAD_I(1, g, angle);
            int i3 = PRJ_CONS_RAD_I(2, g, angle);
            double J1 = u[i1];
            double J2 = u[i2];
            double J3 = u[i3];
            double J1a = P * J1 + (1.0 - P) * J3 / 4.0;
            double J2a = P * J2 + (1.0 - P) * J3 / 4.0;
            double J3a = (1.0 - P) * J1 + (1.0 + P) / 4.0 * J3
                + (1.0 - P) * J2 + (1.0 + P) / 4.0 * J3;

            /* BGK relaxation toward the mixed state J_a with the FFC rate. */
            u[i1] = J1a + (J1 - J1a) * decay;
            u[i2] = J2a + (J2 - J2a) * decay;
            u[i3] = J3a + (J3 - J3a) * decay;
        }
    }
}
#endif

/* Koren slope-limiter function φ(r) = max(0, min(2r, (2+r)/3, 2)). */
static double prj_rad_koren_phi(double r)
{
    double phi = 2.0 * r;
    double t = (2.0 + r) / 3.0;

    if (t < phi) {
        phi = t;
    }
    if (2.0 < phi) {
        phi = 2.0;
    }
    if (phi < 0.0) {
        phi = 0.0;
    }
    return phi;
}

/* Reconstruct cell-value array q[] (one entry per energy group) at the right
 * (side=+1) or left (side=-1) face of cell `gcell` using a Koren-limited linear
 * stencil.  Energy groups are uniformly spaced in log ν, so equal-spaced
 * samples are valid.  The two outermost cells (gcell == 0 or NEGROUP-1) fall
 * back to piecewise constant — at those edges the outer face values are also
 * forced to zero by the caller, so the choice has no effect on the update. */
static double prj_rad_recon_face(const double q[PRJ_NEGROUP], int gcell, int side)
{
    (void)side;
    return q[gcell];
}

#if PRJ_DYNAMIC_GR && PRJ_USE_RADIATION_M1
static double prj_rad_gr_m1_cell_opacity(const prj_block *block,
    int i, int j, int k, int field, int group)
{
    const size_t stride = (size_t)PRJ_NRAD * (size_t)PRJ_NEGROUP;
    const int op_idx = field * PRJ_NEGROUP + group;
    const size_t op_off = (size_t)IDX(i, j, k) * stride + (size_t)op_idx;
    double opacity = 0.0;

    if (block != 0 && block->kappa_cell != 0 && block->sigma_cell != 0) {
        double kappa = block->kappa_cell[op_off];
        double sigma = block->sigma_cell[op_off];

        if (isfinite(kappa) && kappa > 0.0) {
            opacity += kappa;
        }
        if (isfinite(sigma) && sigma > 0.0) {
            opacity += sigma;
        }
    }
    return opacity;
}

static void prj_rad_gr_m1_fill_closure_ctx(const prj_z4c_hydro_geom *geom,
    const double vcon[3], const double dvcon_dx[3][3], int have_shear,
    double opacity, prj_rad_gr_m1_closure_ctx *ctx)
{
    int a;
    int b;
    int d;

    memset(ctx, 0, sizeof(*ctx));
    for (a = 0; a < 3; ++a) {
        ctx->vcon[a] = vcon != 0 ? vcon[a] : 0.0;
        for (b = 0; b < 3; ++b) {
            ctx->gamma[a][b] = geom->gamma[a][b];
            ctx->gamma_inv[a][b] = geom->gamma_inv[a][b];
            ctx->K_dd[a][b] = geom->K_dd[a][b];
            for (d = 0; d < 3; ++d) {
                ctx->dgamma[d][a][b] = geom->dgamma[d][a][b];
            }
        }
    }
    if (dvcon_dx != 0) {
        for (d = 0; d < 3; ++d) {
            for (a = 0; a < 3; ++a) {
                ctx->dvdx[d][a] = dvcon_dx[d][a];
            }
        }
    }
    ctx->opacity = opacity;
    ctx->have_shear = have_shear;
}

static void prj_rad_gr_m1_cell_dvcon_dx(const prj_block *block,
    int ic, int jc, int kc, double dvcon_dx[3][3])
{
    double inv_dx[3];
    int jdir;
    int icomp;

    memset(dvcon_dx, 0, 9 * sizeof(double));
    if (block == 0 || block->v_riemann[0] == 0 ||
        block->v_riemann[1] == 0 || block->v_riemann[2] == 0) {
        return;
    }
    inv_dx[0] = 1.0 / block->dx[0];
    inv_dx[1] = 1.0 / block->dx[1];
    inv_dx[2] = 1.0 / block->dx[2];
    for (jdir = 0; jdir < 3; ++jdir) {
        for (icomp = 0; icomp < 3; ++icomp) {
            int il = ic;
            int jl = jc;
            int kl = kc;
            int ir = ic;
            int jr = jc;
            int kr = kc;
            double vL;
            double vR;

            if (jdir == X1DIR) {
                ir = ic + 1;
            } else if (jdir == X2DIR) {
                jr = jc + 1;
            } else {
                kr = kc + 1;
            }
            vL = block->v_riemann[jdir][VRIDX(icomp, il, jl, kl)];
            vR = block->v_riemann[jdir][VRIDX(icomp, ir, jr, kr)];
            dvcon_dx[jdir][icomp] =
                (vR - vL) * inv_dx[jdir] / PRJ_CLIGHT;
        }
    }
}

static void prj_rad_gr_m1_freq_base_closure_ctx(const prj_z4c_hydro_geom *geom,
    const double *W_mhd, const double dvcon_dx[3][3], int i, int j, int k,
    prj_rad_gr_m1_closure_ctx *ctx)
{
    double vcon[3];
    int a;

    for (a = 0; a < 3; ++a) {
        vcon[a] = W_mhd != 0 ? W_mhd[WIDX(PRJ_PRIM_V1 + a, i, j, k)] /
            PRJ_CLIGHT : 0.0;
    }
    prj_rad_gr_m1_fill_closure_ctx(geom, vcon, dvcon_dx, 1, 0.0, ctx);
}

void prj_rad_freq_flux_apply_gr_m1(const prj_rad *rad, const prj_mesh *mesh,
    const prj_block *block, int z4c_stage, const double *W_state, double *u,
    int ic, int jc, int kc, double dt)
{
    double dvcon_dx[3][3];
    prj_z4c_hydro_geom geom;
    prj_rad_gr_m1_closure_ctx base_closure;
    prj_rad_gr_m1_side_data pside;
    const double *W_mhd;
    double dt_geom;
    int field;

    if (rad == 0 || mesh == 0 || block == 0 || W_state == 0 || u == 0) {
        return;
    }
    if (block->v_riemann[0] == 0 || block->v_riemann[1] == 0 ||
        block->v_riemann[2] == 0) {
        return;
    }
    if (!prj_z4c_load_hydro_geom(mesh, block, z4c_stage, ic, jc, kc, &geom)) {
        fprintf(stderr,
            "prj_rad_freq_flux_apply_gr_m1: failed to load Z4c geometry at cell (%d,%d,%d)\n",
            ic, jc, kc);
        exit(1);
    }
    W_mhd = prj_block_mhd_stage_const(block, z4c_stage);
    dt_geom = dt * geom.alpha * geom.sqrt_gamma;
    if (dt_geom == 0.0) {
        return;
    }

    prj_rad_gr_m1_cell_dvcon_dx(block, ic, jc, kc, dvcon_dx);
    prj_rad_gr_m1_freq_base_closure_ctx(&geom, W_mhd, dvcon_dx, ic, jc, kc,
        &base_closure);
    prj_rad_gr_m1_prepare_side(&base_closure, &pside);

    for (field = 0; field < PRJ_NRAD; ++field) {
        double Eg[PRJ_NEGROUP];
        double Pg[PRJ_NEGROUP][3][3];
        double Mq_cov[PRJ_NEGROUP][3];
        double Acon_spec[PRJ_NEGROUP];
        double Mq_spec[PRJ_NEGROUP][3];
        double inv_dnu[PRJ_NEGROUP];
        double energy_face[PRJ_NEGROUP + 1] = {0.0};
        double momentum_face[PRJ_NEGROUP + 1][PRJ_NDIM] = {{0.0}};
        double energy_available[PRJ_NEGROUP];
        const double *nu_face = rad->eedge[field];
        int g;
        int gf;
        int a;

        if (nu_face == 0) {
            fprintf(stderr, "prj_rad_freq_flux_apply_gr_m1: missing eedge for field %d\n",
                field);
            exit(1);
        }

        for (g = 0; g < PRJ_NEGROUP; ++g) {
            double dnu = nu_face[g + 1] - nu_face[g];
            double Fcov[3];
            prj_rad_gr_m1_closure_ctx closure;
            prj_rad_gr_m1_pressure_data pdata;
            double fbar;
            double Acon = 0.0;

            if (dnu <= 0.0) {
                fprintf(stderr,
                    "prj_rad_freq_flux_apply_gr_m1: non-positive eedge width for field %d group %d\n",
                    field, g);
                exit(1);
            }
            inv_dnu[g] = 1.0 / dnu;
            Eg[g] = W_state[WIDX(PRJ_RAD_PRIM_E(field, g), ic, jc, kc)];
            Fcov[0] = W_state[WIDX(PRJ_RAD_PRIM_F1(field, g), ic, jc, kc)];
            Fcov[1] = W_state[WIDX(PRJ_RAD_PRIM_F2(field, g), ic, jc, kc)];
            Fcov[2] = W_state[WIDX(PRJ_RAD_PRIM_F3(field, g), ic, jc, kc)];

            closure = base_closure;
            closure.opacity = prj_rad_gr_m1_cell_opacity(block, ic, jc, kc,
                field, g);
            prj_rad_gr_m1_prepare_pressure(rad, &closure, &pside, Eg[g], Fcov,
                &pdata);
            /* frequency_drifts (via decompose_m2) needs the full pdata, so
             * prepare_pressure must run; but the geometric source term already
             * solved the identical closure this stage, so reuse its (P, fbar)
             * and skip the expensive implicit root-find on a cache hit. */
            if (!prj_rad_gr_m1_closure_cache_get(block, W_state, z4c_stage,
                    (ic * PRJ_BLOCK_SIZE + jc) * PRJ_BLOCK_SIZE + kc,
                    field * PRJ_NEGROUP + g, Pg[g], &fbar)) {
                prj_rad_gr_m1_pressure_implicit(&pdata, Pg[g], &fbar);
            }
            prj_rad_gr_m1_frequency_drifts(&pdata, Pg[g], fbar, &Acon,
                Mq_cov[g]);
            Acon_spec[g] = Acon * inv_dnu[g];
            for (a = 0; a < 3; ++a) {
                Mq_spec[g][a] = Mq_cov[g][a] * inv_dnu[g];
            }
            energy_available[g] = u[PRJ_CONS_RAD_E(field, g)];
        }

        for (gf = 1; gf < PRJ_NEGROUP; ++gf) {
            double face_drift = Acon_spec[gf - 1] + Acon_spec[gf];
            int gu = face_drift >= 0.0 ? gf : gf - 1;
            double nu = nu_face[gf];

            energy_face[gf] += nu * Acon_spec[gu];
            for (a = 0; a < 3; ++a) {
                momentum_face[gf][a] += nu * Mq_spec[gu][a];
            }
        }

        {
            double outgoing[PRJ_NEGROUP] = {0.0};
            double theta[PRJ_NEGROUP];

            for (gf = 1; gf < PRJ_NEGROUP; ++gf) {
                if (energy_face[gf] > 0.0) {
                    outgoing[gf] += energy_face[gf];
                } else if (energy_face[gf] < 0.0) {
                    outgoing[gf - 1] -= energy_face[gf];
                }
            }
            for (g = 0; g < PRJ_NEGROUP; ++g) {
                double drain = dt_geom * outgoing[g];

                theta[g] = 1.0;
                if (drain > energy_available[g]) {
                    theta[g] = energy_available[g] > 0.0 && drain > 0.0
                        ? nextafter(energy_available[g] / drain, 0.0)
                        : 0.0;
                }
            }
            for (gf = 1; gf < PRJ_NEGROUP; ++gf) {
                int donor = energy_face[gf] > 0.0 ? gf : gf - 1;
                double factor = theta[donor];

                energy_face[gf] *= factor;
                for (a = 0; a < PRJ_NDIM; ++a) {
                    momentum_face[gf][a] *= factor;
                }
            }
        }

        for (g = 0; g < PRJ_NEGROUP; ++g) {
            u[PRJ_CONS_RAD_E(field, g)] += dt_geom *
                (energy_face[g + 1] - energy_face[g]);
            u[PRJ_CONS_RAD_F1(field, g)] += dt_geom *
                (momentum_face[g + 1][0] - momentum_face[g][0]);
            u[PRJ_CONS_RAD_F2(field, g)] += dt_geom *
                (momentum_face[g + 1][1] - momentum_face[g][1]);
            u[PRJ_CONS_RAD_F3(field, g)] += dt_geom *
                (momentum_face[g + 1][2] - momentum_face[g][2]);
        }
    }
}

static void prj_rad_gr_m1_clamp_ef(const prj_z4c_hydro_geom *geom,
    double *E, double Fcov[3])
{
    double Fcon[3];
    double F2 = 0.0;
    double Fmag;
    double cE;
    int a;
    int b;

    if (E == 0 || Fcov == 0) {
        return;
    }
    if (!isfinite(*E) || *E <= 0.0) {
        *E = 0.0;
        Fcov[0] = 0.0;
        Fcov[1] = 0.0;
        Fcov[2] = 0.0;
        return;
    }
    for (a = 0; a < 3; ++a) {
        if (!isfinite(Fcov[a])) {
            Fcov[a] = 0.0;
        }
        Fcon[a] = 0.0;
        for (b = 0; b < 3; ++b) {
            Fcon[a] += geom->gamma_inv[a][b] * Fcov[b];
        }
        F2 += Fcov[a] * Fcon[a];
    }
    if (!isfinite(F2) || F2 <= 0.0) {
        return;
    }
    Fmag = sqrt(F2);
    cE = PRJ_CLIGHT * *E;
    if (Fmag > cE && cE >= 0.0) {
        double scale = cE / Fmag;

        for (a = 0; a < 3; ++a) {
            Fcov[a] *= scale;
        }
    }
}

static void prj_rad_gr_m1_clamp_conserved_cell(
    const prj_z4c_hydro_geom *geom, double *u)
{
    int field;
    int group;

    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            double E = u[PRJ_CONS_RAD_E(field, group)] / geom->sqrt_gamma;
            double Fcov[3];

            Fcov[0] = u[PRJ_CONS_RAD_F1(field, group)] / geom->sqrt_gamma;
            Fcov[1] = u[PRJ_CONS_RAD_F2(field, group)] / geom->sqrt_gamma;
            Fcov[2] = u[PRJ_CONS_RAD_F3(field, group)] / geom->sqrt_gamma;
            prj_rad_gr_m1_clamp_ef(geom, &E, Fcov);
            u[PRJ_CONS_RAD_E(field, group)] = geom->sqrt_gamma * E;
            u[PRJ_CONS_RAD_F1(field, group)] = geom->sqrt_gamma * Fcov[0];
            u[PRJ_CONS_RAD_F2(field, group)] = geom->sqrt_gamma * Fcov[1];
            u[PRJ_CONS_RAD_F3(field, group)] = geom->sqrt_gamma * Fcov[2];
        }
    }
}

static void prj_rad_gr_m1_source_base_closure_ctx(const prj_z4c_hydro_geom *geom,
    const double *W, const double dvcon_dx[3][3],
    prj_rad_gr_m1_closure_ctx *ctx)
{
    double vcon[3];
    int a;

    for (a = 0; a < 3; ++a) {
        vcon[a] = W != 0 ? W[PRJ_PRIM_V1 + a] / PRJ_CLIGHT : 0.0;
    }
    prj_rad_gr_m1_fill_closure_ctx(geom, vcon, dvcon_dx, 1, 0.0, ctx);
}

void prj_rad_gr_m1_matter_update(prj_rad *rad, prj_eos *eos,
    const prj_mesh *mesh, const prj_block *block, int z4c_stage, double *u,
    int i, int j, int k, double dt, double *final_temperature)
{
    prj_z4c_hydro_geom geom;
    double W[PRJ_NVAR_PRIM];
    double eos_q[PRJ_EOS_NQUANT];
    double tmp[PRJ_NVAR_CONS];
    double J_old[PRJ_NRAD * PRJ_NEGROUP];
    double H_old[PRJ_NRAD * PRJ_NEGROUP][4];
    double ucon_src[PRJ_NRAD * PRJ_NEGROUP][4];
    double kappa[PRJ_NRAD * PRJ_NEGROUP];
    double sigma[PRJ_NRAD * PRJ_NEGROUP];
    double delta[PRJ_NRAD * PRJ_NEGROUP];
    double dvcon_dx[3][3];
    prj_rad_gr_m1_closure_ctx base_closure;
    prj_rad_gr_m1_side_data pside;
    double rho;
    double eint;
    double T;
    double Ye_final;
    double dt_lapse;
    double total_dE = 0.0;
    double total_dF[3] = {0.0, 0.0, 0.0};
    int field;
    int group;
    int v;

    if (rad == 0 || eos == 0 || mesh == 0 || block == 0 || u == 0 ||
        dt == 0.0 || !prj_eos_full_dynamic_gr_enabled(mesh)) {
        return;
    }
    if (!prj_z4c_load_hydro_geom(mesh, block, z4c_stage, i, j, k, &geom)) {
        fprintf(stderr,
            "prj_rad_gr_m1_matter_update: failed to load Z4c geometry at cell (%d,%d,%d)\n",
            i, j, k);
        exit(1);
    }
    dt_lapse = dt * geom.alpha;
    if (dt_lapse == 0.0) {
        return;
    }

    prj_rad_gr_m1_clamp_conserved_cell(&geom, u);
    prj_eos_cell_cons2prim(eos, mesh, block, z4c_stage, i, j, k, u, W,
        PRJ_EOS_CTX_MAIN);
    rho = W[PRJ_PRIM_RHO];
    eint = W[PRJ_PRIM_EINT];
    prj_eos_rey(eos, rho, eint, W[PRJ_PRIM_YE], eos_q, PRJ_EOS_CTX_MAIN);
    T = eos_q[PRJ_EOS_TEMPERATURE];
    prj_rad_gr_m1_cell_dvcon_dx(block, i, j, k, dvcon_dx);
    prj_rad_gr_m1_source_base_closure_ctx(&geom, W, dvcon_dx, &base_closure);
    prj_rad_gr_m1_prepare_side(&base_closure, &pside);

    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        tmp[v] = 0.0;
    }
    tmp[PRJ_CONS_RHO] = rho;
    tmp[PRJ_CONS_ETOT] = rho * eint;
    tmp[PRJ_CONS_YE] = rho * W[PRJ_PRIM_YE];

    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            int idx = field * PRJ_NEGROUP + group;
            double E;
            double Fcov[3];
            double P[3][3];
            double J;
            double H[4];
            double L[4][4];
            double ucon[4];
            double ucov[4];
            double hcon[4][4];
            double hcov[4][4];
            prj_rad_gr_m1_closure_ctx closure;
            prj_rad_gr_m1_pressure_data pdata;

            E = u[PRJ_CONS_RAD_E(field, group)] / geom.sqrt_gamma;
            Fcov[0] = u[PRJ_CONS_RAD_F1(field, group)] / geom.sqrt_gamma;
            Fcov[1] = u[PRJ_CONS_RAD_F2(field, group)] / geom.sqrt_gamma;
            Fcov[2] = u[PRJ_CONS_RAD_F3(field, group)] / geom.sqrt_gamma;
            prj_rad_gr_m1_clamp_ef(&geom, &E, Fcov);

            closure = base_closure;
            closure.opacity = prj_rad_gr_m1_cell_opacity(block, i, j, k,
                field, group);
            prj_rad_gr_m1_prepare_pressure(rad, &closure, &pside, E, Fcov,
                &pdata);
            prj_rad_gr_m1_pressure_implicit(&pdata, P, 0);
            prj_rad_gr_m1_decompose_m2(&pdata, P, &J, H, L, ucon, ucov,
                hcon, hcov);

            J_old[idx] = J;
            for (v = 0; v < 4; ++v) {
                H_old[idx][v] = H[v];
                ucon_src[idx][v] = ucon[v];
            }
            tmp[PRJ_CONS_RAD_E(field, group)] = J;
            /* The GR closure uses M^{0i}=F^i/c, so comoving H has energy-density
             * units.  The temporary non-GR M1 source state expects the evolved
             * flux variable, hence cH here.  The inelastic M1 routines still
             * leave these flux slots unchanged, matching the non-GR behavior. */
            tmp[PRJ_CONS_RAD_F1(field, group)] = PRJ_CLIGHT * H[1];
            tmp[PRJ_CONS_RAD_F2(field, group)] = PRJ_CLIGHT * H[2];
            tmp[PRJ_CONS_RAD_F3(field, group)] = PRJ_CLIGHT * H[3];
        }
    }

    prj_rad_eleinel_step(rad, eos, tmp, dt_lapse, T);
    prj_rad_nucinel_step(rad, eos, tmp, dt_lapse, T);
    prj_rad_energy_update(rad, eos, tmp, dt_lapse, 1.0, &T, kappa);
    Ye_final = tmp[PRJ_CONS_RHO] > 0.0 ?
        tmp[PRJ_CONS_YE] / tmp[PRJ_CONS_RHO] : W[PRJ_PRIM_YE];
    prj_rad3_opac_lookup(rad, rho, T, Ye_final, 0, sigma, delta, 0);

    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            int idx = field * PRJ_NEGROUP + group;
            double chi = kappa[idx] + sigma[idx] * (1.0 - delta[idx] / 3.0);
            double hscale = 1.0 / (1.0 + dt_lapse * PRJ_CLIGHT * chi);
            double dJ = tmp[PRJ_CONS_RAD_E(field, group)] - J_old[idx];
            double dH[4];
            double Scon[4];
            double S_cov_spatial[3];
            double dE_cons;
            double dF_cons[3];
            int a;
            int b;

            for (a = 0; a < 4; ++a) {
                double H_new = H_old[idx][a] * hscale;

                dH[a] = H_new - H_old[idx][a];
                Scon[a] = dJ * ucon_src[idx][a] + dH[a];
            }
            dE_cons = geom.sqrt_gamma * Scon[0];
            for (a = 0; a < 3; ++a) {
                S_cov_spatial[a] = 0.0;
                for (b = 0; b < 3; ++b) {
                    S_cov_spatial[a] += geom.gamma[a][b] * Scon[b + 1];
                }
                dF_cons[a] = geom.sqrt_gamma * PRJ_CLIGHT * S_cov_spatial[a];
            }

            u[PRJ_CONS_RAD_E(field, group)] += dE_cons;
            u[PRJ_CONS_RAD_F1(field, group)] += dF_cons[0];
            u[PRJ_CONS_RAD_F2(field, group)] += dF_cons[1];
            u[PRJ_CONS_RAD_F3(field, group)] += dF_cons[2];
            total_dE += dE_cons;
            total_dF[0] += dF_cons[0];
            total_dF[1] += dF_cons[1];
            total_dF[2] += dF_cons[2];
        }
    }

    prj_rad_gr_m1_clamp_conserved_cell(&geom, u);
    u[PRJ_CONS_ETOT] -= total_dE * RAD_SCALE;
    u[PRJ_CONS_MOM1] -= total_dF[0] * RAD_SCALE /
        (PRJ_CLIGHT * PRJ_CLIGHT);
    u[PRJ_CONS_MOM2] -= total_dF[1] * RAD_SCALE /
        (PRJ_CLIGHT * PRJ_CLIGHT);
    u[PRJ_CONS_MOM3] -= total_dF[2] * RAD_SCALE /
        (PRJ_CLIGHT * PRJ_CLIGHT);
    u[PRJ_CONS_YE] = u[PRJ_CONS_RHO] * Ye_final;

    /* This source solve freezes the frame velocity while coupling radiation
     * and matter.  A future fully implicit GR-M1 source solve should include
     * the velocity in the nonlinear system for strongly momentum-coupled cells. */
    if (final_temperature != 0) {
        *final_temperature = T;
    }
}
#endif

#if PRJ_USE_RADIATION_FSA
static void prj_rad_fsa_omega_faces(const prj_rad *rad, int field,
    double omega_face[PRJ_NEGROUP + 1])
{
    int gf;

    if (rad->eedge[field] != 0) {
        for (gf = 0; gf <= PRJ_NEGROUP; ++gf) {
            omega_face[gf] = rad->eedge[field][gf];
        }
        return;
    }
    if (rad->emin[field] <= 0.0 || rad->emax[field] <= rad->emin[field]) {
        fprintf(stderr, "prj_rad_freq_flux_apply: invalid FSA emin/emax for field %d\n",
            field);
        exit(1);
    }
    {
        double log_min = log(rad->emin[field]);
        double log_max = log(rad->emax[field]);
        double dlog = (log_max - log_min) / (double)PRJ_NEGROUP;

        for (gf = 0; gf <= PRJ_NEGROUP; ++gf) {
            omega_face[gf] = exp(log_min + (double)gf * dlog);
        }
        omega_face[0] = rad->emin[field];
        omega_face[PRJ_NEGROUP] = rad->emax[field];
    }
}
#endif

/* Apply the per-cell energy-space flux terms.  The M1 branch applies the
 * energy-space-flux part of the SR redshift terms (Eqs. 21a/21b of the
 * comoving-frame mixed-frame moment equations):
 *
 *   ∂_t E_g    -= - v^i_{;j} [ (ν P^j_{νi})_{g+1/2} - (ν P^j_{νi})_{g-1/2} ]
 *   ∂_t F_{gj} -= - v^i_{;k} [ (ν Q^k_{νji})_{g+1/2} - (ν Q^k_{νji})_{g-1/2} ]
 *
 * (Equivalently dE_g/dt += v^i_{;j} ΔνP^{ji} and dF_{gj}/dt += v^i_{;k} ΔνQ^{kji}.)
 *
 * The FSA branch applies the angular-cell-integrated form of the third LHS term
 * of the intensity equation after moving it to the update side:
 *
 *   ∂_t J_g,a += α A_n ∂_ω(ω J_g,a),  A_n = n · ∇'v · n + a · n / c.
 *
 * The face values are picked by upwinding in frequency space.  With the update
 * written as dU_g/dt = face[g+1] - face[g], a positive drift drains the upper
 * group and uses that group as the donor; a negative drift drains the lower
 * group.
 *
 * M1 closure choices:
 *   - P^{ij}_g from the standard M1 Eddington tensor on cell-centred (E_g, F_g).
 *   - Q^{kij}_g from the Levermore/Vaytet third-moment closure
 *       Q^{ijk} = c E H^{ijk},
 *     with H^{ijk} built from f = |F|/(cE), n = F/|F|, and q(f).
 *
 * Reconstruction: linear-in-log-ν with Koren limiter (groups are uniform in
 * log ν → equal spacing).  Outermost group faces (g = -1/2 and g = NEGROUP-1/2)
 * are set to zero (outflow).  The cell-centred state used for the closure comes
 * from W_state (primitive stage 0 in stage1, primitive stage 1 in stage2), per
 * the user-specified ordering. */
void prj_rad_freq_flux_apply(const prj_rad *rad, const prj_block *block,
    const double *W_state, double *u, int ic, int jc, int kc, double lapse, double dt)
{
#if PRJ_USE_RADIATION_M1
    double dvdx[3][3];
    double inv_dx[3];
    int jdir;
    int icomp;
    int field;

    if (rad == 0 || block == 0 || W_state == 0 || u == 0) {
        return;
    }
    if (block->v_riemann[0] == 0 || block->v_riemann[1] == 0 || block->v_riemann[2] == 0) {
        return;
    }

    inv_dx[0] = 1.0 / block->dx[0];
    inv_dx[1] = 1.0 / block->dx[1];
    inv_dx[2] = 1.0 / block->dx[2];

    /* Cell-centred ∂_jdir v_icomp from the two normal-direction Riemann faces. */
    for (jdir = 0; jdir < 3; ++jdir) {
        for (icomp = 0; icomp < 3; ++icomp) {
            int il = ic;
            int jl = jc;
            int kl = kc;
            int ir = ic;
            int jr = jc;
            int kr = kc;
            double vL;
            double vR;

            if (jdir == X1DIR) {
                ir = ic + 1;
            } else if (jdir == X2DIR) {
                jr = jc + 1;
            } else {
                kr = kc + 1;
            }
            vL = block->v_riemann[jdir][VRIDX(icomp, il, jl, kl)];
            vR = block->v_riemann[jdir][VRIDX(icomp, ir, jr, kr)];
            dvdx[jdir][icomp] = (vR - vL) * inv_dx[jdir];
        }
    }

    int cell_idx = IDX(ic, jc, kc);
    double grad_phi[3];
    double inv_c2 = 1.0 / (PRJ_CLIGHT * PRJ_CLIGHT);

    grad_phi[0] = 0.0;
    grad_phi[1] = 0.0;
    grad_phi[2] = 0.0;
    if (block->grav[0] != 0 && block->grav[1] != 0 && block->grav[2] != 0) {
        grad_phi[0] = -block->grav[0][cell_idx];
        grad_phi[1] = -block->grav[1][cell_idx];
        grad_phi[2] = -block->grav[2][cell_idx];
    }

    for (field = 0; field < PRJ_NRAD; ++field) {
        double dt_lapse = lapse * dt;
        double Eg[PRJ_NEGROUP];
        double Fg[PRJ_NEGROUP][3];
        double Pg[PRJ_NEGROUP][3][3];
        double Mq[PRJ_NEGROUP][3]; /* Q_g : dvdx, the only way Q is ever used */
        double Acon[PRJ_NEGROUP];  /* P_g : dvdx, the energy-space drift per group */
        double inv_dnu[PRJ_NEGROUP];
        double Mq_spec[PRJ_NEGROUP][3];
        double Acon_spec[PRJ_NEGROUP];
        double energy_face[PRJ_NEGROUP + 1] = {0.0};
        double momentum_face[PRJ_NEGROUP + 1][PRJ_NDIM] = {{0.0}};
        double energy_available[PRJ_NEGROUP];
        const double *nu_face = rad->eedge[field];
        int g;
        int ii;
        int jj;

        if (nu_face == 0) {
            fprintf(stderr, "prj_rad_freq_flux_apply: missing eedge for field %d\n",
                field);
            exit(1);
        }

        /* Per-group cell-centred state and closure tensors.  P and Q are built
         * once here and shared: P by the GR redshift terms below, and both P, Q
         * by the SR frequency flux (reconstructed to the frequency faces). */
        for (g = 0; g < PRJ_NEGROUP; ++g) {
            double dnu = nu_face[g + 1] - nu_face[g];

            if (dnu <= 0.0) {
                fprintf(stderr,
                    "prj_rad_freq_flux_apply: non-positive eedge width for field %d group %d\n",
                    field, g);
                exit(1);
            }
            inv_dnu[g] = 1.0 / dnu;

            Eg[g] = W_state[WIDX(PRJ_RAD_PRIM_E(field, g), ic, jc, kc)];
            Fg[g][0] = W_state[WIDX(PRJ_RAD_PRIM_F1(field, g), ic, jc, kc)];
            Fg[g][1] = W_state[WIDX(PRJ_RAD_PRIM_F2(field, g), ic, jc, kc)];
            Fg[g][2] = W_state[WIDX(PRJ_RAD_PRIM_F3(field, g), ic, jc, kc)];
            prj_rad_m1_pressure(rad, Eg[g], Fg[g][0], Fg[g][1], Fg[g][2], Pg[g]);
            prj_rad_m1_third_moment_contract(rad, Eg[g], Fg[g][0], Fg[g][1], Fg[g][2],
                dvdx, Mq[g]);
            Acon[g] = 0.0;
            for (jj = 0; jj < 3; ++jj) {
                for (ii = 0; ii < 3; ++ii) {
                    Acon[g] += Pg[g][jj][ii] * dvdx[jj][ii];
                }
            }
            Acon_spec[g] = Acon[g] * inv_dnu[g];
            for (jj = 0; jj < 3; ++jj) {
                Mq_spec[g][jj] = Mq[g][jj] * inv_dnu[g];
            }
            energy_available[g] = u[PRJ_CONS_RAD_E(field, g)];
        }

        /* SR velocity-gradient energy-space flux (Eqs. 21a/21b).  Each interior
         * frequency face is upwinded by the sign of the actual energy-space
         * drift P:dvdx (the contraction that the flux itself transports),
         * estimated from the two groups adjacent to the face (left+right).  This
         * reduces exactly to the bulk div(v) rule for isotropic radiation
         * (P = E/3·I  =>  P:dvdx = E/3·div v) but stays upwind-correct when the
         * anisotropic/shearing part of ∂_j v_i dominates the trace, where a
         * single trace-based side could disagree with the true drift.  The
         * stored M1 moments are group-integrated.  rad->eedge[field] stores
         * energy-bin edges in MeV, so nu_face[gf] and dnu are both in MeV; the
         * nu_face/dnu factor is dimensionless and needs no PRJ_MEV_TO_ERG.
         * Thus the energy flux at face gf is
         * ν_face[gf]·(P_gu:dvdx)/dν_gu = ν_face[gf]·Acon_spec[gu],
         * and the momentum flux is
         * ν_face[gf]·(Q_gu:dvdx)/dν_gu = ν_face[gf]·Mq_spec[gu].  This is the
         * upper (g+1/2) face of group gf-1 and the lower (g-1/2) face of group
         * gf, so it is scattered into both with opposite signs. */
        {
            int gf;

            for (gf = 1; gf < PRJ_NEGROUP; ++gf) {
                double face_drift = Acon_spec[gf - 1] + Acon_spec[gf];
                int gu = (face_drift >= 0.0) ? gf : gf - 1;
                double nu = nu_face[gf];

                energy_face[gf] += nu * Acon_spec[gu];
                for (jj = 0; jj < 3; ++jj) {
                    momentum_face[gf][jj] += nu * Mq_spec[gu][jj];
                }
            }
        }

        /* GR ε-flux pieces of G^e and G^m_j (the ∂_ε terms only):
         *   G^e          = -F_{sε}·∇φ/c² + (∇φ/c²) · ∂_ε(ε F_{sε})
         *   G^m_j        = -E_{sε} ∇_jφ + ∇_iφ · ∂_ε(ε P^i_{sεj})
         * G is on the source-side (RHS) of the moment equations, so when added
         * to ∂_t E_g / ∂_t F_{gj} it carries an overall minus sign:
         *   ∂_t E_g    -= (∇_i φ / c²) · [(εF^i)_{g+1/2} - (εF^i)_{g-1/2}]
         *   ∂_t F_{gj} -= (∇_i φ)      · [(εP^{ij})_{g+1/2} - (εP^{ij})_{g-1/2}]
         * (sums over i; the −F_{sε}·∇φ/c² and −E_{sε}∇_jφ pieces are NOT done
         * here per user request).
         *
         * Upwind in ε-space follows the same Eq. 22 rule, using the per-i
         * coefficient that multiplies the ε-flux divergence as the "speed":
         *   pick L if coef <  0,  pick R if coef >= 0.
         *
         * ∇_i φ at the cell centre comes from the active monopole gravity:
         * gravitational acceleration a_i = accel(r) · x_i/r, and ∇φ = −a. */
        {
            /* Energy: per i, scalar q[g] = F^i_g, coef = −(∇_i φ)/c². */
            for (ii = 0; ii < 3; ++ii) {
                double coef = -grad_phi[ii] * inv_c2;
                double q[PRJ_NEGROUP];
                double face_val[PRJ_NEGROUP + 1];
                int gf;

                if (coef == 0.0) {
                    continue;
                }
                for (g = 0; g < PRJ_NEGROUP; ++g) {
                    q[g] = Fg[g][ii] * inv_dnu[g];
                }
                face_val[0] = 0.0;
                face_val[PRJ_NEGROUP] = 0.0;
                for (gf = 1; gf < PRJ_NEGROUP; ++gf) {
                    double pick = (coef >= 0.0)
                        ? prj_rad_recon_face(q, gf, -1)
                        : prj_rad_recon_face(q, gf - 1, +1);
                    face_val[gf] = nu_face[gf] * pick;
                }
                for (gf = 1; gf < PRJ_NEGROUP; ++gf) {
                    energy_face[gf] += coef * face_val[gf];
                }
            }

            /* Flux j: per (i, j), scalar q[g] = P^{ij}_g, coef = −∇_i φ. */
            for (jj = 0; jj < 3; ++jj) {
                for (ii = 0; ii < 3; ++ii) {
                    double coef = -grad_phi[ii];
                    double q[PRJ_NEGROUP];
                    double face_val[PRJ_NEGROUP + 1];
                    int gf;

                    if (coef == 0.0) {
                        continue;
                    }
                    for (g = 0; g < PRJ_NEGROUP; ++g) {
                        q[g] = Pg[g][ii][jj] * inv_dnu[g];
                    }
                    face_val[0] = 0.0;
                    face_val[PRJ_NEGROUP] = 0.0;
                    for (gf = 1; gf < PRJ_NEGROUP; ++gf) {
                        double pick = (coef >= 0.0)
                            ? prj_rad_recon_face(q, gf, -1)
                            : prj_rad_recon_face(q, gf - 1, +1);
                        face_val[gf] = nu_face[gf] * pick;
                    }
                    for (gf = 1; gf < PRJ_NEGROUP; ++gf) {
                        momentum_face[gf][jj] += coef * face_val[gf];
                    }
                }
            }
        }

        /* Limit the combined SR+GR frequency-space flux with one factor per
         * donor group.  For dE_g/dt = face[g+1] - face[g], a positive face
         * drains its upper group and a negative face drains its lower group.
         * Both outgoing faces of a group share the same factor. */
        {
            double outgoing[PRJ_NEGROUP] = {0.0};
            double theta[PRJ_NEGROUP];
            int gf;

            for (gf = 1; gf < PRJ_NEGROUP; ++gf) {
                if (energy_face[gf] > 0.0) {
                    outgoing[gf] += energy_face[gf];
                } else if (energy_face[gf] < 0.0) {
                    outgoing[gf - 1] -= energy_face[gf];
                }
            }
            for (g = 0; g < PRJ_NEGROUP; ++g) {
                double drain = dt_lapse * outgoing[g];

                theta[g] = 1.0;
                if (drain > energy_available[g]) {
                    theta[g] = energy_available[g] > 0.0 && drain > 0.0
                        ? nextafter(energy_available[g] / drain, 0.0)
                        : 0.0;
                }
            }
            for (gf = 1; gf < PRJ_NEGROUP; ++gf) {
                int donor = energy_face[gf] > 0.0 ? gf : gf - 1;
                double factor = theta[donor];

                energy_face[gf] *= factor;
                for (ii = 0; ii < PRJ_NDIM; ++ii) {
                    momentum_face[gf][ii] *= factor;
                }
            }
        }

        /* Apply to the stored group-integrated E_g/F_g.  Only the energy-space
         * face states above are spectral; the finite-volume update remains the
         * face difference for each group-integrated conserved variable.  dt is
         * the effective stage weight (full dt in stage1, 0.5·dt in stage2 to
         * match the RK2-Heun mixing of the explicit RHS.  The lapse factor α(r) accounts
         * for the GR proper-time slowdown in the gravitational well, consistent
         * with the lapse multipliers already on the spatial radiation flux and
         * on the gravity source. */
        for (g = 0; g < PRJ_NEGROUP; ++g) {
            u[PRJ_CONS_RAD_E(field, g)] += dt_lapse *
                (energy_face[g + 1] - energy_face[g]);
            u[PRJ_CONS_RAD_F1(field, g)] += dt_lapse *
                (momentum_face[g + 1][0] - momentum_face[g][0]);
            u[PRJ_CONS_RAD_F2(field, g)] += dt_lapse *
                (momentum_face[g + 1][1] - momentum_face[g][1]);
            u[PRJ_CONS_RAD_F3(field, g)] += dt_lapse *
                (momentum_face[g + 1][2] - momentum_face[g][2]);
        }
    }
#elif PRJ_USE_RADIATION_FSA
    double dvdx[3][3] = {{0.0}};
    double a[3] = {0.0, 0.0, 0.0};
    double dt_lapse;
    double inv_c = 1.0 / PRJ_CLIGHT;
    int have_v_riemann;
    int cell_idx;
    int field;

    if (rad == 0 || block == 0 || W_state == 0 || u == 0) {
        return;
    }

    have_v_riemann = block->v_riemann[0] != 0 &&
        block->v_riemann[1] != 0 && block->v_riemann[2] != 0;
    if (have_v_riemann) {
        double inv_dx[3];
        int jdir;
        int icomp;

        inv_dx[0] = 1.0 / block->dx[0];
        inv_dx[1] = 1.0 / block->dx[1];
        inv_dx[2] = 1.0 / block->dx[2];

        for (jdir = 0; jdir < 3; ++jdir) {
            for (icomp = 0; icomp < 3; ++icomp) {
                int il = ic;
                int jl = jc;
                int kl = kc;
                int ir = ic;
                int jr = jc;
                int kr = kc;
                double vL;
                double vR;

                if (jdir == X1DIR) {
                    ir = ic + 1;
                } else if (jdir == X2DIR) {
                    jr = jc + 1;
                } else {
                    kr = kc + 1;
                }
                vL = block->v_riemann[jdir][VRIDX(icomp, il, jl, kl)];
                vR = block->v_riemann[jdir][VRIDX(icomp, ir, jr, kr)];
                dvdx[jdir][icomp] = (vR - vL) * inv_dx[jdir];
            }
        }
    }

    cell_idx = IDX(ic, jc, kc);
    if (block->grav[0] != 0 && block->grav[1] != 0 && block->grav[2] != 0) {
        a[0] = block->grav[0][cell_idx];
        a[1] = block->grav[1][cell_idx];
        a[2] = block->grav[2][cell_idx];
    }

    dt_lapse = lapse * dt;
    if (dt_lapse == 0.0) {
        return;
    }

    for (field = 0; field < PRJ_NRAD; ++field) {
        double omega_face[PRJ_NEGROUP + 1];
        double dnu[PRJ_NEGROUP];
        double inv_dnu[PRJ_NEGROUP];
        int angle;
        int g;

        prj_rad_fsa_omega_faces(rad, field, omega_face);
        /* omega_face is the energy-face coordinate for the FSA frequency flux.
         * It is copied from rad->eedge[field] when present, otherwise rebuilt
         * from rad->emin/emax; all of these radiation energy coordinates are in
         * MeV.  Since the numerator face energy and the denominator group width
         * below both use MeV, omega_face / dnu is dimensionless and no
         * PRJ_MEV_TO_ERG factor appears here. */
        for (g = 0; g < PRJ_NEGROUP; ++g) {
            dnu[g] = omega_face[g + 1] - omega_face[g];
            if (dnu[g] <= 0.0) {
                fprintf(stderr,
                    "prj_rad_freq_flux_apply: non-positive FSA energy width for field %d group %d\n",
                    field, g);
                exit(1);
            }
            inv_dnu[g] = 1.0 / dnu[g];
        }

        for (angle = 0; angle < PRJ_NANGLE; ++angle) {
            double n[3];
            double ndvdxn = 0.0;
            double a_dot_n;
            double drift;
            double J_group[PRJ_NEGROUP];
            double J_spec[PRJ_NEGROUP];
            double freq_face[PRJ_NEGROUP + 1] = {0.0};
            double energy_available[PRJ_NEGROUP];
            double outgoing[PRJ_NEGROUP] = {0.0};
            double theta[PRJ_NEGROUP];
            int ii;
            int jj;
            int gf;

            prj_rad_fsa_rotated_angle_dir(rad, block, angle, ic, jc, kc, n);
            a_dot_n = a[0] * n[0] + a[1] * n[1] + a[2] * n[2];
            for (jj = 0; jj < 3; ++jj) {
                for (ii = 0; ii < 3; ++ii) {
                    ndvdxn += n[jj] * dvdx[jj][ii] * n[ii];
                }
            }
            drift = ndvdxn + a_dot_n * inv_c;
            if (drift == 0.0) {
                continue;
            }

            for (g = 0; g < PRJ_NEGROUP; ++g) {
                int v = PRJ_CONS_RAD_I(field, g, angle);

                J_group[g] = W_state[WIDX(PRJ_RAD_PRIM_I(field, g, angle), ic, jc, kc)];
                J_spec[g] = J_group[g] * inv_dnu[g];
                energy_available[g] = u[v];
            }

            for (gf = 1; gf < PRJ_NEGROUP; ++gf) {
                int gu;
                double J_face;

                if (drift >= 0.0) {
                    gu = gf;
                    J_face = prj_rad_recon_face(J_spec, gu, -1);
                } else {
                    gu = gf - 1;
                    J_face = prj_rad_recon_face(J_spec, gu, +1);
                }
                freq_face[gf] = omega_face[gf] * drift * J_face;
            }

            /* Same donor-group positivity limiter as the M1 frequency flux:
             * dJ_g/dt = face[g+1] - face[g], so a positive face drains its
             * upper group and a negative face drains its lower group. */
            for (gf = 1; gf < PRJ_NEGROUP; ++gf) {
                if (freq_face[gf] > 0.0) {
                    outgoing[gf] += freq_face[gf];
                } else if (freq_face[gf] < 0.0) {
                    outgoing[gf - 1] -= freq_face[gf];
                }
            }
            for (g = 0; g < PRJ_NEGROUP; ++g) {
                double drain = dt_lapse * outgoing[g];

                theta[g] = 1.0;
                if (drain > energy_available[g]) {
                    theta[g] = energy_available[g] > 0.0 && drain > 0.0
                        ? nextafter(energy_available[g] / drain, 0.0)
                        : 0.0;
                }
            }
            for (gf = 1; gf < PRJ_NEGROUP; ++gf) {
                int donor = freq_face[gf] > 0.0 ? gf : gf - 1;
                double factor = theta[donor];

                freq_face[gf] *= factor;
            }

            /* Stored FSA angular variables are group-integrated energy density
             * in each angular cell.  Only the face state above is spectral, so the
             * finite-volume update still applies face[g+1] - face[g] directly
             * to each group-integrated angular variable. */
            for (g = 0; g < PRJ_NEGROUP; ++g) {
                u[PRJ_CONS_RAD_I(field, g, angle)] += dt_lapse *
                    (freq_face[g + 1] - freq_face[g]);
            }
        }
    }
#else
    (void)rad;
    (void)block;
    (void)W_state;
    (void)u;
    (void)ic;
    (void)jc;
    (void)kc;
    (void)lapse;
    (void)dt;
#endif
}

void prj_rad_ang_flux_apply(const prj_rad *rad, const prj_block *block,
    const double *W_state, double *u, int ic, int jc, int kc, double lapse, double dt)
{
#if PRJ_USE_RADIATION_FSA
    double dvdx[3][3] = {{0.0}};
    double a[3] = {0.0, 0.0, 0.0};
    double dt_lapse;
    double inv_c = 1.0 / PRJ_CLIGHT;
    double arc_factor[PRJ_NARC];
    int arc_donor[PRJ_NARC];
    int have_v_riemann;
    int cell_idx;
    int field;
    int group;
    int arc;
    int have_arc_flux;

    if (rad == 0 || block == 0 || W_state == 0 || u == 0) {
        return;
    }
    if (rad->arc_angle == 0 || rad->arc_vec == 0 ||
#if !PRJ_USE_RADIAL_FRAME_FSA
        rad->arc_nface == 0 ||
#endif
        rad->arc_neighbor == 0 || rad->cell_neighbor == 0) {
        return;
    }

    have_v_riemann = block->v_riemann[0] != 0 &&
        block->v_riemann[1] != 0 && block->v_riemann[2] != 0;
    if (have_v_riemann) {
        double inv_dx[3];
        int jdir;
        int icomp;

        inv_dx[0] = 1.0 / block->dx[0];
        inv_dx[1] = 1.0 / block->dx[1];
        inv_dx[2] = 1.0 / block->dx[2];

        for (jdir = 0; jdir < 3; ++jdir) {
            for (icomp = 0; icomp < 3; ++icomp) {
                int il = ic;
                int jl = jc;
                int kl = kc;
                int ir = ic;
                int jr = jc;
                int kr = kc;
                double vL;
                double vR;

                if (jdir == X1DIR) {
                    ir = ic + 1;
                } else if (jdir == X2DIR) {
                    jr = jc + 1;
                } else {
                    kr = kc + 1;
                }
                vL = block->v_riemann[jdir][VRIDX(icomp, il, jl, kl)];
                vR = block->v_riemann[jdir][VRIDX(icomp, ir, jr, kr)];
                dvdx[jdir][icomp] = (vR - vL) * inv_dx[jdir];
            }
        }
    }

    cell_idx = IDX(ic, jc, kc);
    if (block->grav[0] != 0 && block->grav[1] != 0 && block->grav[2] != 0) {
        a[0] = block->grav[0][cell_idx];
        a[1] = block->grav[1][cell_idx];
        a[2] = block->grav[2][cell_idx];
    }

    dt_lapse = lapse * dt;
    if (dt_lapse == 0.0) {
        return;
    }

    have_arc_flux = 0;
    for (arc = 0; arc < PRJ_NARC; ++arc) {
        int c0 = rad->arc_neighbor[2 * arc];
        int c1 = rad->arc_neighbor[2 * arc + 1];
#if PRJ_USE_RADIAL_FRAME_FSA
        const double *arc_vec_ref = &rad->arc_vec[3 * arc];
        double n0[3];
        double n1[3];
        double n_arc[3];
        double arc_vec[3];
        double mag;
        double geom[3] = {0.0, 0.0, 0.0};
#else
        const double *nface = &rad->arc_nface[3 * arc];
        const double *arc_vec = &rad->arc_vec[3 * arc];
#endif
        double b[3];
        double speed;
        int d;
        int jj;

        arc_factor[arc] = 0.0;
        arc_donor[arc] = -1;
        if (c0 < 0 || c0 >= PRJ_NANGLE || c1 < 0 || c1 >= PRJ_NANGLE) {
            continue;
        }
#if PRJ_USE_RADIAL_FRAME_FSA
        prj_rad_fsa_rotated_angle_dir(rad, block, c0, ic, jc, kc, n0);
        prj_rad_fsa_rotated_angle_dir(rad, block, c1, ic, jc, kc, n1);
        for (d = 0; d < 3; ++d) {
            n_arc[d] = n0[d] + n1[d];
        }
        mag = sqrt(prj_rad_fsa_dot(n_arc, n_arc));
        if (mag <= 0.0) {
            continue;
        }
        for (d = 0; d < 3; ++d) {
            n_arc[d] /= mag;
        }
        prj_rad_fsa_rotated_dir(block, ic, jc, kc, arc_vec_ref, arc_vec);
        if (block->ang_geom_fsa != 0) {
            for (d = 0; d < 3; ++d) {
                geom[d] = block->ang_geom_fsa[PRJ_FSA_ANG_GEOM_IDX(arc, d, ic, jc, kc)];
            }
        }
        for (d = 0; d < 3; ++d) {
            b[d] = a[d] * inv_c;
        }
        for (d = 0; d < 3; ++d) {
            for (jj = 0; jj < 3; ++jj) {
                b[d] += n_arc[jj] * dvdx[jj][d];
            }
        }

        speed = b[0] * arc_vec[0] + b[1] * arc_vec[1] + b[2] * arc_vec[2] -
            (geom[0] * arc_vec[0] + geom[1] * arc_vec[1] + geom[2] * arc_vec[2]);
#else
        for (d = 0; d < 3; ++d) {
            b[d] = a[d] * inv_c;
        }
        for (d = 0; d < 3; ++d) {
            for (jj = 0; jj < 3; ++jj) {
                b[d] += nface[jj] * dvdx[jj][d];
            }
        }

        speed = b[0] * arc_vec[0] + b[1] * arc_vec[1] + b[2] * arc_vec[2];
#endif
        if (speed == 0.0) {
            continue;
        }

        /* The angular face speed and donor depend only on this spatial cell,
         * not on radiation species or energy group, so reuse them below. */
        arc_factor[arc] = speed * rad->arc_angle[arc];
        arc_donor[arc] = speed > 0.0 ? c1 : c0;
        have_arc_flux = 1;
    }
    if (!have_arc_flux) {
        return;
    }

    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            double arc_flux[PRJ_NARC];
            double outgoing[PRJ_NANGLE];
            double theta[PRJ_NANGLE];
            int angle;

            for (arc = 0; arc < PRJ_NARC; ++arc) {
                int donor = arc_donor[arc];
                double I_face;

                arc_flux[arc] = 0.0;
                if (donor < 0) {
                    continue;
                }

                /* The stored angular variable is J = ΔΩ I, but the angular
                 * face flux still needs the per-solid-angle intensity I.  The
                 * LHS contains -alpha div_n(I b), so after moving it to the
                 * update side a positive b·arc_vec drains the second angular
                 * cell and fills the first. */
                I_face = W_state[WIDX(PRJ_RAD_PRIM_I(field, group, donor), ic, jc, kc)] /
                    rad->solid_angle[donor];
                arc_flux[arc] = arc_factor[arc] * I_face;
            }

            for (angle = 0; angle < PRJ_NANGLE; ++angle) {
                outgoing[angle] = 0.0;
            }
            for (arc = 0; arc < PRJ_NARC; ++arc) {
                int c0 = rad->arc_neighbor[2 * arc];
                int c1 = rad->arc_neighbor[2 * arc + 1];
                double flux = arc_flux[arc];

                if (flux > 0.0) {
                    outgoing[c1] += flux;
                } else if (flux < 0.0) {
                    outgoing[c0] -= flux;
                }
            }
            for (angle = 0; angle < PRJ_NANGLE; ++angle) {
                int v = PRJ_CONS_RAD_I(field, group, angle);
                double drain = dt_lapse * outgoing[angle];

                theta[angle] = 1.0;
                if (drain > u[v]) {
                    theta[angle] = u[v] > 0.0 && drain > 0.0
                        ? nextafter(u[v] / drain, 0.0)
                        : 0.0;
                }
            }
            for (arc = 0; arc < PRJ_NARC; ++arc) {
                double flux = arc_flux[arc];

                if (flux > 0.0) {
                    int donor = rad->arc_neighbor[2 * arc + 1];

                    arc_flux[arc] *= theta[donor];
                } else if (flux < 0.0) {
                    int donor = rad->arc_neighbor[2 * arc];

                    arc_flux[arc] *= theta[donor];
                }
            }

            for (arc = 0; arc < PRJ_NARC; ++arc) {
                int c0 = rad->arc_neighbor[2 * arc];
                int c1 = rad->arc_neighbor[2 * arc + 1];
                double flux = arc_flux[arc];

                u[PRJ_CONS_RAD_I(field, group, c0)] += dt_lapse * flux;
                u[PRJ_CONS_RAD_I(field, group, c1)] -= dt_lapse * flux;
            }
        }
    }
#else
    (void)rad;
    (void)block;
    (void)W_state;
    (void)u;
    (void)ic;
    (void)jc;
    (void)kc;
    (void)lapse;
    (void)dt;
#endif
}

#else
void prj_rad_energy_update(prj_rad *rad, prj_eos *eos, double *u, double dt, double lapse, double *final_temperature, double *kappa_out)
{
    (void)rad;
    (void)eos;
    (void)u;
    (void)dt;
    (void)final_temperature;
    (void)lapse;
    (void)kappa_out;
}

void prj_rad_momentum_update(prj_rad *rad, prj_eos *eos, double *u, double dt, double lapse, double temperature, const double *kappa_in)
{
    (void)rad;
    (void)eos;
    (void)u;
    (void)dt;
    (void)lapse;
    (void)temperature;
    (void)kappa_in;
}

void prj_rad_freq_flux_apply(const prj_rad *rad, const prj_block *block,
    const double *W_state, double *u, int ic, int jc, int kc, double lapse, double dt)
{
    (void)rad;
    (void)block;
    (void)W_state;
    (void)u;
    (void)ic;
    (void)jc;
    (void)kc;
    (void)lapse;
    (void)dt;
}

void prj_rad_ang_flux_apply(const prj_rad *rad, const prj_block *block,
    const double *W_state, double *u, int ic, int jc, int kc, double lapse, double dt)
{
    (void)rad;
    (void)block;
    (void)W_state;
    (void)u;
    (void)ic;
    (void)jc;
    (void)kc;
    (void)lapse;
    (void)dt;
}
#endif
