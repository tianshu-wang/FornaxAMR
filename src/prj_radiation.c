#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(PRJ_ENABLE_MPI)
#include <mpi.h>
#endif

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
/* Spatial M1 pressure tensor in the Eulerian (lab) frame, algebraic Levermore
 * closure:  P^{ij} = E a(f) gamma^{ij} + b(f) F^i F^j / (c^2 E),
 *   a(f) = (sqrt(4-3 f^2) - 1)/3,  b(f) = 3 / (2 + sqrt(4-3 f^2)),
 *   f = |F| / (c E)  (lab-frame flux factor).
 * Callers guarantee E > 0 and f <= 1; f is clamped only as cheap insurance
 * against non-finite inputs.  gamma_inv is the spatial inverse 3-metric
 * gamma^{ij}; Fcon is the contravariant spatial flux; Fmag = sqrt(F_i F^i).
 * Reduces to E/3 gamma^{ij} at f=0 and E n^i n^j at f=1. */
void prj_rad_grm1_pressure_lab(double E, const double Fcon[3],
    double Fmag, const double gamma_inv[3][3], double P[3][3])
{
    prj_rad_grm1_closure_coeffs coeffs;
    int i;
    int j;

    prj_rad_grm1_closure_coeffs_from_F2(E, Fmag * Fmag, &coeffs);
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            P[i][j] = coeffs.a_coef * gamma_inv[i][j] +
                coeffs.b_coef * Fcon[i] * Fcon[j];
        }
    }
}

int prj_rad_grm1_build_R(const double g_cov[4][4], const double g_con[4][4],
    double alpha, double E, const double Fcon[3], double Rcon[4][4])
{
    double ncon[4];
    double gamma_inv[3][3];
    double P[3][3];
    double F2 = 0.0;
    double Fmag;
    double invc = 1.0 / PRJ_CLIGHT;
    int a;
    int b;

    if (Rcon != 0) {
        memset(Rcon, 0, 16 * sizeof(double));
    }
    if (Rcon == 0 || g_cov == 0 || g_con == 0 || Fcon == 0 ||
        !isfinite(alpha) || alpha <= 0.0 || !isfinite(E) || E < 0.0) {
        return 0;
    }
    for (a = 0; a < 4; ++a) {
        for (b = 0; b < 4; ++b) {
            if (!isfinite(g_cov[a][b]) || !isfinite(g_con[a][b])) {
                return 0;
            }
        }
    }
    for (a = 0; a < 3; ++a) {
        if (!isfinite(Fcon[a])) {
            return 0;
        }
    }

    ncon[0] = 1.0 / alpha;
    for (a = 0; a < 3; ++a) {
        ncon[a + 1] = -alpha * g_con[0][a + 1];
    }

    /* Radiation vanishes for E == 0 (callers guarantee E > 0 for transport;
     * this branch keeps the null seed exact and avoids the f = 0/0 pressure). */
    if (E == 0.0) {
        for (a = 0; a < 3; ++a) {
            if (Fcon[a] != 0.0) {
                return 0;
            }
        }
        return 1;
    }

    /* Spatial inverse 3-metric gamma^{ij} and lab-frame flux magnitude.  The
     * 4-inverse spatial block g^{ij} = gamma^{ij} - beta^i beta^j/alpha^2, and
     * ncon^i = -beta^i/alpha, so gamma^{ij} = g^{ij} + n^i n^j. */
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            gamma_inv[a][b] = g_con[a + 1][b + 1] + ncon[a + 1] * ncon[b + 1];
            F2 += g_cov[a + 1][b + 1] * Fcon[a] * Fcon[b];
        }
    }
    if (!isfinite(F2) || F2 < 0.0) {
        F2 = 0.0;
    }
    Fmag = sqrt(F2);

    /* Algebraic Eulerian-frame Levermore closure:
     *   R^{ab} = E n^a n^b + (n^a F^b + n^b F^a)/c + P^{ab},
     * with F^a purely spatial (F^0 = 0) and P^{ab} spatial (P^{0a} = 0).  The
     * time-space rows R^{0a} = E n^0 n^a + n^0 F^a/c are identical to the
     * previous closure's Eq. (2); only the spatial block R^{ij} changes. */
    prj_rad_grm1_pressure_lab(E, Fcon, Fmag, gamma_inv, P);

    Rcon[0][0] = E * ncon[0] * ncon[0];
    for (a = 0; a < 3; ++a) {
        Rcon[0][a + 1] = E * ncon[0] * ncon[a + 1] + ncon[0] * Fcon[a] * invc;
        Rcon[a + 1][0] = Rcon[0][a + 1];
    }
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            Rcon[a + 1][b + 1] = E * ncon[a + 1] * ncon[b + 1] +
                (ncon[a + 1] * Fcon[b] + ncon[b + 1] * Fcon[a]) * invc +
                P[a][b];
        }
    }
    for (a = 0; a < 4; ++a) {
        for (b = 0; b < 4; ++b) {
            if (!isfinite(Rcon[a][b])) {
                memset(Rcon, 0, 16 * sizeof(double));
                return 0;
            }
        }
    }
    return 1;
}

typedef struct prj_rad_grm1_R_jac_geom {
    double ncon[4];
    double T[4][4];
    double dFcon[4][3];
} prj_rad_grm1_R_jac_geom;

static void prj_rad_grm1_R_jac_geom_init(double alpha,
    const double g_con[4][4], const double gamma_inv[3][3],
    prj_rad_grm1_R_jac_geom *jac_geom)
{
    int a;
    int b;
    int cc;

    memset(jac_geom, 0, sizeof(*jac_geom));
    jac_geom->ncon[0] = 1.0 / alpha;
    for (a = 0; a < 3; ++a) {
        jac_geom->ncon[a + 1] = -alpha * g_con[0][a + 1];
    }

    jac_geom->T[0][0] = jac_geom->ncon[0] * jac_geom->ncon[0];
    for (a = 0; a < 3; ++a) {
        jac_geom->T[0][a + 1] =
            jac_geom->ncon[0] * jac_geom->ncon[a + 1];
    }
    for (cc = 1; cc < 4; ++cc) {
        int jc = cc - 1;

        for (a = 0; a < 3; ++a) {
            jac_geom->T[cc][a + 1] =
                jac_geom->ncon[0] * gamma_inv[a][jc] / PRJ_CLIGHT;
        }
    }
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            jac_geom->dFcon[a + 1][b] = gamma_inv[a][b];
        }
    }
}

/* Analytic Jacobian of the algebraic M1 closure with respect to the lab-frame
 * (Eulerian) primitives p = (E, F_1, F_2, F_3), where F_i is the covariant
 * flux (raised to F^i with the spatial inverse metric `gamma_inv`).  The tensor
 * is reconstructed exactly as prj_rad_grm1_build_R; this routine additionally
 * propagates the derivatives through the closure scalars (see the "lab-frame
 * (E,F_i) reparametrization" note, Section 6).  Outputs the tensor in `Rcon`
 * and dR^{ab}/dp_c in `dRcon[c][a][b]` (c=0 is E, c=1..3 is F_1..F_3).  The R
 * part reproduces build_R bit-for-bit so the residual and Jacobian linearize
 * about the same stress tensor. */
static int prj_rad_grm1_build_R_jac(const double g_cov[4][4],
    const double g_con[4][4], const prj_rad_grm1_R_jac_geom *jac_geom,
    double E, const double Fcov[3], double Rcon[4][4],
    double dRcon[4][4][4])
{
    double Fcon[3];
    double gamma_inv[3][3];
    const double *ncon;
    const double (*T)[4];
    double F2 = 0.0;
    double Fmag;
    double cE;
    double f;
    double f2;
    double s;
    double a_coef;
    double b_coef;
    double D;
    double da_dE;
    double db_dE;
    double invc = 1.0 / PRJ_CLIGHT;
    double c2 = PRJ_CLIGHT * PRJ_CLIGHT;
    int a;
    int b;
    int cc;
    int i;
    int j;
    int m;

    if (Rcon == 0 || dRcon == 0) {
        if (Rcon != 0) {
            memset(Rcon, 0, 16 * sizeof(double));
        }
        if (dRcon != 0) {
            memset(dRcon, 0, 64 * sizeof(double));
        }
        return 0;
    }
    if (g_cov == 0 || g_con == 0 || jac_geom == 0 || Fcov == 0 ||
        !isfinite(E) || E < 0.0) {
        goto fail;
    }
    ncon = jac_geom->ncon;
    T = jac_geom->T;
    for (a = 0; a < 3; ++a) {
        if (!isfinite(Fcov[a])) {
            goto fail;
        }
        Fcon[a] = jac_geom->dFcon[a + 1][0] * Fcov[0] +
            jac_geom->dFcon[a + 1][1] * Fcov[1] +
            jac_geom->dFcon[a + 1][2] * Fcov[2];
        if (!isfinite(Fcon[a])) {
            goto fail;
        }
    }
    /* True inverse 3-metric gamma^{ij} = g^{ij} + n^i n^j (see build_R). */
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            gamma_inv[i][j] = g_con[i + 1][j + 1] + ncon[i + 1] * ncon[j + 1];
        }
    }

    if (E == 0.0) {
        double dFcon[4][4] = {{0.0}};

        memset(Rcon, 0, 16 * sizeof(double));
        for (a = 0; a < 3; ++a) {
            if (Fcon[a] != 0.0) {
                goto fail;
            }
            for (b = 0; b < 3; ++b) {
                dFcon[b + 1][a] = jac_geom->dFcon[b + 1][a];
            }
        }
        /* R = 0 at E = 0.  Seed the Jacobian with the isotropic normal-frame
         * limit, which is the exact F -> 0 limit of the algebraic closure:
         *   dR/dE = (4 n^a n^b + g^{ab})/3,
         *   dR/dF_j = (dF^a_j n^b + n^a dF^b_j)/c.
         */
        for (a = 0; a < 4; ++a) {
            for (b = 0; b < 4; ++b) {
                dRcon[0][a][b] = ((4.0 / 3.0) * ncon[a] * ncon[b]) +
                    (g_con[a][b] / 3.0);
            }
        }
        for (cc = 1; cc < 4; ++cc) {
            int jc = cc - 1;

            for (a = 0; a < 4; ++a) {
                for (b = 0; b < 4; ++b) {
                    dRcon[cc][a][b] =
                        (dFcon[a][jc] * ncon[b] +
                         ncon[a] * dFcon[b][jc]) / PRJ_CLIGHT;
                }
            }
        }
        for (cc = 0; cc < 4; ++cc) {
            for (a = 0; a < 4; ++a) {
                for (b = 0; b < 4; ++b) {
                    if (!isfinite(dRcon[cc][a][b])) {
                        goto fail;
                    }
                }
            }
        }
        return 1;
    }

    /* Lab-frame flux factor and Levermore closure scalars (must match the
     * arithmetic in prj_rad_grm1_pressure_lab so Rcon reproduces build_R). */
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            F2 += g_cov[i + 1][j + 1] * Fcon[i] * Fcon[j];
        }
    }
    if (!isfinite(F2) || F2 < 0.0) {
        F2 = 0.0;
    }
    Fmag = sqrt(F2);
    cE = PRJ_CLIGHT * E;
    f = Fmag / cE;
    if (!isfinite(f) || f < 0.0) {
        f = 0.0;
    } else if (f > 1.0) {
        f = 1.0;
    }
    f2 = f * f;
    s = sqrt(4.0 - 3.0 * f2);
    a_coef = E * (s - 1.0) / 3.0;
    b_coef = 3.0 / ((2.0 + s) * PRJ_CLIGHT * PRJ_CLIGHT * E);

    /* R^{ab} = E n^a n^b + (n^a F^b + n^b F^a)/c + P^{ab}. */
    Rcon[0][0] = E * ncon[0] * ncon[0];
    for (a = 0; a < 3; ++a) {
        Rcon[0][a + 1] = E * ncon[0] * ncon[a + 1] + ncon[0] * Fcon[a] * invc;
        Rcon[a + 1][0] = Rcon[0][a + 1];
    }
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            double P = a_coef * gamma_inv[i][j] + b_coef * Fcon[i] * Fcon[j];

            Rcon[i + 1][j + 1] = E * ncon[i + 1] * ncon[j + 1] +
                (ncon[i + 1] * Fcon[j] + ncon[j + 1] * Fcon[i]) * invc + P;
        }
    }
    for (a = 0; a < 4; ++a) {
        for (b = 0; b < 4; ++b) {
            if (!isfinite(Rcon[a][b])) {
                goto fail;
            }
        }
    }

    /* Closure-scalar derivatives w.r.t. E (f = Fmag/(cE), ds/df = -3f/s):
     *   ds/dE = 3 f^2/(s E),  D = (2+s) E,  b_coef = 3/(c^2 D). */
    D = (2.0 + s) * E;
    da_dE = (s - 1.0) / 3.0 + f2 / s;
    db_dE = -3.0 * ((2.0 + s) + 3.0 * f2 / s) / (c2 * D * D);

    /* dR^{0a}/dp_c are the geometry rows in T (unchanged from build_R). */
    for (cc = 0; cc < 4; ++cc) {
        dRcon[cc][0][0] = T[cc][0];
        for (a = 0; a < 3; ++a) {
            dRcon[cc][0][a + 1] = T[cc][a + 1];
            dRcon[cc][a + 1][0] = T[cc][a + 1];
        }
    }

    /* Spatial block dR^{ij}/dE. */
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            dRcon[0][i + 1][j + 1] = ncon[i + 1] * ncon[j + 1] +
                da_dE * gamma_inv[i][j] + db_dE * Fcon[i] * Fcon[j];
        }
    }

    /* Spatial block dR^{ij}/dF_m (covariant flux component m). */
    for (m = 0; m < 3; ++m) {
        double gm = Fcon[m];
        double da_dFm = -gm / (s * c2 * E);
        double dD_dFm = -3.0 * gm / (s * c2 * E);
        double db_dFm = -3.0 * dD_dFm / (c2 * D * D);

        for (i = 0; i < 3; ++i) {
            for (j = 0; j < 3; ++j) {
                dRcon[m + 1][i + 1][j + 1] =
                    (ncon[i + 1] * gamma_inv[j][m] +
                     ncon[j + 1] * gamma_inv[i][m]) * invc +
                    da_dFm * gamma_inv[i][j] +
                    db_dFm * Fcon[i] * Fcon[j] +
                    b_coef * (gamma_inv[i][m] * Fcon[j] +
                        Fcon[i] * gamma_inv[j][m]);
            }
        }
    }

    for (cc = 0; cc < 4; ++cc) {
        for (a = 0; a < 4; ++a) {
            for (b = 0; b < 4; ++b) {
                if (!isfinite(dRcon[cc][a][b])) {
                    goto fail;
                }
            }
        }
    }
    return 1;

fail:
    memset(Rcon, 0, 16 * sizeof(double));
    memset(dRcon, 0, 64 * sizeof(double));
    return 0;
}

typedef struct prj_rad_grm1_R_projector {
    double ncon[4];
    double gamma_inv[3][3];
    double Fcon[3];
    double E;
    double a_coef;
    double b_coef;
    double da_dE;
    double db_dE;
    double s;
    double D;
    int zero_E;
} prj_rad_grm1_R_projector;

static int prj_rad_grm1_R_projector_init(const double g_cov[4][4],
    const double g_con[4][4], const prj_rad_grm1_R_jac_geom *jac_geom,
    double E, const double Fcov[3], prj_rad_grm1_R_projector *proj)
{
    double F2 = 0.0;
    double Fmag;
    double cE;
    double f;
    double f2;
    double c2 = PRJ_CLIGHT * PRJ_CLIGHT;
    int a;
    int i;
    int j;

    if (proj != 0) {
        memset(proj, 0, sizeof(*proj));
    }
    if (g_cov == 0 || g_con == 0 || jac_geom == 0 || Fcov == 0 ||
        proj == 0 || !isfinite(E) || E < 0.0) {
        return 0;
    }
    proj->E = E;
    for (a = 0; a < 4; ++a) {
        proj->ncon[a] = jac_geom->ncon[a];
    }
    for (a = 0; a < 3; ++a) {
        if (!isfinite(Fcov[a])) {
            return 0;
        }
        proj->Fcon[a] = jac_geom->dFcon[a + 1][0] * Fcov[0] +
            jac_geom->dFcon[a + 1][1] * Fcov[1] +
            jac_geom->dFcon[a + 1][2] * Fcov[2];
        if (!isfinite(proj->Fcon[a])) {
            return 0;
        }
    }
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            proj->gamma_inv[i][j] = g_con[i + 1][j + 1] +
                proj->ncon[i + 1] * proj->ncon[j + 1];
        }
    }
    if (E == 0.0) {
        for (a = 0; a < 3; ++a) {
            if (proj->Fcon[a] != 0.0) {
                return 0;
            }
        }
        proj->zero_E = 1;
        return 1;
    }

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            F2 += g_cov[i + 1][j + 1] * proj->Fcon[i] * proj->Fcon[j];
        }
    }
    if (!isfinite(F2) || F2 < 0.0) {
        F2 = 0.0;
    }
    Fmag = sqrt(F2);
    cE = PRJ_CLIGHT * E;
    f = Fmag / cE;
    if (!isfinite(f) || f < 0.0) {
        f = 0.0;
    } else if (f > 1.0) {
        f = 1.0;
    }
    f2 = f * f;
    proj->s = sqrt(4.0 - 3.0 * f2);
    proj->a_coef = E * (proj->s - 1.0) / 3.0;
    proj->b_coef = 3.0 / ((2.0 + proj->s) * c2 * E);
    proj->D = (2.0 + proj->s) * E;
    proj->da_dE = (proj->s - 1.0) / 3.0 + f2 / proj->s;
    proj->db_dE = -3.0 * ((2.0 + proj->s) + 3.0 * f2 / proj->s) /
        (c2 * proj->D * proj->D);
    return isfinite(proj->a_coef) && isfinite(proj->b_coef) &&
        isfinite(proj->da_dE) && isfinite(proj->db_dE);
}

static void prj_rad_grm1_R_projector_apply(
    const prj_rad_grm1_R_projector *proj, const double cov[4], double out[4])
{
    double invc = 1.0 / PRJ_CLIGHT;
    int i;
    int j;

    for (i = 0; i < 4; ++i) {
        out[i] = 0.0;
    }
    if (proj == 0 || cov == 0 || proj->zero_E) {
        return;
    }
    out[0] = proj->E * proj->ncon[0] * proj->ncon[0] * cov[0];
    for (j = 0; j < 3; ++j) {
        double R0j = proj->E * proj->ncon[0] * proj->ncon[j + 1] +
            proj->ncon[0] * proj->Fcon[j] * invc;

        out[0] += R0j * cov[j + 1];
        out[j + 1] += R0j * cov[0];
    }
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            double Pij = proj->a_coef * proj->gamma_inv[i][j] +
                proj->b_coef * proj->Fcon[i] * proj->Fcon[j];
            double Rij = proj->E * proj->ncon[i + 1] *
                proj->ncon[j + 1] +
                (proj->ncon[i + 1] * proj->Fcon[j] +
                 proj->ncon[j + 1] * proj->Fcon[i]) * invc + Pij;

            out[i + 1] += Rij * cov[j + 1];
        }
    }
}

static void prj_rad_grm1_R_projector_dR_u(
    const prj_rad_grm1_R_projector *proj, const double g_con[4][4],
    const prj_rad_grm1_R_jac_geom *jac_geom, const double ucov[4],
    double dR_u[4][4], double dRuu[4])
{
    double invc = 1.0 / PRJ_CLIGHT;
    double c2 = PRJ_CLIGHT * PRJ_CLIGHT;
    int a;
    int b;
    int cc;
    int i;
    int j;
    int m;

    memset(dR_u, 0, 16 * sizeof(double));
    for (cc = 0; cc < 4; ++cc) {
        dRuu[cc] = 0.0;
    }
    if (proj == 0 || jac_geom == 0 || ucov == 0) {
        return;
    }

    if (proj->zero_E) {
        double dFcon[4][3] = {{0.0}};

        for (a = 0; a < 3; ++a) {
            for (b = 0; b < 3; ++b) {
                dFcon[b + 1][a] = jac_geom->dFcon[b + 1][a];
            }
        }
        for (a = 0; a < 4; ++a) {
            for (b = 0; b < 4; ++b) {
                double dR = ((4.0 / 3.0) * proj->ncon[a] *
                    proj->ncon[b]) + (g_con[a][b] / 3.0);

                dR_u[0][a] += dR * ucov[b];
            }
        }
        for (cc = 1; cc < 4; ++cc) {
            int jc = cc - 1;

            for (a = 0; a < 4; ++a) {
                for (b = 0; b < 4; ++b) {
                    double dR = (dFcon[a][jc] * proj->ncon[b] +
                        proj->ncon[a] * dFcon[b][jc]) * invc;

                    dR_u[cc][a] += dR * ucov[b];
                }
            }
        }
        for (cc = 0; cc < 4; ++cc) {
            for (a = 0; a < 4; ++a) {
                dRuu[cc] += dR_u[cc][a] * ucov[a];
            }
        }
        return;
    }

    for (cc = 0; cc < 4; ++cc) {
        dR_u[cc][0] = jac_geom->T[cc][0] * ucov[0];
        for (a = 0; a < 3; ++a) {
            dR_u[cc][0] += jac_geom->T[cc][a + 1] * ucov[a + 1];
            dR_u[cc][a + 1] = jac_geom->T[cc][a + 1] * ucov[0];
        }
    }

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            double dR = proj->ncon[i + 1] * proj->ncon[j + 1] +
                proj->da_dE * proj->gamma_inv[i][j] +
                proj->db_dE * proj->Fcon[i] * proj->Fcon[j];

            dR_u[0][i + 1] += dR * ucov[j + 1];
        }
    }

    for (m = 0; m < 3; ++m) {
        double gm = proj->Fcon[m];
        double da_dFm = -gm / (proj->s * c2 * proj->E);
        double dD_dFm = -3.0 * gm / (proj->s * c2 * proj->E);
        double db_dFm = -3.0 * dD_dFm / (c2 * proj->D * proj->D);

        for (i = 0; i < 3; ++i) {
            for (j = 0; j < 3; ++j) {
                double dR =
                    (proj->ncon[i + 1] * proj->gamma_inv[j][m] +
                     proj->ncon[j + 1] * proj->gamma_inv[i][m]) * invc +
                    da_dFm * proj->gamma_inv[i][j] +
                    db_dFm * proj->Fcon[i] * proj->Fcon[j] +
                    proj->b_coef * (proj->gamma_inv[i][m] *
                        proj->Fcon[j] + proj->Fcon[i] *
                        proj->gamma_inv[j][m]);

                dR_u[m + 1][i + 1] += dR * ucov[j + 1];
            }
        }
    }
    for (cc = 0; cc < 4; ++cc) {
        for (a = 0; a < 4; ++a) {
            dRuu[cc] += dR_u[cc][a] * ucov[a];
        }
    }
}

typedef struct prj_rad_grm1_m3_data {
    double ucon[4];
    double ucov[4];
    double hcon[4][4];
    double hcov[4][4];
    double H[4];
    double L[4][4];
    double J;
    double Hnorm2;
    double inv_Hnorm3;
    double thin_w;
    double thick_w;
} prj_rad_grm1_m3_data;

static int prj_rad_grm1_prepare_m3_data(const double g_cov[4][4],
    const double g_con[4][4], const double ucon[4],
    const double Rcon[4][4], prj_rad_grm1_m3_data *m3)
{
    double J = 0.0;
    double unorm = 0.0;
    double Hnorm2 = 0.0;
    double xi;
    double chi;
    int a;
    int b;

    if (m3 != 0) {
        memset(m3, 0, sizeof(*m3));
    }
    if (m3 == 0 || g_cov == 0 || g_con == 0 || ucon == 0 || Rcon == 0) {
        return 0;
    }
    for (a = 0; a < 4; ++a) {
        if (!isfinite(ucon[a])) {
            return 0;
        }
        for (b = 0; b < 4; ++b) {
            if (!isfinite(g_cov[a][b]) || !isfinite(g_con[a][b]) ||
                !isfinite(Rcon[a][b])) {
                return 0;
            }
        }
    }

    for (a = 0; a < 4; ++a) {
        m3->ucon[a] = ucon[a];
        m3->ucov[a] = 0.0;
        for (b = 0; b < 4; ++b) {
            m3->ucov[a] += g_cov[a][b] * ucon[b];
        }
        unorm += m3->ucov[a] * ucon[a];
    }
    if (!isfinite(unorm) || fabs(unorm + 1.0) > 1.0e-8) {
        return 0;
    }

    for (a = 0; a < 4; ++a) {
        for (b = 0; b < 4; ++b) {
            J += Rcon[a][b] * m3->ucov[a] * m3->ucov[b];
        }
    }
    if (!isfinite(J) || J < 0.0) {
        return 0;
    }
    if (J == 0.0) {
        for (a = 0; a < 4; ++a) {
            for (b = 0; b < 4; ++b) {
                if (Rcon[a][b] != 0.0) {
                    return 0;
                }
            }
        }
        return 1;
    }
    m3->J = J;

    for (a = 0; a < 4; ++a) {
        double Ru = 0.0;

        for (b = 0; b < 4; ++b) {
            Ru += Rcon[b][a] * m3->ucov[b];
        }
        m3->H[a] = -Ru - J * ucon[a];
        if (!isfinite(m3->H[a])) {
            return 0;
        }
    }
    for (a = 0; a < 4; ++a) {
        for (b = 0; b < 4; ++b) {
            m3->hcon[a][b] = g_con[a][b] + ucon[a] * ucon[b];
            m3->hcov[a][b] = g_cov[a][b] + m3->ucov[a] * m3->ucov[b];
            m3->L[a][b] = Rcon[a][b] - J * ucon[a] * ucon[b] -
                m3->H[a] * ucon[b] - m3->H[b] * ucon[a];
            Hnorm2 += m3->hcov[a][b] * m3->H[a] * m3->H[b];
            if (!isfinite(m3->hcon[a][b]) || !isfinite(m3->hcov[a][b]) ||
                !isfinite(m3->L[a][b])) {
                return 0;
            }
        }
    }

    if (!isfinite(Hnorm2) || Hnorm2 < -1.0e-12 * J * J) {
        return 0;
    }
    if (Hnorm2 < 0.0) {
        Hnorm2 = 0.0;
    }
    m3->Hnorm2 = Hnorm2;
    if (Hnorm2 > 0.0) {
        m3->inv_Hnorm3 = 1.0 / (Hnorm2 * sqrt(Hnorm2));
    }
    xi = sqrt(Hnorm2) / J;
    if (!isfinite(xi) || xi < 0.0) {
        return 0;
    }
    if (xi > 1.0) {
        if (xi > 1.0 + 1.0e-10) {
            return 0;
        }
        xi = 1.0;
    }
    chi = prj_rad_m1_chi_exact(xi);
    m3->thin_w = 0.5 * (3.0 * chi - 1.0);
    m3->thick_w = 1.5 * (1.0 - chi);
    return 1;
}

static double prj_rad_grm1_m3_component(
    const prj_rad_grm1_m3_data *m3, int a, int b, int c)
{
    double Nthin = 0.0;
    double Nthick;

    if (m3->inv_Hnorm3 > 0.0) {
        Nthin = m3->J * m3->H[a] * m3->H[b] * m3->H[c] *
            m3->inv_Hnorm3;
    }
    Nthick = 0.2 * (m3->H[a] * m3->hcon[b][c] +
        m3->H[b] * m3->hcon[a][c] + m3->H[c] * m3->hcon[a][b]);
    return
        m3->J * m3->ucon[a] * m3->ucon[b] * m3->ucon[c] +
        m3->H[a] * m3->ucon[b] * m3->ucon[c] +
        m3->H[b] * m3->ucon[a] * m3->ucon[c] +
        m3->H[c] * m3->ucon[a] * m3->ucon[b] +
        m3->L[a][b] * m3->ucon[c] + m3->L[a][c] * m3->ucon[b] +
        m3->L[b][c] * m3->ucon[a] +
        m3->thin_w * Nthin + m3->thick_w * Nthick;
}

int prj_rad_grm1_freq_drift(const double g_cov[4][4],
    const double g_con[4][4], const double ucon[4],
    const double Rcon[4][4], const double ducov[4][4], double drift[4])
{
    prj_rad_grm1_m3_data m3;
    double D_u[4] = {0.0, 0.0, 0.0, 0.0};
    double u_D[4] = {0.0, 0.0, 0.0, 0.0};
    double D_H[4] = {0.0, 0.0, 0.0, 0.0};
    double H_D[4] = {0.0, 0.0, 0.0, 0.0};
    double Duu = 0.0;
    double DHu = 0.0;
    double DuH = 0.0;
    double DL = 0.0;
    double DHH = 0.0;
    double Dh = 0.0;
    int a;
    int b;
    int c;

    if (drift != 0) {
        memset(drift, 0, 4 * sizeof(double));
    }
    if (drift == 0 || ducov == 0 ||
        !prj_rad_grm1_prepare_m3_data(g_cov, g_con, ucon, Rcon, &m3)) {
        return 0;
    }
    for (a = 0; a < 4; ++a) {
        for (b = 0; b < 4; ++b) {
            if (!isfinite(ducov[a][b])) {
                return 0;
            }
        }
    }
    if (m3.J == 0.0) {
        return 1;
    }

    for (b = 0; b < 4; ++b) {
        for (c = 0; c < 4; ++c) {
            double D = ducov[b][c];

            D_u[b] += D * m3.ucon[c];
            u_D[c] += m3.ucon[b] * D;
            D_H[b] += D * m3.H[c];
            H_D[c] += m3.H[b] * D;
            Duu += m3.ucon[b] * m3.ucon[c] * D;
            DHu += m3.H[b] * m3.ucon[c] * D;
            DuH += m3.ucon[b] * m3.H[c] * D;
            DL += m3.L[b][c] * D;
            DHH += m3.H[b] * m3.H[c] * D;
            Dh += m3.hcon[b][c] * D;
        }
    }

    for (a = 0; a < 4; ++a) {
        drift[a] = (m3.J * m3.ucon[a] + m3.H[a]) * Duu +
            m3.ucon[a] * (DHu + DuH + DL);
        for (b = 0; b < 4; ++b) {
            drift[a] += m3.L[a][b] * (D_u[b] + u_D[b]);
            drift[a] += 0.2 * m3.thick_w * m3.hcon[a][b] *
                (D_H[b] + H_D[b]);
        }
        drift[a] += m3.thin_w * m3.J * m3.H[a] * DHH *
            m3.inv_Hnorm3;
        drift[a] += 0.2 * m3.thick_w * m3.H[a] * Dh;
        if (!isfinite(drift[a])) {
            memset(drift, 0, 4 * sizeof(double));
            return 0;
        }
    }
    return 1;
}

/* Cell-local pieces of the GR M1 frequency drift.  These depend on the metric,
 * fluid velocity, velocity gradients, and observer acceleration, but not on the
 * radiation energy group. */
typedef struct prj_rad_gr_m1_freq_drift_cache {
    double g_cov[4][4];
    double g_con[4][4];
    double ucon[4];
    double ucov[4];
    double gcon_ucov[4];
    double hcon[4][4];
    double hcov[4][4];
    double g_ucov_ucov;
    double energy_Z;
    double energy_Y[3];
    double energy_X[3][3];
    double mom_Y[3][3];
    double mom_X[3][3][3];
    double mom_W[3][3][3][3];
    double q_uuu[4];
    double q_H[4][4];
    double q_L[4][4][4];
    double q_h[4][4];
} prj_rad_gr_m1_freq_drift_cache;

static void prj_rad_gr_m1_freq_precontract_add(
    prj_rad_gr_m1_freq_drift_cache *cache, int out,
    int a, int b, int c, double value)
{
    cache->q_uuu[out] += value * cache->ucon[a] * cache->ucon[b] *
        cache->ucon[c];
    cache->q_H[out][a] += value * cache->ucon[b] * cache->ucon[c];
    cache->q_H[out][b] += value * cache->ucon[a] * cache->ucon[c];
    cache->q_H[out][c] += value * cache->ucon[a] * cache->ucon[b];
    cache->q_L[out][a][b] += value * cache->ucon[c];
    cache->q_L[out][a][c] += value * cache->ucon[b];
    cache->q_L[out][b][c] += value * cache->ucon[a];
    cache->q_h[out][a] += value * cache->hcon[b][c];
    cache->q_h[out][b] += value * cache->hcon[a][c];
    cache->q_h[out][c] += value * cache->hcon[a][b];
}

static void prj_rad_gr_m1_freq_drift_precontract(
    prj_rad_gr_m1_freq_drift_cache *cache)
{
    int a;
    int x;
    int i;
    int j;
    int k;

    memset(cache->q_uuu, 0, sizeof(cache->q_uuu));
    memset(cache->q_H, 0, sizeof(cache->q_H));
    memset(cache->q_L, 0, sizeof(cache->q_L));
    memset(cache->q_h, 0, sizeof(cache->q_h));

    prj_rad_gr_m1_freq_precontract_add(cache, 0, 0, 0, 0,
        cache->energy_Z);
    for (x = 0; x < 3; ++x) {
        prj_rad_gr_m1_freq_precontract_add(cache, 0, 0, 0, x + 1,
            cache->energy_Y[x]);
        for (i = 0; i < 3; ++i) {
            prj_rad_gr_m1_freq_precontract_add(cache, 0, 0, x + 1, i + 1,
                cache->energy_X[x][i]);
        }
    }
    for (j = 0; j < 3; ++j) {
        int out_mom = j + 1;

        for (x = 0; x < 3; ++x) {
            prj_rad_gr_m1_freq_precontract_add(cache, out_mom, 0, 0, x + 1,
                cache->mom_Y[j][x]);
            for (i = 0; i < 3; ++i) {
                prj_rad_gr_m1_freq_precontract_add(cache, out_mom, 0, x + 1,
                    i + 1, cache->mom_X[j][x][i]);
            }
        }
        for (a = 0; a < 3; ++a) {
            for (k = 0; k < 3; ++k) {
                for (i = 0; i < 3; ++i) {
                    prj_rad_gr_m1_freq_precontract_add(cache, out_mom, a + 1,
                        k + 1, i + 1, cache->mom_W[j][a][k][i]);
                }
            }
        }
    }
}

static void prj_rad_gr_m1_freq_m3_contract_all(
    const prj_rad_gr_m1_freq_drift_cache *cache,
    const prj_rad_grm1_m3_data *m3, double value[4])
{
    double cubic[4] = {0.0, 0.0, 0.0, 0.0};
    int out;
    int a;
    int b;
    int c;
    int i;
    int j;
    int k;
    int x;

    for (out = 0; out < 4; ++out) {
        value[out] = m3->J * cache->q_uuu[out];
        for (a = 0; a < 4; ++a) {
            value[out] += m3->H[a] * cache->q_H[out][a];
        }
        for (a = 0; a < 4; ++a) {
            for (b = 0; b < 4; ++b) {
                value[out] += m3->L[a][b] * cache->q_L[out][a][b];
            }
        }
        for (a = 0; a < 4; ++a) {
            value[out] += 0.2 * m3->thick_w * m3->H[a] *
                cache->q_h[out][a];
        }
    }
    if (m3->inv_Hnorm3 > 0.0) {
        double hhh[4][4][4];

        for (a = 0; a < 4; ++a) {
            for (b = 0; b < 4; ++b) {
                for (c = 0; c < 4; ++c) {
                    hhh[a][b][c] = m3->H[a] * m3->H[b] * m3->H[c];
                }
            }
        }
        cubic[0] = cache->energy_Z * hhh[0][0][0];
        for (x = 0; x < 3; ++x) {
            cubic[0] += cache->energy_Y[x] * hhh[0][0][x + 1];
            for (i = 0; i < 3; ++i) {
                cubic[0] += cache->energy_X[x][i] *
                    hhh[0][x + 1][i + 1];
            }
        }
        for (j = 0; j < 3; ++j) {
            int out_mom = j + 1;

            for (x = 0; x < 3; ++x) {
                cubic[out_mom] += cache->mom_Y[j][x] *
                    hhh[0][0][x + 1];
                for (i = 0; i < 3; ++i) {
                    cubic[out_mom] += cache->mom_X[j][x][i] *
                        hhh[0][x + 1][i + 1];
                }
            }
            for (a = 0; a < 3; ++a) {
                for (k = 0; k < 3; ++k) {
                    for (i = 0; i < 3; ++i) {
                        cubic[out_mom] += cache->mom_W[j][a][k][i] *
                            hhh[a + 1][k + 1][i + 1];
                    }
                }
            }
        }
        for (out = 0; out < 4; ++out) {
            value[out] += m3->thin_w * m3->J * m3->inv_Hnorm3 *
                cubic[out];
        }
    }
}

static void prj_rad_gr_m1_freq_drift_cache_init(
    const prj_rad_gr_m1_closure_ctx *ctx,
    const prj_rad_gr_m1_side_data *side, const prj_z4c_hydro_geom *geom,
    const double observer_time_derivative[4],
    prj_rad_gr_m1_freq_drift_cache *cache)
{
    double dLam[3];
    double dLv[3][3];
    double dLv_cov[3][3];
    double nDLam;
    double nDLv[3];
    double da[3];
    double vdb[3] = {0.0, 0.0, 0.0};
    double vdb_cov[3] = {0.0, 0.0, 0.0};
    double nDLv_cov[3] = {0.0, 0.0, 0.0};
    double vdg[3][3] = {{0.0}};
    double sum_v_da = 0.0;
    double inv_alpha;
    int a;
    int b;
    int d;
    int i;
    int j;
    int k;
    int m;
    int x;

    memset(cache, 0, sizeof(*cache));

    cache->g_cov[0][0] = -1.0;
    cache->g_con[0][0] = -1.0;
    for (a = 0; a < 3; ++a) {
        cache->ucon[a + 1] = side->wlor * side->vcon[a];
        for (b = 0; b < 3; ++b) {
            cache->g_cov[a + 1][b + 1] = ctx->gamma[a][b];
            cache->g_con[a + 1][b + 1] = ctx->gamma_inv[a][b];
        }
    }
    cache->ucon[0] = side->wlor;

    for (a = 0; a < 4; ++a) {
        cache->ucov[a] = 0.0;
        for (b = 0; b < 4; ++b) {
            cache->ucov[a] += cache->g_cov[a][b] * cache->ucon[b];
        }
    }
    for (a = 0; a < 4; ++a) {
        cache->gcon_ucov[a] = 0.0;
        for (b = 0; b < 4; ++b) {
            cache->gcon_ucov[a] += cache->g_con[b][a] * cache->ucov[b];
        }
        cache->g_ucov_ucov += cache->gcon_ucov[a] * cache->ucov[a];
    }
    for (a = 0; a < 4; ++a) {
        for (b = 0; b < 4; ++b) {
            cache->hcon[a][b] = cache->g_con[a][b] +
                cache->ucon[a] * cache->ucon[b];
            cache->hcov[a][b] = cache->g_cov[a][b] +
                cache->ucov[a] * cache->ucov[b];
        }
    }

    for (d = 0; d < 3; ++d) {
        double dv2 = 0.0;

        for (a = 0; a < 3; ++a) {
            for (b = 0; b < 3; ++b) {
                dv2 += ctx->dgamma[d][a][b] * side->vcon[a] *
                    side->vcon[b] +
                    2.0 * ctx->gamma[a][b] * side->vcon[a] *
                    ctx->dvdx[d][b];
            }
        }
        dLam[d] = 0.5 * side->wlor * side->wlor * side->wlor * dv2;
        for (k = 0; k < 3; ++k) {
            dLv[d][k] = dLam[d] * side->vcon[k] +
                side->wlor * ctx->dvdx[d][k];
        }
    }
    for (d = 0; d < 3; ++d) {
        for (b = 0; b < 3; ++b) {
            dLv_cov[d][b] = 0.0;
            for (k = 0; k < 3; ++k) {
                dLv_cov[d][b] += ctx->gamma[b][k] * dLv[d][k];
            }
        }
    }

    inv_alpha = 1.0 / geom->alpha;
    nDLam = observer_time_derivative[0] / PRJ_CLIGHT;
    for (i = 0; i < 3; ++i) {
        nDLam -= geom->beta[i] * dLam[i];
    }
    nDLam *= inv_alpha;
    for (k = 0; k < 3; ++k) {
        nDLv[k] = observer_time_derivative[1 + k] / PRJ_CLIGHT;
        for (i = 0; i < 3; ++i) {
            nDLv[k] -= geom->beta[i] * dLv[i][k];
        }
        nDLv[k] *= inv_alpha;
    }

    for (i = 0; i < 3; ++i) {
        da[i] = geom->dalpha[i] * inv_alpha;
        sum_v_da += side->vcon[i] * da[i];
        for (k = 0; k < 3; ++k) {
            vdb[k] += side->vcon[i] * geom->dbeta[i][k] * inv_alpha;
        }
    }
    for (m = 0; m < 3; ++m) {
        for (k = 0; k < 3; ++k) {
            for (i = 0; i < 3; ++i) {
                vdg[k][i] += side->vcon[m] * ctx->dgamma[m][k][i];
            }
        }
    }
    for (x = 0; x < 3; ++x) {
        for (k = 0; k < 3; ++k) {
            vdb_cov[x] += ctx->gamma[k][x] * vdb[k];
            nDLv_cov[x] += ctx->gamma[k][x] * nDLv[k];
        }
    }

    cache->energy_Z = -(side->wlor * sum_v_da + nDLam);
    for (x = 0; x < 3; ++x) {
        double Rn_Y = -da[x] - vdb_cov[x];
        double On_Y = dLam[x] - nDLv_cov[x];

        cache->energy_Y[x] = -(side->wlor * Rn_Y + On_Y);
        for (i = 0; i < 3; ++i) {
            double Rn_X = ctx->K_dd[x][i] - 0.5 * vdg[x][i];
            double On_X = 0.0;

            for (k = 0; k < 3; ++k) {
                On_X -= ctx->gamma[k][x] * dLv[i][k];
            }
            cache->energy_X[x][i] = -(side->wlor * Rn_X + On_X);
        }
    }

    for (j = 0; j < 3; ++j) {
        for (x = 0; x < 3; ++x) {
            cache->mom_Y[j][x] =
                -(side->wlor * ctx->gamma[j][x] * sum_v_da +
                    ctx->gamma[j][x] * nDLam);
            for (i = 0; i < 3; ++i) {
                double Rg_X;
                double Og_X;

                Rg_X = ctx->gamma[j][x] * (-da[i] - vdb_cov[i]);
                Og_X = ctx->gamma[j][x] * (dLam[i] - nDLv_cov[i]);
                cache->mom_X[j][x][i] = -(side->wlor * Rg_X + Og_X);
            }
            for (k = 0; k < 3; ++k) {
                for (i = 0; i < 3; ++i) {
                    double Rg_W = ctx->gamma[j][x] *
                        (ctx->K_dd[k][i] - 0.5 * vdg[k][i]);
                    double Og_W = -ctx->gamma[j][x] * dLv_cov[i][k];

                    cache->mom_W[j][x][k][i] =
                        -(side->wlor * Rg_W + Og_W);
                }
            }
        }
    }
    prj_rad_gr_m1_freq_drift_precontract(cache);
}

static void prj_rad_grm1_prepare_m3_data_freq(
    const prj_rad_gr_m1_freq_drift_cache *cache, double E,
    const double Fcov_in[3], prj_rad_grm1_m3_data *m3)
{
    double Fcov[3];
    double Fcon[3] = {0.0, 0.0, 0.0};
    double gamma_inv[3][3];
    double P[3][3];
    double Rcon[4][4];
    double F2 = 0.0;
    double Fmag;
    double Fmax;
    double Hnorm2 = 0.0;
    double xi;
    double chi;
    double J;
    double invc = 1.0 / PRJ_CLIGHT;
    int a;
    int b;
    int i;
    int j;

    memset(m3, 0, sizeof(*m3));
    if (E < PRJ_RAD_GR_M1_E_FLOOR) {
        E = PRJ_RAD_GR_M1_E_FLOOR;
    }
    for (a = 0; a < 3; ++a) {
        Fcov[a] = Fcov_in[a];
    }
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            Fcon[a] += cache->g_con[a + 1][b + 1] * Fcov[b];
        }
    }
    for (a = 0; a < 3; ++a) {
        F2 += Fcov[a] * Fcon[a];
    }
    Fmag = sqrt(F2);
    Fmax = (1.0 - PRJ_RAD_GR_M1_F_MARGIN) * PRJ_CLIGHT * E;
    if (Fmag > Fmax) {
        double scale = Fmax / Fmag;

        for (a = 0; a < 3; ++a) {
            Fcov[a] *= scale;
            Fcon[a] *= scale;
        }
        Fmag = Fmax;
    }

    /* Algebraic lab-frame Levermore closure in the cache's local frame
     * (alpha = 1, beta = 0), where n^a = (1,0,0,0): R^{00}=E, R^{0i}=F^i/c,
     * and the spatial block is the pressure tensor P^{ij}. */
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            gamma_inv[i][j] = cache->g_con[i + 1][j + 1];
        }
    }
    prj_rad_grm1_pressure_lab(E, Fcon, Fmag, gamma_inv, P);
    Rcon[0][0] = E;
    for (a = 0; a < 3; ++a) {
        Rcon[0][a + 1] = Fcon[a] * invc;
        Rcon[a + 1][0] = Rcon[0][a + 1];
    }
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            Rcon[i + 1][j + 1] = P[i][j];
        }
    }

    /* Fluid-frame decomposition (uses the cached u and projectors). */
    for (a = 0; a < 4; ++a) {
        m3->ucon[a] = cache->ucon[a];
        m3->ucov[a] = cache->ucov[a];
    }
    J = 0.0;
    for (a = 0; a < 4; ++a) {
        for (b = 0; b < 4; ++b) {
            J += Rcon[a][b] * cache->ucov[a] * cache->ucov[b];
        }
    }
    m3->J = J;
    for (a = 0; a < 4; ++a) {
        double Ru = 0.0;

        for (b = 0; b < 4; ++b) {
            Ru += Rcon[b][a] * cache->ucov[b];
        }
        m3->H[a] = -Ru - J * cache->ucon[a];
    }
    for (a = 0; a < 4; ++a) {
        for (b = 0; b < 4; ++b) {
            m3->hcon[a][b] = cache->hcon[a][b];
            m3->hcov[a][b] = cache->hcov[a][b];
            m3->L[a][b] = Rcon[a][b] - J * cache->ucon[a] * cache->ucon[b] -
                m3->H[a] * cache->ucon[b] - m3->H[b] * cache->ucon[a];
            Hnorm2 += m3->hcov[a][b] * m3->H[a] * m3->H[b];
        }
    }
    m3->Hnorm2 = Hnorm2;
    if (Hnorm2 > 0.0) {
        m3->inv_Hnorm3 = 1.0 / (Hnorm2 * sqrt(Hnorm2));
    }
    xi = sqrt(Hnorm2) / J;
    chi = prj_rad_m1_chi_exact(xi);
    m3->thin_w = 0.5 * (3.0 * chi - 1.0);
    m3->thick_w = 1.5 * (1.0 - chi);
}

static void prj_rad_gr_m1_frequency_drifts_component_ref(
    const prj_rad_gr_m1_freq_drift_cache *cache,
    const prj_rad_grm1_m3_data *m3, double *energy_drift,
    double momentum_drift_cov[3])
{
    double Ycon[3];
    double Xcon[3][3];
    double energy;
    double momentum[3];
    int a;
    int i;
    int j;
    int k;
    int x;

    *energy_drift = 0.0;
    for (j = 0; j < 3; ++j) {
        momentum_drift_cov[j] = 0.0;
        momentum[j] = 0.0;
    }
    if (m3 == 0 || m3->J == 0.0) {
        return;
    }

    energy = cache->energy_Z * prj_rad_grm1_m3_component(m3, 0, 0, 0);
    for (x = 0; x < 3; ++x) {
        Ycon[x] = prj_rad_grm1_m3_component(m3, 0, 0, x + 1);
        energy += cache->energy_Y[x] * Ycon[x];
        for (i = 0; i < 3; ++i) {
            Xcon[x][i] = prj_rad_grm1_m3_component(m3, 0, x + 1, i + 1);
            energy += cache->energy_X[x][i] * Xcon[x][i];
        }
    }
    for (j = 0; j < 3; ++j) {
        for (x = 0; x < 3; ++x) {
            momentum[j] += cache->mom_Y[j][x] * Ycon[x];
            for (i = 0; i < 3; ++i) {
                momentum[j] += cache->mom_X[j][x][i] * Xcon[x][i];
            }
        }
    }
    for (a = 0; a < 3; ++a) {
        for (k = 0; k < 3; ++k) {
            for (i = 0; i < 3; ++i) {
                double Wcon = prj_rad_grm1_m3_component(m3, a + 1, k + 1,
                    i + 1);

                for (j = 0; j < 3; ++j) {
                    momentum[j] += cache->mom_W[j][a][k][i] * Wcon;
                }
            }
        }
    }
    *energy_drift = PRJ_CLIGHT * energy;
    for (j = 0; j < 3; ++j) {
        momentum_drift_cov[j] = PRJ_CLIGHT * PRJ_CLIGHT * momentum[j];
    }
}

static void prj_rad_gr_m1_frequency_drifts_cached(
    const prj_rad_gr_m1_freq_drift_cache *cache, double E_in,
    const double Fcov_in[3], double *energy_drift,
    double momentum_drift_cov[3])
{
    prj_rad_grm1_m3_data m3;
    double drift[4];
    int j;

    *energy_drift = 0.0;
    for (j = 0; j < 3; ++j) {
        momentum_drift_cov[j] = 0.0;
    }
    prj_rad_grm1_prepare_m3_data_freq(cache, E_in, Fcov_in, &m3);
    if (m3.J == 0.0) {
        return;
    }

    prj_rad_gr_m1_freq_m3_contract_all(cache, &m3, drift);
    *energy_drift = PRJ_CLIGHT * drift[0];
    for (j = 0; j < 3; ++j) {
        momentum_drift_cov[j] = PRJ_CLIGHT * PRJ_CLIGHT *
            drift[j + 1];
    }
}

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

static void prj_rad_grm1_zero_pressure(double P[3][3])
{
    int a;
    int b;

    if (P == 0) {
        return;
    }
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            P[a][b] = 0.0;
        }
    }
}

static void prj_rad_grm1_normal_metric_from_ctx(
    const prj_rad_gr_m1_closure_ctx *ctx, double g_cov[4][4],
    double g_con[4][4])
{
    int a;
    int b;

    memset(g_cov, 0, 16 * sizeof(double));
    memset(g_con, 0, 16 * sizeof(double));
    g_cov[0][0] = -1.0;
    g_con[0][0] = -1.0;
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            g_cov[a + 1][b + 1] = ctx->gamma[a][b];
            g_con[a + 1][b + 1] = ctx->gamma_inv[a][b];
        }
    }
}

static void prj_rad_grm1_ucon_from_side(
    const prj_rad_gr_m1_side_data *side, double ucon[4])
{
    int a;

    ucon[0] = side->wlor;
    for (a = 0; a < 3; ++a) {
        ucon[a + 1] = side->wlor * side->vcon[a];
    }
}

static int prj_rad_grm1_build_R_from_ctx(
    const prj_rad_gr_m1_closure_ctx *ctx, double E,
    const double Fcov_in[3], double g_cov[4][4], double g_con[4][4],
    double Rcon[4][4], double Fcon_out[3])
{
    double Fcov[3];
    double Fcon[3];
    double F2 = 0.0;
    double Fmag;
    double cE;
    int a;

    if (Fcon_out != 0) {
        for (a = 0; a < 3; ++a) {
            Fcon_out[a] = 0.0;
        }
    }
    if (ctx == 0 || Fcov_in == 0 || g_cov == 0 || g_con == 0 ||
        Rcon == 0) {
        return 0;
    }
    if (!isfinite(E) || E < 0.0) {
        E = 0.0;
    }
    for (a = 0; a < 3; ++a) {
        Fcov[a] = isfinite(Fcov_in[a]) ? Fcov_in[a] : 0.0;
    }
    prj_rad_gr_m1_raise_vec_ctx(ctx, Fcov, Fcon);
    for (a = 0; a < 3; ++a) {
        F2 += Fcov[a] * Fcon[a];
    }
    if (!isfinite(F2) || F2 < 0.0) {
        F2 = 0.0;
    }
    Fmag = sqrt(F2);
    cE = PRJ_CLIGHT * E;
    if (Fmag > cE && Fmag > 0.0) {
        double scale = cE / Fmag;

        for (a = 0; a < 3; ++a) {
            Fcon[a] *= scale;
        }
    }
    if (Fcon_out != 0) {
        for (a = 0; a < 3; ++a) {
            Fcon_out[a] = Fcon[a];
        }
    }
    prj_rad_grm1_normal_metric_from_ctx(ctx, g_cov, g_con);
    return prj_rad_grm1_build_R(g_cov, g_con, 1.0, E, Fcon, Rcon);
}

void prj_rad_gr_m1_prepare_side(const prj_rad_gr_m1_closure_ctx *ctx,
    prj_rad_gr_m1_side_data *side)
{
    double beta2;
    int a;

    if (side == 0) {
        return;
    }
    memset(side, 0, sizeof(*side));
    if (ctx == 0) {
        side->wlor = 1.0;
        return;
    }
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
 * prj_rad_gr_m1_prepare_side and reused, so only the E/F-dependent tensor work
 * runs per group. */
void prj_rad_gr_m1_pressure_cached(const prj_rad *rad,
    const prj_rad_gr_m1_closure_ctx *ctx,
    const prj_rad_gr_m1_side_data *side, double E, const double Fcov_in[3],
    double P[3][3])
{
    prj_rad_gr_m1_side_data local_side;

    (void)rad;
    if (P == 0) {
        return;
    }
    prj_rad_grm1_zero_pressure(P);
    if (ctx == 0 || Fcov_in == 0) {
        return;
    }
    if (side == 0) {
        prj_rad_gr_m1_prepare_side(ctx, &local_side);
        side = &local_side;
    }
    prj_rad_gr_m1_pressure_fbar_cached(rad, ctx, side, E, Fcov_in, P,
        0, 0, 0);
}

void prj_rad_gr_m1_pressure_fbar_cached(const prj_rad *rad,
    const prj_rad_gr_m1_closure_ctx *ctx,
    const prj_rad_gr_m1_side_data *side, double E, const double Fcov_in[3],
    double P[3][3], double *fbar_out, double *J0_out, double H0_out[3])
{
    prj_rad_gr_m1_side_data local_side;
    prj_rad_grm1_m3_data m3;
    double g_cov[4][4];
    double g_con[4][4];
    double Rcon[4][4];
    double ucon[4];
    int a;
    int b;

    (void)rad;
    if (fbar_out != 0) {
        *fbar_out = 0.0;
    }
    if (J0_out != 0) {
        *J0_out = 0.0;
    }
    if (H0_out != 0) {
        for (a = 0; a < 3; ++a) {
            H0_out[a] = 0.0;
        }
    }
    if (P == 0) {
        return;
    }
    prj_rad_grm1_zero_pressure(P);
    if (ctx == 0 || Fcov_in == 0) {
        return;
    }
    if (side == 0) {
        prj_rad_gr_m1_prepare_side(ctx, &local_side);
        side = &local_side;
    }
    if (!prj_rad_grm1_build_R_from_ctx(ctx, E, Fcov_in, g_cov, g_con,
            Rcon, 0)) {
        return;
    }
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            P[a][b] = Rcon[a + 1][b + 1];
        }
    }

    if (fbar_out == 0 && J0_out == 0 && H0_out == 0) {
        return;
    }
    prj_rad_grm1_ucon_from_side(side, ucon);
    if (!prj_rad_grm1_prepare_m3_data(g_cov, g_con, ucon, Rcon, &m3)) {
        return;
    }
    if (fbar_out != 0 && m3.J > 0.0) {
        double fbar = sqrt(m3.Hnorm2) / m3.J;

        if (!isfinite(fbar) || fbar < 0.0) {
            fbar = 0.0;
        } else if (fbar > 1.0) {
            fbar = 1.0;
        }
        *fbar_out = fbar;
    }
    if (J0_out != 0) {
        *J0_out = m3.J;
    }
    if (H0_out != 0) {
        for (a = 0; a < 3; ++a) {
            H0_out[a] = m3.H[a + 1];
        }
    }
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

        /* The analytic T^3 benchmark EOS starts cold and is heated
         * impulsively, so its Newton solve must cross temperature decades.
         * Preserve the production limiter for every other EOS. */
        step_scale = 1.0;
#if PRJ_EOS_T3
        const double max_temp_fraction = 1.0;
#else
        const double max_temp_fraction = 0.03;
#endif
        if (fabs(dT) > max_temp_fraction * T) {
            step_scale = max_temp_fraction * T / fabs(dT);
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
#if PRJ_PRESCRIBED_MATTER
    double eos_q[PRJ_EOS_NQUANT];
    double kappa[PRJ_NRAD * PRJ_NEGROUP];
    double eta[PRJ_NRAD * PRJ_NEGROUP];
    double rho = u[PRJ_CONS_RHO];
    double ye = u[PRJ_CONS_YE] / rho;
    double v1 = u[PRJ_CONS_MOM1] / rho;
    double v2 = u[PRJ_CONS_MOM2] / rho;
    double v3 = u[PRJ_CONS_MOM3] / rho;
    double eint = u[PRJ_CONS_ETOT] / rho -
        0.5 * (v1 * v1 + v2 * v2 + v3 * v3);
    int field, group;

    prj_eos_rey(eos, rho, eint, ye, eos_q, PRJ_EOS_CTX_MAIN);
    prj_rad3_opac_lookup(rad, rho, eos_q[PRJ_EOS_TEMPERATURE], ye,
        kappa, 0, 0, eta);
    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            int q = field * PRJ_NEGROUP + group;
            int e = PRJ_CONS_RAD_E(field, group);
            u[e] = (u[e] + dt * lapse * eta[q]) /
                (1.0 + dt * lapse * PRJ_CLIGHT * kappa[q]);
        }
    }
    if (final_temperature) *final_temperature = eos_q[PRJ_EOS_TEMPERATURE];
    if (kappa_out) memcpy(kappa_out, kappa, sizeof(kappa));
#else
    prj_rad_energy_update_impl(rad, eos, u, dt, lapse, final_temperature, 0, kappa_out, 0);
#endif
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
#if !PRJ_STATIC_MATTER && !PRJ_PRESCRIBED_MATTER
    double inv_c2;
    double dmom[3];
    double e_unchanged;
#endif
    int nu;
    int g;
    int d;

    (void)eos;

    rho = u[PRJ_CONS_RHO];
    Ye = u[PRJ_CONS_YE] / rho;
    prj_rad3_opac_lookup(rad, rho, temperature, Ye, 0, sigma, delta, 0);

#if !PRJ_STATIC_MATTER && !PRJ_PRESCRIBED_MATTER
    inv_c2 = 1.0 / (PRJ_CLIGHT * PRJ_CLIGHT);
    dmom[0] = 0.0;
    dmom[1] = 0.0;
    dmom[2] = 0.0;
    e_unchanged = u[PRJ_CONS_ETOT] - 0.5 * (u[PRJ_CONS_MOM1] * u[PRJ_CONS_MOM1] +
                                            u[PRJ_CONS_MOM2] * u[PRJ_CONS_MOM2] +
                                            u[PRJ_CONS_MOM3] * u[PRJ_CONS_MOM3]) / rho;
#endif

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
#if !PRJ_STATIC_MATTER && !PRJ_PRESCRIBED_MATTER
                dmom[d] += dF * inv_c2;
#endif
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

#if !PRJ_STATIC_MATTER && !PRJ_PRESCRIBED_MATTER
    /* dmom/detot accumulate radiation-flux changes in RAD_SCALE*erg units;
       multiply back to erg for the gas momentum/energy back-reaction. */
    u[PRJ_CONS_MOM1] -= dmom[0] * RAD_SCALE;
    u[PRJ_CONS_MOM2] -= dmom[1] * RAD_SCALE;
    u[PRJ_CONS_MOM3] -= dmom[2] * RAD_SCALE;
    u[PRJ_CONS_ETOT] = e_unchanged + 0.5 * (u[PRJ_CONS_MOM1] * u[PRJ_CONS_MOM1] +
                                            u[PRJ_CONS_MOM2] * u[PRJ_CONS_MOM2] +
                                            u[PRJ_CONS_MOM3] * u[PRJ_CONS_MOM3]) / rho;
#endif
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
static void prj_rad_gr_m1_fill_closure_ctx(const prj_z4c_hydro_geom *geom,
    const double vcon[3], const double dvcon_dx[3][3],
    prj_rad_gr_m1_closure_ctx *ctx)
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
    prj_rad_gr_m1_fill_closure_ctx(geom, vcon, dvcon_dx, ctx);
}

void prj_rad_freq_flux_apply_gr_m1_geom(const prj_rad *rad,
    const prj_mesh *mesh, const prj_block *block, int z4c_stage,
    const prj_z4c_hydro_geom *geom, const double *W_state, double *u,
    int ic, int jc, int kc, double dt,
    const double observer_time_derivative[4])
{
    double dvcon_dx[3][3];
    prj_rad_gr_m1_closure_ctx base_closure;
    prj_rad_gr_m1_side_data pside;
    prj_rad_gr_m1_freq_drift_cache drift_cache;
    const double *W_mhd;
    double dt_geom;
    int field;

    if (rad == 0 || mesh == 0 || block == 0 || geom == 0 || W_state == 0 ||
        u == 0) {
        return;
    }
    if (block->v_riemann[0] == 0 || block->v_riemann[1] == 0 ||
        block->v_riemann[2] == 0) {
        return;
    }
    W_mhd = prj_block_mhd_stage_const(block, z4c_stage);
    dt_geom = dt * geom->alpha * geom->sqrt_gamma;
    if (dt_geom == 0.0) {
        return;
    }

    prj_rad_gr_m1_cell_dvcon_dx(block, ic, jc, kc, dvcon_dx);
    prj_rad_gr_m1_freq_base_closure_ctx(geom, W_mhd, dvcon_dx, ic, jc, kc,
        &base_closure);
    prj_rad_gr_m1_prepare_side(&base_closure, &pside);
    prj_rad_gr_m1_freq_drift_cache_init(&base_closure, &pside, geom,
        observer_time_derivative, &drift_cache);

    for (field = 0; field < PRJ_NRAD; ++field) {
        double Acon_spec[PRJ_NEGROUP];
        double Mq_spec[PRJ_NEGROUP][3];
        double energy_face[PRJ_NEGROUP + 1] = {0.0};
        double momentum_face[PRJ_NEGROUP + 1][PRJ_NDIM] = {{0.0}};
        const double *nu_face = rad->eedge[field];
        int cached_grid = rad->freq_grid_valid[field] != 0;
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
            double Mq_cov[3];
            double Acon = 0.0;
            double inv_dnu;

            if (!cached_grid && dnu <= 0.0) {
                fprintf(stderr,
                    "prj_rad_freq_flux_apply_gr_m1: non-positive eedge width for field %d group %d\n",
                    field, g);
                exit(1);
            }
            inv_dnu = cached_grid ? rad->freq_inv_dnu[field][g] : 1.0 / dnu;
            Fcov[0] = W_state[WIDX(PRJ_RAD_PRIM_F1(field, g), ic, jc, kc)];
            Fcov[1] = W_state[WIDX(PRJ_RAD_PRIM_F2(field, g), ic, jc, kc)];
            Fcov[2] = W_state[WIDX(PRJ_RAD_PRIM_F3(field, g), ic, jc, kc)];

            prj_rad_gr_m1_frequency_drifts_cached(&drift_cache,
                W_state[WIDX(PRJ_RAD_PRIM_E(field, g), ic, jc, kc)],
                Fcov, &Acon, Mq_cov);
            Acon_spec[g] = Acon * inv_dnu;
            for (a = 0; a < 3; ++a) {
                Mq_spec[g][a] = Mq_cov[a] * inv_dnu;
            }
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
                double energy_available = u[PRJ_CONS_RAD_E(field, g)];

                if (drain > energy_available) {
                    theta[g] = energy_available > 0.0 && drain > 0.0
                        ? nextafter(energy_available / drain, 0.0)
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

void prj_rad_freq_flux_apply_gr_m1(const prj_rad *rad, const prj_mesh *mesh,
    const prj_block *block, int z4c_stage, const double *W_state, double *u,
    int ic, int jc, int kc, double dt,
    const double observer_time_derivative[4])
{
    prj_z4c_hydro_geom geom;

    if (rad == 0 || mesh == 0 || block == 0 || W_state == 0 || u == 0) {
        return;
    }
    if (!prj_z4c_load_hydro_geom(mesh, block, z4c_stage, ic, jc, kc, &geom)) {
        fprintf(stderr,
            "prj_rad_freq_flux_apply_gr_m1: failed to load Z4c geometry at cell (%d,%d,%d)\n",
            ic, jc, kc);
        exit(1);
    }
    prj_rad_freq_flux_apply_gr_m1_geom(rad, mesh, block, z4c_stage, &geom,
        W_state, u, ic, jc, kc, dt, observer_time_derivative);
}

#if PRJ_NRAD > 0
#define PRJ_RAD_GR_M1_NP (6 + 4 * PRJ_NRAD * PRJ_NEGROUP)
#define PRJ_RAD_GR_M1_NFLUID 6
#define PRJ_RAD_GR_M1_NRAD_BLOCK 4
#define PRJ_RAD_GR_M1_NGROUPS (PRJ_NRAD * PRJ_NEGROUP)
#define PRJ_RAD_GR_M1_NRAD_RHS (1 + PRJ_RAD_GR_M1_NFLUID)
#define PRJ_RAD_GR_M1_SOLVE_TOL_DEFAULT 1.0e-6
#define PRJ_RAD_GR_M1_SOLVE_MAXITER_DEFAULT 50
#define PRJ_RAD_GR_M1_SOLVE_VMIN 1.0e5
/* Radiation energy floor in stored code units; with RAD_SCALE this is
 * 1.0e-25 erg cm^-3 for the default scale. */
#define PRJ_RAD_GR_M1_SOLVE_EMIN PRJ_RAD_GR_M1_E_FLOOR
#define PRJ_RAD_GR_M1_LINESEARCH_MAX 16
#define PRJ_RAD_GR_M1_APPROX_MAXITER 8
/* Velocity gate for the Eulerian (u = n) fast candidate: cells with
 * (v/c)^2 below this threshold build the fast candidate in the Eulerian normal
 * frame, dropping the O(v/c) fluid boost.  The exact residual check that gates
 * acceptance still uses the true four-velocity, so this only affects candidate
 * quality, never correctness.  Constant for now (v/c < 0.01); to be replaced by
 * a tol-derived per-cell estimate once realistic residual data is available. */
#define PRJ_RAD_GR_M1_EULERIAN_V2_GATE 1.0e-4

typedef struct prj_rad_gr_m1_jac_blocks {
    double fluid[PRJ_RAD_GR_M1_NFLUID][PRJ_RAD_GR_M1_NFLUID];
    double fluid_rad[PRJ_RAD_GR_M1_NGROUPS][PRJ_RAD_GR_M1_NFLUID][PRJ_RAD_GR_M1_NRAD_BLOCK];
    double rad_fluid[PRJ_RAD_GR_M1_NGROUPS][PRJ_RAD_GR_M1_NRAD_BLOCK][PRJ_RAD_GR_M1_NFLUID];
    double rad_rad[PRJ_RAD_GR_M1_NGROUPS][PRJ_RAD_GR_M1_NRAD_BLOCK][PRJ_RAD_GR_M1_NRAD_BLOCK];
} prj_rad_gr_m1_jac_blocks;

typedef struct prj_rad_gr_m1_frozen_group {
    double dRcon[PRJ_RAD_GR_M1_NRAD_BLOCK][4][4];
} prj_rad_gr_m1_frozen_group;

typedef prj_rad3_opac_interp_result prj_rad_gr_m1_opac_cache;

typedef struct prj_rad_gr_m1_cell_context {
    double g_cov[4][4];
    double g_con[4][4];
    double n_cov[4];
    double ucon[4];
    double ucov[4];
    prj_eos_gr_geom eos_geom;
} prj_rad_gr_m1_cell_context;

typedef struct prj_rad_gr_m1_solve_diag {
    int iter;
    int maxiter;
    int line_search_trials;
    int invalid_trials;
    int best_ls;
    int block_solve_ok;
    int dense_solve_ok;
    int directional_check_ok;
    int directional_row;
    int fd_check_ok;
    int fd_skipped_cols;
    int fd_max_abs_row;
    int fd_max_abs_col;
    int fd_max_rel_row;
    int fd_max_rel_col;
    int dP_max_abs_col;
    int dP_max_rel_col;
    int linear_min_pivot_group;
    int linear_min_pivot_k;
    double threshold;
    double norm;
    double first_trial_norm;
    double best_trial_norm;
    double last_trial_norm;
    double best_lambda;
    double dP_max_abs;
    double dP_max_rel;
    double linear_min_pivot;
    double linear_min_pivot_rel;
    double directional_eps;
    double directional_max_abs;
    double directional_max_rel;
    double directional_jdp_norm;
    double directional_fd_norm;
    double fd_max_abs;
    double fd_max_abs_got;
    double fd_max_abs_expected;
    double fd_max_abs_h;
    double fd_max_rel;
    double fd_max_rel_got;
    double fd_max_rel_expected;
    double fd_max_rel_h;
} prj_rad_gr_m1_solve_diag;

static prj_rad_gr_m1_solve_diag prj_rad_gr_m1_last_solve_diag;

static void prj_rad_gr_m1_solve_diag_reset(double threshold, int maxiter)
{
    memset(&prj_rad_gr_m1_last_solve_diag, 0,
        sizeof(prj_rad_gr_m1_last_solve_diag));
    prj_rad_gr_m1_last_solve_diag.iter = -1;
    prj_rad_gr_m1_last_solve_diag.maxiter = maxiter;
    prj_rad_gr_m1_last_solve_diag.best_ls = -1;
    prj_rad_gr_m1_last_solve_diag.directional_row = -1;
    prj_rad_gr_m1_last_solve_diag.fd_max_abs_row = -1;
    prj_rad_gr_m1_last_solve_diag.fd_max_abs_col = -1;
    prj_rad_gr_m1_last_solve_diag.fd_max_rel_row = -1;
    prj_rad_gr_m1_last_solve_diag.fd_max_rel_col = -1;
    prj_rad_gr_m1_last_solve_diag.dP_max_abs_col = -1;
    prj_rad_gr_m1_last_solve_diag.dP_max_rel_col = -1;
    prj_rad_gr_m1_last_solve_diag.linear_min_pivot_group = -1;
    prj_rad_gr_m1_last_solve_diag.linear_min_pivot_k = -1;
    prj_rad_gr_m1_last_solve_diag.threshold = threshold;
    prj_rad_gr_m1_last_solve_diag.norm = HUGE_VAL;
    prj_rad_gr_m1_last_solve_diag.first_trial_norm = HUGE_VAL;
    prj_rad_gr_m1_last_solve_diag.best_trial_norm = HUGE_VAL;
    prj_rad_gr_m1_last_solve_diag.last_trial_norm = HUGE_VAL;
    prj_rad_gr_m1_last_solve_diag.best_lambda = 0.0;
    prj_rad_gr_m1_last_solve_diag.directional_eps = 0.0;
    prj_rad_gr_m1_last_solve_diag.directional_max_abs = HUGE_VAL;
    prj_rad_gr_m1_last_solve_diag.directional_max_rel = HUGE_VAL;
    prj_rad_gr_m1_last_solve_diag.directional_jdp_norm = HUGE_VAL;
    prj_rad_gr_m1_last_solve_diag.directional_fd_norm = HUGE_VAL;
    prj_rad_gr_m1_last_solve_diag.fd_max_abs = 0.0;
    prj_rad_gr_m1_last_solve_diag.fd_max_abs_got = 0.0;
    prj_rad_gr_m1_last_solve_diag.fd_max_abs_expected = 0.0;
    prj_rad_gr_m1_last_solve_diag.fd_max_abs_h = 0.0;
    prj_rad_gr_m1_last_solve_diag.fd_max_rel = 0.0;
    prj_rad_gr_m1_last_solve_diag.fd_max_rel_got = 0.0;
    prj_rad_gr_m1_last_solve_diag.fd_max_rel_expected = 0.0;
    prj_rad_gr_m1_last_solve_diag.fd_max_rel_h = 0.0;
}

static void prj_rad_gr_m1_metric4_from_geom(const prj_z4c_hydro_geom *geom,
    double g_cov[4][4], double g_con[4][4])
{
    double beta_cov[3] = {0.0, 0.0, 0.0};
    double beta2 = 0.0;
    double inv_alpha2;
    int a;
    int b;

    memset(g_cov, 0, 16 * sizeof(double));
    memset(g_con, 0, 16 * sizeof(double));
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            beta_cov[a] += geom->gamma[a][b] * geom->beta[b];
        }
        beta2 += beta_cov[a] * geom->beta[a];
    }

    g_cov[0][0] = -geom->alpha * geom->alpha + beta2;
    for (a = 0; a < 3; ++a) {
        g_cov[0][a + 1] = beta_cov[a];
        g_cov[a + 1][0] = beta_cov[a];
        for (b = 0; b < 3; ++b) {
            g_cov[a + 1][b + 1] = geom->gamma[a][b];
        }
    }

    inv_alpha2 = 1.0 / (geom->alpha * geom->alpha);
    g_con[0][0] = -inv_alpha2;
    for (a = 0; a < 3; ++a) {
        g_con[0][a + 1] = geom->beta[a] * inv_alpha2;
        g_con[a + 1][0] = g_con[0][a + 1];
        for (b = 0; b < 3; ++b) {
            g_con[a + 1][b + 1] = geom->gamma_inv[a][b] -
                geom->beta[a] * geom->beta[b] * inv_alpha2;
        }
    }
}

static void prj_rad_gr_m1_lower4(const double g_cov[4][4],
    const double ucon[4], double ucov[4])
{
    int a;
    int b;

    for (a = 0; a < 4; ++a) {
        ucov[a] = 0.0;
        for (b = 0; b < 4; ++b) {
            ucov[a] += g_cov[a][b] * ucon[b];
        }
    }
}

static int prj_rad_gr_m1_fluid_four_velocity(
    const prj_z4c_hydro_geom *geom, const double g_cov[4][4],
    const double *P, double ucon[4], double ucov[4])
{
    double vhat[3];
    double v2 = 0.0;
    double wlor;
    int a;
    int b;

    if (geom->alpha <= 0.0 || !isfinite(geom->alpha)) {
        return 0;
    }
    for (a = 0; a < 3; ++a) {
        vhat[a] = P[1 + a] / PRJ_CLIGHT;
        if (!isfinite(vhat[a])) {
            return 0;
        }
    }
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            v2 += geom->gamma[a][b] * vhat[a] * vhat[b];
        }
    }
    if (!isfinite(v2) || v2 < 0.0 || v2 >= 1.0) {
        return 0;
    }
    wlor = 1.0 / sqrt(1.0 - v2);
    ucon[0] = wlor / geom->alpha;
    for (a = 0; a < 3; ++a) {
        ucon[a + 1] = wlor * (vhat[a] - geom->beta[a] / geom->alpha);
    }
    prj_rad_gr_m1_lower4(g_cov, ucon, ucov);
    return 1;
}

static int prj_rad_gr_m1_cell_context_init(
    const prj_z4c_hydro_geom *geom, const double *P,
    prj_rad_gr_m1_cell_context *ctx)
{
    int a;
    int b;

    if (geom == 0 || P == 0 || ctx == 0 ||
        !isfinite(geom->alpha) || geom->alpha <= 0.0 ||
        !isfinite(geom->sqrt_gamma) || geom->sqrt_gamma <= 0.0) {
        return 0;
    }
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            if (!isfinite(geom->gamma[a][b]) ||
                !isfinite(geom->gamma_inv[a][b])) {
                return 0;
            }
            ctx->eos_geom.gamma[a][b] = geom->gamma[a][b];
        }
    }
    prj_rad_gr_m1_metric4_from_geom(geom, ctx->g_cov, ctx->g_con);
    ctx->n_cov[0] = -geom->alpha;
    ctx->n_cov[1] = 0.0;
    ctx->n_cov[2] = 0.0;
    ctx->n_cov[3] = 0.0;
    return prj_rad_gr_m1_fluid_four_velocity(geom, ctx->g_cov, P,
        ctx->ucon, ctx->ucov);
}

static int prj_rad_gr_m1_hydro_cons_cached(
    const prj_z4c_hydro_geom *geom, const double *u_old, const double *P,
    double eint, double pressure, double *u_new)
{
    double beta_con[3];
    double beta_cov[3] = {0.0, 0.0, 0.0};
    double beta2 = 0.0;
    double sqrt_one_minus_beta2;
    double wlor;
    double wlor2;
    double wlor_m1;
    double Bcon[3] = {0.0, 0.0, 0.0};
    double Bcov[3] = {0.0, 0.0, 0.0};
    double Bsq = 0.0;
    double Bbeta = 0.0;
    double w;
    double sqrtg;
    double c = PRJ_CLIGHT;
    double c2 = PRJ_CLIGHT * PRJ_CLIGHT;
    int a;
    int b;

    if (geom == 0 || u_old == 0 || P == 0 || u_new == 0 ||
        !isfinite(geom->sqrt_gamma) || geom->sqrt_gamma <= 0.0 ||
        !isfinite(P[0]) || P[0] <= 0.0 || !isfinite(eint) || eint < 0.0 ||
        !isfinite(pressure) || !isfinite(P[5])) {
        return 0;
    }
    sqrtg = geom->sqrt_gamma;
    for (a = 0; a < 3; ++a) {
        beta_con[a] = P[1 + a] / c;
        if (!isfinite(beta_con[a])) {
            return 0;
        }
    }
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            beta_cov[a] += geom->gamma[a][b] * beta_con[b];
        }
        beta2 += beta_cov[a] * beta_con[a];
    }
    if (!isfinite(beta2) || beta2 < 0.0 || beta2 >= 1.0) {
        return 0;
    }
    sqrt_one_minus_beta2 = sqrt(1.0 - beta2);
    wlor = 1.0 / sqrt_one_minus_beta2;
    wlor2 = wlor * wlor;
    wlor_m1 = beta2 /
        (sqrt_one_minus_beta2 * (1.0 + sqrt_one_minus_beta2));
#if PRJ_MHD
    Bcon[0] = u_old[PRJ_CONS_B1] / sqrtg;
    Bcon[1] = u_old[PRJ_CONS_B2] / sqrtg;
    Bcon[2] = u_old[PRJ_CONS_B3] / sqrtg;
#endif
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            Bcov[a] += geom->gamma[a][b] * Bcon[b];
        }
        Bsq += Bcov[a] * Bcon[a];
        Bbeta += Bcov[a] * beta_con[a];
    }
    if (!isfinite(Bsq) || !isfinite(Bbeta) || Bsq < 0.0) {
        return 0;
    }
    w = P[0] * c2 + P[0] * eint + pressure;
    u_new[PRJ_CONS_RHO] = sqrtg * P[0] * wlor;
    for (a = 0; a < 3; ++a) {
        u_new[PRJ_CONS_MOM1 + a] = sqrtg *
            ((w * wlor2 + Bsq) * beta_cov[a] - Bbeta * Bcov[a]) / c;
    }
    u_new[PRJ_CONS_ETOT] = sqrtg *
        ((P[0] * eint + pressure) * wlor2 +
         P[0] * c2 * wlor * wlor_m1 - pressure + Bsq -
         0.5 * (Bbeta * Bbeta + Bsq / wlor2));
    u_new[PRJ_CONS_YE] = u_new[PRJ_CONS_RHO] * P[5];
#if PRJ_MHD
    u_new[PRJ_CONS_B1] = sqrtg * Bcon[0];
    u_new[PRJ_CONS_B2] = sqrtg * Bcon[1];
    u_new[PRJ_CONS_B3] = sqrtg * Bcon[2];
#endif
    return isfinite(u_new[PRJ_CONS_RHO]) &&
        isfinite(u_new[PRJ_CONS_MOM1]) &&
        isfinite(u_new[PRJ_CONS_MOM2]) &&
        isfinite(u_new[PRJ_CONS_MOM3]) &&
        isfinite(u_new[PRJ_CONS_ETOT]) && isfinite(u_new[PRJ_CONS_YE]);
}

static int prj_rad_gr_m1_p_to_prim(prj_eos *eos,
    const prj_z4c_hydro_geom *geom, const double g_cov[4][4],
    const double g_con[4][4], const double *u_old, const double *P,
    double *W, double Rcon_out[][4][4])
{
    double eos_q[PRJ_EOS_NQUANT];
    int field;
    int group;
    int a;
    int b;
    int v;

    if (eos == 0 || geom == 0 || g_cov == 0 || g_con == 0 ||
        u_old == 0 || P == 0 || W == 0 ||
        !isfinite(geom->alpha) || geom->alpha <= 0.0 ||
        !isfinite(geom->sqrt_gamma) || geom->sqrt_gamma <= 0.0 ||
        !isfinite(P[0]) || P[0] <= 0.0 ||
        !isfinite(P[4]) || P[4] <= 0.0 || !isfinite(P[5])) {
        return 0;
    }

    for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
        W[v] = 0.0;
    }
    W[PRJ_PRIM_RHO] = P[0];
    W[PRJ_PRIM_V1] = P[1];
    W[PRJ_PRIM_V2] = P[2];
    W[PRJ_PRIM_V3] = P[3];
    W[PRJ_PRIM_YE] = P[5];
    prj_eos_rty(eos, P[0], P[4], P[5], eos_q, PRJ_EOS_CTX_MAIN);
    W[PRJ_PRIM_EINT] = eos_q[PRJ_EOS_EINT];
    if (!isfinite(W[PRJ_PRIM_EINT]) || W[PRJ_PRIM_EINT] < 0.0) {
        return 0;
    }
#if PRJ_MHD
    W[PRJ_PRIM_B1] = u_old[PRJ_CONS_B1] / geom->sqrt_gamma;
    W[PRJ_PRIM_B2] = u_old[PRJ_CONS_B2] / geom->sqrt_gamma;
    W[PRJ_PRIM_B3] = u_old[PRJ_CONS_B3] / geom->sqrt_gamma;
#endif

    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            int idx = field * PRJ_NEGROUP + group;
            int pidx = 6 + 4 * idx;
            double E = P[pidx];
            double Fcov[3];
            double Fcon[3];
            double Rcon[4][4];

            /* Newton variables are the lab-frame (Eulerian) moments
             * p = (E, F_i); F_i is the covariant flux stored in the primitive
             * list.  The conserved rad rows collapse to sqrt(gamma) (E, F_i),
             * so W carries E and F_i directly.  R^{ab} is still needed for the
             * radiation four-force, so build it from the exact algebraic M1
             * closure and export it for the caller (the residual). */
            if (!isfinite(E) || E < 0.0) {
                return 0;
            }
            Fcov[0] = P[pidx + 1];
            Fcov[1] = P[pidx + 2];
            Fcov[2] = P[pidx + 3];
            for (a = 0; a < 3; ++a) {
                if (!isfinite(Fcov[a])) {
                    return 0;
                }
                Fcon[a] = 0.0;
                for (b = 0; b < 3; ++b) {
                    Fcon[a] += geom->gamma_inv[a][b] * Fcov[b];
                }
            }
            if (!prj_rad_grm1_build_R(g_cov, g_con, geom->alpha, E, Fcon,
                    Rcon)) {
                return 0;
            }
            if (Rcon_out != 0) {
                for (a = 0; a < 4; ++a) {
                    for (b = 0; b < 4; ++b) {
                        Rcon_out[idx][a][b] = Rcon[a][b];
                    }
                }
            }

            W[PRJ_PRIM_RAD_E(field, group)] = E;
            W[PRJ_PRIM_RAD_F1(field, group)] = Fcov[0];
            W[PRJ_PRIM_RAD_F2(field, group)] = Fcov[1];
            W[PRJ_PRIM_RAD_F3(field, group)] = Fcov[2];
        }
    }
    for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
        if (!isfinite(W[v])) {
            return 0;
        }
    }
    return 1;
}

static int prj_rad_gr_m1_fluid_row_from_index(int idx)
{
    if (idx == 0) {
        return PRJ_CONS_RHO;
    }
    if (idx >= 1 && idx <= 3) {
        return PRJ_CONS_MOM1 + (idx - 1);
    }
    if (idx == 4) {
        return PRJ_CONS_ETOT;
    }
    return PRJ_CONS_YE;
}

static void prj_rad_gr_m1_jac_blocks_zero(prj_rad_gr_m1_jac_blocks *blocks)
{
    if (blocks != 0) {
        memset(blocks, 0, sizeof(*blocks));
    }
}

static void prj_rad_gr_m1_jac_blocks_to_dense(
    const prj_rad_gr_m1_jac_blocks *blocks, double *dense)
{
    const int np = PRJ_RAD_GR_M1_NP;
    int frow;
    int fcol;
    int gidx;
    int r;
    int c;
    int v;

    if (blocks == 0 || dense == 0) {
        return;
    }
    for (v = 0; v < PRJ_NVAR_CONS * PRJ_RAD_GR_M1_NP; ++v) {
        dense[v] = 0.0;
    }
    for (frow = 0; frow < PRJ_RAD_GR_M1_NFLUID; ++frow) {
        int row = prj_rad_gr_m1_fluid_row_from_index(frow);

        for (fcol = 0; fcol < PRJ_RAD_GR_M1_NFLUID; ++fcol) {
            dense[(size_t)row * (size_t)np + (size_t)fcol] =
                blocks->fluid[frow][fcol];
        }
    }
    for (gidx = 0; gidx < PRJ_RAD_GR_M1_NGROUPS; ++gidx) {
        int field = gidx / PRJ_NEGROUP;
        int group = gidx % PRJ_NEGROUP;
        int pidx = PRJ_RAD_GR_M1_NFLUID +
            PRJ_RAD_GR_M1_NRAD_BLOCK * gidx;
        int rad_row[PRJ_RAD_GR_M1_NRAD_BLOCK];

        rad_row[0] = PRJ_CONS_RAD_E(field, group);
        rad_row[1] = PRJ_CONS_RAD_F1(field, group);
        rad_row[2] = PRJ_CONS_RAD_F2(field, group);
        rad_row[3] = PRJ_CONS_RAD_F3(field, group);

        for (frow = 0; frow < PRJ_RAD_GR_M1_NFLUID; ++frow) {
            int row = prj_rad_gr_m1_fluid_row_from_index(frow);

            for (r = 0; r < PRJ_RAD_GR_M1_NRAD_BLOCK; ++r) {
                dense[(size_t)row * (size_t)np + (size_t)(pidx + r)] =
                    blocks->fluid_rad[gidx][frow][r];
            }
        }
        for (r = 0; r < PRJ_RAD_GR_M1_NRAD_BLOCK; ++r) {
            for (fcol = 0; fcol < PRJ_RAD_GR_M1_NFLUID; ++fcol) {
                dense[(size_t)rad_row[r] * (size_t)np + (size_t)fcol] =
                    blocks->rad_fluid[gidx][r][fcol];
            }
            for (c = 0; c < PRJ_RAD_GR_M1_NRAD_BLOCK; ++c) {
                dense[(size_t)rad_row[r] * (size_t)np +
                    (size_t)(pidx + c)] = blocks->rad_rad[gidx][r][c];
            }
        }
    }
}

static int prj_rad_gr_m1_freeze_group(const double g_cov[4][4],
    const double g_con[4][4], const prj_rad_grm1_R_jac_geom *R_jac_geom,
    const double p[PRJ_RAD_GR_M1_NRAD_BLOCK],
    prj_rad_gr_m1_frozen_group *frozen)
{
    double Rcon[4][4];
    double Fcov[3];
    int a;
    int b;
    int cc;

    if (g_cov == 0 || g_con == 0 || R_jac_geom == 0 || p == 0 ||
        frozen == 0 || !isfinite(p[0]) || p[0] < 0.0) {
        return 0;
    }
    Fcov[0] = p[1];
    Fcov[1] = p[2];
    Fcov[2] = p[3];
    for (a = 0; a < 3; ++a) {
        if (!isfinite(Fcov[a])) {
            return 0;
        }
    }
    if (!prj_rad_grm1_build_R_jac(g_cov, g_con, R_jac_geom, p[0], Fcov,
            Rcon, frozen->dRcon)) {
        return 0;
    }
    for (cc = 0; cc < PRJ_RAD_GR_M1_NRAD_BLOCK; ++cc) {
        for (a = 0; a < 4; ++a) {
            for (b = 0; b < 4; ++b) {
                if (!isfinite(frozen->dRcon[cc][a][b])) {
                    return 0;
                }
            }
        }
    }
    return 1;
}

static int prj_rad_gr_m1_freeze_closure(const prj_z4c_hydro_geom *geom,
    const double *P, prj_rad_gr_m1_frozen_group frozen[PRJ_RAD_GR_M1_NGROUPS])
{
    double g_cov[4][4];
    double g_con[4][4];
    prj_rad_grm1_R_jac_geom R_jac_geom;
    int gidx;

    if (geom == 0 || P == 0 || frozen == 0 ||
        !isfinite(geom->alpha) || geom->alpha <= 0.0) {
        return 0;
    }
    prj_rad_gr_m1_metric4_from_geom(geom, g_cov, g_con);
    prj_rad_grm1_R_jac_geom_init(geom->alpha, g_con, geom->gamma_inv,
        &R_jac_geom);
    for (gidx = 0; gidx < PRJ_RAD_GR_M1_NGROUPS; ++gidx) {
        int pidx = PRJ_RAD_GR_M1_NFLUID +
            PRJ_RAD_GR_M1_NRAD_BLOCK * gidx;

        if (!prj_rad_gr_m1_freeze_group(g_cov, g_con, &R_jac_geom,
                &P[pidx], &frozen[gidx])) {
            return 0;
        }
    }
    return 1;
}

static int prj_rad_gr_m1_frozen_R(
    const prj_rad_gr_m1_frozen_group *frozen, const double p[4],
    double Rcon[4][4])
{
    int a;
    int b;
    int cc;

    if (frozen == 0 || p == 0 || Rcon == 0 || !isfinite(p[0]) ||
        p[0] < 0.0) {
        return 0;
    }
    for (cc = 0; cc < PRJ_RAD_GR_M1_NRAD_BLOCK; ++cc) {
        if (!isfinite(p[cc])) {
            return 0;
        }
    }
    for (a = 0; a < 4; ++a) {
        for (b = 0; b < 4; ++b) {
            Rcon[a][b] = 0.0;
            for (cc = 0; cc < PRJ_RAD_GR_M1_NRAD_BLOCK; ++cc) {
                Rcon[a][b] += frozen->dRcon[cc][a][b] * p[cc];
            }
            if (!isfinite(Rcon[a][b])) {
                return 0;
            }
        }
    }
    return 1;
}

static int prj_rad_gr_m1_residual(const prj_rad *rad, prj_eos *eos,
    const prj_z4c_hydro_geom *geom, const double *u_old, const double *P,
    double dt, double *resid, double *u_new_out, double *W_out)
{
    double W[PRJ_NVAR_PRIM];
    double u_new[PRJ_NVAR_CONS];
    double resid_tmp[PRJ_NVAR_CONS];
    double Rcon_store[PRJ_NRAD * PRJ_NEGROUP][4][4];
    prj_eos_gr_geom eos_geom;
    double g_cov[4][4];
    double g_con[4][4];
    double n_cov[4];
    double ucon[4];
    double ucov[4];
    double kappa[PRJ_NRAD * PRJ_NEGROUP];
    double sigma[PRJ_NRAD * PRJ_NEGROUP];
    double delta[PRJ_NRAD * PRJ_NEGROUP];
    double eta[PRJ_NRAD * PRJ_NEGROUP];
    double sum_Gn = 0.0;
    double sum_Gu_xe = 0.0;
    double sum_Ggamma[3] = {0.0, 0.0, 0.0};
    double group_Gn[PRJ_NRAD * PRJ_NEGROUP];
    double group_Ggamma[PRJ_NRAD * PRJ_NEGROUP][3];
    int field;
    int group;
    int a;
    int b;
    int d;
    int v;
    int ok;

    if (resid != 0) {
        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
            resid[v] = 0.0;
        }
    }
    if (u_new_out != 0) {
        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
            u_new_out[v] = 0.0;
        }
    }
    if (W_out != 0) {
        for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
            W_out[v] = 0.0;
        }
    }
    if (rad == 0 || eos == 0 || geom == 0 || u_old == 0 || P == 0 ||
        resid == 0 || !isfinite(dt) || !isfinite(geom->sqrt_gamma) ||
        geom->sqrt_gamma <= 0.0 || !isfinite(geom->alpha) ||
        geom->alpha <= 0.0) {
        return 0;
    }
    if (!isfinite(P[0]) || P[0] <= 0.0 || !isfinite(P[4]) || P[4] <= 0.0 ||
        !isfinite(P[5])) {
        return 0;
    }
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            if (!isfinite(geom->gamma[a][b]) ||
                !isfinite(geom->gamma_inv[a][b])) {
                return 0;
            }
            eos_geom.gamma[a][b] = geom->gamma[a][b];
        }
    }

    prj_rad_gr_m1_metric4_from_geom(geom, g_cov, g_con);
    n_cov[0] = -geom->alpha;
    n_cov[1] = 0.0;
    n_cov[2] = 0.0;
    n_cov[3] = 0.0;
    if (!prj_rad_gr_m1_fluid_four_velocity(geom, g_cov, P, ucon, ucov)) {
        return 0;
    }
    /* TEMP TIMER: remove after rad-matter coupling profiling. */
    PRJ_TIMER_CURRENT_START("rad_matter_temp_resid_p_to_prim");
    ok = prj_rad_gr_m1_p_to_prim(eos, geom, g_cov, g_con, u_old, P, W,
        Rcon_store);
    PRJ_TIMER_CURRENT_STOP("rad_matter_temp_resid_p_to_prim");
    if (!ok) {
        return 0;
    }

    /* TEMP TIMER: remove after rad-matter coupling profiling. */
    PRJ_TIMER_CURRENT_START("rad_matter_temp_resid_opac");
    prj_rad3_opac_lookup(rad, P[0], P[4], P[5], kappa, sigma, delta, eta);
    PRJ_TIMER_CURRENT_STOP("rad_matter_temp_resid_opac");
    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            int idx = field * PRJ_NEGROUP + group;
            /* Reuse the R^{ab} that p_to_prim already built from the same
             * lab-frame moments and geometry. */
            const double (*Rcon)[4] = Rcon_store[idx];
            double R_u[4];
            double Ruu = 0.0;
            double Gcon[4];
            double kappa_eff;
            double eta_eff;
            double sigma_eff;
            double Gdotn = 0.0;
            double Gdotu = 0.0;

            for (a = 0; a < 4; ++a) {
                R_u[a] = 0.0;
                for (b = 0; b < 4; ++b) {
                    R_u[a] += Rcon[a][b] * ucov[b];
                }
                Ruu += R_u[a] * ucov[a];
            }
            kappa_eff = kappa[idx] + rad->gr_m1_eleinel_dkappa[idx];
            eta_eff = eta[idx] + rad->gr_m1_eleinel_deta[idx];
            sigma_eff = sigma[idx] * (1.0 - delta[idx] / 3.0);
            if (!isfinite(kappa_eff) || !isfinite(sigma_eff) ||
                !isfinite(eta_eff)) {
                return 0;
            }
            for (a = 0; a < 4; ++a) {
                Gcon[a] = -(kappa_eff + sigma_eff) * R_u[a] -
                    (sigma_eff * Ruu + eta_eff / PRJ_CLIGHT) * ucon[a];
                if (!isfinite(Gcon[a])) {
                    return 0;
                }
                Gdotn += Gcon[a] * n_cov[a];
                Gdotu += Gcon[a] * ucov[a];
            }

            group_Gn[idx] = Gdotn;
            sum_Gn += Gdotn;
            sum_Gu_xe += Gdotu * rad->x_e[field][group];
            for (d = 0; d < 3; ++d) {
                double Ggamma = 0.0;

                /* n_cov[d+1] == 0 (spatial), so the n_cov n_cov term drops. */
                for (a = 0; a < 4; ++a) {
                    Ggamma += Gcon[a] * g_cov[d + 1][a];
                }
                group_Ggamma[idx][d] = Ggamma;
                sum_Ggamma[d] += Ggamma;
            }
        }
    }

    /* TEMP TIMER: remove after rad-matter coupling profiling. */
    PRJ_TIMER_CURRENT_START("rad_matter_temp_resid_prim2cons");
    ok = prj_eos_gr_prim2cons(eos, &eos_geom, W, u_new,
        PRJ_EOS_CTX_MAIN) == PRJ_EOS_GR_OK;
    PRJ_TIMER_CURRENT_STOP("rad_matter_temp_resid_prim2cons");
    if (!ok) {
        return 0;
    }

    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        resid_tmp[v] = u_new[v] - u_old[v];
    }
    resid_tmp[PRJ_CONS_ETOT] += geom->sqrt_gamma * dt *
        RAD_SCALE * PRJ_CLIGHT * sum_Gn;
    resid_tmp[PRJ_CONS_MOM1] -= geom->alpha * geom->sqrt_gamma * dt *
        RAD_SCALE * sum_Ggamma[0];
    resid_tmp[PRJ_CONS_MOM2] -= geom->alpha * geom->sqrt_gamma * dt *
        RAD_SCALE * sum_Ggamma[1];
    resid_tmp[PRJ_CONS_MOM3] -= geom->alpha * geom->sqrt_gamma * dt *
        RAD_SCALE * sum_Ggamma[2];
    resid_tmp[PRJ_CONS_YE] += geom->alpha * geom->sqrt_gamma * dt *
        PRJ_CLIGHT * sum_Gu_xe;
    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            int idx = field * PRJ_NEGROUP + group;

            resid_tmp[PRJ_CONS_RAD_E(field, group)] -=
                geom->alpha * geom->sqrt_gamma * dt *
                PRJ_CLIGHT * group_Gn[idx];
            resid_tmp[PRJ_CONS_RAD_F1(field, group)] +=
                geom->alpha * geom->sqrt_gamma * dt *
                PRJ_CLIGHT * PRJ_CLIGHT * group_Ggamma[idx][0];
            resid_tmp[PRJ_CONS_RAD_F2(field, group)] +=
                geom->alpha * geom->sqrt_gamma * dt *
                PRJ_CLIGHT * PRJ_CLIGHT * group_Ggamma[idx][1];
            resid_tmp[PRJ_CONS_RAD_F3(field, group)] +=
                geom->alpha * geom->sqrt_gamma * dt *
                PRJ_CLIGHT * PRJ_CLIGHT * group_Ggamma[idx][2];
        }
    }

    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        resid[v] = resid_tmp[v];
        if (u_new_out != 0) {
            u_new_out[v] = u_new[v];
        }
    }
    if (W_out != 0) {
        memcpy(W_out, W, sizeof(W));
    }
    return 1;
}

static int prj_rad_gr_m1_residual_jacobian(const prj_rad *rad, prj_eos *eos,
    const prj_z4c_hydro_geom *geom, const double *u_old, const double *P,
    double dt, double *resid, prj_rad_gr_m1_jac_blocks *blocks, double *jac,
    double *u_new_out, double *W_out)
{
    double resid_local[PRJ_NVAR_CONS];
    double *resid_use = resid != 0 ? resid : resid_local;
    double W[PRJ_NVAR_PRIM];
    double u_new[PRJ_NVAR_CONS];
    double resid_tmp[PRJ_NVAR_CONS];
    prj_eos_gr_geom eos_geom;
    double g_cov[4][4];
    double g_con[4][4];
    prj_rad_grm1_R_jac_geom R_jac_geom;
    double n_cov[4];
    double ucon[4];
    double ucov[4];
    double ducon[4][PRJ_RAD_GR_M1_NFLUID];
    double ducov[4][PRJ_RAD_GR_M1_NFLUID];
    double kappa[PRJ_NRAD * PRJ_NEGROUP];
    double sigma[PRJ_NRAD * PRJ_NEGROUP];
    double delta[PRJ_NRAD * PRJ_NEGROUP];
    double eta[PRJ_NRAD * PRJ_NEGROUP];
    double dkappa_drho[PRJ_NRAD * PRJ_NEGROUP];
    double dkappa_dT[PRJ_NRAD * PRJ_NEGROUP];
    double dkappa_dYe[PRJ_NRAD * PRJ_NEGROUP];
    double dsigma_drho[PRJ_NRAD * PRJ_NEGROUP];
    double dsigma_dT[PRJ_NRAD * PRJ_NEGROUP];
    double dsigma_dYe[PRJ_NRAD * PRJ_NEGROUP];
    double ddelta_drho[PRJ_NRAD * PRJ_NEGROUP];
    double ddelta_dT[PRJ_NRAD * PRJ_NEGROUP];
    double ddelta_dYe[PRJ_NRAD * PRJ_NEGROUP];
    double deta_drho[PRJ_NRAD * PRJ_NEGROUP];
    double deta_dT[PRJ_NRAD * PRJ_NEGROUP];
    double deta_dYe[PRJ_NRAD * PRJ_NEGROUP];
    double eint;
    double pressure;
    double deint_drho;
    double deint_dT;
    double deint_dYe;
    double dpressure_drho;
    double dpressure_dT;
    double dpressure_dYe;
    double beta_con[3];
    double beta_cov[3] = {0.0, 0.0, 0.0};
    double beta2 = 0.0;
    double wlor;
    double wlor2;
    double wlor_m1;
    double Bcon[3] = {0.0, 0.0, 0.0};
    double Bcov[3] = {0.0, 0.0, 0.0};
    double Bsq = 0.0;
    double Bbeta = 0.0;
    double w;
    double D;
    double c = PRJ_CLIGHT;
    double c2 = PRJ_CLIGHT * PRJ_CLIGHT;
    double sqrtg;
    double alpha;
    /* Loop-invariant scale factors for the block accumulations, hoisted out of
     * the per-group/per-column loops (evaluation order preserved so the result
     * is bit-identical).  Each is a product of only function-constant terms. */
    double s_gn;
    double s_gg_f;
    double s_gu;
    double s_rgg;
    double sum_Gn = 0.0;
    double sum_Gu_xe = 0.0;
    double sum_Ggamma[3] = {0.0, 0.0, 0.0};
    double group_Gn[PRJ_RAD_GR_M1_NGROUPS];
    double group_Ggamma[PRJ_RAD_GR_M1_NGROUPS][3];
    int a;
    int b;
    int d;
    int field;
    int group;
    int block_col;
    int n;
    int v;
    int ok;

    if (resid != 0) {
        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
            resid[v] = 0.0;
        }
    }
    if (u_new_out != 0) {
        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
            u_new_out[v] = 0.0;
        }
    }
    if (W_out != 0) {
        for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
            W_out[v] = 0.0;
        }
    }
    if (jac != 0) {
        for (v = 0; v < PRJ_NVAR_CONS * PRJ_RAD_GR_M1_NP; ++v) {
            jac[v] = 0.0;
        }
    }
    prj_rad_gr_m1_jac_blocks_zero(blocks);

    if (rad == 0 || eos == 0 || geom == 0 || u_old == 0 || P == 0 ||
        blocks == 0 || !isfinite(dt) ||
        !isfinite(geom->sqrt_gamma) || geom->sqrt_gamma <= 0.0 ||
        !isfinite(geom->alpha) || geom->alpha <= 0.0 ||
        !isfinite(P[0]) || P[0] <= 0.0 ||
        !isfinite(P[4]) || P[4] <= 0.0 || !isfinite(P[5])) {
        return 0;
    }
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            if (!isfinite(geom->gamma[a][b]) ||
                !isfinite(geom->gamma_inv[a][b])) {
                return 0;
            }
            eos_geom.gamma[a][b] = geom->gamma[a][b];
        }
    }
    sqrtg = geom->sqrt_gamma;
    alpha = geom->alpha;
    /* Preserve the left-to-right evaluation order of the original inline
     * products so these hoisted factors are bit-identical to recomputing. */
    s_gn = sqrtg * dt * RAD_SCALE * c;
    s_gg_f = alpha * sqrtg * dt * RAD_SCALE;
    s_gu = alpha * sqrtg * dt * c;
    s_rgg = alpha * sqrtg * dt * c * c;
    prj_rad_gr_m1_metric4_from_geom(geom, g_cov, g_con);
    prj_rad_grm1_R_jac_geom_init(alpha, g_con, geom->gamma_inv,
        &R_jac_geom);
    n_cov[0] = -geom->alpha;
    n_cov[1] = 0.0;
    n_cov[2] = 0.0;
    n_cov[3] = 0.0;
    /* TEMP TIMER: remove after rad-matter coupling profiling. */
    PRJ_TIMER_CURRENT_START("rad_matter_temp_resjac_eos_derivs");
    ok = prj_rad_gr_m1_fluid_four_velocity(geom, g_cov, P, ucon, ucov) &&
        prj_eos_rty_derivs(eos, P[0], P[4], P[5], &eint, &pressure,
            &deint_drho, &deint_dT, &deint_dYe, &dpressure_drho,
            &dpressure_dT, &dpressure_dYe, PRJ_EOS_CTX_MAIN);
    PRJ_TIMER_CURRENT_STOP("rad_matter_temp_resjac_eos_derivs");
    if (!ok) {
        return 0;
    }
    if (!isfinite(eint) || eint < 0.0 || !isfinite(pressure)) {
        return 0;
    }
    /* TEMP TIMER: remove after rad-matter coupling profiling. */
    PRJ_TIMER_CURRENT_START("rad_matter_temp_resjac_opac_derivs");
    prj_rad3_opac_lookup_derivs(rad, P[0], P[4], P[5], kappa, sigma, delta,
        eta, dkappa_drho, dkappa_dT, dkappa_dYe, dsigma_drho, dsigma_dT,
        dsigma_dYe, ddelta_drho, ddelta_dT, ddelta_dYe, deta_drho, deta_dT,
        deta_dYe);
    PRJ_TIMER_CURRENT_STOP("rad_matter_temp_resjac_opac_derivs");

    for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
        W[v] = 0.0;
    }
    W[PRJ_PRIM_RHO] = P[0];
    W[PRJ_PRIM_V1] = P[1];
    W[PRJ_PRIM_V2] = P[2];
    W[PRJ_PRIM_V3] = P[3];
    W[PRJ_PRIM_EINT] = eint;
    W[PRJ_PRIM_YE] = P[5];
#if PRJ_MHD
    W[PRJ_PRIM_B1] = u_old[PRJ_CONS_B1] / sqrtg;
    W[PRJ_PRIM_B2] = u_old[PRJ_CONS_B2] / sqrtg;
    W[PRJ_PRIM_B3] = u_old[PRJ_CONS_B3] / sqrtg;
#endif

    memset(ducon, 0, sizeof(ducon));
    memset(ducov, 0, sizeof(ducov));
    for (d = 0; d < 3; ++d) {
        beta_con[d] = P[1 + d] / c;
    }
    for (d = 0; d < 3; ++d) {
        for (b = 0; b < 3; ++b) {
            beta_cov[d] += geom->gamma[d][b] * beta_con[b];
        }
        beta2 += beta_cov[d] * beta_con[d];
    }
    if (!isfinite(beta2) || beta2 < 0.0 || beta2 >= 1.0) {
        return 0;
    }
    wlor = 1.0 / sqrt(1.0 - beta2);
    wlor2 = wlor * wlor;
    wlor_m1 = beta2 / (sqrt(1.0 - beta2) * (1.0 + sqrt(1.0 - beta2)));
#if PRJ_MHD
    Bcon[0] = u_old[PRJ_CONS_B1] / sqrtg;
    Bcon[1] = u_old[PRJ_CONS_B2] / sqrtg;
    Bcon[2] = u_old[PRJ_CONS_B3] / sqrtg;
#endif
    for (d = 0; d < 3; ++d) {
        for (b = 0; b < 3; ++b) {
            Bcov[d] += geom->gamma[d][b] * Bcon[b];
        }
        Bsq += Bcov[d] * Bcon[d];
        Bbeta += Bcov[d] * beta_con[d];
    }
    w = P[0] * c2 + P[0] * eint + pressure;
    D = P[0] * wlor;

    for (n = 0; n < PRJ_RAD_GR_M1_NFLUID; ++n) {
        double drho = n == 0 ? 1.0 : 0.0;
        double dYe = n == 5 ? 1.0 : 0.0;
        double dbeta_con[3] = {0.0, 0.0, 0.0};
        double dbeta_cov[3] = {0.0, 0.0, 0.0};
        double dbeta2 = 0.0;
        double dwlor;
        double dwlor2;
        double dBbeta = 0.0;
        double deint = 0.0;
        double dpressure = 0.0;
        double dw;
        double dD;
        double dA;
        double dUtmp;

        if (n >= 1 && n <= 3) {
            dbeta_con[n - 1] = 1.0 / c;
        }
        for (d = 0; d < 3; ++d) {
            int m;

            for (m = 0; m < 3; ++m) {
                dbeta_cov[d] += geom->gamma[d][m] * dbeta_con[m];
            }
            dbeta2 += dbeta_cov[d] * beta_con[d] +
                beta_cov[d] * dbeta_con[d];
            dBbeta += Bcov[d] * dbeta_con[d];
        }
        dwlor = 0.5 * wlor * wlor * wlor * dbeta2;
        dwlor2 = 2.0 * wlor * dwlor;
        if (n == 0) {
            deint = deint_drho;
            dpressure = dpressure_drho;
        } else if (n == 4) {
            deint = deint_dT;
            dpressure = dpressure_dT;
        } else if (n == 5) {
            deint = deint_dYe;
            dpressure = dpressure_dYe;
        }
        dw = drho * (c2 + eint) + P[0] * deint + dpressure;

        ducon[0][n] = dwlor / geom->alpha;
        for (d = 0; d < 3; ++d) {
            ducon[d + 1][n] = dwlor *
                (beta_con[d] - geom->beta[d] / geom->alpha) +
                wlor * dbeta_con[d];
        }
        for (a = 0; a < 4; ++a) {
            for (b = 0; b < 4; ++b) {
                ducov[a][n] += g_cov[a][b] * ducon[b][n];
            }
        }

        dD = drho * wlor + P[0] * dwlor;
        blocks->fluid[0][n] += sqrtg * dD;
        for (d = 0; d < 3; ++d) {
            dUtmp = ((dw * wlor2 + w * dwlor2) * beta_cov[d] +
                (w * wlor2 + Bsq) * dbeta_cov[d] -
                dBbeta * Bcov[d]) / c;
            blocks->fluid[1 + d][n] += sqrtg * dUtmp;
        }
        dA = drho * eint + P[0] * deint + dpressure;
        dUtmp = dA * wlor2 + (P[0] * eint + pressure) * dwlor2 +
            c2 * (drho * wlor * wlor_m1 +
                P[0] * (dwlor * wlor_m1 + wlor * dwlor)) -
            dpressure - Bbeta * dBbeta;
        if (wlor2 > 0.0) {
            dUtmp += 0.5 * Bsq * dwlor2 / (wlor2 * wlor2);
        }
        blocks->fluid[4][n] += sqrtg * dUtmp;
        blocks->fluid[5][n] += sqrtg * (dD * P[5] + D * dYe);
    }

    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            int idx = field * PRJ_NEGROUP + group;
            int pidx = 6 + 4 * idx;
            double E = P[pidx];
            double xe = rad->x_e[field][group];
            double Fcov[3];
            prj_rad_grm1_R_projector Rproj;
            double R_u[4];
            double Ruu = 0.0;
            double dR_u_rad[PRJ_RAD_GR_M1_NRAD_BLOCK][4];
            double dRuu_rad[PRJ_RAD_GR_M1_NRAD_BLOCK];
            double Gcon[4];
            double kappa_eff = kappa[idx] + rad->gr_m1_eleinel_dkappa[idx];
            double eta_eff = eta[idx] + rad->gr_m1_eleinel_deta[idx];
            double sigma_eff = sigma[idx] * (1.0 - delta[idx] / 3.0);
            double kt = kappa_eff + sigma_eff;
            double scalar;
            double Gdotn = 0.0;
            double Gdotu = 0.0;

            /* Newton variables are the lab-frame moments p = (E, F_i).  The
             * conserved rad rows are sqrt(gamma) (E, F_i), so the transport
             * part of the block is the constant, invertible sqrt(gamma) I
             * (written below), the rad_fluid transport is exactly zero, and the
             * closure derivative dR^{ab}/dp enters only the opacity source. */
            if (!isfinite(E) || E < 0.0) {
                return 0;
            }
            Fcov[0] = P[pidx + 1];
            Fcov[1] = P[pidx + 2];
            Fcov[2] = P[pidx + 3];
            for (d = 0; d < 3; ++d) {
                if (!isfinite(Fcov[d])) {
                    return 0;
                }
            }
            if (!prj_rad_grm1_R_projector_init(g_cov, g_con, &R_jac_geom,
                    E, Fcov, &Rproj)) {
                return 0;
            }
            prj_rad_grm1_R_projector_apply(&Rproj, ucov, R_u);
            for (a = 0; a < 4; ++a) {
                Ruu += R_u[a] * ucov[a];
            }
            prj_rad_grm1_R_projector_dR_u(&Rproj, g_con, &R_jac_geom,
                ucov, dR_u_rad, dRuu_rad);
            for (a = 0; a < 4; ++a) {
                if (!isfinite(R_u[a]) || !isfinite(dRuu_rad[a])) {
                    return 0;
                }
                for (b = 0; b < 4; ++b) {
                    if (!isfinite(dR_u_rad[a][b])) {
                        return 0;
                    }
                }
            }
            scalar = sigma_eff * Ruu + eta_eff / c;
            if (!isfinite(kappa_eff) || !isfinite(sigma_eff) ||
                !isfinite(eta_eff) || !isfinite(scalar)) {
                return 0;
            }
            for (a = 0; a < 4; ++a) {
                Gcon[a] = -kt * R_u[a] - scalar * ucon[a];
                if (!isfinite(Gcon[a])) {
                    return 0;
                }
                Gdotn += Gcon[a] * n_cov[a];
                Gdotu += Gcon[a] * ucov[a];
            }

            W[PRJ_PRIM_RAD_E(field, group)] = E;
            W[PRJ_PRIM_RAD_F1(field, group)] = Fcov[0];
            W[PRJ_PRIM_RAD_F2(field, group)] = Fcov[1];
            W[PRJ_PRIM_RAD_F3(field, group)] = Fcov[2];

            /* Transport block: d(sqrt(gamma) E)/dE = sqrt(gamma) and
             * d(sqrt(gamma) F_i)/dF_j = sqrt(gamma) delta_ij; no dependence on
             * the fluid columns or the other radiation moments. */
            blocks->rad_rad[idx][0][0] += sqrtg;
            blocks->rad_rad[idx][1][1] += sqrtg;
            blocks->rad_rad[idx][2][2] += sqrtg;
            blocks->rad_rad[idx][3][3] += sqrtg;

            group_Gn[idx] = Gdotn;
            sum_Gn += Gdotn;
            sum_Gu_xe += Gdotu * rad->x_e[field][group];
            for (d = 0; d < 3; ++d) {
                double Ggamma = 0.0;

                /* n_cov[d+1] == 0 (spatial), so the n_cov n_cov term drops. */
                for (a = 0; a < 4; ++a) {
                    Ggamma += Gcon[a] * g_cov[d + 1][a];
                }
                group_Ggamma[idx][d] = Ggamma;
                sum_Ggamma[d] += Ggamma;
            }

            {
                double delta_fac = 1.0 - delta[idx] / 3.0;
                double dsigma_eff_col[PRJ_RAD_GR_M1_NFLUID] =
                    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
                double dkt_col[PRJ_RAD_GR_M1_NFLUID] =
                    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
                double deta_col[PRJ_RAD_GR_M1_NFLUID] =
                    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
                dsigma_eff_col[0] = dsigma_drho[idx] * delta_fac -
                    sigma[idx] * ddelta_drho[idx] / 3.0;
                dkt_col[0] = dkappa_drho[idx] + dsigma_eff_col[0];
                deta_col[0] = deta_drho[idx];
                dsigma_eff_col[4] = dsigma_dT[idx] * delta_fac -
                    sigma[idx] * ddelta_dT[idx] / 3.0;
                dkt_col[4] = dkappa_dT[idx] + dsigma_eff_col[4];
                deta_col[4] = deta_dT[idx];
                dsigma_eff_col[5] = dsigma_dYe[idx] * delta_fac -
                    sigma[idx] * ddelta_dYe[idx] / 3.0;
                dkt_col[5] = dkappa_dYe[idx] + dsigma_eff_col[5];
                deta_col[5] = deta_dYe[idx];

                /* Fluid columns: R is fixed; only u^a/u_a and opacities move. */
                for (block_col = 0; block_col < PRJ_RAD_GR_M1_NFLUID;
                     ++block_col) {
                    double dR_u[4] = {0.0, 0.0, 0.0, 0.0};
                    double dRuu = 0.0;
                    double ducov_col[4];
                    double dscalar;
                    double dGcon[4];
                    double dGn = 0.0;
                    double dGu = 0.0;
                    double dGgamma[3] = {0.0, 0.0, 0.0};

                    for (a = 0; a < 4; ++a) {
                        ducov_col[a] = ducov[a][block_col];
                    }
                    prj_rad_grm1_R_projector_apply(&Rproj, ducov_col,
                        dR_u);
                    for (a = 0; a < 4; ++a) {
                        double ducov_an = ducov[a][block_col];

                        dRuu += dR_u[a] * ucov[a] + R_u[a] * ducov_an;
                    }
                    dscalar = dsigma_eff_col[block_col] * Ruu +
                        sigma_eff * dRuu + deta_col[block_col] / c;
                    for (a = 0; a < 4; ++a) {
                        double ducon_an = ducon[a][block_col];
                        double ducov_an = ducov[a][block_col];

                        dGcon[a] = -dkt_col[block_col] * R_u[a] -
                            kt * dR_u[a] - dscalar * ucon[a] -
                            scalar * ducon_an;
                        dGn += dGcon[a] * n_cov[a];
                        dGu += dGcon[a] * ucov[a] + Gcon[a] * ducov_an;
                        for (d = 0; d < 3; ++d) {
                            dGgamma[d] += dGcon[a] * g_cov[d + 1][a];
                        }
                    }

                    blocks->fluid[4][block_col] += s_gn * dGn;
                    for (d = 0; d < 3; ++d) {
                        blocks->fluid[1 + d][block_col] -=
                            s_gg_f * dGgamma[d];
                    }
                    blocks->fluid[5][block_col] += s_gu * dGu * xe;
                    blocks->rad_fluid[idx][0][block_col] -= s_gu * dGn;
                    for (d = 0; d < 3; ++d) {
                        blocks->rad_fluid[idx][1 + d][block_col] +=
                            s_rgg * dGgamma[d];
                    }
                }

                /* Radiation columns: u^a/u_a and opacities are fixed; only the
                 * exact algebraic M1 closure derivative contributes. */
                for (block_col = 0; block_col < PRJ_RAD_GR_M1_NRAD_BLOCK;
                     ++block_col) {
                    const double *dR_u = dR_u_rad[block_col];
                    double dRuu = dRuu_rad[block_col];
                    double dscalar;
                    double dGcon[4];
                    double dGn = 0.0;
                    double dGu = 0.0;
                    double dGgamma[3] = {0.0, 0.0, 0.0};

                    dscalar = sigma_eff * dRuu;
                    for (a = 0; a < 4; ++a) {
                        dGcon[a] = -kt * dR_u[a] - dscalar * ucon[a];
                        dGn += dGcon[a] * n_cov[a];
                        dGu += dGcon[a] * ucov[a];
                        for (d = 0; d < 3; ++d) {
                            dGgamma[d] += dGcon[a] * g_cov[d + 1][a];
                        }
                    }

                    blocks->fluid_rad[idx][4][block_col] += s_gn * dGn;
                    for (d = 0; d < 3; ++d) {
                        blocks->fluid_rad[idx][1 + d][block_col] -=
                            s_gg_f * dGgamma[d];
                    }
                    blocks->fluid_rad[idx][5][block_col] +=
                        s_gu * dGu * xe;
                    blocks->rad_rad[idx][0][block_col] -= s_gu * dGn;
                    for (d = 0; d < 3; ++d) {
                        blocks->rad_rad[idx][1 + d][block_col] +=
                            s_rgg * dGgamma[d];
                    }
                }
            }
        }
    }

    /* TEMP TIMER: remove after rad-matter coupling profiling. */
    PRJ_TIMER_CURRENT_START("rad_matter_temp_resjac_prim2cons");
    ok = prj_eos_gr_prim2cons(eos, &eos_geom, W, u_new,
        PRJ_EOS_CTX_MAIN) == PRJ_EOS_GR_OK;
    PRJ_TIMER_CURRENT_STOP("rad_matter_temp_resjac_prim2cons");
    if (!ok) {
        return 0;
    }

    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        resid_tmp[v] = u_new[v] - u_old[v];
    }
    resid_tmp[PRJ_CONS_ETOT] += sqrtg * dt * RAD_SCALE * c * sum_Gn;
    resid_tmp[PRJ_CONS_MOM1] -= geom->alpha * sqrtg * dt *
        RAD_SCALE * sum_Ggamma[0];
    resid_tmp[PRJ_CONS_MOM2] -= geom->alpha * sqrtg * dt *
        RAD_SCALE * sum_Ggamma[1];
    resid_tmp[PRJ_CONS_MOM3] -= geom->alpha * sqrtg * dt *
        RAD_SCALE * sum_Ggamma[2];
    resid_tmp[PRJ_CONS_YE] += geom->alpha * sqrtg * dt * c * sum_Gu_xe;
    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            int idx = field * PRJ_NEGROUP + group;

            resid_tmp[PRJ_CONS_RAD_E(field, group)] -=
                geom->alpha * sqrtg * dt * c * group_Gn[idx];
            resid_tmp[PRJ_CONS_RAD_F1(field, group)] +=
                geom->alpha * sqrtg * dt * c * c * group_Ggamma[idx][0];
            resid_tmp[PRJ_CONS_RAD_F2(field, group)] +=
                geom->alpha * sqrtg * dt * c * c * group_Ggamma[idx][1];
            resid_tmp[PRJ_CONS_RAD_F3(field, group)] +=
                geom->alpha * sqrtg * dt * c * c * group_Ggamma[idx][2];
        }
    }

    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        resid_use[v] = resid_tmp[v];
        if (u_new_out != 0) {
            u_new_out[v] = u_new[v];
        }
    }
    if (W_out != 0) {
        memcpy(W_out, W, sizeof(W));
    }
    if (jac != 0) {
        /* TEMP TIMER: remove after rad-matter coupling profiling. */
        PRJ_TIMER_CURRENT_START("rad_matter_temp_resjac_dense");
        prj_rad_gr_m1_jac_blocks_to_dense(blocks, jac);
        PRJ_TIMER_CURRENT_STOP("rad_matter_temp_resjac_dense");
    }
    return 1;
}

static int prj_rad_gr_m1_jacobian(const prj_rad *rad, prj_eos *eos,
    const prj_z4c_hydro_geom *geom, const double *u_old, const double *P,
    double dt, double *resid, double *jac, double *u_new_out)
{
    prj_rad_gr_m1_jac_blocks blocks;

    if (jac == 0) {
        return 0;
    }
    return prj_rad_gr_m1_residual_jacobian(rad, eos, geom, u_old, P, dt,
        resid, &blocks, jac, u_new_out, 0);
}

static int prj_rad_gr_m1_active_row(int eq)
{
    if (eq == 0) {
        return PRJ_CONS_RHO;
    }
    if (eq >= 1 && eq <= 3) {
        return PRJ_CONS_MOM1 + (eq - 1);
    }
    if (eq == 4) {
        return PRJ_CONS_ETOT;
    }
    if (eq == 5) {
        return PRJ_CONS_YE;
    }
    {
        int idx = (eq - 6) / 4;
        int comp = (eq - 6) % 4;
        int field = idx / PRJ_NEGROUP;
        int group = idx % PRJ_NEGROUP;

        return PRJ_CONS_RAD_E(field, group) + comp;
    }
}

static double prj_rad_gr_m1_residual_norm(const double *u_old,
    const double *resid, double threshold)
{
    double max_norm = 0.0;
    double rho_scale;
    double etot_scale;
    double mom2;
    double mom_norm;
    int field;
    int group;
    int d;

    if (u_old == 0 || resid == 0 || !isfinite(threshold) ||
        threshold <= 0.0 || !isfinite(u_old[PRJ_CONS_YE]) ||
        u_old[PRJ_CONS_YE] <= 0.0) {
        return HUGE_VAL;
    }

    rho_scale = fabs(u_old[PRJ_CONS_RHO]);
    etot_scale = fabs(u_old[PRJ_CONS_ETOT]);
    if (!isfinite(rho_scale) || rho_scale <= 0.0 ||
        !isfinite(etot_scale) || etot_scale <= 0.0) {
        return HUGE_VAL;
    }
    if (!isfinite(resid[PRJ_CONS_RHO]) ||
        !isfinite(resid[PRJ_CONS_ETOT]) ||
        !isfinite(resid[PRJ_CONS_YE])) {
        return HUGE_VAL;
    }

    max_norm = fmax(max_norm, fabs(resid[PRJ_CONS_RHO]) / rho_scale);
    max_norm = fmax(max_norm, fabs(resid[PRJ_CONS_ETOT]) / etot_scale);
    max_norm = fmax(max_norm,
        fabs(resid[PRJ_CONS_YE]) / u_old[PRJ_CONS_YE]);

    mom2 = 0.0;
    for (d = 0; d < 3; ++d) {
        double mom = u_old[PRJ_CONS_MOM1 + d];

        if (!isfinite(mom) || !isfinite(resid[PRJ_CONS_MOM1 + d])) {
            return HUGE_VAL;
        }
        mom2 += mom * mom;
    }
    if (!isfinite(mom2)) {
        return HUGE_VAL;
    }
    if (mom2 == 0.0) {
        mom2 = rho_scale * rho_scale * PRJ_RAD_GR_M1_SOLVE_VMIN *
            PRJ_RAD_GR_M1_SOLVE_VMIN;
    }
    mom_norm = sqrt(mom2);
    if (!isfinite(mom_norm) || mom_norm <= 0.0) {
        return HUGE_VAL;
    }
    {
        double rmom2 = 0.0;

        for (d = 0; d < 3; ++d) {
            double r = resid[PRJ_CONS_MOM1 + d];

            rmom2 += r * r;
        }
        max_norm = fmax(max_norm, sqrt(rmom2) / mom_norm);
    }

    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            int eidx = PRJ_CONS_RAD_E(field, group);
            int fidx = PRJ_CONS_RAD_F1(field, group);
            /* Scale the radiation residual by a per-group share of the hydro
             * energy/momentum, hydro/(NEGROUP*NRAD), rather than by
             * threshold*hydro.  The old floor made the effective raw tolerance
             * threshold^2*hydro (comparing to threshold after the divide),
             * which was far too strict for the radiation sector. */
            double scale_e = fmax(etot_scale /
                    (RAD_SCALE * PRJ_NEGROUP * PRJ_NRAD),
                fabs(u_old[eidx]));
            double flux2 = 0.0;
            double rflux2 = 0.0;
            double scale_f;

            if (!isfinite(u_old[eidx]) || !isfinite(resid[eidx]) ||
                !isfinite(scale_e) || scale_e <= 0.0) {
                return HUGE_VAL;
            }
            max_norm = fmax(max_norm, fabs(resid[eidx]) / scale_e);

            for (d = 0; d < 3; ++d) {
                double f = u_old[fidx + d];
                double r = resid[fidx + d];

                if (!isfinite(f) || !isfinite(r)) {
                    return HUGE_VAL;
                }
                flux2 += f * f;
                rflux2 += r * r;
            }
            if (!isfinite(flux2) || !isfinite(rflux2)) {
                return HUGE_VAL;
            }
            scale_f = fmax(mom_norm * PRJ_CLIGHT * PRJ_CLIGHT /
                    (RAD_SCALE * PRJ_NEGROUP * PRJ_NRAD), sqrt(flux2));
            if (!isfinite(scale_f) || scale_f <= 0.0) {
                return HUGE_VAL;
            }
            max_norm = fmax(max_norm, sqrt(rflux2) / scale_f);
        }
    }

    return max_norm;
}

typedef struct prj_rad_gr_m1_residual_norm_state {
    const double *u_old;
    double max_norm;
    double rho_scale;
    double etot_scale;
    double mom_norm;
} prj_rad_gr_m1_residual_norm_state;

static int prj_rad_gr_m1_residual_norm_state_init(const double *u_old,
    double threshold, prj_rad_gr_m1_residual_norm_state *state)
{
    double rho_scale;
    double etot_scale;
    double mom2 = 0.0;
    double mom_norm;
    int d;

    if (u_old == 0 || state == 0 || !isfinite(threshold) ||
        threshold <= 0.0 || !isfinite(u_old[PRJ_CONS_YE]) ||
        u_old[PRJ_CONS_YE] <= 0.0) {
        return 0;
    }
    rho_scale = fabs(u_old[PRJ_CONS_RHO]);
    etot_scale = fabs(u_old[PRJ_CONS_ETOT]);
    if (!isfinite(rho_scale) || rho_scale <= 0.0 ||
        !isfinite(etot_scale) || etot_scale <= 0.0) {
        return 0;
    }
    for (d = 0; d < 3; ++d) {
        double mom = u_old[PRJ_CONS_MOM1 + d];

        if (!isfinite(mom)) {
            return 0;
        }
        mom2 += mom * mom;
    }
    if (!isfinite(mom2)) {
        return 0;
    }
    if (mom2 == 0.0) {
        mom2 = rho_scale * rho_scale * PRJ_RAD_GR_M1_SOLVE_VMIN *
            PRJ_RAD_GR_M1_SOLVE_VMIN;
    }
    mom_norm = sqrt(mom2);
    if (!isfinite(mom_norm) || mom_norm <= 0.0) {
        return 0;
    }
    state->u_old = u_old;
    state->max_norm = 0.0;
    state->rho_scale = rho_scale;
    state->etot_scale = etot_scale;
    state->mom_norm = mom_norm;
    return 1;
}

static int prj_rad_gr_m1_residual_norm_add_hydro(
    prj_rad_gr_m1_residual_norm_state *state, double r_rho,
    const double r_mom[3], double r_etot, double r_ye)
{
    double rmom2 = 0.0;
    int d;

    if (state == 0 || r_mom == 0 || !isfinite(r_rho) ||
        !isfinite(r_etot) || !isfinite(r_ye)) {
        return 0;
    }
    state->max_norm = fmax(state->max_norm, fabs(r_rho) / state->rho_scale);
    state->max_norm = fmax(state->max_norm,
        fabs(r_etot) / state->etot_scale);
    state->max_norm = fmax(state->max_norm,
        fabs(r_ye) / state->u_old[PRJ_CONS_YE]);
    for (d = 0; d < 3; ++d) {
        if (!isfinite(r_mom[d])) {
            return 0;
        }
        rmom2 += r_mom[d] * r_mom[d];
    }
    if (!isfinite(rmom2)) {
        return 0;
    }
    state->max_norm = fmax(state->max_norm, sqrt(rmom2) / state->mom_norm);
    return isfinite(state->max_norm);
}

static int prj_rad_gr_m1_residual_norm_add_rad_group(
    prj_rad_gr_m1_residual_norm_state *state, int field, int group,
    double r_e, const double r_f[3])
{
    int eidx;
    int fidx;
    double scale_e;
    double flux2 = 0.0;
    double rflux2 = 0.0;
    double scale_f;
    int d;

    if (state == 0 || r_f == 0 || field < 0 || field >= PRJ_NRAD ||
        group < 0 || group >= PRJ_NEGROUP) {
        return 0;
    }
    eidx = PRJ_CONS_RAD_E(field, group);
    fidx = PRJ_CONS_RAD_F1(field, group);
    scale_e = fmax(state->etot_scale / (RAD_SCALE * PRJ_NEGROUP * PRJ_NRAD),
        fabs(state->u_old[eidx]));
    if (!isfinite(state->u_old[eidx]) || !isfinite(r_e) ||
        !isfinite(scale_e) || scale_e <= 0.0) {
        return 0;
    }
    state->max_norm = fmax(state->max_norm, fabs(r_e) / scale_e);
    for (d = 0; d < 3; ++d) {
        double f = state->u_old[fidx + d];
        double r = r_f[d];

        if (!isfinite(f) || !isfinite(r)) {
            return 0;
        }
        flux2 += f * f;
        rflux2 += r * r;
    }
    if (!isfinite(flux2) || !isfinite(rflux2)) {
        return 0;
    }
    scale_f = fmax(state->mom_norm * PRJ_CLIGHT * PRJ_CLIGHT /
            (RAD_SCALE * PRJ_NEGROUP * PRJ_NRAD), sqrt(flux2));
    if (!isfinite(scale_f) || scale_f <= 0.0) {
        return 0;
    }
    state->max_norm = fmax(state->max_norm, sqrt(rflux2) / scale_f);
    return isfinite(state->max_norm);
}

static int prj_rad_gr_m1_exact_residual_check(const prj_rad *rad,
    prj_eos *eos, const prj_z4c_hydro_geom *geom, const double *u_old,
    const double *P, double dt, const prj_eos_rty_interp_result *eos_in,
    const prj_rad_gr_m1_opac_cache *opac_in,
    const prj_rad_gr_m1_cell_context *ctx_in, double threshold,
    double *norm_out, double *resid_out, double *u_new_out, double *W_out)
{
    double W[PRJ_NVAR_PRIM];
    double u_new[PRJ_NVAR_CONS];
    double resid_tmp[PRJ_NVAR_CONS];
    double eos_q[PRJ_EOS_NQUANT];
    prj_rad_gr_m1_cell_context ctx_local;
    const prj_rad_gr_m1_cell_context *ctx;
    prj_rad_gr_m1_residual_norm_state norm_state;
    const double (*g_cov)[4];
    const double (*g_con)[4];
    const double *n_cov;
    const double *ucon;
    const double *ucov;
    prj_rad_grm1_R_jac_geom R_jac_geom;
    double kappa_local[PRJ_RAD_GR_M1_NGROUPS];
    double sigma_local[PRJ_RAD_GR_M1_NGROUPS];
    double delta_local[PRJ_RAD_GR_M1_NGROUPS];
    double eta_local[PRJ_RAD_GR_M1_NGROUPS];
    const double *kappa;
    const double *sigma;
    const double *delta;
    const double *eta;
    double pressure;
    double sum_Gn = 0.0;
    double sum_Gu_xe = 0.0;
    double sum_Ggamma[3] = {0.0, 0.0, 0.0};
    double group_Gn[PRJ_RAD_GR_M1_NGROUPS];
    double group_Ggamma[PRJ_RAD_GR_M1_NGROUPS][3];
    int field;
    int group;
    int a;
    int d;
    int v;

    if (resid_out != 0) {
        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
            resid_out[v] = 0.0;
        }
    }
    if (u_new_out != 0) {
        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
            u_new_out[v] = 0.0;
        }
    }
    if (W_out != 0) {
        for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
            W_out[v] = 0.0;
        }
    }
    if (norm_out != 0) {
        *norm_out = HUGE_VAL;
    }
    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        u_new[v] = 0.0;
    }
    if (rad == 0 || eos == 0 || geom == 0 || u_old == 0 || P == 0 ||
        norm_out == 0 || !isfinite(dt) ||
        !isfinite(geom->sqrt_gamma) || geom->sqrt_gamma <= 0.0 ||
        !isfinite(geom->alpha) || geom->alpha <= 0.0 ||
        !isfinite(P[0]) || P[0] <= 0.0 ||
        !isfinite(P[4]) || P[4] <= 0.0 || !isfinite(P[5]) ||
        !prj_rad_gr_m1_residual_norm_state_init(u_old, threshold,
            &norm_state)) {
        return 0;
    }
    if (ctx_in != 0) {
        ctx = ctx_in;
    } else {
        if (!prj_rad_gr_m1_cell_context_init(geom, P, &ctx_local)) {
            return 0;
        }
        ctx = &ctx_local;
    }
    g_cov = ctx->g_cov;
    g_con = ctx->g_con;
    n_cov = ctx->n_cov;
    ucon = ctx->ucon;
    ucov = ctx->ucov;

    for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
        W[v] = 0.0;
    }
    W[PRJ_PRIM_RHO] = P[0];
    W[PRJ_PRIM_V1] = P[1];
    W[PRJ_PRIM_V2] = P[2];
    W[PRJ_PRIM_V3] = P[3];
    W[PRJ_PRIM_YE] = P[5];
    if (eos_in != 0) {
        W[PRJ_PRIM_EINT] = eos_in->eint;
        pressure = eos_in->pressure;
    } else {
        prj_eos_rty(eos, P[0], P[4], P[5], eos_q, PRJ_EOS_CTX_MAIN);
        W[PRJ_PRIM_EINT] = eos_q[PRJ_EOS_EINT];
        pressure = eos_q[PRJ_EOS_PRESSURE];
    }
    if (!isfinite(W[PRJ_PRIM_EINT]) || W[PRJ_PRIM_EINT] < 0.0 ||
        !isfinite(pressure)) {
        return 0;
    }
#if PRJ_MHD
    W[PRJ_PRIM_B1] = u_old[PRJ_CONS_B1] / geom->sqrt_gamma;
    W[PRJ_PRIM_B2] = u_old[PRJ_CONS_B2] / geom->sqrt_gamma;
    W[PRJ_PRIM_B3] = u_old[PRJ_CONS_B3] / geom->sqrt_gamma;
#endif
    if (!prj_rad_gr_m1_hydro_cons_cached(geom, u_old, P,
            W[PRJ_PRIM_EINT], pressure, u_new)) {
        return 0;
    }
    prj_rad_grm1_R_jac_geom_init(geom->alpha, g_con, geom->gamma_inv,
        &R_jac_geom);

    if (opac_in != 0) {
        kappa = opac_in->kappa;
        sigma = opac_in->sigma;
        delta = opac_in->delta;
        eta = opac_in->eta;
    } else {
        PRJ_TIMER_CURRENT_START("rad_matter_temp_exact_check_opac");
        prj_rad3_opac_lookup(rad, P[0], P[4], P[5], kappa_local,
            sigma_local, delta_local, eta_local);
        PRJ_TIMER_CURRENT_STOP("rad_matter_temp_exact_check_opac");
        kappa = kappa_local;
        sigma = sigma_local;
        delta = delta_local;
        eta = eta_local;
    }

    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            int idx = field * PRJ_NEGROUP + group;
            int pidx = PRJ_RAD_GR_M1_NFLUID +
                PRJ_RAD_GR_M1_NRAD_BLOCK * idx;
            double E = P[pidx];
            double Fcov[3];
            double R_u[4];
            double Ruu = 0.0;
            double Gcon[4];
            prj_rad_grm1_R_projector Rproj;
            double kappa_eff = kappa[idx] + rad->gr_m1_eleinel_dkappa[idx];
            double eta_eff = eta[idx] + rad->gr_m1_eleinel_deta[idx];
            double sigma_eff = sigma[idx] * (1.0 - delta[idx] / 3.0);
            double scalar;
            double Gdotn = 0.0;
            double Gdotu = 0.0;

            if (!isfinite(E) || E < 0.0 || !isfinite(kappa_eff) ||
                !isfinite(sigma_eff) || !isfinite(eta_eff)) {
                return 0;
            }
            Fcov[0] = P[pidx + 1];
            Fcov[1] = P[pidx + 2];
            Fcov[2] = P[pidx + 3];
            for (a = 0; a < 3; ++a) {
                if (!isfinite(Fcov[a])) {
                    return 0;
                }
            }
            if (!prj_rad_grm1_R_projector_init(g_cov, g_con, &R_jac_geom,
                    E, Fcov, &Rproj)) {
                return 0;
            }
            prj_rad_grm1_R_projector_apply(&Rproj, ucov, R_u);
            W[PRJ_PRIM_RAD_E(field, group)] = E;
            W[PRJ_PRIM_RAD_F1(field, group)] = Fcov[0];
            W[PRJ_PRIM_RAD_F2(field, group)] = Fcov[1];
            W[PRJ_PRIM_RAD_F3(field, group)] = Fcov[2];
            u_new[PRJ_CONS_RAD_E(field, group)] = geom->sqrt_gamma * E;
            u_new[PRJ_CONS_RAD_F1(field, group)] =
                geom->sqrt_gamma * Fcov[0];
            u_new[PRJ_CONS_RAD_F2(field, group)] =
                geom->sqrt_gamma * Fcov[1];
            u_new[PRJ_CONS_RAD_F3(field, group)] =
                geom->sqrt_gamma * Fcov[2];

            for (a = 0; a < 4; ++a) {
                Ruu += R_u[a] * ucov[a];
            }
            scalar = sigma_eff * Ruu + eta_eff / PRJ_CLIGHT;
            if (!isfinite(scalar)) {
                return 0;
            }
            for (a = 0; a < 4; ++a) {
                Gcon[a] = -(kappa_eff + sigma_eff) * R_u[a] -
                    scalar * ucon[a];
                if (!isfinite(Gcon[a])) {
                    return 0;
                }
                Gdotn += Gcon[a] * n_cov[a];
                Gdotu += Gcon[a] * ucov[a];
            }
            group_Gn[idx] = Gdotn;
            sum_Gn += Gdotn;
            sum_Gu_xe += Gdotu * rad->x_e[field][group];
            for (d = 0; d < 3; ++d) {
                double Ggamma = 0.0;

                for (a = 0; a < 4; ++a) {
                    Ggamma += Gcon[a] * g_cov[d + 1][a];
                }
                group_Ggamma[idx][d] = Ggamma;
                sum_Ggamma[d] += Ggamma;
            }
        }
    }

    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        resid_tmp[v] = u_new[v] - u_old[v];
        if (!isfinite(resid_tmp[v])) {
            return 0;
        }
    }
    resid_tmp[PRJ_CONS_ETOT] += geom->sqrt_gamma * dt *
        RAD_SCALE * PRJ_CLIGHT * sum_Gn;
    resid_tmp[PRJ_CONS_MOM1] -= geom->alpha * geom->sqrt_gamma * dt *
        RAD_SCALE * sum_Ggamma[0];
    resid_tmp[PRJ_CONS_MOM2] -= geom->alpha * geom->sqrt_gamma * dt *
        RAD_SCALE * sum_Ggamma[1];
    resid_tmp[PRJ_CONS_MOM3] -= geom->alpha * geom->sqrt_gamma * dt *
        RAD_SCALE * sum_Ggamma[2];
    resid_tmp[PRJ_CONS_YE] += geom->alpha * geom->sqrt_gamma * dt *
        PRJ_CLIGHT * sum_Gu_xe;
    {
        double r_mom[3];

        for (d = 0; d < 3; ++d) {
            r_mom[d] = resid_tmp[PRJ_CONS_MOM1 + d];
        }
        if (!prj_rad_gr_m1_residual_norm_add_hydro(&norm_state,
                resid_tmp[PRJ_CONS_RHO], r_mom, resid_tmp[PRJ_CONS_ETOT],
                resid_tmp[PRJ_CONS_YE])) {
            return 0;
        }
    }
    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            int idx = field * PRJ_NEGROUP + group;
            double r_f[3];

            resid_tmp[PRJ_CONS_RAD_E(field, group)] -=
                geom->alpha * geom->sqrt_gamma * dt *
                PRJ_CLIGHT * group_Gn[idx];
            resid_tmp[PRJ_CONS_RAD_F1(field, group)] +=
                geom->alpha * geom->sqrt_gamma * dt *
                PRJ_CLIGHT * PRJ_CLIGHT * group_Ggamma[idx][0];
            resid_tmp[PRJ_CONS_RAD_F2(field, group)] +=
                geom->alpha * geom->sqrt_gamma * dt *
                PRJ_CLIGHT * PRJ_CLIGHT * group_Ggamma[idx][1];
            resid_tmp[PRJ_CONS_RAD_F3(field, group)] +=
                geom->alpha * geom->sqrt_gamma * dt *
                PRJ_CLIGHT * PRJ_CLIGHT * group_Ggamma[idx][2];
            for (d = 0; d < 3; ++d) {
                r_f[d] = resid_tmp[PRJ_CONS_RAD_F1(field, group) + d];
            }
            if (!prj_rad_gr_m1_residual_norm_add_rad_group(&norm_state,
                    field, group, resid_tmp[PRJ_CONS_RAD_E(field, group)],
                    r_f)) {
                return 0;
            }
        }
    }

    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        if (!isfinite(resid_tmp[v])) {
            return 0;
        }
        if (resid_out != 0) {
            resid_out[v] = resid_tmp[v];
        }
        if (u_new_out != 0) {
            u_new_out[v] = u_new[v];
        }
    }
    if (W_out != 0) {
        memcpy(W_out, W, sizeof(W));
    }
    *norm_out = norm_state.max_norm;
    return isfinite(*norm_out);
}

static int prj_rad_gr_m1_solve_dense_n(int n, double *A, double *b, double *x)
{
    int i;
    int j;
    int k;

    if (n <= 0 || A == 0 || b == 0 || x == 0) {
        return 0;
    }

    for (k = 0; k < n; ++k) {
        int pivot = k;
        double piv_abs = fabs(A[(size_t)k * (size_t)n + (size_t)k]);

        for (i = k + 1; i < n; ++i) {
            double cand = fabs(A[(size_t)i * (size_t)n + (size_t)k]);

            if (cand > piv_abs) {
                piv_abs = cand;
                pivot = i;
            }
        }
        if (!isfinite(piv_abs) || piv_abs <= 1.0e-300) {
            return 0;
        }
        if (pivot != k) {
            double tmp;

            for (j = k; j < n; ++j) {
                tmp = A[(size_t)k * (size_t)n + (size_t)j];
                A[(size_t)k * (size_t)n + (size_t)j] =
                    A[(size_t)pivot * (size_t)n + (size_t)j];
                A[(size_t)pivot * (size_t)n + (size_t)j] = tmp;
            }
            tmp = b[k];
            b[k] = b[pivot];
            b[pivot] = tmp;
        }
        for (i = k + 1; i < n; ++i) {
            double factor = A[(size_t)i * (size_t)n + (size_t)k] /
                A[(size_t)k * (size_t)n + (size_t)k];

            if (!isfinite(factor)) {
                return 0;
            }
            A[(size_t)i * (size_t)n + (size_t)k] = 0.0;
            for (j = k + 1; j < n; ++j) {
                A[(size_t)i * (size_t)n + (size_t)j] -=
                    factor * A[(size_t)k * (size_t)n + (size_t)j];
            }
            b[i] -= factor * b[k];
        }
    }

    for (i = n - 1; i >= 0; --i) {
        double sum = b[i];

        for (j = i + 1; j < n; ++j) {
            sum -= A[(size_t)i * (size_t)n + (size_t)j] * x[j];
        }
        if (!isfinite(sum) ||
            !isfinite(A[(size_t)i * (size_t)n + (size_t)i]) ||
            fabs(A[(size_t)i * (size_t)n + (size_t)i]) <= 1.0e-300) {
            return 0;
        }
        x[i] = sum / A[(size_t)i * (size_t)n + (size_t)i];
        if (!isfinite(x[i])) {
            return 0;
        }
    }

    return 1;
}

static int prj_rad_gr_m1_factor4(
    double A[PRJ_RAD_GR_M1_NRAD_BLOCK][PRJ_RAD_GR_M1_NRAD_BLOCK],
    int pivot[PRJ_RAD_GR_M1_NRAD_BLOCK], int gidx)
{
    int i;
    int j;
    int k;
    double block_max_abs = 0.0;

    if (A == 0 || pivot == 0) {
        return 0;
    }
    for (i = 0; i < PRJ_RAD_GR_M1_NRAD_BLOCK; ++i) {
        for (j = 0; j < PRJ_RAD_GR_M1_NRAD_BLOCK; ++j) {
            block_max_abs = fmax(block_max_abs, fabs(A[i][j]));
        }
    }
    if (block_max_abs <= 0.0 || !isfinite(block_max_abs)) {
        block_max_abs = 1.0;
    }

    for (k = 0; k < PRJ_RAD_GR_M1_NRAD_BLOCK; ++k) {
        int p = k;
        double piv_abs = fabs(A[k][k]);

        for (i = k + 1; i < PRJ_RAD_GR_M1_NRAD_BLOCK; ++i) {
            double cand = fabs(A[i][k]);

            if (cand > piv_abs) {
                piv_abs = cand;
                p = i;
            }
        }
        if (piv_abs < prj_rad_gr_m1_last_solve_diag.linear_min_pivot) {
            prj_rad_gr_m1_last_solve_diag.linear_min_pivot = piv_abs;
            prj_rad_gr_m1_last_solve_diag.linear_min_pivot_rel =
                piv_abs / block_max_abs;
            prj_rad_gr_m1_last_solve_diag.linear_min_pivot_group = gidx;
            prj_rad_gr_m1_last_solve_diag.linear_min_pivot_k = k;
        }
        if (!isfinite(piv_abs) || piv_abs <= 1.0e-300) {
            return 0;
        }
        pivot[k] = p;
        if (p != k) {
            for (j = 0; j < PRJ_RAD_GR_M1_NRAD_BLOCK; ++j) {
                double tmp = A[k][j];

                A[k][j] = A[p][j];
                A[p][j] = tmp;
            }
        }
        for (i = k + 1; i < PRJ_RAD_GR_M1_NRAD_BLOCK; ++i) {
            double factor = A[i][k] / A[k][k];

            if (!isfinite(factor)) {
                return 0;
            }
            A[i][k] = factor;
            for (j = k + 1; j < PRJ_RAD_GR_M1_NRAD_BLOCK; ++j) {
                A[i][j] -= factor * A[k][j];
            }
        }
    }
    return 1;
}

static int prj_rad_gr_m1_solve4_many(
    const double A[PRJ_RAD_GR_M1_NRAD_BLOCK][PRJ_RAD_GR_M1_NRAD_BLOCK],
    const int pivot[PRJ_RAD_GR_M1_NRAD_BLOCK], int nrhs,
    double rhs[PRJ_RAD_GR_M1_NRAD_BLOCK][PRJ_RAD_GR_M1_NRAD_RHS])
{
    int col;
    int i;
    int j;
    int k;

    if (A == 0 || pivot == 0 || rhs == 0 || nrhs <= 0 ||
        nrhs > PRJ_RAD_GR_M1_NRAD_RHS) {
        return 0;
    }

    for (k = 0; k < PRJ_RAD_GR_M1_NRAD_BLOCK; ++k) {
        if (pivot[k] != k) {
            for (col = 0; col < nrhs; ++col) {
                double tmp = rhs[k][col];

                rhs[k][col] = rhs[pivot[k]][col];
                rhs[pivot[k]][col] = tmp;
            }
        }
    }

    for (i = 0; i < PRJ_RAD_GR_M1_NRAD_BLOCK; ++i) {
        for (col = 0; col < nrhs; ++col) {
            double sum = rhs[i][col];

            for (j = 0; j < i; ++j) {
                sum -= A[i][j] * rhs[j][col];
            }
            if (!isfinite(sum)) {
                return 0;
            }
            rhs[i][col] = sum;
        }
    }

    for (i = PRJ_RAD_GR_M1_NRAD_BLOCK - 1; i >= 0; --i) {
        double diag = A[i][i];

        if (!isfinite(diag) || fabs(diag) <= 1.0e-300) {
            return 0;
        }
        for (col = 0; col < nrhs; ++col) {
            double sum = rhs[i][col];

            for (j = i + 1; j < PRJ_RAD_GR_M1_NRAD_BLOCK; ++j) {
                sum -= A[i][j] * rhs[j][col];
            }
            rhs[i][col] = sum / diag;
            if (!isfinite(rhs[i][col])) {
                return 0;
            }
        }
    }
    return 1;
}

static int prj_rad_gr_m1_approx_reduced_resjac(const prj_rad *rad,
    prj_eos *eos, const prj_z4c_hydro_geom *geom, const double *u_old,
    const prj_rad_gr_m1_frozen_group frozen[PRJ_RAD_GR_M1_NGROUPS],
    const prj_rad_gr_m1_cell_context *ctx_in, const double *P_hydro,
    double dt, double *P_full, double *resid,
    double red_jac[PRJ_RAD_GR_M1_NFLUID * PRJ_RAD_GR_M1_NFLUID],
    double *u_new_out, const prj_eos_rty_interp_result *eos_in,
    prj_eos_rty_interp_result *eos_out,
    const prj_rad_gr_m1_opac_cache *opac_in,
    prj_rad_gr_m1_opac_cache *opac_out, double *norm_out,
    double norm_threshold, int candidate_only)
{
    double W[PRJ_NVAR_PRIM];
    double u_new[PRJ_NVAR_CONS];
    double resid_tmp[PRJ_NVAR_CONS];
    prj_rad_gr_m1_residual_norm_state norm_state;
    prj_rad_gr_m1_cell_context ctx_local;
    const prj_rad_gr_m1_cell_context *ctx;
    prj_rad_grm1_R_jac_geom streaming_R_jac_geom;
    const double (*g_cov)[4];
    const double (*g_con)[4];
    const double *n_cov;
    const double *ucon;
    const double *ucov;
    double ducon[4][PRJ_RAD_GR_M1_NFLUID];
    double ducov[4][PRJ_RAD_GR_M1_NFLUID];
    prj_eos_rty_interp_result eos_local;
    const prj_eos_rty_interp_result *eos_used;
    prj_eos_rty_interp_result *eos_fill;
    prj_rad_gr_m1_opac_cache opac_local;
    const prj_rad_gr_m1_opac_cache *opac_used;
    prj_rad_gr_m1_opac_cache *opac_fill;
    const double *kappa;
    const double *sigma;
    const double *delta;
    const double *eta;
    double eint;
    double pressure;
    double deint_drho;
    double deint_dT;
    double deint_dYe;
    double dpressure_drho;
    double dpressure_dT;
    double dpressure_dYe;
    double beta_con[3];
    double beta_cov[3] = {0.0, 0.0, 0.0};
    double beta2 = 0.0;
    double wlor;
    double wlor2;
    double wlor_m1;
    double Bcon[3] = {0.0, 0.0, 0.0};
    double Bcov[3] = {0.0, 0.0, 0.0};
    double Bsq = 0.0;
    double Bbeta = 0.0;
    double w;
    double D;
    double c_light = PRJ_CLIGHT;
    double c2 = PRJ_CLIGHT * PRJ_CLIGHT;
    double sqrtg;
    double alpha;
    double s_gn;
    double s_gg_f;
    double s_gu;
    double s_rgg;
    double sum_Gn = 0.0;
    double sum_Gu_xe = 0.0;
    double sum_Ggamma[3] = {0.0, 0.0, 0.0};
    double fluid_jac[PRJ_RAD_GR_M1_NFLUID][PRJ_RAD_GR_M1_NFLUID];
    double schur_sub[PRJ_RAD_GR_M1_NFLUID][PRJ_RAD_GR_M1_NFLUID];
    double group_Gn[PRJ_RAD_GR_M1_NGROUPS];
    double group_Ggamma[PRJ_RAD_GR_M1_NGROUPS][3];
    int a;
    int b;
    int d;
    int field;
    int group;
    int block_col;
    int frow;
    int fcol;
    int n;
    int v;
    int ok;
    int need_jac;
    int norm_only;

    if (resid != 0) {
        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
            resid[v] = 0.0;
        }
    }
    if (u_new_out != 0) {
        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
            u_new_out[v] = 0.0;
        }
    }
    if (norm_out != 0) {
        *norm_out = HUGE_VAL;
    }
    need_jac = red_jac != 0;
    norm_only = resid == 0 && !candidate_only;
    if (rad == 0 || eos == 0 || geom == 0 || u_old == 0 ||
        (frozen == 0 && !candidate_only) || P_hydro == 0 || P_full == 0 ||
        (!candidate_only && resid == 0 && norm_out == 0) ||
        (candidate_only && (resid != 0 || need_jac || u_new_out != 0 ||
            norm_out != 0)) || (norm_only && need_jac) ||
        (norm_only && u_new_out != 0) ||
        !isfinite(dt) ||
        !isfinite(geom->sqrt_gamma) || geom->sqrt_gamma <= 0.0 ||
        !isfinite(geom->alpha) || geom->alpha <= 0.0 ||
        !isfinite(P_hydro[0]) || P_hydro[0] <= 0.0 ||
        !isfinite(P_hydro[4]) || P_hydro[4] <= 0.0 ||
        !isfinite(P_hydro[5])) {
        return 0;
    }
    if (norm_only &&
        !prj_rad_gr_m1_residual_norm_state_init(u_old, norm_threshold,
            &norm_state)) {
        return 0;
    }
    for (v = 0; v < PRJ_RAD_GR_M1_NFLUID; ++v) {
        P_full[v] = P_hydro[v];
    }
    sqrtg = geom->sqrt_gamma;
    alpha = geom->alpha;
    s_gn = sqrtg * dt * RAD_SCALE * c_light;
    s_gg_f = alpha * sqrtg * dt * RAD_SCALE;
    s_gu = alpha * sqrtg * dt * c_light;
    s_rgg = alpha * sqrtg * dt * c_light * c_light;
    if (ctx_in != 0) {
        ctx = ctx_in;
    } else {
        if (!prj_rad_gr_m1_cell_context_init(geom, P_hydro, &ctx_local)) {
            return 0;
        }
        ctx = &ctx_local;
    }
    g_cov = ctx->g_cov;
    g_con = ctx->g_con;
    n_cov = ctx->n_cov;
    ucon = ctx->ucon;
    ucov = ctx->ucov;
    if (candidate_only && frozen == 0) {
        prj_rad_grm1_R_jac_geom_init(geom->alpha, g_con,
            geom->gamma_inv, &streaming_R_jac_geom);
    }

    if (eos_in != 0) {
        eos_used = eos_in;
        ok = 1;
    } else {
        eos_fill = eos_out != 0 ? eos_out : &eos_local;
        PRJ_TIMER_CURRENT_START("rad_matter_temp_approx_eos_table");
        ok = prj_eos_rty_interp(eos, P_hydro[0], P_hydro[4],
            P_hydro[5], eos_fill, PRJ_EOS_CTX_MAIN);
        PRJ_TIMER_CURRENT_STOP("rad_matter_temp_approx_eos_table");
        eos_used = eos_fill;
    }
    eint = eos_used->eint;
    pressure = eos_used->pressure;
    if (!ok || !isfinite(eint) || eint < 0.0 || !isfinite(pressure)) {
        return 0;
    }
    if (need_jac) {
        PRJ_TIMER_CURRENT_START("rad_matter_temp_approx_eos_deriv_convert");
        ok = prj_eos_rty_interp_derivs(eos_used, &deint_drho,
            &deint_dT, &deint_dYe, &dpressure_drho, &dpressure_dT,
            &dpressure_dYe);
        PRJ_TIMER_CURRENT_STOP("rad_matter_temp_approx_eos_deriv_convert");
        if (!ok) {
            return 0;
        }
    }
    if (opac_in != 0) {
        opac_used = opac_in;
    } else {
        opac_fill = opac_out != 0 ? opac_out : &opac_local;
        PRJ_TIMER_CURRENT_START("rad_matter_temp_approx_opac_table");
        prj_rad3_opac_lookup_interp(rad, P_hydro[0], P_hydro[4],
            P_hydro[5], opac_fill);
        PRJ_TIMER_CURRENT_STOP("rad_matter_temp_approx_opac_table");
        opac_used = opac_fill;
    }
    kappa = opac_used->kappa;
    sigma = opac_used->sigma;
    delta = opac_used->delta;
    eta = opac_used->eta;
    if (need_jac) {
        memset(fluid_jac, 0, sizeof(fluid_jac));
        memset(schur_sub, 0, sizeof(schur_sub));
    }
    if (!candidate_only) {
        for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
            W[v] = 0.0;
        }
        W[PRJ_PRIM_RHO] = P_hydro[0];
        W[PRJ_PRIM_V1] = P_hydro[1];
        W[PRJ_PRIM_V2] = P_hydro[2];
        W[PRJ_PRIM_V3] = P_hydro[3];
        W[PRJ_PRIM_EINT] = eint;
        W[PRJ_PRIM_YE] = P_hydro[5];
#if PRJ_MHD
        W[PRJ_PRIM_B1] = u_old[PRJ_CONS_B1] / sqrtg;
        W[PRJ_PRIM_B2] = u_old[PRJ_CONS_B2] / sqrtg;
        W[PRJ_PRIM_B3] = u_old[PRJ_CONS_B3] / sqrtg;
#endif
    }

    if (need_jac) {
        memset(ducon, 0, sizeof(ducon));
        memset(ducov, 0, sizeof(ducov));
        for (d = 0; d < 3; ++d) {
            beta_con[d] = P_hydro[1 + d] / c_light;
        }
        for (d = 0; d < 3; ++d) {
            for (b = 0; b < 3; ++b) {
                beta_cov[d] += geom->gamma[d][b] * beta_con[b];
            }
            beta2 += beta_cov[d] * beta_con[d];
        }
        if (!isfinite(beta2) || beta2 < 0.0 || beta2 >= 1.0) {
            return 0;
        }
        wlor = 1.0 / sqrt(1.0 - beta2);
        wlor2 = wlor * wlor;
        wlor_m1 = beta2 / (sqrt(1.0 - beta2) *
            (1.0 + sqrt(1.0 - beta2)));
#if PRJ_MHD
        Bcon[0] = u_old[PRJ_CONS_B1] / sqrtg;
        Bcon[1] = u_old[PRJ_CONS_B2] / sqrtg;
        Bcon[2] = u_old[PRJ_CONS_B3] / sqrtg;
#endif
        for (d = 0; d < 3; ++d) {
            for (b = 0; b < 3; ++b) {
                Bcov[d] += geom->gamma[d][b] * Bcon[b];
            }
            Bsq += Bcov[d] * Bcon[d];
            Bbeta += Bcov[d] * beta_con[d];
        }
        w = P_hydro[0] * c2 + P_hydro[0] * eint + pressure;
        D = P_hydro[0] * wlor;

        for (n = 0; n < PRJ_RAD_GR_M1_NFLUID; ++n) {
            double drho = n == 0 ? 1.0 : 0.0;
            double dYe = n == 5 ? 1.0 : 0.0;
            double dbeta_con[3] = {0.0, 0.0, 0.0};
            double dbeta_cov[3] = {0.0, 0.0, 0.0};
            double dbeta2 = 0.0;
            double dwlor;
            double dwlor2;
            double dBbeta = 0.0;
            double deint = 0.0;
            double dpressure = 0.0;
            double dw;
            double dD;
            double dA;
            double dUtmp;

            if (n >= 1 && n <= 3) {
                dbeta_con[n - 1] = 1.0 / c_light;
            }
            for (d = 0; d < 3; ++d) {
                int m;

                for (m = 0; m < 3; ++m) {
                    dbeta_cov[d] += geom->gamma[d][m] * dbeta_con[m];
                }
                dbeta2 += dbeta_cov[d] * beta_con[d] +
                    beta_cov[d] * dbeta_con[d];
                dBbeta += Bcov[d] * dbeta_con[d];
            }
            dwlor = 0.5 * wlor * wlor * wlor * dbeta2;
            dwlor2 = 2.0 * wlor * dwlor;
            if (n == 0) {
                deint = deint_drho;
                dpressure = dpressure_drho;
            } else if (n == 4) {
                deint = deint_dT;
                dpressure = dpressure_dT;
            } else if (n == 5) {
                deint = deint_dYe;
                dpressure = dpressure_dYe;
            }
            dw = drho * (c2 + eint) + P_hydro[0] * deint + dpressure;

            ducon[0][n] = dwlor / alpha;
            for (d = 0; d < 3; ++d) {
                ducon[d + 1][n] = dwlor *
                    (beta_con[d] - geom->beta[d] / alpha) +
                    wlor * dbeta_con[d];
            }
            for (a = 0; a < 4; ++a) {
                for (b = 0; b < 4; ++b) {
                    ducov[a][n] += g_cov[a][b] * ducon[b][n];
                }
            }

            dD = drho * wlor + P_hydro[0] * dwlor;
            fluid_jac[0][n] += sqrtg * dD;
            for (d = 0; d < 3; ++d) {
                dUtmp = ((dw * wlor2 + w * dwlor2) * beta_cov[d] +
                    (w * wlor2 + Bsq) * dbeta_cov[d] -
                    dBbeta * Bcov[d]) / c_light;
                fluid_jac[1 + d][n] += sqrtg * dUtmp;
            }
            dA = drho * eint + P_hydro[0] * deint + dpressure;
            dUtmp = dA * wlor2 +
                (P_hydro[0] * eint + pressure) * dwlor2 +
                c2 * (drho * wlor * wlor_m1 +
                    P_hydro[0] * (dwlor * wlor_m1 + wlor * dwlor)) -
                dpressure - Bbeta * dBbeta;
            if (wlor2 > 0.0) {
                dUtmp += 0.5 * Bsq * dwlor2 / (wlor2 * wlor2);
            }
            fluid_jac[4][n] += sqrtg * dUtmp;
            fluid_jac[5][n] += sqrtg * (dD * P_hydro[5] + D * dYe);
        }
    }

    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            int idx = field * PRJ_NEGROUP + group;
            int pidx = PRJ_RAD_GR_M1_NFLUID +
                PRJ_RAD_GR_M1_NRAD_BLOCK * idx;
            double xe = rad->x_e[field][group];
            const prj_rad_gr_m1_frozen_group *frozen_group;
            const double (*L)[4][4];
            prj_rad_grm1_R_projector streamed_proj;
            double streamed_dR_u[PRJ_RAD_GR_M1_NRAD_BLOCK][4];
            double streamed_dRuu[PRJ_RAD_GR_M1_NRAD_BLOCK];
            double A[PRJ_RAD_GR_M1_NRAD_BLOCK][PRJ_RAD_GR_M1_NRAD_BLOCK];
            double Bmat[PRJ_RAD_GR_M1_NRAD_BLOCK][PRJ_RAD_GR_M1_NRAD_BLOCK];
            double (*candidate_mat)[PRJ_RAD_GR_M1_NRAD_BLOCK] =
                candidate_only ? Bmat : A;
            double rhs[PRJ_RAD_GR_M1_NRAD_BLOCK][PRJ_RAD_GR_M1_NRAD_RHS] =
                {{0.0}};
            int pivot[PRJ_RAD_GR_M1_NRAD_BLOCK];
            double p[PRJ_RAD_GR_M1_NRAD_BLOCK];
            double Rcon[4][4];
            double R_u[4];
            double Ruu = 0.0;
            double Gcon[4];
            double Gcoeff[PRJ_RAD_GR_M1_NRAD_BLOCK][4];
            double Gn_coeff[PRJ_RAD_GR_M1_NRAD_BLOCK] =
                {0.0, 0.0, 0.0, 0.0};
            double Gu_coeff[PRJ_RAD_GR_M1_NRAD_BLOCK] =
                {0.0, 0.0, 0.0, 0.0};
            double Ggamma_coeff[3][PRJ_RAD_GR_M1_NRAD_BLOCK] = {{0.0}};
            double Gconst[4];
            double Gn_const = 0.0;
            double Gu_const = 0.0;
            double Ggamma_const[3] = {0.0, 0.0, 0.0};
            double de_kappa = rad->gr_m1_eleinel_dkappa[idx];
            double de_eta = rad->gr_m1_eleinel_deta[idx];
            double kappa_eff = kappa[idx] + de_kappa;
            double eta_eff = eta[idx] + de_eta;
            double sigma_eff = sigma[idx] * (1.0 - delta[idx] / 3.0);
            double kt = kappa_eff + sigma_eff;
            double scalar;
            double Gdotn = 0.0;
            double Gdotu = 0.0;
            double Ggamma[3] = {0.0, 0.0, 0.0};
            int cc;
            int r;

            if (frozen != 0) {
                frozen_group = &frozen[idx];
                L = frozen_group->dRcon;
            } else {
                double Fcov_stream[3];

                Fcov_stream[0] = P_hydro[pidx + 1];
                Fcov_stream[1] = P_hydro[pidx + 2];
                Fcov_stream[2] = P_hydro[pidx + 3];
                if (!prj_rad_grm1_R_projector_init(g_cov, g_con,
                        &streaming_R_jac_geom, P_hydro[pidx], Fcov_stream,
                        &streamed_proj)) {
                    return 0;
                }
                prj_rad_grm1_R_projector_dR_u(&streamed_proj, g_con,
                    &streaming_R_jac_geom, ucov, streamed_dR_u,
                    streamed_dRuu);
                frozen_group = 0;
                L = 0;
            }
            if (!isfinite(kappa_eff) || !isfinite(sigma_eff) ||
                !isfinite(kt) || !isfinite(eta_eff)) {
                return 0;
            }
            for (cc = 0; cc < PRJ_RAD_GR_M1_NRAD_BLOCK; ++cc) {
                double dR_u[4] = {0.0, 0.0, 0.0, 0.0};
                double dRuu = 0.0;

                if (L != 0) {
                    for (a = 0; a < 4; ++a) {
                        for (b = 0; b < 4; ++b) {
                            dR_u[a] += L[cc][a][b] * ucov[b];
                        }
                        dRuu += dR_u[a] * ucov[a];
                    }
                } else {
                    for (a = 0; a < 4; ++a) {
                        dR_u[a] = streamed_dR_u[cc][a];
                    }
                    dRuu = streamed_dRuu[cc];
                }
                for (a = 0; a < 4; ++a) {
                    Gcoeff[cc][a] = -kt * dR_u[a] -
                        sigma_eff * dRuu * ucon[a];
                    if (!isfinite(Gcoeff[cc][a])) {
                        return 0;
                    }
                    Gn_coeff[cc] += Gcoeff[cc][a] * n_cov[a];
                    Gu_coeff[cc] += Gcoeff[cc][a] * ucov[a];
                    for (d = 0; d < 3; ++d) {
                        Ggamma_coeff[d][cc] +=
                            Gcoeff[cc][a] * g_cov[d + 1][a];
                    }
                }
            }
            for (a = 0; a < 4; ++a) {
                Gconst[a] = -(eta_eff / c_light) * ucon[a];
                if (!isfinite(Gconst[a])) {
                    return 0;
                }
                Gn_const += Gconst[a] * n_cov[a];
                Gu_const += Gconst[a] * ucov[a];
                for (d = 0; d < 3; ++d) {
                    Ggamma_const[d] += Gconst[a] * g_cov[d + 1][a];
                }
            }

            for (r = 0; r < PRJ_RAD_GR_M1_NRAD_BLOCK; ++r) {
                for (cc = 0; cc < PRJ_RAD_GR_M1_NRAD_BLOCK; ++cc) {
                    candidate_mat[r][cc] = 0.0;
                }
            }
            for (cc = 0; cc < PRJ_RAD_GR_M1_NRAD_BLOCK; ++cc) {
                candidate_mat[0][cc] = (cc == 0 ? sqrtg : 0.0) -
                    s_gu * Gn_coeff[cc];
                for (d = 0; d < 3; ++d) {
                    candidate_mat[1 + d][cc] =
                        (cc == 1 + d ? sqrtg : 0.0) +
                        s_rgg * Ggamma_coeff[d][cc];
                }
            }
            rhs[0][0] = u_old[PRJ_CONS_RAD_E(field, group)] +
                s_gu * Gn_const;
            for (d = 0; d < 3; ++d) {
                rhs[1 + d][0] =
                    u_old[PRJ_CONS_RAD_F1(field, group) + d] -
                    s_rgg * Ggamma_const[d];
            }
            if (!candidate_only) {
                for (r = 0; r < PRJ_RAD_GR_M1_NRAD_BLOCK; ++r) {
                    for (cc = 0; cc < PRJ_RAD_GR_M1_NRAD_BLOCK; ++cc) {
                        Bmat[r][cc] = A[r][cc];
                    }
                }
            }
            if (!prj_rad_gr_m1_factor4(Bmat, pivot, idx) ||
                !prj_rad_gr_m1_solve4_many(Bmat, pivot, 1, rhs)) {
                return 0;
            }
            for (cc = 0; cc < PRJ_RAD_GR_M1_NRAD_BLOCK; ++cc) {
                p[cc] = rhs[cc][0];
                if (!isfinite(p[cc]) || (cc == 0 && p[cc] < 0.0)) {
                    return 0;
                }
                P_full[pidx + cc] = p[cc];
            }
            if (!candidate_only) {
                W[PRJ_PRIM_RAD_E(field, group)] = p[0];
                W[PRJ_PRIM_RAD_F1(field, group)] = p[1];
                W[PRJ_PRIM_RAD_F2(field, group)] = p[2];
                W[PRJ_PRIM_RAD_F3(field, group)] = p[3];
            } else {
                continue;
            }

            if (need_jac) {
                if (!prj_rad_gr_m1_frozen_R(frozen_group, p, Rcon)) {
                    return 0;
                }
                for (a = 0; a < 4; ++a) {
                    R_u[a] = 0.0;
                    for (b = 0; b < 4; ++b) {
                        R_u[a] += Rcon[a][b] * ucov[b];
                    }
                    Ruu += R_u[a] * ucov[a];
                }
                scalar = sigma_eff * Ruu + eta_eff / c_light;
                if (!isfinite(scalar)) {
                    return 0;
                }
                for (a = 0; a < 4; ++a) {
                    Gcon[a] = -kt * R_u[a] - scalar * ucon[a];
                    if (!isfinite(Gcon[a])) {
                        return 0;
                    }
                    Gdotn += Gcon[a] * n_cov[a];
                    Gdotu += Gcon[a] * ucov[a];
                }
                for (d = 0; d < 3; ++d) {
                    for (a = 0; a < 4; ++a) {
                        Ggamma[d] += Gcon[a] * g_cov[d + 1][a];
                    }
                }
            } else {
                Gdotn = Gn_const;
                Gdotu = Gu_const;
                for (d = 0; d < 3; ++d) {
                    Ggamma[d] = Ggamma_const[d];
                }
                for (cc = 0; cc < PRJ_RAD_GR_M1_NRAD_BLOCK; ++cc) {
                    Gdotn += p[cc] * Gn_coeff[cc];
                    Gdotu += p[cc] * Gu_coeff[cc];
                    for (d = 0; d < 3; ++d) {
                        Ggamma[d] += p[cc] * Ggamma_coeff[d][cc];
                    }
                }
                if (!isfinite(Gdotn) || !isfinite(Gdotu) ||
                    !isfinite(Ggamma[0]) || !isfinite(Ggamma[1]) ||
                    !isfinite(Ggamma[2])) {
                    return 0;
                }
            }
            group_Gn[idx] = Gdotn;
            sum_Gn += Gdotn;
            sum_Gu_xe += Gdotu * xe;
            for (d = 0; d < 3; ++d) {
                group_Ggamma[idx][d] = Ggamma[d];
                sum_Ggamma[d] += Ggamma[d];
            }

            if (need_jac) {
                double rad_fluid[PRJ_RAD_GR_M1_NRAD_BLOCK]
                    [PRJ_RAD_GR_M1_NFLUID] = {{0.0}};
                double fluid_rad[PRJ_RAD_GR_M1_NFLUID]
                    [PRJ_RAD_GR_M1_NRAD_BLOCK] = {{0.0}};
                double delta_fac = 1.0 - delta[idx] / 3.0;
                double dsigma_eff_col[PRJ_RAD_GR_M1_NFLUID] =
                    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
                double dkt_col[PRJ_RAD_GR_M1_NFLUID] =
                    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
                double deta_col[PRJ_RAD_GR_M1_NFLUID] =
                    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
                double dkappa_phys[3];
                double dsigma_phys[3];
                double ddelta_phys[3];
                double deta_phys[3];

                PRJ_TIMER_CURRENT_START(
                    "rad_matter_temp_approx_opac_deriv_convert");
                ok = prj_rad3_opac_interp_group_derivs(opac_used, idx,
                    dkappa_phys, dsigma_phys, ddelta_phys, deta_phys);
                PRJ_TIMER_CURRENT_STOP(
                    "rad_matter_temp_approx_opac_deriv_convert");
                if (!ok) {
                    return 0;
                }
                dsigma_eff_col[0] = dsigma_phys[0] * delta_fac -
                    sigma[idx] * ddelta_phys[0] / 3.0;
                dkt_col[0] = dkappa_phys[0] + dsigma_eff_col[0];
                deta_col[0] = deta_phys[0];
                dsigma_eff_col[4] = dsigma_phys[1] * delta_fac -
                    sigma[idx] * ddelta_phys[1] / 3.0;
                dkt_col[4] = dkappa_phys[1] + dsigma_eff_col[4];
                deta_col[4] = deta_phys[1];
                dsigma_eff_col[5] = dsigma_phys[2] * delta_fac -
                    sigma[idx] * ddelta_phys[2] / 3.0;
                dkt_col[5] = dkappa_phys[2] + dsigma_eff_col[5];
                deta_col[5] = deta_phys[2];

                for (block_col = 0; block_col < PRJ_RAD_GR_M1_NFLUID;
                     ++block_col) {
                    double dR_u[4] = {0.0, 0.0, 0.0, 0.0};
                    double dRuu = 0.0;
                    double dscalar;
                    double dGcon[4];
                    double dGn = 0.0;
                    double dGu = 0.0;
                    double dGgamma[3] = {0.0, 0.0, 0.0};

                    for (a = 0; a < 4; ++a) {
                        double ducov_an = ducov[a][block_col];

                        for (b = 0; b < 4; ++b) {
                            dR_u[a] += Rcon[a][b] *
                                ducov[b][block_col];
                        }
                        dRuu += dR_u[a] * ucov[a] + R_u[a] * ducov_an;
                    }
                    dscalar = dsigma_eff_col[block_col] * Ruu +
                        sigma_eff * dRuu + deta_col[block_col] / c_light;
                    for (a = 0; a < 4; ++a) {
                        double ducon_an = ducon[a][block_col];
                        double ducov_an = ducov[a][block_col];

                        dGcon[a] = -dkt_col[block_col] * R_u[a] -
                            kt * dR_u[a] - dscalar * ucon[a] -
                            scalar * ducon_an;
                        dGn += dGcon[a] * n_cov[a];
                        dGu += dGcon[a] * ucov[a] + Gcon[a] * ducov_an;
                        for (d = 0; d < 3; ++d) {
                            dGgamma[d] += dGcon[a] * g_cov[d + 1][a];
                        }
                    }

                    fluid_jac[4][block_col] += s_gn * dGn;
                    for (d = 0; d < 3; ++d) {
                        fluid_jac[1 + d][block_col] -=
                            s_gg_f * dGgamma[d];
                    }
                    fluid_jac[5][block_col] += s_gu * dGu * xe;
                    rad_fluid[0][block_col] -= s_gu * dGn;
                    for (d = 0; d < 3; ++d) {
                        rad_fluid[1 + d][block_col] += s_rgg * dGgamma[d];
                    }
                }

                for (block_col = 0; block_col < PRJ_RAD_GR_M1_NRAD_BLOCK;
                     ++block_col) {
                    double dGn = 0.0;
                    double dGu = 0.0;
                    double dGgamma[3] = {0.0, 0.0, 0.0};

                    for (a = 0; a < 4; ++a) {
                        double dGcon_a = Gcoeff[block_col][a];

                        dGn += dGcon_a * n_cov[a];
                        dGu += dGcon_a * ucov[a];
                        for (d = 0; d < 3; ++d) {
                            dGgamma[d] += dGcon_a * g_cov[d + 1][a];
                        }
                    }

                    fluid_rad[4][block_col] += s_gn * dGn;
                    for (d = 0; d < 3; ++d) {
                        fluid_rad[1 + d][block_col] -= s_gg_f * dGgamma[d];
                    }
                    fluid_rad[5][block_col] += s_gu * dGu * xe;
                }

                {
                    double rad_fluid_sol[PRJ_RAD_GR_M1_NRAD_BLOCK]
                        [PRJ_RAD_GR_M1_NRAD_RHS] = {{0.0}};

                    for (r = 0; r < PRJ_RAD_GR_M1_NRAD_BLOCK; ++r) {
                        for (fcol = 0; fcol < PRJ_RAD_GR_M1_NFLUID; ++fcol) {
                            rad_fluid_sol[r][fcol] = rad_fluid[r][fcol];
                        }
                    }
                    if (!prj_rad_gr_m1_solve4_many(Bmat, pivot,
                            PRJ_RAD_GR_M1_NFLUID, rad_fluid_sol)) {
                        return 0;
                    }
                    for (frow = 0; frow < PRJ_RAD_GR_M1_NFLUID; ++frow) {
                        for (r = 0; r < PRJ_RAD_GR_M1_NRAD_BLOCK; ++r) {
                            double fluid_to_rad = fluid_rad[frow][r];

                            for (fcol = 0; fcol < PRJ_RAD_GR_M1_NFLUID;
                                 ++fcol) {
                                schur_sub[frow][fcol] += fluid_to_rad *
                                    rad_fluid_sol[r][fcol];
                            }
                        }
                    }
                }
            }
        }
    }

    if (candidate_only) {
        return 1;
    }
    ok = prj_eos_gr_prim2cons(eos, &ctx->eos_geom, W, u_new,
        PRJ_EOS_CTX_MAIN) == PRJ_EOS_GR_OK;
    if (!ok) {
        return 0;
    }
    if (norm_only) {
        double r_rho;
        double r_mom[3];
        double r_etot;
        double r_ye;

        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
            double r = u_new[v] - u_old[v];

            if (!isfinite(r)) {
                return 0;
            }
        }
        r_rho = u_new[PRJ_CONS_RHO] - u_old[PRJ_CONS_RHO];
        r_etot = u_new[PRJ_CONS_ETOT] - u_old[PRJ_CONS_ETOT];
        r_etot += sqrtg * dt * RAD_SCALE * c_light * sum_Gn;
        for (d = 0; d < 3; ++d) {
            r_mom[d] = u_new[PRJ_CONS_MOM1 + d] -
                u_old[PRJ_CONS_MOM1 + d];
            r_mom[d] -= alpha * sqrtg * dt * RAD_SCALE * sum_Ggamma[d];
        }
        r_ye = u_new[PRJ_CONS_YE] - u_old[PRJ_CONS_YE];
        r_ye += alpha * sqrtg * dt * c_light * sum_Gu_xe;
        if (!prj_rad_gr_m1_residual_norm_add_hydro(&norm_state, r_rho,
                r_mom, r_etot, r_ye)) {
            return 0;
        }
        for (field = 0; field < PRJ_NRAD; ++field) {
            for (group = 0; group < PRJ_NEGROUP; ++group) {
                int idx = field * PRJ_NEGROUP + group;
                double r_e = u_new[PRJ_CONS_RAD_E(field, group)] -
                    u_old[PRJ_CONS_RAD_E(field, group)];
                double r_f[3];

                r_e -= alpha * sqrtg * dt * c_light * group_Gn[idx];
                for (d = 0; d < 3; ++d) {
                    r_f[d] = u_new[PRJ_CONS_RAD_F1(field, group) + d] -
                        u_old[PRJ_CONS_RAD_F1(field, group) + d];
                    r_f[d] += alpha * sqrtg * dt * c_light * c_light *
                        group_Ggamma[idx][d];
                }
                if (!prj_rad_gr_m1_residual_norm_add_rad_group(&norm_state,
                        field, group, r_e, r_f)) {
                    return 0;
                }
            }
        }
        *norm_out = norm_state.max_norm;
        return 1;
    }
    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        resid_tmp[v] = u_new[v] - u_old[v];
    }
    resid_tmp[PRJ_CONS_ETOT] += sqrtg * dt * RAD_SCALE * c_light * sum_Gn;
    resid_tmp[PRJ_CONS_MOM1] -= alpha * sqrtg * dt *
        RAD_SCALE * sum_Ggamma[0];
    resid_tmp[PRJ_CONS_MOM2] -= alpha * sqrtg * dt *
        RAD_SCALE * sum_Ggamma[1];
    resid_tmp[PRJ_CONS_MOM3] -= alpha * sqrtg * dt *
        RAD_SCALE * sum_Ggamma[2];
    resid_tmp[PRJ_CONS_YE] += alpha * sqrtg * dt * c_light * sum_Gu_xe;
    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            int idx = field * PRJ_NEGROUP + group;

            resid_tmp[PRJ_CONS_RAD_E(field, group)] -=
                alpha * sqrtg * dt * c_light * group_Gn[idx];
            resid_tmp[PRJ_CONS_RAD_F1(field, group)] +=
                alpha * sqrtg * dt * c_light * c_light *
                group_Ggamma[idx][0];
            resid_tmp[PRJ_CONS_RAD_F2(field, group)] +=
                alpha * sqrtg * dt * c_light * c_light *
                group_Ggamma[idx][1];
            resid_tmp[PRJ_CONS_RAD_F3(field, group)] +=
                alpha * sqrtg * dt * c_light * c_light *
                group_Ggamma[idx][2];
        }
    }

    if (need_jac) {
        for (frow = 0; frow < PRJ_RAD_GR_M1_NFLUID; ++frow) {
            for (fcol = 0; fcol < PRJ_RAD_GR_M1_NFLUID; ++fcol) {
                double value = fluid_jac[frow][fcol] -
                    schur_sub[frow][fcol];

                if (!isfinite(value)) {
                    return 0;
                }
                red_jac[(size_t)frow * PRJ_RAD_GR_M1_NFLUID +
                    (size_t)fcol] = value;
            }
        }
    }
    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        if (!isfinite(resid_tmp[v])) {
            return 0;
        }
        resid[v] = resid_tmp[v];
        if (u_new_out != 0) {
            u_new_out[v] = u_new[v];
        }
    }
    return 1;
}

static int prj_rad_gr_m1_solve_blocks(
    const prj_rad_gr_m1_jac_blocks *restrict blocks,
    const double *restrict resid, double *restrict x)
{
    /* The Newton matrix couples all groups through the six fluid columns, while
     * each radiation group owns an independent 4x4 block.  Eliminate those
     * group blocks and solve the resulting 6x6 Schur complement. */
    double schur[PRJ_RAD_GR_M1_NFLUID * PRJ_RAD_GR_M1_NFLUID];
    double rhs_s[PRJ_RAD_GR_M1_NFLUID];
    double fluid_x[PRJ_RAD_GR_M1_NFLUID];
    double rad_rhs_sol[PRJ_RAD_GR_M1_NGROUPS][PRJ_RAD_GR_M1_NRAD_BLOCK];
    double rad_fluid_sol[PRJ_RAD_GR_M1_NGROUPS][PRJ_RAD_GR_M1_NRAD_BLOCK][PRJ_RAD_GR_M1_NFLUID];
    int frow;
    int fcol;
    int gidx;
    int r;
    int c;

    if (blocks == 0 || resid == 0 || x == 0) {
        return 0;
    }

    for (frow = 0; frow < PRJ_RAD_GR_M1_NFLUID; ++frow) {
        int row = prj_rad_gr_m1_active_row(frow);

        rhs_s[frow] = -resid[row];
        for (fcol = 0; fcol < PRJ_RAD_GR_M1_NFLUID; ++fcol) {
            schur[(size_t)frow * PRJ_RAD_GR_M1_NFLUID + (size_t)fcol] =
                blocks->fluid[frow][fcol];
        }
    }

    for (gidx = 0; gidx < PRJ_RAD_GR_M1_NGROUPS; ++gidx) {
        int field = gidx / PRJ_NEGROUP;
        int group = gidx % PRJ_NEGROUP;
        int rad_row[PRJ_RAD_GR_M1_NRAD_BLOCK];
        double B[PRJ_RAD_GR_M1_NRAD_BLOCK][PRJ_RAD_GR_M1_NRAD_BLOCK];
        double group_sol[PRJ_RAD_GR_M1_NRAD_BLOCK][PRJ_RAD_GR_M1_NRAD_RHS];
        int pivot[PRJ_RAD_GR_M1_NRAD_BLOCK];

        rad_row[0] = PRJ_CONS_RAD_E(field, group);
        rad_row[1] = PRJ_CONS_RAD_F1(field, group);
        rad_row[2] = PRJ_CONS_RAD_F2(field, group);
        rad_row[3] = PRJ_CONS_RAD_F3(field, group);

        for (r = 0; r < PRJ_RAD_GR_M1_NRAD_BLOCK; ++r) {
            group_sol[r][0] = -resid[rad_row[r]];
            for (fcol = 0; fcol < PRJ_RAD_GR_M1_NFLUID; ++fcol) {
                group_sol[r][1 + fcol] = blocks->rad_fluid[gidx][r][fcol];
            }
            for (c = 0; c < PRJ_RAD_GR_M1_NRAD_BLOCK; ++c) {
                B[r][c] = blocks->rad_rad[gidx][r][c];
            }
        }
        if (!prj_rad_gr_m1_factor4(B, pivot, gidx) ||
            !prj_rad_gr_m1_solve4_many(B, pivot, PRJ_RAD_GR_M1_NRAD_RHS,
                group_sol)) {
            return 0;
        }
        for (r = 0; r < PRJ_RAD_GR_M1_NRAD_BLOCK; ++r) {
            rad_rhs_sol[gidx][r] = group_sol[r][0];
            for (fcol = 0; fcol < PRJ_RAD_GR_M1_NFLUID; ++fcol) {
                rad_fluid_sol[gidx][r][fcol] = group_sol[r][1 + fcol];
            }
        }

        for (frow = 0; frow < PRJ_RAD_GR_M1_NFLUID; ++frow) {
            for (r = 0; r < PRJ_RAD_GR_M1_NRAD_BLOCK; ++r) {
                double fluid_to_rad = blocks->fluid_rad[gidx][frow][r];

                rhs_s[frow] -= fluid_to_rad * rad_rhs_sol[gidx][r];
                for (fcol = 0; fcol < PRJ_RAD_GR_M1_NFLUID; ++fcol) {
                    schur[(size_t)frow * PRJ_RAD_GR_M1_NFLUID +
                        (size_t)fcol] -= fluid_to_rad *
                        rad_fluid_sol[gidx][r][fcol];
                }
            }
        }
    }

    if (!prj_rad_gr_m1_solve_dense_n(PRJ_RAD_GR_M1_NFLUID, schur, rhs_s,
            fluid_x)) {
        return 0;
    }
    for (fcol = 0; fcol < PRJ_RAD_GR_M1_NFLUID; ++fcol) {
        x[fcol] = fluid_x[fcol];
    }
    for (gidx = 0; gidx < PRJ_RAD_GR_M1_NGROUPS; ++gidx) {
        int pidx = PRJ_RAD_GR_M1_NFLUID +
            PRJ_RAD_GR_M1_NRAD_BLOCK * gidx;

        for (r = 0; r < PRJ_RAD_GR_M1_NRAD_BLOCK; ++r) {
            double value = rad_rhs_sol[gidx][r];

            for (fcol = 0; fcol < PRJ_RAD_GR_M1_NFLUID; ++fcol) {
                value -= rad_fluid_sol[gidx][r][fcol] * fluid_x[fcol];
            }
            if (!isfinite(value)) {
                return 0;
            }
            x[pidx + r] = value;
        }
    }
    return 1;
}

static int prj_rad_gr_m1_solve_blocks_active_dense(
    const prj_rad_gr_m1_jac_blocks *restrict blocks,
    const double *restrict resid, double *restrict x)
{
    const int np = PRJ_RAD_GR_M1_NP;
    double A[PRJ_RAD_GR_M1_NP * PRJ_RAD_GR_M1_NP];
    double rhs[PRJ_RAD_GR_M1_NP];
    int eq;
    int col;
    int gidx;
    int r;
    int c;

    if (blocks == 0 || resid == 0 || x == 0) {
        return 0;
    }
    for (eq = 0; eq < np; ++eq) {
        rhs[eq] = -resid[prj_rad_gr_m1_active_row(eq)];
        for (col = 0; col < np; ++col) {
            A[(size_t)eq * (size_t)np + (size_t)col] = 0.0;
        }
    }
    for (eq = 0; eq < PRJ_RAD_GR_M1_NFLUID; ++eq) {
        for (col = 0; col < PRJ_RAD_GR_M1_NFLUID; ++col) {
            A[(size_t)eq * (size_t)np + (size_t)col] =
                blocks->fluid[eq][col];
        }
    }
    for (gidx = 0; gidx < PRJ_RAD_GR_M1_NGROUPS; ++gidx) {
        int pidx = PRJ_RAD_GR_M1_NFLUID +
            PRJ_RAD_GR_M1_NRAD_BLOCK * gidx;

        for (eq = 0; eq < PRJ_RAD_GR_M1_NFLUID; ++eq) {
            for (r = 0; r < PRJ_RAD_GR_M1_NRAD_BLOCK; ++r) {
                A[(size_t)eq * (size_t)np + (size_t)(pidx + r)] =
                    blocks->fluid_rad[gidx][eq][r];
            }
        }
        for (r = 0; r < PRJ_RAD_GR_M1_NRAD_BLOCK; ++r) {
            int arow = pidx + r;

            for (col = 0; col < PRJ_RAD_GR_M1_NFLUID; ++col) {
                A[(size_t)arow * (size_t)np + (size_t)col] =
                    blocks->rad_fluid[gidx][r][col];
            }
            for (c = 0; c < PRJ_RAD_GR_M1_NRAD_BLOCK; ++c) {
                A[(size_t)arow * (size_t)np + (size_t)(pidx + c)] =
                    blocks->rad_rad[gidx][r][c];
            }
        }
    }
    return prj_rad_gr_m1_solve_dense_n(np, A, rhs, x);
}

static void prj_rad_gr_m1_diagnose_directional_jacobian(
    const prj_rad *rad, prj_eos *eos, const prj_z4c_hydro_geom *geom,
    const double *u_old, const double *P, double dt, const double *resid,
    const prj_rad_gr_m1_jac_blocks *blocks, const double *dP)
{
    const int np = PRJ_RAD_GR_M1_NP;
    double jac[PRJ_NVAR_CONS * PRJ_RAD_GR_M1_NP];
    double P_eps[PRJ_RAD_GR_M1_NP];
    double resid_eps[PRJ_NVAR_CONS];
    double u_eps[PRJ_NVAR_CONS];
    double eps;
    int attempt;
    int col;
    int eq;
    int ok = 0;

    if (rad == 0 || eos == 0 || geom == 0 || u_old == 0 || P == 0 ||
        resid == 0 || blocks == 0 || dP == 0) {
        return;
    }
    prj_rad_gr_m1_jac_blocks_to_dense(blocks, jac);
    eps = 1.0e-6 / fmax(1.0, prj_rad_gr_m1_last_solve_diag.dP_max_rel);
    eps = fmin(eps, 1.0e-7);
    eps = fmax(eps, 1.0e-24);
    for (attempt = 0; attempt < 8; ++attempt) {
        for (col = 0; col < np; ++col) {
            P_eps[col] = P[col] + eps * dP[col];
        }
        ok = prj_rad_gr_m1_residual(rad, eos, geom, u_old, P_eps, dt,
            resid_eps, u_eps, 0);
        if (ok) {
            break;
        }
        eps *= 0.1;
    }
    if (!ok) {
        prj_rad_gr_m1_last_solve_diag.directional_check_ok = 0;
        prj_rad_gr_m1_last_solve_diag.directional_eps = eps;
        return;
    }

    prj_rad_gr_m1_last_solve_diag.directional_check_ok = 1;
    prj_rad_gr_m1_last_solve_diag.directional_eps = eps;
    prj_rad_gr_m1_last_solve_diag.directional_max_abs = 0.0;
    prj_rad_gr_m1_last_solve_diag.directional_max_rel = 0.0;
    prj_rad_gr_m1_last_solve_diag.directional_jdp_norm = 0.0;
    prj_rad_gr_m1_last_solve_diag.directional_fd_norm = 0.0;
    for (eq = 0; eq < np; ++eq) {
        int row = prj_rad_gr_m1_active_row(eq);
        double jdp = 0.0;
        double fd = (resid_eps[row] - resid[row]) / eps;
        double diff;
        double scale;

        for (col = 0; col < np; ++col) {
            jdp += jac[(size_t)row * (size_t)np + (size_t)col] * dP[col];
        }
        diff = fabs(fd - jdp);
        scale = fmax(1.0, fmax(fabs(fd), fabs(jdp)));
        prj_rad_gr_m1_last_solve_diag.directional_jdp_norm += jdp * jdp;
        prj_rad_gr_m1_last_solve_diag.directional_fd_norm += fd * fd;
        if (diff > prj_rad_gr_m1_last_solve_diag.directional_max_abs) {
            prj_rad_gr_m1_last_solve_diag.directional_max_abs = diff;
            prj_rad_gr_m1_last_solve_diag.directional_row = row;
        }
        prj_rad_gr_m1_last_solve_diag.directional_max_rel =
            fmax(prj_rad_gr_m1_last_solve_diag.directional_max_rel,
                diff / scale);
    }
    prj_rad_gr_m1_last_solve_diag.directional_jdp_norm =
        sqrt(prj_rad_gr_m1_last_solve_diag.directional_jdp_norm);
    prj_rad_gr_m1_last_solve_diag.directional_fd_norm =
        sqrt(prj_rad_gr_m1_last_solve_diag.directional_fd_norm);
}

static double prj_rad_gr_m1_fd_step_for_col(const double *P, int col)
{
    double h;

    if (P == 0 || col < 0 || col >= PRJ_RAD_GR_M1_NP) {
        return 0.0;
    }
    h = 1.0e-6 * fmax(1.0, fabs(P[col]));
    if (col >= PRJ_RAD_GR_M1_NFLUID) {
        int local = col - PRJ_RAD_GR_M1_NFLUID;
        int comp = local % PRJ_RAD_GR_M1_NRAD_BLOCK;

        h = 1.0e-5 * fmax(1.0, fabs(P[col]));
        if (comp == 0) {
            h = fmin(h, 0.25 * P[col]);
        }
        /* comp == 0 is the radiation energy E (kept positive by the cap
         * above); comp > 0 is a covariant flux component F_i for which the
         * plain relative step is appropriate. */
    }
    if (!isfinite(h) || h <= 0.0) {
        return 0.0;
    }
    return h;
}

static void prj_rad_gr_m1_diagnose_elementwise_jacobian(
    const prj_rad *rad, prj_eos *eos, const prj_z4c_hydro_geom *geom,
    const double *u_old, const double *P, double dt,
    const prj_rad_gr_m1_jac_blocks *blocks)
{
    const int np = PRJ_RAD_GR_M1_NP;
    double jac[PRJ_NVAR_CONS * PRJ_RAD_GR_M1_NP];
    double Pp[PRJ_RAD_GR_M1_NP];
    double Pm[PRJ_RAD_GR_M1_NP];
    double resp[PRJ_NVAR_CONS];
    double resm[PRJ_NVAR_CONS];
    double u_tmp[PRJ_NVAR_CONS];
    int col;
    int eq;

    /* TEMP DIAGNOSTIC: remove after rad-matter Jacobian debugging. */
    if (rad == 0 || eos == 0 || geom == 0 || u_old == 0 || P == 0 ||
        blocks == 0) {
        return;
    }
    prj_rad_gr_m1_jac_blocks_to_dense(blocks, jac);
    prj_rad_gr_m1_last_solve_diag.fd_check_ok = 1;
    prj_rad_gr_m1_last_solve_diag.fd_skipped_cols = 0;
    prj_rad_gr_m1_last_solve_diag.fd_max_abs = 0.0;
    prj_rad_gr_m1_last_solve_diag.fd_max_rel = 0.0;
    prj_rad_gr_m1_last_solve_diag.fd_max_abs_row = -1;
    prj_rad_gr_m1_last_solve_diag.fd_max_abs_col = -1;
    prj_rad_gr_m1_last_solve_diag.fd_max_rel_row = -1;
    prj_rad_gr_m1_last_solve_diag.fd_max_rel_col = -1;

    for (col = 0; col < np; ++col) {
        double h = prj_rad_gr_m1_fd_step_for_col(P, col);
        int ok = 0;
        int attempt;

        for (attempt = 0; attempt < 12 && h > 0.0; ++attempt) {
            memcpy(Pp, P, sizeof(Pp));
            memcpy(Pm, P, sizeof(Pm));
            Pp[col] += h;
            Pm[col] -= h;
            ok = prj_rad_gr_m1_residual(rad, eos, geom, u_old, Pp, dt,
                    resp, u_tmp, 0) &&
                prj_rad_gr_m1_residual(rad, eos, geom, u_old, Pm, dt,
                    resm, u_tmp, 0);
            if (ok) {
                break;
            }
            h *= 0.1;
        }
        if (!ok) {
            prj_rad_gr_m1_last_solve_diag.fd_check_ok = 0;
            prj_rad_gr_m1_last_solve_diag.fd_skipped_cols += 1;
            continue;
        }

        for (eq = 0; eq < np; ++eq) {
            int row = prj_rad_gr_m1_active_row(eq);
            double fd = (resp[row] - resm[row]) / (2.0 * h);
            double got = jac[(size_t)row * (size_t)np + (size_t)col];
            double diff = fabs(fd - got);
            double scale = fmax(1.0, fmax(fabs(fd), fabs(got)));
            double rel = diff / scale;

            if (!isfinite(fd) || !isfinite(got)) {
                prj_rad_gr_m1_last_solve_diag.fd_check_ok = 0;
                continue;
            }
            if (diff > prj_rad_gr_m1_last_solve_diag.fd_max_abs) {
                prj_rad_gr_m1_last_solve_diag.fd_max_abs = diff;
                prj_rad_gr_m1_last_solve_diag.fd_max_abs_row = row;
                prj_rad_gr_m1_last_solve_diag.fd_max_abs_col = col;
                prj_rad_gr_m1_last_solve_diag.fd_max_abs_got = got;
                prj_rad_gr_m1_last_solve_diag.fd_max_abs_expected = fd;
                prj_rad_gr_m1_last_solve_diag.fd_max_abs_h = h;
            }
            if (rel > prj_rad_gr_m1_last_solve_diag.fd_max_rel) {
                prj_rad_gr_m1_last_solve_diag.fd_max_rel = rel;
                prj_rad_gr_m1_last_solve_diag.fd_max_rel_row = row;
                prj_rad_gr_m1_last_solve_diag.fd_max_rel_col = col;
                prj_rad_gr_m1_last_solve_diag.fd_max_rel_got = got;
                prj_rad_gr_m1_last_solve_diag.fd_max_rel_expected = fd;
                prj_rad_gr_m1_last_solve_diag.fd_max_rel_h = h;
            }
        }
    }
}

/* Closed-form Eulerian (u = n) radiation candidate for the matter coupling.
 * In the Eulerian normal frame the per-group 4x4 implicit system built by
 * prj_rad_gr_m1_approx_reduced_resjac becomes exactly diagonal and closure-
 * independent: every M1 closure-derivative term is contracted with the spatial
 * part of n_cov (which is zero) and the flux columns of the geometry tensor T
 * vanish, leaving a scalar absorption/emission update for E and an independent
 * scalar drag on each flux component,
 *
 *   E_new   = (E_old + s_gu * eta/c) / (sqrtg + s_gu * kappa_eff)
 *   F_new_d =  F_old_d              / (sqrtg + s_rgg * kt/c)
 *
 * with kappa_eff = kappa, sigma_eff = sigma*(1 - delta/3), kt = kappa_eff +
 * sigma_eff, and the same s_gu = alpha*sqrtg*dt*c, s_rgg = alpha*sqrtg*dt*c^2
 * as the general assembly.  This reproduces the general candidate at u = n up
 * to LU-vs-division rounding, so the exact residual check (which always uses
 * the true four-velocity) still gates acceptance.  Used only below the velocity
 * gate; it skips the four-velocity boost, the closure projector, and the 4x4
 * factor/solve entirely, but still fills the eos/opac caches the exact check
 * reuses.  Returns 0 (fall through to the general path) on any invalid input. */
static int prj_rad_gr_m1_eulerian_candidate(const prj_rad *rad,
    prj_eos *eos, const prj_z4c_hydro_geom *geom, const double *u_old,
    double dt, const double *P_hydro, double *P_full,
    prj_eos_rty_interp_result *eos_out,
    prj_rad_gr_m1_opac_cache *opac_out)
{
    double c_light = PRJ_CLIGHT;
    double sqrtg;
    double alpha;
    double s_gu;
    double s_rgg;
    double eint;
    double pressure;
    int field;
    int group;
    int v;
    int ok;

    if (rad == 0 || eos == 0 || geom == 0 || u_old == 0 || P_hydro == 0 ||
        P_full == 0 || eos_out == 0 || opac_out == 0 || !isfinite(dt) ||
        !isfinite(geom->sqrt_gamma) || geom->sqrt_gamma <= 0.0 ||
        !isfinite(geom->alpha) || geom->alpha <= 0.0 ||
        !isfinite(P_hydro[0]) || P_hydro[0] <= 0.0 ||
        !isfinite(P_hydro[4]) || P_hydro[4] <= 0.0 ||
        !isfinite(P_hydro[5])) {
        return 0;
    }
    sqrtg = geom->sqrt_gamma;
    alpha = geom->alpha;
    s_gu = alpha * sqrtg * dt * c_light;
    s_rgg = alpha * sqrtg * dt * c_light * c_light;

    for (v = 0; v < PRJ_RAD_GR_M1_NFLUID; ++v) {
        P_full[v] = P_hydro[v];
    }

    ok = prj_eos_rty_interp(eos, P_hydro[0], P_hydro[4], P_hydro[5],
        eos_out, PRJ_EOS_CTX_MAIN);
    eint = eos_out->eint;
    pressure = eos_out->pressure;
    if (!ok || !isfinite(eint) || eint < 0.0 || !isfinite(pressure)) {
        return 0;
    }
    prj_rad3_opac_lookup_interp(rad, P_hydro[0], P_hydro[4], P_hydro[5],
        opac_out);

    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            int idx = field * PRJ_NEGROUP + group;
            int pidx = PRJ_RAD_GR_M1_NFLUID +
                PRJ_RAD_GR_M1_NRAD_BLOCK * idx;
            double kappa_eff = opac_out->kappa[idx] +
                rad->gr_m1_eleinel_dkappa[idx];
            double sigma_eff = opac_out->sigma[idx] *
                (1.0 - opac_out->delta[idx] / 3.0);
            double kt = kappa_eff + sigma_eff;
            double eta = opac_out->eta[idx] + rad->gr_m1_eleinel_deta[idx];
            double E_old = u_old[PRJ_CONS_RAD_E(field, group)];
            double E_den = sqrtg + s_gu * kappa_eff;
            double F_den = sqrtg + s_rgg * kt / c_light;
            double E_new;
            int d;

            if (!isfinite(kappa_eff) || !isfinite(sigma_eff) ||
                !isfinite(kt) || !isfinite(eta) || !isfinite(E_old) ||
                !(E_den > 0.0) || !(F_den > 0.0)) {
                return 0;
            }
            E_new = (E_old + s_gu * eta / c_light) / E_den;
            if (!isfinite(E_new) || E_new < 0.0) {
                return 0;
            }
            P_full[pidx] = E_new;
            for (d = 0; d < 3; ++d) {
                double F_old = u_old[PRJ_CONS_RAD_F1(field, group) + d];
                double F_new = F_old / F_den;

                if (!isfinite(F_new)) {
                    return 0;
                }
                P_full[pidx + 1 + d] = F_new;
            }
        }
    }
    return 1;
}

static int prj_rad_gr_m1_implicit_solve_approx(const prj_rad *rad,
    prj_eos *eos, const prj_z4c_hydro_geom *geom, const double *u_old,
    double dt, double *P, double *resid_out, double *u_new_out,
    double *W_out)
{
    const int np = PRJ_RAD_GR_M1_NP;
    double threshold = PRJ_RAD_GR_M1_SOLVE_TOL_DEFAULT;
    int maxiter = PRJ_RAD_GR_M1_SOLVE_MAXITER_DEFAULT;
    prj_rad_gr_m1_frozen_group frozen[PRJ_RAD_GR_M1_NGROUPS];
    prj_rad_gr_m1_cell_context fast_ctx;
    int use_eulerian = 0;
    prj_eos_rty_interp_result eos_cache;
    prj_eos_rty_interp_result eos_trial_cache;
    prj_rad_gr_m1_opac_cache opac_cache;
    prj_rad_gr_m1_opac_cache opac_trial_cache;
    double P_current[PRJ_RAD_GR_M1_NP];
    double P_eval[PRJ_RAD_GR_M1_NP];
    double P_trial_in[PRJ_RAD_GR_M1_NP];
    double P_trial[PRJ_RAD_GR_M1_NP];
    double resid[PRJ_NVAR_CONS];
    double resid_exact[PRJ_NVAR_CONS];
    double u_new[PRJ_NVAR_CONS];
    double u_exact[PRJ_NVAR_CONS];
    double W_exact[PRJ_NVAR_PRIM];
    double red_jac[PRJ_RAD_GR_M1_NFLUID * PRJ_RAD_GR_M1_NFLUID];
    double rhs[PRJ_RAD_GR_M1_NFLUID];
    double dfluid[PRJ_RAD_GR_M1_NFLUID];
    int iter;
    int col;
    int ok;
    int have_cache = 0;
    int force_jac = 0;

    if (rad == 0 || eos == 0 || geom == 0 || u_old == 0 || P == 0 ||
        !isfinite(dt) || !isfinite(u_old[PRJ_CONS_YE]) ||
        u_old[PRJ_CONS_YE] <= 0.0) {
        return 0;
    }
    if (isfinite(rad->implicit_err_tol) && rad->implicit_err_tol > 0.0) {
        threshold = rad->implicit_err_tol;
    }
    if (rad->maxiter > 0) {
        maxiter = rad->maxiter;
    }
    if (maxiter > PRJ_RAD_GR_M1_APPROX_MAXITER) {
        maxiter = PRJ_RAD_GR_M1_APPROX_MAXITER;
    }
    if (maxiter <= 0) {
        return 0;
    }
    prj_rad_gr_m1_solve_diag_reset(threshold, maxiter);
    memcpy(P_current, P, sizeof(P_current));
    if (!prj_rad_gr_m1_cell_context_init(geom, P_current, &fast_ctx)) {
        return 0;
    }

    /* Velocity gate: for slow cells solve the fast candidate in the Eulerian
     * normal frame (u = n) via the closed-form scalar update, which drops the
     * O(v/c) boost, the M1 closure projector, and the per-group 4x4 solve.  The
     * exact residual check below always uses the true four-velocity (fast_ctx),
     * so a mis-gated cell merely fails the check and falls through to the full
     * solve -- correctness is intact. */
    {
        double vhat[3];
        double v2 = 0.0;
        int gi;
        int gj;

        for (gi = 0; gi < 3; ++gi) {
            vhat[gi] = P_current[1 + gi] / PRJ_CLIGHT;
        }
        for (gi = 0; gi < 3; ++gi) {
            for (gj = 0; gj < 3; ++gj) {
                v2 += geom->gamma[gi][gj] * vhat[gi] * vhat[gj];
            }
        }
        use_eulerian = isfinite(v2) && v2 >= 0.0 &&
            v2 < PRJ_RAD_GR_M1_EULERIAN_V2_GATE;
    }

    PRJ_TIMER_CURRENT_START("rad_matter_temp_approx_fast_candidate");
    if (use_eulerian) {
        ok = prj_rad_gr_m1_eulerian_candidate(rad, eos, geom, u_old, dt,
            P_current, P_eval, &eos_cache, &opac_cache);
    } else {
        ok = prj_rad_gr_m1_approx_reduced_resjac(rad, eos, geom, u_old, 0,
            &fast_ctx, P_current, dt, P_eval, 0, 0, 0, 0, &eos_cache, 0,
            &opac_cache, 0, threshold, 1);
    }
    PRJ_TIMER_CURRENT_STOP("rad_matter_temp_approx_fast_candidate");
    if (!ok) {
        return 0;
    }
    have_cache = 1;
    {
        double exact_norm;

        PRJ_TIMER_CURRENT_START("rad_matter_temp_approx_exact_check");
        ok = prj_rad_gr_m1_exact_residual_check(rad, eos, geom, u_old,
            P_eval, dt, &eos_cache, &opac_cache, &fast_ctx, threshold,
            &exact_norm, resid_exact, u_exact, W_exact);
        PRJ_TIMER_CURRENT_STOP("rad_matter_temp_approx_exact_check");
        if (ok && isfinite(exact_norm) && exact_norm < threshold) {
            memcpy(P, P_eval, (size_t)np * sizeof(double));
            if (resid_out != 0) {
                memcpy(resid_out, resid_exact, sizeof(resid_exact));
            }
            if (u_new_out != 0) {
                memcpy(u_new_out, u_exact, sizeof(u_exact));
            }
            if (W_out != 0) {
                memcpy(W_out, W_exact, sizeof(W_exact));
            }
            prj_rad_gr_m1_last_solve_diag.norm = exact_norm;
            return 1;
        }
    }

    PRJ_TIMER_CURRENT_START("rad_matter_temp_approx_materialize_closure");
    ok = prj_rad_gr_m1_freeze_closure(geom, P_current, frozen);
    PRJ_TIMER_CURRENT_STOP("rad_matter_temp_approx_materialize_closure");
    if (!ok) {
        return 0;
    }
    force_jac = 1;

    for (iter = 0; iter < maxiter; ++iter) {
        double norm;
        int accepted = 0;
        int ls;

        prj_rad_gr_m1_last_solve_diag.iter = iter;
        if (!force_jac) {
            PRJ_TIMER_CURRENT_START("rad_matter_temp_approx_resid");
            ok = prj_rad_gr_m1_approx_reduced_resjac(rad, eos, geom, u_old,
                frozen, 0, P_current, dt, P_eval, 0, 0, 0,
                have_cache ? &eos_cache : 0,
                have_cache ? 0 : &eos_cache,
                have_cache ? &opac_cache : 0,
                have_cache ? 0 : &opac_cache, &norm, threshold, 0);
            PRJ_TIMER_CURRENT_STOP("rad_matter_temp_approx_resid");
            if (!ok) {
                return 0;
            }
            prj_rad_gr_m1_last_solve_diag.norm = norm;
            if (!isfinite(norm)) {
                return 0;
            }
            if (norm < threshold) {
                double exact_norm;

                PRJ_TIMER_CURRENT_START(
                    "rad_matter_temp_approx_exact_check");
                ok = prj_rad_gr_m1_exact_residual_check(rad, eos, geom,
                    u_old, P_eval, dt, &eos_cache, &opac_cache, 0,
                    threshold, &exact_norm, resid_exact, u_exact, W_exact);
                PRJ_TIMER_CURRENT_STOP(
                    "rad_matter_temp_approx_exact_check");
                if (ok && isfinite(exact_norm) && exact_norm < threshold) {
                    memcpy(P, P_eval, (size_t)np * sizeof(double));
                    if (resid_out != 0) {
                        memcpy(resid_out, resid_exact, sizeof(resid_exact));
                    }
                    if (u_new_out != 0) {
                        memcpy(u_new_out, u_exact, sizeof(u_exact));
                    }
                    if (W_out != 0) {
                        memcpy(W_out, W_exact, sizeof(W_exact));
                    }
                    prj_rad_gr_m1_last_solve_diag.norm = exact_norm;
                    return 1;
                }
                return 0;
            }
        }
        force_jac = 0;

        PRJ_TIMER_CURRENT_START("rad_matter_temp_approx_resjac");
        ok = prj_rad_gr_m1_approx_reduced_resjac(rad, eos, geom, u_old,
            frozen, 0, P_current, dt, P_eval, resid, red_jac, u_new,
            &eos_cache, 0, &opac_cache, 0, 0, threshold, 0);
        PRJ_TIMER_CURRENT_STOP("rad_matter_temp_approx_resjac");
        if (!ok) {
            return 0;
        }
        PRJ_TIMER_CURRENT_START("rad_matter_temp_approx_norm");
        norm = prj_rad_gr_m1_residual_norm(u_old, resid, threshold);
        PRJ_TIMER_CURRENT_STOP("rad_matter_temp_approx_norm");
        prj_rad_gr_m1_last_solve_diag.norm = norm;
        if (!isfinite(norm)) {
            return 0;
        }
        if (norm < threshold) {
            double exact_norm;

            PRJ_TIMER_CURRENT_START("rad_matter_temp_approx_exact_check");
            ok = prj_rad_gr_m1_exact_residual_check(rad, eos, geom, u_old,
                P_eval, dt, &eos_cache, &opac_cache, 0, threshold,
                &exact_norm, resid_exact, u_exact, W_exact);
            PRJ_TIMER_CURRENT_STOP("rad_matter_temp_approx_exact_check");
            if (ok && isfinite(exact_norm) && exact_norm < threshold) {
                memcpy(P, P_eval, (size_t)np * sizeof(double));
                if (resid_out != 0) {
                    memcpy(resid_out, resid_exact, sizeof(resid_exact));
                }
                if (u_new_out != 0) {
                    memcpy(u_new_out, u_exact, sizeof(u_exact));
                }
                if (W_out != 0) {
                    memcpy(W_out, W_exact, sizeof(W_exact));
                }
                prj_rad_gr_m1_last_solve_diag.norm = exact_norm;
                return 1;
            }
            return 0;
        }
        for (col = 0; col < PRJ_RAD_GR_M1_NFLUID; ++col) {
            rhs[col] = -resid[prj_rad_gr_m1_active_row(col)];
        }
        PRJ_TIMER_CURRENT_START("rad_matter_temp_approx_linear");
        ok = prj_rad_gr_m1_solve_dense_n(PRJ_RAD_GR_M1_NFLUID, red_jac,
            rhs, dfluid);
        PRJ_TIMER_CURRENT_STOP("rad_matter_temp_approx_linear");
        prj_rad_gr_m1_last_solve_diag.dense_solve_ok = ok;
        if (!ok) {
            return 0;
        }
        prj_rad_gr_m1_last_solve_diag.dP_max_abs = 0.0;
        prj_rad_gr_m1_last_solve_diag.dP_max_rel = 0.0;
        prj_rad_gr_m1_last_solve_diag.dP_max_abs_col = -1;
        prj_rad_gr_m1_last_solve_diag.dP_max_rel_col = -1;
        for (col = 0; col < PRJ_RAD_GR_M1_NFLUID; ++col) {
            double abs_dP = fabs(dfluid[col]);
            double rel_dP = abs_dP / fmax(1.0, fabs(P_eval[col]));

            if (abs_dP > prj_rad_gr_m1_last_solve_diag.dP_max_abs) {
                prj_rad_gr_m1_last_solve_diag.dP_max_abs = abs_dP;
                prj_rad_gr_m1_last_solve_diag.dP_max_abs_col = col;
            }
            if (rel_dP > prj_rad_gr_m1_last_solve_diag.dP_max_rel) {
                prj_rad_gr_m1_last_solve_diag.dP_max_rel = rel_dP;
                prj_rad_gr_m1_last_solve_diag.dP_max_rel_col = col;
            }
        }

        PRJ_TIMER_CURRENT_START("rad_matter_temp_approx_linesearch");
        prj_rad_gr_m1_last_solve_diag.line_search_trials = 0;
        prj_rad_gr_m1_last_solve_diag.invalid_trials = 0;
        prj_rad_gr_m1_last_solve_diag.best_ls = -1;
        prj_rad_gr_m1_last_solve_diag.first_trial_norm = HUGE_VAL;
        prj_rad_gr_m1_last_solve_diag.best_trial_norm = HUGE_VAL;
        prj_rad_gr_m1_last_solve_diag.last_trial_norm = HUGE_VAL;
        prj_rad_gr_m1_last_solve_diag.best_lambda = 0.0;
        for (ls = 0; ls < PRJ_RAD_GR_M1_LINESEARCH_MAX; ++ls) {
            double lambda = ldexp(1.0, -ls);
            double trial_norm;

            prj_rad_gr_m1_last_solve_diag.line_search_trials += 1;
            memcpy(P_trial_in, P_eval, sizeof(P_trial_in));
            for (col = 0; col < PRJ_RAD_GR_M1_NFLUID; ++col) {
                P_trial_in[col] = P_eval[col] + lambda * dfluid[col];
            }
            ok = prj_rad_gr_m1_approx_reduced_resjac(rad, eos, geom, u_old,
                frozen, 0, P_trial_in, dt, P_trial, 0, 0, 0, 0,
                &eos_trial_cache, 0, &opac_trial_cache,
                &trial_norm, threshold, 0);
            if (!ok) {
                prj_rad_gr_m1_last_solve_diag.invalid_trials += 1;
                continue;
            }
            if (prj_rad_gr_m1_last_solve_diag.first_trial_norm == HUGE_VAL) {
                prj_rad_gr_m1_last_solve_diag.first_trial_norm = trial_norm;
            }
            prj_rad_gr_m1_last_solve_diag.last_trial_norm = trial_norm;
            if (isfinite(trial_norm) &&
                trial_norm < prj_rad_gr_m1_last_solve_diag.best_trial_norm) {
                prj_rad_gr_m1_last_solve_diag.best_trial_norm = trial_norm;
                prj_rad_gr_m1_last_solve_diag.best_ls = ls;
                prj_rad_gr_m1_last_solve_diag.best_lambda = lambda;
            }
            if (isfinite(trial_norm) && trial_norm < norm) {
                memcpy(P_current, P_trial, sizeof(P_current));
                eos_cache = eos_trial_cache;
                opac_cache = opac_trial_cache;
                have_cache = 1;
                norm = trial_norm;
                accepted = 1;
                break;
            }
        }
        PRJ_TIMER_CURRENT_STOP("rad_matter_temp_approx_linesearch");
        if (!accepted) {
            return 0;
        }
    }
    return 0;
}

static int prj_rad_gr_m1_implicit_solve_exact(const prj_rad *rad, prj_eos *eos,
    const prj_z4c_hydro_geom *geom, const double *u_old, double dt, double *P,
    double *resid_out, double *u_new_out, double *W_out)
{
    const int np = PRJ_RAD_GR_M1_NP;
    double threshold = PRJ_RAD_GR_M1_SOLVE_TOL_DEFAULT;
    int maxiter = PRJ_RAD_GR_M1_SOLVE_MAXITER_DEFAULT;
    double resid[PRJ_NVAR_CONS];
    double resid_trial[PRJ_NVAR_CONS];
    double u_new[PRJ_NVAR_CONS];
    double u_new_trial[PRJ_NVAR_CONS];
    double W_new[PRJ_NVAR_PRIM];
    double W_trial[PRJ_NVAR_PRIM];
    prj_rad_gr_m1_jac_blocks jac_blocks;
    double dP[PRJ_RAD_GR_M1_NP];
    double P_trial[PRJ_RAD_GR_M1_NP];
    int iter;
    int col;
    int ok;

    if (rad == 0 || eos == 0 || geom == 0 || u_old == 0 || P == 0 ||
        !isfinite(dt) || !isfinite(u_old[PRJ_CONS_YE]) ||
        u_old[PRJ_CONS_YE] <= 0.0) {
        return 0;
    }
    if (isfinite(rad->implicit_err_tol) && rad->implicit_err_tol > 0.0) {
        threshold = rad->implicit_err_tol;
    }
    if (rad->maxiter > 0) {
        maxiter = rad->maxiter;
    }
    prj_rad_gr_m1_solve_diag_reset(threshold, maxiter);

    for (iter = 0; iter < maxiter; ++iter) {
        double norm;
        int accepted = 0;
        int ls;

        prj_rad_gr_m1_last_solve_diag.iter = iter;

        /* TEMP TIMER: remove after rad-matter coupling profiling. */
        PRJ_TIMER_CURRENT_START("rad_matter_temp_newton_resjac");
        ok = prj_rad_gr_m1_residual_jacobian(rad, eos, geom, u_old, P, dt,
            resid, &jac_blocks, 0, u_new, W_new);
        PRJ_TIMER_CURRENT_STOP("rad_matter_temp_newton_resjac");
        if (!ok) {
            return 0;
        }
        /* TEMP TIMER: remove after rad-matter coupling profiling. */
        PRJ_TIMER_CURRENT_START("rad_matter_temp_newton_norm");
        norm = prj_rad_gr_m1_residual_norm(u_old, resid, threshold);
        PRJ_TIMER_CURRENT_STOP("rad_matter_temp_newton_norm");
        prj_rad_gr_m1_last_solve_diag.norm = norm;
        if (isfinite(norm) && norm < threshold) {
            if (resid_out != 0) {
                memcpy(resid_out, resid, sizeof(resid));
            }
            if (u_new_out != 0) {
                memcpy(u_new_out, u_new, sizeof(u_new));
            }
            if (W_out != 0) {
                memcpy(W_out, W_new, sizeof(W_new));
            }
            return 1;
        }
        if (!isfinite(norm)) {
            return 0;
        }

        /* In the lab-frame (E, F_i) parametrization the radiation transport
         * block is the constant sqrt(gamma) I, so each per-group 4x4 block is
         * diagonally dominant and invertible even at near-vacuum -- the old
         * near-singular uR^i diagonal regularization is no longer needed. */

        /* TEMP TIMER: remove after rad-matter coupling profiling. */
        PRJ_TIMER_CURRENT_START("rad_matter_temp_newton_linear");
        prj_rad_gr_m1_last_solve_diag.linear_min_pivot = HUGE_VAL;
        prj_rad_gr_m1_last_solve_diag.linear_min_pivot_rel = HUGE_VAL;
        ok = prj_rad_gr_m1_solve_blocks(&jac_blocks, resid, dP);
        prj_rad_gr_m1_last_solve_diag.block_solve_ok = ok;
        if (!ok) {
            ok = prj_rad_gr_m1_solve_blocks_active_dense(&jac_blocks, resid,
                dP);
            prj_rad_gr_m1_last_solve_diag.dense_solve_ok = ok;
        } else {
            prj_rad_gr_m1_last_solve_diag.dense_solve_ok = 0;
        }
        PRJ_TIMER_CURRENT_STOP("rad_matter_temp_newton_linear");
        if (!ok) {
            return 0;
        }
        prj_rad_gr_m1_last_solve_diag.dP_max_abs = 0.0;
        prj_rad_gr_m1_last_solve_diag.dP_max_rel = 0.0;
        prj_rad_gr_m1_last_solve_diag.dP_max_abs_col = -1;
        prj_rad_gr_m1_last_solve_diag.dP_max_rel_col = -1;
        for (col = 0; col < np; ++col) {
            double abs_dP = fabs(dP[col]);
            double rel_dP = abs_dP / fmax(1.0, fabs(P[col]));

            if (abs_dP > prj_rad_gr_m1_last_solve_diag.dP_max_abs) {
                prj_rad_gr_m1_last_solve_diag.dP_max_abs = abs_dP;
                prj_rad_gr_m1_last_solve_diag.dP_max_abs_col = col;
            }
            if (rel_dP > prj_rad_gr_m1_last_solve_diag.dP_max_rel) {
                prj_rad_gr_m1_last_solve_diag.dP_max_rel = rel_dP;
                prj_rad_gr_m1_last_solve_diag.dP_max_rel_col = col;
            }
        }

        /* TEMP TIMER: remove after rad-matter coupling profiling. */
        PRJ_TIMER_CURRENT_START("rad_matter_temp_newton_linesearch");
        prj_rad_gr_m1_last_solve_diag.line_search_trials = 0;
        prj_rad_gr_m1_last_solve_diag.invalid_trials = 0;
        prj_rad_gr_m1_last_solve_diag.best_ls = -1;
        prj_rad_gr_m1_last_solve_diag.first_trial_norm = HUGE_VAL;
        prj_rad_gr_m1_last_solve_diag.best_trial_norm = HUGE_VAL;
        prj_rad_gr_m1_last_solve_diag.last_trial_norm = HUGE_VAL;
        prj_rad_gr_m1_last_solve_diag.best_lambda = 0.0;
        for (ls = 0; ls < PRJ_RAD_GR_M1_LINESEARCH_MAX; ++ls) {
            double lambda = ldexp(1.0, -ls);
            double trial_norm;

            prj_rad_gr_m1_last_solve_diag.line_search_trials += 1;
            for (col = 0; col < np; ++col) {
                P_trial[col] = P[col] + lambda * dP[col];
            }
            /* TEMP TIMER: remove after rad-matter coupling profiling. */
            PRJ_TIMER_CURRENT_START("rad_matter_temp_newton_trial_resid");
            ok = prj_rad_gr_m1_residual(rad, eos, geom, u_old, P_trial, dt,
                resid_trial, u_new_trial, W_trial);
            PRJ_TIMER_CURRENT_STOP("rad_matter_temp_newton_trial_resid");
            if (!ok) {
                prj_rad_gr_m1_last_solve_diag.invalid_trials += 1;
                continue;
            }
            /* TEMP TIMER: remove after rad-matter coupling profiling. */
            PRJ_TIMER_CURRENT_START("rad_matter_temp_newton_trial_norm");
            trial_norm = prj_rad_gr_m1_residual_norm(u_old, resid_trial,
                threshold);
            PRJ_TIMER_CURRENT_STOP("rad_matter_temp_newton_trial_norm");
            if (prj_rad_gr_m1_last_solve_diag.first_trial_norm == HUGE_VAL) {
                prj_rad_gr_m1_last_solve_diag.first_trial_norm = trial_norm;
            }
            prj_rad_gr_m1_last_solve_diag.last_trial_norm = trial_norm;
            if (isfinite(trial_norm) &&
                trial_norm < prj_rad_gr_m1_last_solve_diag.best_trial_norm) {
                prj_rad_gr_m1_last_solve_diag.best_trial_norm = trial_norm;
                prj_rad_gr_m1_last_solve_diag.best_ls = ls;
                prj_rad_gr_m1_last_solve_diag.best_lambda = lambda;
            }
            if (isfinite(trial_norm) && trial_norm < norm) {
                memcpy(P, P_trial, (size_t)np * sizeof(double));
                memcpy(resid, resid_trial, sizeof(resid));
                memcpy(u_new, u_new_trial, sizeof(u_new));
                memcpy(W_new, W_trial, sizeof(W_new));
                norm = trial_norm;
                accepted = 1;
                break;
            }
        }
        PRJ_TIMER_CURRENT_STOP("rad_matter_temp_newton_linesearch");
        if (!accepted) {
            prj_rad_gr_m1_diagnose_directional_jacobian(rad, eos, geom,
                u_old, P, dt, resid, &jac_blocks, dP);
            prj_rad_gr_m1_diagnose_elementwise_jacobian(rad, eos, geom,
                u_old, P, dt, &jac_blocks);
            if (resid_out != 0) {
                memcpy(resid_out, resid, sizeof(resid));
            }
            if (u_new_out != 0) {
                memcpy(u_new_out, u_new, sizeof(u_new));
            }
            return 0;
        }
        if (norm < threshold) {
            if (resid_out != 0) {
                memcpy(resid_out, resid, sizeof(resid));
            }
            if (u_new_out != 0) {
                memcpy(u_new_out, u_new, sizeof(u_new));
            }
            if (W_out != 0) {
                memcpy(W_out, W_new, sizeof(W_new));
            }
            return 1;
        }
    }

    if (resid_out != 0) {
        /* TEMP TIMER: remove after rad-matter coupling profiling. */
        PRJ_TIMER_CURRENT_START("rad_matter_temp_newton_final_resid");
        ok = prj_rad_gr_m1_residual(rad, eos, geom, u_old, P, dt, resid_out,
            u_new_out, W_out);
        PRJ_TIMER_CURRENT_STOP("rad_matter_temp_newton_final_resid");
    } else {
        ok = 0;
    }
    if (resid_out != 0 && ok) {
        return 0;
    }
    return 0;
}

static int prj_rad_gr_m1_implicit_solve(const prj_rad *rad, prj_eos *eos,
    const prj_z4c_hydro_geom *geom, const double *u_old, double dt, double *P,
    double *resid_out, double *u_new_out, double *W_out)
{
    double P_approx[PRJ_RAD_GR_M1_NP];
    int ok;

    if (P == 0) {
        return 0;
    }
    memcpy(P_approx, P, sizeof(P_approx));
    PRJ_TIMER_CURRENT_START("rad_matter_temp_implicit_approx");
    ok = prj_rad_gr_m1_implicit_solve_approx(rad, eos, geom, u_old, dt,
        P_approx, resid_out, u_new_out, W_out);
    PRJ_TIMER_CURRENT_STOP("rad_matter_temp_implicit_approx");
    if (ok) {
        memcpy(P, P_approx, sizeof(P_approx));
        return 1;
    }

    PRJ_TIMER_CURRENT_START("rad_matter_temp_implicit_exact_fallback");
    ok = prj_rad_gr_m1_implicit_solve_exact(rad, eos, geom, u_old, dt, P,
        resid_out, u_new_out, W_out);
    PRJ_TIMER_CURRENT_STOP("rad_matter_temp_implicit_exact_fallback");
    return ok;
}

static int prj_rad_gr_m1_moments_to_p(const prj_z4c_hydro_geom *geom,
    const double g_cov[4][4], const double g_con[4][4], double E,
    const double Fcov[3], double *P_group)
{
    int a;

    /* The Newton primitives are the lab-frame Eulerian moments (E, F_i), which
     * are exactly the stored radiation primitives, so packing is a direct copy
     * (fluxes are already inside the physical |F| <= c E bound by
     * prj_rad_gr_m1_clamp_fluxes_for_solve). */
    (void)g_cov;
    (void)g_con;
    if (geom == 0 || Fcov == 0 || P_group == 0 || !isfinite(E) || E < 0.0) {
        return 0;
    }
    for (a = 0; a < 3; ++a) {
        if (!isfinite(Fcov[a])) {
            return 0;
        }
    }
    P_group[0] = E;
    P_group[1] = Fcov[0];
    P_group[2] = Fcov[1];
    P_group[3] = Fcov[2];
    return 1;
}

static int prj_rad_gr_m1_clamp_fluxes_for_solve(
    const prj_z4c_hydro_geom *geom, double *u)
{
    int field;
    int group;
    int a;
    int b;

    if (geom == 0 || u == 0) {
        return 0;
    }
    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            int eidx = PRJ_CONS_RAD_E(field, group);
            int fidx = PRJ_CONS_RAD_F1(field, group);
            double E = u[eidx];
            double F2 = 0.0;
            double cE;
            double fmax_bound;

            if (!isfinite(E)) {
                return 0;
            }
            for (a = 0; a < 3; ++a) {
                if (!isfinite(u[fidx + a])) {
                    return 0;
                }
                for (b = 0; b < 3; ++b) {
                    F2 += geom->gamma_inv[a][b] * u[fidx + a] *
                        u[fidx + b];
                }
            }
            if (!isfinite(F2) || F2 < 0.0) {
                return 0;
            }
            if (E < PRJ_RAD_GR_M1_SOLVE_EMIN) {
                E = PRJ_RAD_GR_M1_SOLVE_EMIN;
                u[eidx] = E;
            }
            /* The physical M1 bound is |F|/c <= E.  The evolved F_i stores
             * the explicit c, so the stored conserved bound is |F| <= c E.
             * Keep solve inputs just inside the free-streaming boundary; the
             * lab-frame variables remain nonsingular at vacuum, but the closure
             * derivative is better conditioned away from f = 1. */
            cE = PRJ_CLIGHT * E;
            fmax_bound = (1.0 - PRJ_RAD_GR_M1_F_MARGIN) * cE;
            if (F2 > fmax_bound * fmax_bound && F2 > 0.0) {
                double scale = fmax_bound / sqrt(F2);

                u[fidx] *= scale;
                u[fidx + 1] *= scale;
                u[fidx + 2] *= scale;
            }
        }
    }
    return 1;
}

static int prj_rad_gr_m1_copy_radiation_cons_to_prim(
    const prj_z4c_hydro_geom *geom, const double *u, double *prim)
{
    double inv_sqrtg;
    int field;
    int group;

    if (geom == 0 || u == 0 || prim == 0 ||
        !isfinite(geom->sqrt_gamma) || geom->sqrt_gamma <= 0.0) {
        return 0;
    }
    inv_sqrtg = 1.0 / geom->sqrt_gamma;
    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            prim[PRJ_PRIM_RAD_E(field, group)] =
                inv_sqrtg * u[PRJ_CONS_RAD_E(field, group)];
            prim[PRJ_PRIM_RAD_F1(field, group)] =
                inv_sqrtg * u[PRJ_CONS_RAD_F1(field, group)];
            prim[PRJ_PRIM_RAD_F2(field, group)] =
                inv_sqrtg * u[PRJ_CONS_RAD_F2(field, group)];
            prim[PRJ_PRIM_RAD_F3(field, group)] =
                inv_sqrtg * u[PRJ_CONS_RAD_F3(field, group)];
        }
    }
    return 1;
}

static void prj_rad_gr_m1_copy_m3_test_outputs(
    const prj_rad_grm1_m3_data *m3, double *J_out, double H_out[4],
    double L_out[4][4], double Q_out[4][4][4])
{
    int a;
    int b;
    int c;

    if (J_out != 0) {
        *J_out = m3 != 0 ? m3->J : 0.0;
    }
    if (H_out != 0) {
        for (a = 0; a < 4; ++a) {
            H_out[a] = m3 != 0 ? m3->H[a] : 0.0;
        }
    }
    if (L_out != 0) {
        for (a = 0; a < 4; ++a) {
            for (b = 0; b < 4; ++b) {
                L_out[a][b] = m3 != 0 ? m3->L[a][b] : 0.0;
            }
        }
    }
    if (Q_out != 0) {
        for (a = 0; a < 4; ++a) {
            for (b = 0; b < 4; ++b) {
                for (c = 0; c < 4; ++c) {
                    Q_out[a][b][c] = m3 != 0 ?
                        prj_rad_grm1_m3_component(m3, a, b, c) : 0.0;
                }
            }
        }
    }
}

int prj_rad_gr_m1_m3_from_build_R_test_wrapper(
    const prj_rad_gr_m1_closure_ctx *ctx, double E, const double Fcov[3],
    double *J_out, double H_out[4], double L_out[4][4],
    double Q_out[4][4][4])
{
    prj_rad_gr_m1_side_data side;
    prj_rad_grm1_m3_data m3;
    double g_cov[4][4];
    double g_con[4][4];
    double Rcon[4][4];
    double ucon[4];

    prj_rad_gr_m1_copy_m3_test_outputs(0, J_out, H_out, L_out, Q_out);
    if (ctx == 0 || Fcov == 0) {
        return 0;
    }
    prj_rad_gr_m1_prepare_side(ctx, &side);
    if (!prj_rad_grm1_build_R_from_ctx(ctx, E, Fcov, g_cov, g_con,
            Rcon, 0)) {
        return 0;
    }
    prj_rad_grm1_ucon_from_side(&side, ucon);
    if (!prj_rad_grm1_prepare_m3_data(g_cov, g_con, ucon, Rcon, &m3)) {
        return 0;
    }
    prj_rad_gr_m1_copy_m3_test_outputs(&m3, J_out, H_out, L_out, Q_out);
    return 1;
}

int prj_rad_gr_m1_freq_m3_test_wrapper(
    const prj_rad_gr_m1_closure_ctx *ctx,
    const prj_z4c_hydro_geom *geom, double E, const double Fcov[3],
    double *J_out, double H_out[4], double L_out[4][4],
    double Q_out[4][4][4])
{
    prj_rad_gr_m1_side_data side;
    prj_rad_gr_m1_freq_drift_cache cache;
    prj_rad_grm1_m3_data m3;
    double observer_time_derivative[4] = {0.0, 0.0, 0.0, 0.0};

    prj_rad_gr_m1_copy_m3_test_outputs(0, J_out, H_out, L_out, Q_out);
    if (ctx == 0 || geom == 0 || Fcov == 0) {
        return 0;
    }
    prj_rad_gr_m1_prepare_side(ctx, &side);
    prj_rad_gr_m1_freq_drift_cache_init(ctx, &side, geom,
        observer_time_derivative, &cache);
    prj_rad_grm1_prepare_m3_data_freq(&cache, E, Fcov, &m3);
    prj_rad_gr_m1_copy_m3_test_outputs(&m3, J_out, H_out, L_out, Q_out);
    return 1;
}

int prj_rad_gr_m1_freq_drift_test_wrapper(
    const prj_rad_gr_m1_closure_ctx *ctx,
    const prj_z4c_hydro_geom *geom, double E, const double Fcov[3],
    double fast[4], double ref[4])
{
    prj_rad_gr_m1_side_data side;
    prj_rad_gr_m1_freq_drift_cache cache;
    prj_rad_grm1_m3_data m3;
    double observer_time_derivative[4] = {0.0, 0.0, 0.0, 0.0};
    int a;

    if (fast != 0) {
        for (a = 0; a < 4; ++a) {
            fast[a] = 0.0;
        }
    }
    if (ref != 0) {
        for (a = 0; a < 4; ++a) {
            ref[a] = 0.0;
        }
    }
    if (ctx == 0 || geom == 0 || Fcov == 0 || fast == 0 || ref == 0) {
        return 0;
    }
    prj_rad_gr_m1_prepare_side(ctx, &side);
    prj_rad_gr_m1_freq_drift_cache_init(ctx, &side, geom,
        observer_time_derivative, &cache);
    prj_rad_gr_m1_frequency_drifts_cached(&cache, E, Fcov, &fast[0],
        &fast[1]);
    prj_rad_grm1_prepare_m3_data_freq(&cache, E, Fcov, &m3);
    prj_rad_gr_m1_frequency_drifts_component_ref(&cache, &m3, &ref[0],
        &ref[1]);
    return 1;
}

int prj_rad_gr_m1_build_R_jac_test_wrapper(
    const prj_z4c_hydro_geom *geom, double E, const double Fcov[3],
    double Rcon[4][4])
{
    double g_cov[4][4];
    double g_con[4][4];
    double dRcon[4][4][4];
    prj_rad_grm1_R_jac_geom R_jac_geom;

    if (Rcon != 0) {
        memset(Rcon, 0, 16 * sizeof(double));
    }
    if (geom == 0 || Fcov == 0 || Rcon == 0 ||
        !isfinite(geom->alpha) || geom->alpha <= 0.0) {
        return 0;
    }
    prj_rad_gr_m1_metric4_from_geom(geom, g_cov, g_con);
    prj_rad_grm1_R_jac_geom_init(geom->alpha, g_con, geom->gamma_inv,
        &R_jac_geom);
    return prj_rad_grm1_build_R_jac(g_cov, g_con, &R_jac_geom, E, Fcov,
        Rcon, dRcon);
}

int prj_rad_gr_m1_projected_jac_test_wrapper(
    const prj_z4c_hydro_geom *geom, double E, const double Fcov[3],
    const double cov[4], double R_u[4], double dR_u[4][4],
    double dRuu[4])
{
    double g_cov[4][4];
    double g_con[4][4];
    prj_rad_grm1_R_jac_geom R_jac_geom;
    prj_rad_grm1_R_projector Rproj;

    if (R_u != 0) {
        memset(R_u, 0, 4 * sizeof(double));
    }
    if (dR_u != 0) {
        memset(dR_u, 0, 16 * sizeof(double));
    }
    if (dRuu != 0) {
        memset(dRuu, 0, 4 * sizeof(double));
    }
    if (geom == 0 || Fcov == 0 || cov == 0 || R_u == 0 || dR_u == 0 ||
        dRuu == 0 || !isfinite(geom->alpha) || geom->alpha <= 0.0) {
        return 0;
    }
    prj_rad_gr_m1_metric4_from_geom(geom, g_cov, g_con);
    prj_rad_grm1_R_jac_geom_init(geom->alpha, g_con, geom->gamma_inv,
        &R_jac_geom);
    if (!prj_rad_grm1_R_projector_init(g_cov, g_con, &R_jac_geom, E,
            Fcov, &Rproj)) {
        return 0;
    }
    prj_rad_grm1_R_projector_apply(&Rproj, cov, R_u);
    prj_rad_grm1_R_projector_dR_u(&Rproj, g_con, &R_jac_geom, cov, dR_u,
        dRuu);
    return 1;
}

int prj_rad_gr_m1_full_jac_projection_test_wrapper(
    const prj_z4c_hydro_geom *geom, double E, const double Fcov[3],
    const double cov[4], double R_u[4], double dR_u[4][4],
    double dRuu[4])
{
    double g_cov[4][4];
    double g_con[4][4];
    double Rcon[4][4];
    double dRcon[4][4][4];
    prj_rad_grm1_R_jac_geom R_jac_geom;
    int a;
    int b;
    int cc;

    if (R_u != 0) {
        memset(R_u, 0, 4 * sizeof(double));
    }
    if (dR_u != 0) {
        memset(dR_u, 0, 16 * sizeof(double));
    }
    if (dRuu != 0) {
        memset(dRuu, 0, 4 * sizeof(double));
    }
    if (geom == 0 || Fcov == 0 || cov == 0 || R_u == 0 || dR_u == 0 ||
        dRuu == 0 || !isfinite(geom->alpha) || geom->alpha <= 0.0) {
        return 0;
    }
    prj_rad_gr_m1_metric4_from_geom(geom, g_cov, g_con);
    prj_rad_grm1_R_jac_geom_init(geom->alpha, g_con, geom->gamma_inv,
        &R_jac_geom);
    if (!prj_rad_grm1_build_R_jac(g_cov, g_con, &R_jac_geom, E, Fcov,
            Rcon, dRcon)) {
        return 0;
    }
    for (a = 0; a < 4; ++a) {
        for (b = 0; b < 4; ++b) {
            R_u[a] += Rcon[a][b] * cov[b];
        }
    }
    for (cc = 0; cc < 4; ++cc) {
        for (a = 0; a < 4; ++a) {
            for (b = 0; b < 4; ++b) {
                dR_u[cc][a] += dRcon[cc][a][b] * cov[b];
            }
            dRuu[cc] += dR_u[cc][a] * cov[a];
        }
    }
    return 1;
}

int prj_rad_gr_m1_residual_test_wrapper(const prj_rad *rad, prj_eos *eos,
    const prj_z4c_hydro_geom *geom, const double *u_old, const double *P,
    double dt, double *resid, double *u_new_out)
{
    return prj_rad_gr_m1_residual(rad, eos, geom, u_old, P, dt, resid,
        u_new_out, 0);
}

int prj_rad_gr_m1_jacobian_test_wrapper(const prj_rad *rad, prj_eos *eos,
    const prj_z4c_hydro_geom *geom, const double *u_old, const double *P,
    double dt, double *resid, double *jac, double *u_new_out)
{
    return prj_rad_gr_m1_jacobian(rad, eos, geom, u_old, P, dt, resid, jac,
        u_new_out);
}

int prj_rad_gr_m1_fast_candidate_compare_test_wrapper(const prj_rad *rad,
    prj_eos *eos, const prj_z4c_hydro_geom *geom, const double *u_old,
    const double *P, double dt, double *projected_out, double *dense_out,
    double *legacy_norm_out, double *projected_exact_norm_out,
    double *dense_exact_norm_out)
{
    const double threshold = PRJ_RAD_GR_M1_SOLVE_TOL_DEFAULT;
    prj_rad_gr_m1_frozen_group frozen[PRJ_RAD_GR_M1_NGROUPS];
    prj_rad_gr_m1_cell_context ctx;
    prj_eos_rty_interp_result eos_cache;
    prj_rad_gr_m1_opac_cache opac_cache;
    double projected[PRJ_RAD_GR_M1_NP];
    double dense[PRJ_RAD_GR_M1_NP];
    double legacy[PRJ_RAD_GR_M1_NP];
    double exact_resid[PRJ_NVAR_CONS];
    double exact_u[PRJ_NVAR_CONS];
    double exact_W[PRJ_NVAR_PRIM];
    double projected_norm;
    double dense_norm;
    double legacy_norm;

    if (rad == 0 || eos == 0 || geom == 0 || u_old == 0 || P == 0 ||
        projected_out == 0 || dense_out == 0 || legacy_norm_out == 0 ||
        projected_exact_norm_out == 0 || dense_exact_norm_out == 0 ||
        !prj_rad_gr_m1_cell_context_init(geom, P, &ctx)) {
        return 0;
    }
    if (!prj_rad_gr_m1_approx_reduced_resjac(rad, eos, geom, u_old, 0,
            &ctx, P, dt, projected, 0, 0, 0, 0, &eos_cache, 0,
            &opac_cache, 0, threshold, 1)) {
        return 0;
    }
    if (!prj_rad_gr_m1_freeze_closure(geom, P, frozen)) {
        return 0;
    }
    if (!prj_rad_gr_m1_approx_reduced_resjac(rad, eos, geom, u_old, frozen,
            &ctx, P, dt, dense, 0, 0, 0, &eos_cache, 0, &opac_cache, 0,
            0, threshold, 1)) {
        return 0;
    }
    if (!prj_rad_gr_m1_approx_reduced_resjac(rad, eos, geom, u_old, frozen,
            &ctx, P, dt, legacy, 0, 0, 0, &eos_cache, 0, &opac_cache, 0,
            &legacy_norm, threshold, 0)) {
        return 0;
    }
    if (!prj_rad_gr_m1_exact_residual_check(rad, eos, geom, u_old,
            projected, dt, &eos_cache, &opac_cache, &ctx, threshold,
            &projected_norm, exact_resid, exact_u, exact_W)) {
        return 0;
    }
    if (!prj_rad_gr_m1_exact_residual_check(rad, eos, geom, u_old, dense,
            dt, &eos_cache, &opac_cache, &ctx, threshold, &dense_norm,
            exact_resid, exact_u, exact_W)) {
        return 0;
    }
    memcpy(projected_out, projected, sizeof(projected));
    memcpy(dense_out, dense, sizeof(dense));
    *legacy_norm_out = legacy_norm;
    *projected_exact_norm_out = projected_norm;
    *dense_exact_norm_out = dense_norm;
    return 1;
}

int prj_rad_gr_m1_implicit_solve_test_wrapper(const prj_rad *rad,
    prj_eos *eos, const prj_z4c_hydro_geom *geom, const double *u_old,
    double dt, double *P, double *resid_out, double *u_new_out)
{
    return prj_rad_gr_m1_implicit_solve(rad, eos, geom, u_old, dt, P,
        resid_out, u_new_out, 0);
}

static void prj_rad_gr_m1_matter_abort(const char *reason,
    const prj_rad *rad, const prj_block *block, int z4c_stage,
    int i, int j, int k, double dt, const double *u_old,
    const double *P, const double *W, const prj_z4c_hydro_geom *geom,
    const double *resid)
{
    double threshold = PRJ_RAD_GR_M1_SOLVE_TOL_DEFAULT;
    double norm = HUGE_VAL;

    if (rad != 0 && isfinite(rad->implicit_err_tol) &&
        rad->implicit_err_tol > 0.0) {
        threshold = rad->implicit_err_tol;
    }
    if (u_old != 0 && resid != 0) {
        norm = prj_rad_gr_m1_residual_norm(u_old, resid, threshold);
    }

    fprintf(stderr,
        "prj_rad_gr_m1_matter_update: %s at block id=%d level=%d "
        "z4c_stage=%d cell=(%d,%d,%d) dt=%.17e residual_norm=%.17e "
        "threshold=%.17e\n",
        reason,
        block != 0 ? block->id : -1,
        block != 0 ? block->level : -1,
        z4c_stage, i, j, k, dt, norm, threshold);
    if (P != 0) {
        fprintf(stderr,
            "  Newton state: rho=%.17e v=(%.17e, %.17e, %.17e) "
            "T=%.17e Ye=%.17e\n",
            P[0], P[1], P[2], P[3], P[4], P[5]);
    }
    if (resid != 0) {
        fprintf(stderr,
            "  residual fluid rows: rho=%.17e mom=(%.17e, %.17e, %.17e) "
            "etot=%.17e Ye=%.17e\n",
            resid[PRJ_CONS_RHO], resid[PRJ_CONS_MOM1],
            resid[PRJ_CONS_MOM2], resid[PRJ_CONS_MOM3],
            resid[PRJ_CONS_ETOT], resid[PRJ_CONS_YE]);
    }
    if (W != 0 && geom != 0) {
        double beta2 = 0.0;
        double min_E = HUGE_VAL;
        double max_flux_factor = 0.0;
        double min_P_E = HUGE_VAL;
        int max_flux_field = -1;
        int max_flux_group = -1;
        int min_E_field = -1;
        int min_E_group = -1;
        int field;
        int group;
        int a;
        int b;

        for (a = 0; a < 3; ++a) {
            for (b = 0; b < 3; ++b) {
                beta2 += geom->gamma[a][b] *
                    (W[PRJ_PRIM_V1 + a] / PRJ_CLIGHT) *
                    (W[PRJ_PRIM_V1 + b] / PRJ_CLIGHT);
            }
        }
        for (field = 0; field < PRJ_NRAD; ++field) {
            for (group = 0; group < PRJ_NEGROUP; ++group) {
                int idx = field * PRJ_NEGROUP + group;
                int pidx = 6 + 4 * idx;
                double E = W[PRJ_PRIM_RAD_E(field, group)];
                double Fcov[3];
                double F2 = 0.0;
                double flux_factor = HUGE_VAL;

                Fcov[0] = W[PRJ_PRIM_RAD_F1(field, group)];
                Fcov[1] = W[PRJ_PRIM_RAD_F2(field, group)];
                Fcov[2] = W[PRJ_PRIM_RAD_F3(field, group)];
                for (a = 0; a < 3; ++a) {
                    for (b = 0; b < 3; ++b) {
                        F2 += geom->gamma_inv[a][b] * Fcov[a] * Fcov[b];
                    }
                }
                if (isfinite(E) && E > 0.0 && isfinite(F2) && F2 >= 0.0) {
                    flux_factor = sqrt(F2) / (PRJ_CLIGHT * E);
                }
                if (E < min_E) {
                    min_E = E;
                    min_E_field = field;
                    min_E_group = group;
                }
                if (flux_factor > max_flux_factor) {
                    max_flux_factor = flux_factor;
                    max_flux_field = field;
                    max_flux_group = group;
                }
                if (P != 0 && P[pidx] < min_P_E) {
                    min_P_E = P[pidx];
                }
            }
        }
        fprintf(stderr,
            "  input primitive: rho=%.17e eint=%.17e Ye=%.17e "
            "v2/c2=%.17e min_E=%.17e(F%d,g%d) max_|F|/cE=%.17e(F%d,g%d) "
            "min_P_E=%.17e\n",
            W[PRJ_PRIM_RHO], W[PRJ_PRIM_EINT], W[PRJ_PRIM_YE],
            beta2, min_E, min_E_field, min_E_group, max_flux_factor,
            max_flux_field, max_flux_group, min_P_E);
    }
    fprintf(stderr,
        "  Newton diagnostics: iter=%d/%d block_solve_ok=%d dense_solve_ok=%d "
        "dP_max_abs=%.17e(col=%d) dP_max_rel=%.17e(col=%d) "
        "min_rad_block_pivot=%.17e rel=%.17e group=%d k=%d "
        "trials=%d invalid_trials=%d first_trial_norm=%.17e "
        "best_trial_norm=%.17e best_ls=%d best_lambda=%.17e "
        "last_trial_norm=%.17e\n",
        prj_rad_gr_m1_last_solve_diag.iter,
        prj_rad_gr_m1_last_solve_diag.maxiter,
        prj_rad_gr_m1_last_solve_diag.block_solve_ok,
        prj_rad_gr_m1_last_solve_diag.dense_solve_ok,
        prj_rad_gr_m1_last_solve_diag.dP_max_abs,
        prj_rad_gr_m1_last_solve_diag.dP_max_abs_col,
        prj_rad_gr_m1_last_solve_diag.dP_max_rel,
        prj_rad_gr_m1_last_solve_diag.dP_max_rel_col,
        prj_rad_gr_m1_last_solve_diag.linear_min_pivot,
        prj_rad_gr_m1_last_solve_diag.linear_min_pivot_rel,
        prj_rad_gr_m1_last_solve_diag.linear_min_pivot_group,
        prj_rad_gr_m1_last_solve_diag.linear_min_pivot_k,
        prj_rad_gr_m1_last_solve_diag.line_search_trials,
        prj_rad_gr_m1_last_solve_diag.invalid_trials,
        prj_rad_gr_m1_last_solve_diag.first_trial_norm,
        prj_rad_gr_m1_last_solve_diag.best_trial_norm,
        prj_rad_gr_m1_last_solve_diag.best_ls,
        prj_rad_gr_m1_last_solve_diag.best_lambda,
        prj_rad_gr_m1_last_solve_diag.last_trial_norm);
    fprintf(stderr,
        "  J*dP check: ok=%d eps=%.17e max_abs=%.17e max_rel=%.17e "
        "row=%d |JdP|=%.17e |fd|=%.17e\n",
        prj_rad_gr_m1_last_solve_diag.directional_check_ok,
        prj_rad_gr_m1_last_solve_diag.directional_eps,
        prj_rad_gr_m1_last_solve_diag.directional_max_abs,
        prj_rad_gr_m1_last_solve_diag.directional_max_rel,
        prj_rad_gr_m1_last_solve_diag.directional_row,
        prj_rad_gr_m1_last_solve_diag.directional_jdp_norm,
        prj_rad_gr_m1_last_solve_diag.directional_fd_norm);
    fprintf(stderr,
        "  FD Jacobian check: ok=%d skipped_cols=%d "
        "max_abs=%.17e row=%d col=%d got=%.17e fd=%.17e h=%.17e "
        "max_rel=%.17e row=%d col=%d got=%.17e fd=%.17e h=%.17e\n",
        prj_rad_gr_m1_last_solve_diag.fd_check_ok,
        prj_rad_gr_m1_last_solve_diag.fd_skipped_cols,
        prj_rad_gr_m1_last_solve_diag.fd_max_abs,
        prj_rad_gr_m1_last_solve_diag.fd_max_abs_row,
        prj_rad_gr_m1_last_solve_diag.fd_max_abs_col,
        prj_rad_gr_m1_last_solve_diag.fd_max_abs_got,
        prj_rad_gr_m1_last_solve_diag.fd_max_abs_expected,
        prj_rad_gr_m1_last_solve_diag.fd_max_abs_h,
        prj_rad_gr_m1_last_solve_diag.fd_max_rel,
        prj_rad_gr_m1_last_solve_diag.fd_max_rel_row,
        prj_rad_gr_m1_last_solve_diag.fd_max_rel_col,
        prj_rad_gr_m1_last_solve_diag.fd_max_rel_got,
        prj_rad_gr_m1_last_solve_diag.fd_max_rel_expected,
        prj_rad_gr_m1_last_solve_diag.fd_max_rel_h);
    fflush(stderr);
#if defined(PRJ_ENABLE_MPI)
    {
        int mpi_ready = 0;

        MPI_Initialized(&mpi_ready);
        if (mpi_ready) {
            int finalized = 0;

            MPI_Finalized(&finalized);
            if (!finalized) {
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            }
        }
    }
#endif
    exit(EXIT_FAILURE);
}
#endif

#if PRJ_NRAD > 0
/* Freeze the electron-inelastic-scattering (eleinel) opacity/emissivity for one
 * cell at the start of the GR M1 matter-coupling implicit step.  Using the
 * comoving radiation energy J = R^{ab} u_a u_b (built from the SAME canonical
 * M1 R^{ab} the residual/jacobian use) and the INITIAL temperature, this calls
 * the exact non-GR eleinel kernel and stores per-group additions to eta/kappa
 * that are then held constant during the Newton iteration.  The source/sink map
 * onto (eta, kappa) because the non-GR normal energy update and the non-GR
 * eleinel step share the identical backward-Euler form
 *   E_new = (E_old + dt*eta) / (1 + dt*c*kappa)   [normal]
 *   E_new = (E_old + dt*source_phys) / (1 + dt*c*sink)   [eleinel].
 * When eleinel is absent (table not loaded, or rho below threshold) the arrays
 * stay zero, so the assembly-site additions are a no-op (bit-identical). */
static void prj_rad_gr_m1_freeze_eleinel(prj_rad *rad, prj_eos *eos,
    const prj_z4c_hydro_geom *geom, const double g_cov[4][4],
    const double g_con[4][4], const double *u_old, const double *P)
{
    const int ngroups = PRJ_NRAD * PRJ_NEGROUP;
    const double eta_factor = 1.0 / (4.0 * M_PI);
    double je[PRJ_NRAD * PRJ_NEGROUP];
    double he[PRJ_NRAD * PRJ_NEGROUP * PRJ_NDIM];
    double source_arr[PRJ_NRAD * PRJ_NEGROUP];
    double sink_arr[PRJ_NRAD * PRJ_NEGROUP];
    double ucon[4];
    double ucov[4];
    double rho;
    double T;
    double Ye;
    double etael;
    int field;
    int group;
    int d;
    int idx;

    (void)u_old;
    for (idx = 0; idx < ngroups; ++idx) {
        rad->gr_m1_eleinel_dkappa[idx] = 0.0;
        rad->gr_m1_eleinel_deta[idx] = 0.0;
        je[idx] = 0.0;
        for (d = 0; d < PRJ_NDIM; ++d) {
            he[idx * PRJ_NDIM + d] = 0.0;
        }
    }
    if (!rad->eleinel_table_loaded) {
        return;
    }
    rho = P[0];
    T = P[4];
    Ye = P[5];
    if (rho <= rad->min_inel_density) {
        return;
    }
    if (!prj_rad_gr_m1_fluid_four_velocity(geom, g_cov, P, ucon, ucov)) {
        return;
    }
    etael = prj_eos_rty_geteta(eos, rho, T, Ye, PRJ_EOS_CTX_MAIN);
    if (etael < -20.0) {
        etael = -20.0;
    }

    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            int gidx = field * PRJ_NEGROUP + group;
            int pidx = 6 + 4 * gidx;
            double E = P[pidx];
            double Fcov[3];
            double Fcon[3];
            double Rcon[4][4];
            double R_u[4];
            double Ruu = 0.0;
            double Hcon[4];
            double he_mag = 0.0;
            double conv = RAD_SCALE * eta_factor * PRJ_MEV_TO_ERG /
                rad->degroup_erg[field][group];
            int a;
            int b;

            Fcov[0] = P[pidx + 1];
            Fcov[1] = P[pidx + 2];
            Fcov[2] = P[pidx + 3];
            for (a = 0; a < 3; ++a) {
                Fcon[a] = 0.0;
                for (b = 0; b < 3; ++b) {
                    Fcon[a] += geom->gamma_inv[a][b] * Fcov[b];
                }
            }
            if (E < 0.0 || !prj_rad_grm1_build_R(g_cov, g_con, geom->alpha, E,
                    Fcon, Rcon)) {
                continue; /* leave this group's je/he at zero */
            }
            for (a = 0; a < 4; ++a) {
                R_u[a] = 0.0;
                for (b = 0; b < 4; ++b) {
                    R_u[a] += Rcon[a][b] * ucov[b];
                }
                Ruu += R_u[a] * ucov[a];
            }
            if (Ruu < 0.0) {
                Ruu = 0.0;
            }
            /* Comoving energy J = Ruu; comoving flux H^a = -(R^{ab}u_b + u^a J)
             * (projector uses u.u = -1 in this geometric metric).  Flat/static
             * rest limit: J -> E, H_i -> F_i, matching the non-GR je/he. */
            je[gidx] = PRJ_CLIGHT * Ruu * conv;
            for (a = 0; a < 4; ++a) {
                Hcon[a] = -(R_u[a] + ucon[a] * Ruu);
            }
            for (d = 0; d < PRJ_NDIM; ++d) {
                double Hcov_d = 0.0;
                double ratio;

                for (b = 0; b < 3; ++b) {
                    Hcov_d += geom->gamma[d][b] * Hcon[b + 1];
                }
                he[gidx * PRJ_NDIM + d] = Hcov_d * conv;
                ratio = he[gidx * PRJ_NDIM + d] / (je[gidx] + 1.0e-15);
                he_mag += ratio * ratio;
            }
            he_mag = sqrt(he_mag);
            if (he_mag > 1.0) {
                for (d = 0; d < PRJ_NDIM; ++d) {
                    he[gidx * PRJ_NDIM + d] /= he_mag;
                }
            }
        }
    }

    prj_rad_eleinel_lookup(rad, rho, T, Ye, etael, je, he,
        source_arr, sink_arr, 0);

    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            int gidx = field * PRJ_NEGROUP + group;
            double source_phys = source_arr[gidx] *
                rad->degroup_erg[field][group] /
                (PRJ_MEV_TO_ERG * eta_factor);

            rad->gr_m1_eleinel_deta[gidx] = source_phys;
            rad->gr_m1_eleinel_dkappa[gidx] = sink_arr[gidx];
        }
    }
}
#endif

void prj_rad_gr_m1_matter_update_geom(prj_rad *rad, prj_eos *eos,
    const prj_mesh *mesh, const prj_block *block, int z4c_stage,
    const prj_z4c_hydro_geom *geom, double *u, double *prim,
    int i, int j, int k, double dt, double *final_temperature)
{
#if PRJ_NRAD > 0
    prj_eos_gr_geom eos_geom;
    double g_cov[4][4];
    double g_con[4][4];
    double u_old[PRJ_NVAR_CONS];
    double u_new[PRJ_NVAR_CONS];
    double W[PRJ_NVAR_PRIM];
    double W_new[PRJ_NVAR_PRIM];
    double P[PRJ_RAD_GR_M1_NP];
    double resid[PRJ_NVAR_CONS];
    double eos_q[PRJ_EOS_NQUANT];
    int field;
    int group;
    int a;
    int b;
    int v;
    int ok;

    ok = rad != 0 && eos != 0 && mesh != 0 && block != 0 && u != 0 &&
        geom != 0 && isfinite(dt);
    if (!ok) {
        return;
    }
    /* TEMP TIMER: remove after rad-matter coupling profiling. */
    PRJ_TIMER_CURRENT_START("rad_matter_temp_clamp");
    ok = prj_rad_gr_m1_clamp_fluxes_for_solve(geom, u);
    PRJ_TIMER_CURRENT_STOP("rad_matter_temp_clamp");
    if (!ok) {
        return;
    }
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            eos_geom.gamma[a][b] = geom->gamma[a][b];
        }
    }
    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        u_old[v] = u[v];
    }
    /* TEMP TIMER: remove after rad-matter coupling profiling. */
    PRJ_TIMER_CURRENT_START("rad_matter_temp_cons2prim");
    ok = prj_eos_gr_cons2prim(eos, &eos_geom, u_old, W,
        PRJ_EOS_CTX_MAIN) == PRJ_EOS_GR_OK;
    PRJ_TIMER_CURRENT_STOP("rad_matter_temp_cons2prim");
    if (!ok) {
        return;
    }
    if (!isfinite(W[PRJ_PRIM_RHO]) || W[PRJ_PRIM_RHO] <= 0.0 ||
        !isfinite(W[PRJ_PRIM_EINT]) || W[PRJ_PRIM_EINT] < 0.0 ||
        !isfinite(W[PRJ_PRIM_YE])) {
        return;
    }
    /* TEMP TIMER: remove after rad-matter coupling profiling. */
    PRJ_TIMER_CURRENT_START("rad_matter_temp_eos_rey");
    prj_eos_rey(eos, W[PRJ_PRIM_RHO], W[PRJ_PRIM_EINT],
        W[PRJ_PRIM_YE], eos_q, PRJ_EOS_CTX_MAIN);
    PRJ_TIMER_CURRENT_STOP("rad_matter_temp_eos_rey");
    if (!isfinite(eos_q[PRJ_EOS_TEMPERATURE]) ||
        eos_q[PRJ_EOS_TEMPERATURE] <= 0.0) {
        return;
    }
    if (final_temperature != 0) {
        *final_temperature = eos_q[PRJ_EOS_TEMPERATURE];
    }

    P[0] = W[PRJ_PRIM_RHO];
    P[1] = W[PRJ_PRIM_V1];
    P[2] = W[PRJ_PRIM_V2];
    P[3] = W[PRJ_PRIM_V3];
    P[4] = eos_q[PRJ_EOS_TEMPERATURE];
    P[5] = W[PRJ_PRIM_YE];

    /* TEMP TIMER: remove after rad-matter coupling profiling. */
    PRJ_TIMER_CURRENT_START("rad_matter_temp_moment_pack");
    prj_rad_gr_m1_metric4_from_geom(geom, g_cov, g_con);
    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            int idx = field * PRJ_NEGROUP + group;
            int pidx = 6 + 4 * idx;
            double Fcov[3];

            Fcov[0] = W[PRJ_PRIM_RAD_F1(field, group)];
            Fcov[1] = W[PRJ_PRIM_RAD_F2(field, group)];
            Fcov[2] = W[PRJ_PRIM_RAD_F3(field, group)];
            if (!prj_rad_gr_m1_moments_to_p(geom, g_cov, g_con,
                    W[PRJ_PRIM_RAD_E(field, group)], Fcov, &P[pidx])) {
                PRJ_TIMER_CURRENT_STOP("rad_matter_temp_moment_pack");
                return;
            }
        }
    }
    PRJ_TIMER_CURRENT_STOP("rad_matter_temp_moment_pack");

    /* Freeze the eleinel opacity/emissivity from the start-of-step comoving J
     * and the initial temperature; added, held constant, to eta/kappa inside
     * the implicit solve below. */
    PRJ_TIMER_CURRENT_START("rad_matter_temp_eleinel_freeze");
    prj_rad_gr_m1_freeze_eleinel(rad, eos, geom, g_cov, g_con, u_old, P);
    PRJ_TIMER_CURRENT_STOP("rad_matter_temp_eleinel_freeze");

    /* TEMP TIMER: remove after rad-matter coupling profiling. */
    PRJ_TIMER_CURRENT_START("rad_matter_temp_implicit");
    ok = prj_rad_gr_m1_implicit_solve(rad, eos, geom, u_old, dt, P, resid,
        u_new, W_new);
    PRJ_TIMER_CURRENT_STOP("rad_matter_temp_implicit");
    if (!ok) {
        prj_rad_gr_m1_matter_abort("implicit solve failed", rad, block,
            z4c_stage, i, j, k, dt, u_old, P, W, geom, resid);
    }
    ok = prj_rad_gr_m1_clamp_fluxes_for_solve(geom, u_new);
    if (!ok) {
        prj_rad_gr_m1_matter_abort("post-solve radiation projection failed",
            rad, block, z4c_stage, i, j, k, dt, u_old, P, W, geom, resid);
    }
    ok = prj_rad_gr_m1_copy_radiation_cons_to_prim(geom, u_new, W_new);
    if (!ok) {
        prj_rad_gr_m1_matter_abort("post-solve primitive projection failed",
            rad, block, z4c_stage, i, j, k, dt, u_old, P, W, geom, resid);
    }
    /* TEMP TIMER: remove after rad-matter coupling profiling. */
    PRJ_TIMER_CURRENT_START("rad_matter_temp_copy_out");
    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        u[v] = u_new[v];
    }
    if (prim != 0) {
        for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
            prim[v] = W_new[v];
        }
    }
    if (final_temperature != 0) {
        *final_temperature = P[4];
    }
    PRJ_TIMER_CURRENT_STOP("rad_matter_temp_copy_out");
#else
    (void)rad;
    (void)eos;
    (void)mesh;
    (void)block;
    (void)z4c_stage;
    (void)geom;
    (void)u;
    (void)prim;
    (void)i;
    (void)j;
    (void)k;
    (void)dt;
    (void)final_temperature;
#endif
}

void prj_rad_gr_m1_matter_update(prj_rad *rad, prj_eos *eos,
    const prj_mesh *mesh, const prj_block *block, int z4c_stage, double *u,
    double *prim, int i, int j, int k, double dt,
    double *final_temperature)
{
    prj_z4c_hydro_geom geom;
    int ok;

    /* TEMP TIMER: remove after rad-matter coupling profiling. */
    PRJ_TIMER_CURRENT_START("rad_matter_temp_geom");
    ok = rad != 0 && eos != 0 && mesh != 0 && block != 0 && u != 0 &&
        isfinite(dt) &&
        prj_z4c_load_hydro_geom(mesh, block, z4c_stage, i, j, k, &geom);
    PRJ_TIMER_CURRENT_STOP("rad_matter_temp_geom");
    if (!ok) {
        return;
    }
    prj_rad_gr_m1_matter_update_geom(rad, eos, mesh, block, z4c_stage, &geom,
        u, prim, i, j, k, dt, final_temperature);
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
