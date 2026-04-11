#include "prj.h"

static double prj_min_double(double a, double b)
{
    return a < b ? a : b;
}

static double prj_max_double(double a, double b)
{
    return a > b ? a : b;
}

double prj_reconstruct_slope(double stencil[3], double dx)
{
    double sl;
    double sr;
    double v;
    double phi;

    if (stencil == 0 || dx <= 0.0) {
        return 0.0;
    }

    sl = (stencil[1] - stencil[0]) / dx;
    sr = (stencil[2] - stencil[1]) / dx;

    if (sl == 0.0 || sr == 0.0 || sl * sr <= 0.0) {
        return 0.0;
    }

    v = sl / sr;
    phi = prj_max_double(0.0, prj_min_double(0.5 * (1.0 + v), prj_min_double(2.0, 2.0 * v)));
    return sr * phi;
}

int prj_reconstruct_step(void)
{
    return 0;
}
