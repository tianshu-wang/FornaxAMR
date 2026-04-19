#ifndef PRJ_RECONSTRUCT_H
#define PRJ_RECONSTRUCT_H

static inline double prj_reconstruct_slope(double stencil[3], double dx)
{
    double sl;
    double sr;
    double v;
    double a;
    double b;
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
    a = 0.5 * (1.0 + v);
    b = 2.0 * v;
    phi = (a < b) ? a : b;
    if (phi > 2.0) phi = 2.0;
    if (phi < 0.0) phi = 0.0;
    return sr * phi;
}

int prj_reconstruct_step(void);

#endif
