#ifndef PRJ_RECONSTRUCT_H
#define PRJ_RECONSTRUCT_H

static inline double prj_reconstruct_abs(double x)
{
    return x < 0.0 ? -x : x;
}

static inline double prj_reconstruct_minmod_slope(const double stencil[3])
{
    double sl;
    double sr;

    if (stencil == 0) {
        return 0.0;
    }

    sl = stencil[1] - stencil[0];
    sr = stencil[2] - stencil[1];
    if (sl * sr <= 0.0) {
        return 0.0;
    }
    return prj_reconstruct_abs(sl) < prj_reconstruct_abs(sr) ? sl : sr;
}

static inline double prj_reconstruct_mc_slope(const double stencil[3])
{
    double sl;
    double sr;
    double v;
    double a;
    double b;
    double phi;

    if (stencil == 0) {
        return 0.0;
    }

    sl = stencil[1] - stencil[0];
    sr = stencil[2] - stencil[1];

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

/* target_location values are offsets from stencil[1] in stencil-spacing units. */
static inline void prj_reconstruct_for_prolongate(const double stencil[3],
    int ntarget, const double target_location[], double target_value[])
{
    int n;
    double slope;

    if (stencil == 0 || target_location == 0 || target_value == 0 || ntarget <= 0) {
        return;
    }

    slope = prj_reconstruct_minmod_slope(stencil);
    for (n = 0; n < ntarget; ++n) {
        target_value[n] = stencil[1] + slope * target_location[n];
    }
}

static inline void prj_reconstruct_for_riemann(const double stencil[3],
    int ntarget, const double target_location[], double target_value[])
{
    int n;
    double slope;

    if (stencil == 0 || target_location == 0 || target_value == 0 || ntarget <= 0) {
        return;
    }

    slope = prj_reconstruct_mc_slope(stencil);
    for (n = 0; n < ntarget; ++n) {
        target_value[n] = stencil[1] + slope * target_location[n];
    }
}

int prj_reconstruct_step(void);

#endif
