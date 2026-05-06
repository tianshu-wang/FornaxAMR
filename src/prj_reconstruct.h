#ifndef PRJ_RECONSTRUCT_H
#define PRJ_RECONSTRUCT_H

static inline double prj_reconstruct_abs(double x)
{
    return x < 0.0 ? -x : x;
}

static inline double prj_reconstruct_min(double a, double b)
{
    return a < b ? a : b;
}

static inline double prj_reconstruct_max(double a, double b)
{
    return a > b ? a : b;
}

static inline double prj_reconstruct_clamp01(double x)
{
    if (x < 0.0) {
        return 0.0;
    }
    if (x > 1.0) {
        return 1.0;
    }
    return x;
}

static inline int prj_reconstruct_stencil3_index(int di, int dj, int dk)
{
    return (di + 1) * 9 + (dj + 1) * 3 + (dk + 1);
}

static inline int prj_reconstruct_stencil2_index(int di, int dj)
{
    return (di + 1) * 3 + (dj + 1);
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

/* target values are offsets from the coarse-cell center in coarse-cell units. */
static inline double prj_reconstruct_cell_for_prolongate_minmod(const double stencil[27],
    const double target_location[3])
{
    double base;
    double stx[3];
    double sty[3];
    double stz[3];
    double tx[1];
    double ty[1];
    double tz[1];
    double vx[1];
    double vy[1];
    double vz[1];

    base = stencil[prj_reconstruct_stencil3_index(0, 0, 0)];
    stx[0] = stencil[prj_reconstruct_stencil3_index(-1, 0, 0)];
    stx[1] = base;
    stx[2] = stencil[prj_reconstruct_stencil3_index(1, 0, 0)];
    sty[0] = stencil[prj_reconstruct_stencil3_index(0, -1, 0)];
    sty[1] = base;
    sty[2] = stencil[prj_reconstruct_stencil3_index(0, 1, 0)];
    stz[0] = stencil[prj_reconstruct_stencil3_index(0, 0, -1)];
    stz[1] = base;
    stz[2] = stencil[prj_reconstruct_stencil3_index(0, 0, 1)];

    tx[0] = target_location[0];
    ty[0] = target_location[1];
    tz[0] = target_location[2];
    prj_reconstruct_for_prolongate(stx, 1, tx, vx);
    prj_reconstruct_for_prolongate(sty, 1, ty, vy);
    prj_reconstruct_for_prolongate(stz, 1, tz, vz);
    return vx[0] + vy[0] + vz[0] - 2.0 * base;
}

static inline double prj_reconstruct_cell_for_prolongate_bj(const double stencil[27],
    const double target_location[3])
{
    static const double fine_offset[2] = {-0.25, 0.25};
    double base;
    double grad[3];
    double qmin;
    double qmax;
    double alpha = 1.0;
    int di;
    int dj;
    int dk;
    int ii;
    int jj;
    int kk;

    base = stencil[prj_reconstruct_stencil3_index(0, 0, 0)];
    grad[0] = 0.5 * (stencil[prj_reconstruct_stencil3_index(1, 0, 0)] -
        stencil[prj_reconstruct_stencil3_index(-1, 0, 0)]);
    grad[1] = 0.5 * (stencil[prj_reconstruct_stencil3_index(0, 1, 0)] -
        stencil[prj_reconstruct_stencil3_index(0, -1, 0)]);
    grad[2] = 0.5 * (stencil[prj_reconstruct_stencil3_index(0, 0, 1)] -
        stencil[prj_reconstruct_stencil3_index(0, 0, -1)]);

    qmin = stencil[prj_reconstruct_stencil3_index(-1, -1, -1)];
    qmax = qmin;
    for (di = -1; di <= 1; ++di) {
        for (dj = -1; dj <= 1; ++dj) {
            for (dk = -1; dk <= 1; ++dk) {
                double q;

                if (di == 0 && dj == 0 && dk == 0) {
                    continue;
                }
                q = stencil[prj_reconstruct_stencil3_index(di, dj, dk)];
                qmin = prj_reconstruct_min(qmin, q);
                qmax = prj_reconstruct_max(qmax, q);
            }
        }
    }

    for (ii = 0; ii < 2; ++ii) {
        for (jj = 0; jj < 2; ++jj) {
            for (kk = 0; kk < 2; ++kk) {
                double delta = grad[0] * fine_offset[ii] +
                    grad[1] * fine_offset[jj] + grad[2] * fine_offset[kk];
                double q = base + delta;
                double a;

                if (q > qmax) {
                    a = (delta != 0.0) ? (qmax - base) / delta : 0.0;
                    alpha = prj_reconstruct_min(alpha, prj_reconstruct_clamp01(a));
                } else if (q < qmin) {
                    a = (delta != 0.0) ? (qmin - base) / delta : 0.0;
                    alpha = prj_reconstruct_min(alpha, prj_reconstruct_clamp01(a));
                }
            }
        }
    }

    return base + alpha * (grad[0] * target_location[0] +
        grad[1] * target_location[1] + grad[2] * target_location[2]);
}

static inline double prj_reconstruct_cell_for_prolongate(const double stencil[27],
    const double target_location[3], int use_BJ)
{
    if (stencil == 0 || target_location == 0) {
        return 0.0;
    }
    if (use_BJ != 0) {
        return prj_reconstruct_cell_for_prolongate_bj(stencil, target_location);
    }
    return prj_reconstruct_cell_for_prolongate_minmod(stencil, target_location);
}

/* Face-centered Bf prolongation uses offsets in the two tangential directions. */
static inline double prj_reconstruct_face_for_prolongate_minmod(const double stencil[9],
    const double target_location[2])
{
    double base;
    double st0[3];
    double st1[3];
    double tx[1];
    double ty[1];
    double vx[1];
    double vy[1];

    base = stencil[prj_reconstruct_stencil2_index(0, 0)];
    st0[0] = stencil[prj_reconstruct_stencil2_index(-1, 0)];
    st0[1] = base;
    st0[2] = stencil[prj_reconstruct_stencil2_index(1, 0)];
    st1[0] = stencil[prj_reconstruct_stencil2_index(0, -1)];
    st1[1] = base;
    st1[2] = stencil[prj_reconstruct_stencil2_index(0, 1)];

    tx[0] = target_location[0];
    ty[0] = target_location[1];
    prj_reconstruct_for_prolongate(st0, 1, tx, vx);
    prj_reconstruct_for_prolongate(st1, 1, ty, vy);
    return vx[0] + vy[0] - base;
}

static inline double prj_reconstruct_face_for_prolongate_bj(const double stencil[9],
    const double target_location[2])
{
    static const double fine_offset[2] = {-0.25, 0.25};
    double base;
    double grad[2];
    double qmin;
    double qmax;
    double alpha = 1.0;
    int di;
    int dj;
    int ii;
    int jj;

    base = stencil[prj_reconstruct_stencil2_index(0, 0)];
    grad[0] = 0.5 * (stencil[prj_reconstruct_stencil2_index(1, 0)] -
        stencil[prj_reconstruct_stencil2_index(-1, 0)]);
    grad[1] = 0.5 * (stencil[prj_reconstruct_stencil2_index(0, 1)] -
        stencil[prj_reconstruct_stencil2_index(0, -1)]);

    qmin = stencil[prj_reconstruct_stencil2_index(-1, -1)];
    qmax = qmin;
    for (di = -1; di <= 1; ++di) {
        for (dj = -1; dj <= 1; ++dj) {
            double q;

            if (di == 0 && dj == 0) {
                continue;
            }
            q = stencil[prj_reconstruct_stencil2_index(di, dj)];
            qmin = prj_reconstruct_min(qmin, q);
            qmax = prj_reconstruct_max(qmax, q);
        }
    }

    for (ii = 0; ii < 2; ++ii) {
        for (jj = 0; jj < 2; ++jj) {
            double delta = grad[0] * fine_offset[ii] + grad[1] * fine_offset[jj];
            double q = base + delta;
            double a;

            if (q > qmax) {
                a = (delta != 0.0) ? (qmax - base) / delta : 0.0;
                alpha = prj_reconstruct_min(alpha, prj_reconstruct_clamp01(a));
            } else if (q < qmin) {
                a = (delta != 0.0) ? (qmin - base) / delta : 0.0;
                alpha = prj_reconstruct_min(alpha, prj_reconstruct_clamp01(a));
            }
        }
    }

    return base + alpha * (grad[0] * target_location[0] +
        grad[1] * target_location[1]);
}

static inline double prj_reconstruct_face_for_prolongate(const double stencil[9],
    const double target_location[2], int use_BJ)
{
    if (stencil == 0 || target_location == 0) {
        return 0.0;
    }
    if (use_BJ != 0) {
        return prj_reconstruct_face_for_prolongate_bj(stencil, target_location);
    }
    return prj_reconstruct_face_for_prolongate_minmod(stencil, target_location);
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
