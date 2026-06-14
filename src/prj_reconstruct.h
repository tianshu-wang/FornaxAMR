#ifndef PRJ_RECONSTRUCT_H
#define PRJ_RECONSTRUCT_H

#include "prj_defs.h"

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

static inline double prj_reconstruct_fine_offset(int idx)
{
    return idx == 0 ? -0.25 : 0.25;
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
                double delta = grad[0] * prj_reconstruct_fine_offset(ii) +
                    grad[1] * prj_reconstruct_fine_offset(jj) +
                    grad[2] * prj_reconstruct_fine_offset(kk);
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
            double delta = grad[0] * prj_reconstruct_fine_offset(ii) +
                grad[1] * prj_reconstruct_fine_offset(jj);
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

/* ---------- Riemann face-state reconstruction (compile-time selectable) ----------
 *
 * PRJ_RECON_HYDRO selects the reconstruction scheme used to build hydro/MHD
 * primitive face states, while PRJ_RECON_RADIATION selects radiation primitive
 * face states:
 *
 *   PRJ_RECON_* == PRJ_RECON_MC     monotonized-central piecewise linear
 *   PRJ_RECON_* == PRJ_RECON_WENO3  3rd-order WENO-Z from a 3-cell stencil
 *   PRJ_RECON_* == PRJ_RECON_WENO7  7th-order WENO-Z from a 7-cell stencil
 *                                   (requires PRJ_NGHOST >= 4)
 *
 * Dispatchers take a stencil of `ncells` cell-averaged values centered on the
 * source cell (q[ncells / 2] is q_i) and a `target` offset in cell-widths from
 * q_i: +0.5 for the left state of the i+1/2 face, -0.5 for the right state of
 * the i-1/2 face. The hydro/EOS dispatcher uses PRJ_RECON_HYDRO_NCELLS and the
 * radiation dispatcher uses PRJ_RECON_RADIATION_NCELLS, so each physics carries
 * only the stencil width its own scheme requires.
 */

/* MC (monotonized-central) limited slope used by prj_reconstruct_mc_face, from
 * the 3-cell stencil (qm, q0, qp). The face value is q0 + target * slope, so both
 * faces of a cell (target = +/-0.5) share one slope -- callers that need both can
 * compute the slope once. (Distinct from prj_reconstruct_mc_slope above, which is
 * the prolongation limiter.) */
static inline double prj_reconstruct_mc_face_slope(double qm, double q0, double qp)
{
    double sl = q0 - qm;
    double sr = qp - q0;
    double c0 = 2.0 * sl;
    double c1 = sl + 0.5 * (sr - sl);
    double c2 = 2.0 * sr;
    /* Branchless, bit-identical to the if-form: when sl,sr share a sign the
     * limited slope is min(c0,c1,c2) (positive) or max(c0,c1,c2) (negative),
     * otherwise 0. Compute both min and max unconditionally (the ternaries lower
     * to vmin/vmax), select by sign, and mask out the opposite-sign/zero case --
     * no data-dependent branch, so the reconstruction loop can vectorize. */
    double mn = c0 < c1 ? c0 : c1;
    double mx = c0 > c1 ? c0 : c1;
    double slope;

    mn = mn < c2 ? mn : c2;
    mx = mx > c2 ? mx : c2;
    slope = sl > 0.0 ? mn : mx;
    return (sl * sr > 0.0) ? slope : 0.0;
}

static inline double prj_reconstruct_mc_face(double qm, double q0, double qp,
    double target)
{
    return q0 + target * prj_reconstruct_mc_face_slope(qm, q0, qp);
}

/* WENO-Z3 reconstruction (Don & Borges 2013) of a face value from the 3-cell
 * stencil (qm, q0, qp). `target` selects which face of cell i: +0.5 for the
 * i+1/2 face (left state) and -0.5 for the i-1/2 face (right state). The
 * Z-weights use the global smoothness indicator tau3 = |beta0 - beta1| to
 * approach the optimal linear weights faster than Jiang–Shu near smooth and
 * critical-point regions:
 *
 *   alpha_k = gamma_k * (1 + tau3 / (eps + beta_k))
 *
 * (p = 1; the standard choice for the third-order variant.) */
static inline double prj_reconstruct_weno3_face(double qm, double q0, double qp,
    double target)
{
    const double eps = 1.0e-40;
    double p0 = q0 + target * (q0 - qm);  /* left-biased sub-stencil (qm, q0)  */
    double p1 = q0 + target * (qp - q0);  /* right-biased sub-stencil (q0, qp) */
    double dl = q0 - qm;
    double dr = qp - q0;
    double beta0 = dl * dl;
    double beta1 = dr * dr;
    double tau3 = beta0 - beta1;
    double gamma0;
    double gamma1;
    double alpha0;
    double alpha1;
    double sum;

    if (tau3 < 0.0) {
        tau3 = -tau3;
    }
    if (target > 0.0) {
        gamma0 = 1.0 / 3.0;
        gamma1 = 2.0 / 3.0;
    } else {
        gamma0 = 2.0 / 3.0;
        gamma1 = 1.0 / 3.0;
    }
    alpha0 = gamma0 * (1.0 + tau3 / (eps + beta0));
    alpha1 = gamma1 * (1.0 + tau3 / (eps + beta1));
    sum = alpha0 + alpha1;
    return (alpha0 * p0 + alpha1 * p1) / sum;
}

static inline double prj_reconstruct_weno7_beta0(double q0, double q1,
    double q2, double q3)
{
    return (547.0 * q0 * q0 - 3882.0 * q0 * q1 + 4642.0 * q0 * q2 -
        1854.0 * q0 * q3 + 7043.0 * q1 * q1 - 17246.0 * q1 * q2 +
        7042.0 * q1 * q3 + 11003.0 * q2 * q2 - 9402.0 * q2 * q3 +
        2107.0 * q3 * q3) / 240.0;
}

static inline double prj_reconstruct_weno7_beta1(double q0, double q1,
    double q2, double q3)
{
    return (267.0 * q0 * q0 - 1642.0 * q0 * q1 + 1602.0 * q0 * q2 -
        494.0 * q0 * q3 + 2843.0 * q1 * q1 - 5966.0 * q1 * q2 +
        1922.0 * q1 * q3 + 3443.0 * q2 * q2 - 2522.0 * q2 * q3 +
        547.0 * q3 * q3) / 240.0;
}

static inline double prj_reconstruct_weno7_beta2(double q0, double q1,
    double q2, double q3)
{
    return (547.0 * q0 * q0 - 2522.0 * q0 * q1 + 1922.0 * q0 * q2 -
        494.0 * q0 * q3 + 3443.0 * q1 * q1 - 5966.0 * q1 * q2 +
        1602.0 * q1 * q3 + 2843.0 * q2 * q2 - 1642.0 * q2 * q3 +
        267.0 * q3 * q3) / 240.0;
}

static inline double prj_reconstruct_weno7_beta3(double q0, double q1,
    double q2, double q3)
{
    return (2107.0 * q0 * q0 - 9402.0 * q0 * q1 + 7042.0 * q0 * q2 -
        1854.0 * q0 * q3 + 11003.0 * q1 * q1 - 17246.0 * q1 * q2 +
        4642.0 * q1 * q3 + 7043.0 * q2 * q2 - 3882.0 * q2 * q3 +
        547.0 * q3 * q3) / 240.0;
}

static inline double prj_reconstruct_weno7_face_plus(const double q[7])
{
    const double eps = 1.0e-40;
    double p0 = (-3.0 * q[0] + 13.0 * q[1] - 23.0 * q[2] +
        25.0 * q[3]) / 12.0;
    double p1 = (q[1] - 5.0 * q[2] + 13.0 * q[3] + 3.0 * q[4]) / 12.0;
    double p2 = (-q[2] + 7.0 * q[3] + 7.0 * q[4] - q[5]) / 12.0;
    double p3 = (3.0 * q[3] + 13.0 * q[4] - 5.0 * q[5] + q[6]) / 12.0;
    double beta0 = prj_reconstruct_weno7_beta0(q[0], q[1], q[2], q[3]);
    double beta1 = prj_reconstruct_weno7_beta1(q[1], q[2], q[3], q[4]);
    double beta2 = prj_reconstruct_weno7_beta2(q[2], q[3], q[4], q[5]);
    double beta3 = prj_reconstruct_weno7_beta3(q[3], q[4], q[5], q[6]);
    double tau7 = beta0 - beta3 + 3.0 * (beta1 - beta2);
    double r0;
    double r1;
    double r2;
    double r3;
    double alpha0;
    double alpha1;
    double alpha2;
    double alpha3;
    double sum;

    if (tau7 < 0.0) {
        tau7 = -tau7;
    }

    r0 = tau7 / (eps + beta0);
    r1 = tau7 / (eps + beta1);
    r2 = tau7 / (eps + beta2);
    r3 = tau7 / (eps + beta3);
    alpha0 = (1.0 / 35.0) * (1.0 + r0 * r0);
    alpha1 = (12.0 / 35.0) * (1.0 + r1 * r1);
    alpha2 = (18.0 / 35.0) * (1.0 + r2 * r2);
    alpha3 = (4.0 / 35.0) * (1.0 + r3 * r3);
    sum = alpha0 + alpha1 + alpha2 + alpha3;

    return (alpha0 * p0 + alpha1 * p1 + alpha2 * p2 + alpha3 * p3) / sum;
}

/* WENO-Z7 reconstruction of the face value from q[i-3]...q[i+3]. The
 * positive face uses the left-biased seventh-order candidate set; the negative
 * face is the same formula on the mirrored stencil. */
static inline double prj_reconstruct_weno7_face(const double q[7], double target)
{
    double qr[7];

    if (target > 0.0) {
        return prj_reconstruct_weno7_face_plus(q);
    }
    if (target < 0.0) {
        int n;

        for (n = 0; n < 7; ++n) {
            qr[n] = q[6 - n];
        }
        return prj_reconstruct_weno7_face_plus(qr);
    }
    return q[3];
}

/* MC and WENO3 read a 3-cell stencil (q_i at q[1]); WENO7 reads a 7-cell
 * stencil (q_i at q[3]). Each helper owns its stencil width, so callers only
 * have to center q_i at q[ncells / 2] for the matching scheme. */
static inline double prj_reconstruct_mc_face_value(const double q[3],
    double target)
{
    return prj_reconstruct_mc_face(q[0], q[1], q[2], target);
}

static inline double prj_reconstruct_weno3_face_value(const double q[3],
    double target)
{
    return prj_reconstruct_weno3_face(q[0], q[1], q[2], target);
}

static inline double prj_reconstruct_weno7_face_value(const double q[7],
    double target)
{
    return prj_reconstruct_weno7_face(q, target);
}

static inline double prj_reconstruct_hydro_face_value(const double q[PRJ_RECON_HYDRO_NCELLS],
    double target)
{
#if PRJ_RECON_HYDRO == PRJ_RECON_WENO3
    return prj_reconstruct_weno3_face_value(q, target);
#elif PRJ_RECON_HYDRO == PRJ_RECON_WENO7
    return prj_reconstruct_weno7_face_value(q, target);
#else
    return prj_reconstruct_mc_face_value(q, target);
#endif
}

static inline double prj_reconstruct_radiation_face_value(const double q[PRJ_RECON_RADIATION_NCELLS],
    double target)
{
#if PRJ_RECON_RADIATION == PRJ_RECON_WENO3
    return prj_reconstruct_weno3_face_value(q, target);
#elif PRJ_RECON_RADIATION == PRJ_RECON_WENO7
    return prj_reconstruct_weno7_face_value(q, target);
#else
    return prj_reconstruct_mc_face_value(q, target);
#endif
}

int prj_reconstruct_step(void);

#endif
