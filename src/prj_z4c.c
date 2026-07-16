#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(PRJ_ENABLE_MPI)
#include <mpi.h>
#endif

#include "prj.h"

#define PRJ_Z4C_MPI_TAG 122

static const char *const prj_z4c_names[PRJ_NZ4C] = {
    "z4c_chi",
    "z4c_gxx", "z4c_gxy", "z4c_gxz", "z4c_gyy", "z4c_gyz", "z4c_gzz",
    "z4c_Khat",
    "z4c_Axx", "z4c_Axy", "z4c_Axz", "z4c_Ayy", "z4c_Ayz", "z4c_Azz",
    "z4c_Gamx", "z4c_Gamy", "z4c_Gamz",
    "z4c_Theta",
    "z4c_alpha",
    "z4c_betax", "z4c_betay", "z4c_betaz"
};

static const char *const prj_z4c_tmunu_names[PRJ_NTMUNU] = {
    "z4c_tmunu_Sxx", "z4c_tmunu_Sxy", "z4c_tmunu_Sxz",
    "z4c_tmunu_Syy", "z4c_tmunu_Syz", "z4c_tmunu_Szz",
    "z4c_tmunu_E", "z4c_tmunu_Sx", "z4c_tmunu_Sy", "z4c_tmunu_Sz"
};

const char *prj_z4c_var_name(int var)
{
    return var >= 0 && var < PRJ_NZ4C ? prj_z4c_names[var] : "z4c_unknown";
}

const char *prj_z4c_tmunu_var_name(int var)
{
    return var >= 0 && var < PRJ_NTMUNU ? prj_z4c_tmunu_names[var] :
        "z4c_tmunu_unknown";
}

void prj_z4c_init_params(prj_z4c_params *params)
{
    if (params == 0) {
        return;
    }
    params->chi_psi_power = -4.0;
    params->chi_div_floor = -1000.0;
    params->chi_min_floor = 1.0e-12;
    params->floor_chi = 0;
    params->diss = 0.0;
    params->damp_kappa1_inv_cm = 0.0;
    params->damp_kappa2 = 0.0;
    params->lapse_harmonicf = 1.0;
    params->lapse_harmonic = 0.0;
    params->lapse_oplog = 2.0;
    params->lapse_advect = 1.0;
    params->slow_start_lapse = 0;
    params->ssl_damping_amp_inv_cm = 0.0;
    params->ssl_damping_time_cm = 1.0;
    params->ssl_damping_index = 1;
    params->shift_Gamma = 1.0;
    params->shift_advect = 1.0;
    params->shift_alpha2Gamma = 0.0;
    params->shift_H = 0.0;
    params->shift_eta_inv_cm = 0.0;
    params->puncture_mass_cm = 1.0;
    params->puncture_center_cm[0] = 0.0;
    params->puncture_center_cm[1] = 0.0;
    params->puncture_center_cm[2] = 0.0;
    params->puncture_floor_radius_cm = 0.0;
    params->puncture_mass_plus_cm = 0.483;
    params->puncture_mass_minus_cm = 0.483;
    params->puncture_half_separation_cm = 3.257;
    params->puncture_momentum_plus_cm[0] = 0.0;
    params->puncture_momentum_plus_cm[1] = -0.133;
    params->puncture_momentum_plus_cm[2] = 0.0;
    params->puncture_momentum_minus_cm[0] = 0.0;
    params->puncture_momentum_minus_cm[1] = 0.133;
    params->puncture_momentum_minus_cm[2] = 0.0;
    params->use_z4c = 1;
    params->user_Sbc = 0;
}

#if PRJ_DYNAMIC_GR

static void prj_z4c_fail(const char *message)
{
    fprintf(stderr, "%s\n", message);
#if defined(PRJ_ENABLE_MPI)
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#endif
    exit(EXIT_FAILURE);
}

int prj_z4c_runtime_enabled(const prj_mesh *mesh)
{
    return mesh != 0;
}

static int prj_z4c_local_block(const prj_mpi *mpi, const prj_block *block)
{
    return block != 0 && block->id >= 0 && block->active == 1 &&
        block->z4c != 0 && block->z4c_rhs != 0 &&
        (mpi == 0 || block->rank == mpi->rank);
}

static int prj_z4c_g_var(int a, int b)
{
    if (a > b) {
        int t = a;
        a = b;
        b = t;
    }
    if (a == 0 && b == 0) return PRJ_Z4C_GXX;
    if (a == 0 && b == 1) return PRJ_Z4C_GXY;
    if (a == 0 && b == 2) return PRJ_Z4C_GXZ;
    if (a == 1 && b == 1) return PRJ_Z4C_GYY;
    if (a == 1 && b == 2) return PRJ_Z4C_GYZ;
    return PRJ_Z4C_GZZ;
}

static int prj_z4c_A_var(int a, int b)
{
    if (a > b) {
        int t = a;
        a = b;
        b = t;
    }
    if (a == 0 && b == 0) return PRJ_Z4C_AXX;
    if (a == 0 && b == 1) return PRJ_Z4C_AXY;
    if (a == 0 && b == 2) return PRJ_Z4C_AXZ;
    if (a == 1 && b == 1) return PRJ_Z4C_AYY;
    if (a == 1 && b == 2) return PRJ_Z4C_AYZ;
    return PRJ_Z4C_AZZ;
}

static int prj_z4c_Gam_var(int a)
{
    return PRJ_Z4C_GAMX + a;
}

static int prj_z4c_beta_var(int a)
{
    return PRJ_Z4C_BETAX + a;
}

static int prj_z4c_tmunu_stress_var(int a, int b)
{
    static const int vars[3][3] = {
        {PRJ_TMUNU_SXX, PRJ_TMUNU_SXY, PRJ_TMUNU_SXZ},
        {PRJ_TMUNU_SXY, PRJ_TMUNU_SYY, PRJ_TMUNU_SYZ},
        {PRJ_TMUNU_SXZ, PRJ_TMUNU_SYZ, PRJ_TMUNU_SZZ}
    };
    return vars[a][b];
}

static int prj_z4c_tmunu_mom_var(int a)
{
    return PRJ_TMUNU_SX + a;
}

static double prj_z4c_tmunu_cgs_to_geo(double value)
{
    double c2 = PRJ_CLIGHT * PRJ_CLIGHT;
    return value * PRJ_GNEWT / (c2 * c2);
}

static double prj_z4c_get(const double *z, int v, int i, int j, int k)
{
    return z[Z4CIDX(v, i, j, k)];
}

static void prj_z4c_set(double *z, int v, int i, int j, int k, double value)
{
    z[Z4CIDX(v, i, j, k)] = value;
}

static int prj_z4c_cell_in_storage(int i, int j, int k)
{
    return i >= -PRJ_NGHOST_Z4C && i < PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C &&
        j >= -PRJ_NGHOST_Z4C && j < PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C &&
        k >= -PRJ_NGHOST_Z4C && k < PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C;
}

static void prj_z4c_shift_index(int dir, int offset, int *i, int *j, int *k)
{
    if (dir == X1DIR) {
        *i += offset;
    } else if (dir == X2DIR) {
        *j += offset;
    } else {
        *k += offset;
    }
}

static void prj_z4c_fd1_coeff(int *n, int off[6], double c[6])
{
#if PRJ_NGHOST_Z4C == 2
    *n = 2;
    off[0] = -1; c[0] = -0.5;
    off[1] =  1; c[1] =  0.5;
#elif PRJ_NGHOST_Z4C == 3
    *n = 4;
    off[0] = -2; c[0] =  1.0 / 12.0;
    off[1] = -1; c[1] = -2.0 / 3.0;
    off[2] =  1; c[2] =  2.0 / 3.0;
    off[3] =  2; c[3] = -1.0 / 12.0;
#elif PRJ_NGHOST_Z4C == 4
    *n = 6;
    off[0] = -3; c[0] = -1.0 / 60.0;
    off[1] = -2; c[1] =  3.0 / 20.0;
    off[2] = -1; c[2] = -3.0 / 4.0;
    off[3] =  1; c[3] =  3.0 / 4.0;
    off[4] =  2; c[4] = -3.0 / 20.0;
    off[5] =  3; c[5] =  1.0 / 60.0;
#else
#error "Z4c finite differences support PRJ_NGHOST_Z4C == 2, 3, or 4"
#endif
}

static void prj_z4c_fd2_coeff(int *n, int off[7], double c[7])
{
#if PRJ_NGHOST_Z4C == 2
    *n = 3;
    off[0] = -1; c[0] =  1.0;
    off[1] =  0; c[1] = -2.0;
    off[2] =  1; c[2] =  1.0;
#elif PRJ_NGHOST_Z4C == 3
    *n = 5;
    off[0] = -2; c[0] = -1.0 / 12.0;
    off[1] = -1; c[1] =  4.0 / 3.0;
    off[2] =  0; c[2] = -5.0 / 2.0;
    off[3] =  1; c[3] =  4.0 / 3.0;
    off[4] =  2; c[4] = -1.0 / 12.0;
#elif PRJ_NGHOST_Z4C == 4
    *n = 7;
    off[0] = -3; c[0] =  1.0 / 90.0;
    off[1] = -2; c[1] = -3.0 / 20.0;
    off[2] = -1; c[2] =  3.0 / 2.0;
    off[3] =  0; c[3] = -49.0 / 18.0;
    off[4] =  1; c[4] =  3.0 / 2.0;
    off[5] =  2; c[5] = -3.0 / 20.0;
    off[6] =  3; c[6] =  1.0 / 90.0;
#else
#error "Z4c finite differences support PRJ_NGHOST_Z4C == 2, 3, or 4"
#endif
}

static double prj_z4c_Dx(const double *z, int v, int dir, const double idx[3],
    int i, int j, int k)
{
    int off[6];
    double c[6];
    int n;
    double out = 0.0;
    int p;

    prj_z4c_fd1_coeff(&n, off, c);
    for (p = 0; p < n; ++p) {
        int ii = i, jj = j, kk = k;
        prj_z4c_shift_index(dir, off[p], &ii, &jj, &kk);
        out += c[p] * prj_z4c_get(z, v, ii, jj, kk);
    }
    return out * idx[dir];
}

static double prj_z4c_Dx2nd(const double *z, int v, int dir, const double idx[3],
    int i, int j, int k)
{
    int im = i, jm = j, km = k;
    int ip = i, jp = j, kp = k;

    prj_z4c_shift_index(dir, -1, &im, &jm, &km);
    prj_z4c_shift_index(dir, 1, &ip, &jp, &kp);
    return 0.5 * idx[dir] *
        (prj_z4c_get(z, v, ip, jp, kp) - prj_z4c_get(z, v, im, jm, km));
}

static double prj_z4c_Dxx(const double *z, int v, int dir, const double idx[3],
    int i, int j, int k)
{
    int off[7];
    double c[7];
    int n;
    double out = 0.0;
    int p;

    prj_z4c_fd2_coeff(&n, off, c);
    for (p = 0; p < n; ++p) {
        int ii = i, jj = j, kk = k;
        prj_z4c_shift_index(dir, off[p], &ii, &jj, &kk);
        out += c[p] * prj_z4c_get(z, v, ii, jj, kk);
    }
    return out * idx[dir] * idx[dir];
}

static double prj_z4c_Dxy(const double *z, int v, int dirx, int diry,
    const double idx[3], int i, int j, int k)
{
    int offx[6], offy[6];
    double cx[6], cy[6];
    int nx, ny;
    double out = 0.0;
    int px, py;

    prj_z4c_fd1_coeff(&nx, offx, cx);
    prj_z4c_fd1_coeff(&ny, offy, cy);
    for (px = 0; px < nx; ++px) {
        for (py = 0; py < ny; ++py) {
            int ii = i, jj = j, kk = k;
            prj_z4c_shift_index(dirx, offx[px], &ii, &jj, &kk);
            prj_z4c_shift_index(diry, offy[py], &ii, &jj, &kk);
            out += cx[px] * cy[py] * prj_z4c_get(z, v, ii, jj, kk);
        }
    }
    return out * idx[dirx] * idx[diry];
}

static double prj_z4c_Lx(const double *z, int beta_var, int q_var, int dir,
    const double idx[3], int i, int j, int k)
{
    double beta = prj_z4c_get(z, beta_var, i, j, k);
    double q0 = prj_z4c_get(z, q_var, i, j, k);
    double dl;
    double dr;

#if PRJ_NGHOST_Z4C == 2
    {
        int im2 = i, jm2 = j, km2 = k, im1 = i, jm1 = j, km1 = k;
        int ip1 = i, jp1 = j, kp1 = k, ip2 = i, jp2 = j, kp2 = k;
        prj_z4c_shift_index(dir, -2, &im2, &jm2, &km2);
        prj_z4c_shift_index(dir, -1, &im1, &jm1, &km1);
        prj_z4c_shift_index(dir, 1, &ip1, &jp1, &kp1);
        prj_z4c_shift_index(dir, 2, &ip2, &jp2, &kp2);
        dl = 0.5 * prj_z4c_get(z, q_var, im2, jm2, km2) -
            2.0 * prj_z4c_get(z, q_var, im1, jm1, km1) + 1.5 * q0;
        dr = -0.5 * prj_z4c_get(z, q_var, ip2, jp2, kp2) +
            2.0 * prj_z4c_get(z, q_var, ip1, jp1, kp1) - 1.5 * q0;
    }
#elif PRJ_NGHOST_Z4C == 3
    {
        int im3 = i, jm3 = j, km3 = k, im2 = i, jm2 = j, km2 = k;
        int im1 = i, jm1 = j, km1 = k, ip1 = i, jp1 = j, kp1 = k;
        int ip2 = i, jp2 = j, kp2 = k, ip3 = i, jp3 = j, kp3 = k;
        prj_z4c_shift_index(dir, -3, &im3, &jm3, &km3);
        prj_z4c_shift_index(dir, -2, &im2, &jm2, &km2);
        prj_z4c_shift_index(dir, -1, &im1, &jm1, &km1);
        prj_z4c_shift_index(dir, 1, &ip1, &jp1, &kp1);
        prj_z4c_shift_index(dir, 2, &ip2, &jp2, &kp2);
        prj_z4c_shift_index(dir, 3, &ip3, &jp3, &kp3);
        dl = (-prj_z4c_get(z, q_var, im3, jm3, km3) +
            6.0 * prj_z4c_get(z, q_var, im2, jm2, km2) -
            18.0 * prj_z4c_get(z, q_var, im1, jm1, km1) + 10.0 * q0 +
            3.0 * prj_z4c_get(z, q_var, ip1, jp1, kp1)) / 12.0;
        dr = (prj_z4c_get(z, q_var, ip3, jp3, kp3) -
            6.0 * prj_z4c_get(z, q_var, ip2, jp2, kp2) +
            18.0 * prj_z4c_get(z, q_var, ip1, jp1, kp1) - 10.0 * q0 -
            3.0 * prj_z4c_get(z, q_var, im1, jm1, km1)) / 12.0;
    }
#else
    {
        int im4 = i, jm4 = j, km4 = k, im3 = i, jm3 = j, km3 = k;
        int im2 = i, jm2 = j, km2 = k, im1 = i, jm1 = j, km1 = k;
        int ip1 = i, jp1 = j, kp1 = k, ip2 = i, jp2 = j, kp2 = k;
        int ip3 = i, jp3 = j, kp3 = k, ip4 = i, jp4 = j, kp4 = k;
        prj_z4c_shift_index(dir, -4, &im4, &jm4, &km4);
        prj_z4c_shift_index(dir, -3, &im3, &jm3, &km3);
        prj_z4c_shift_index(dir, -2, &im2, &jm2, &km2);
        prj_z4c_shift_index(dir, -1, &im1, &jm1, &km1);
        prj_z4c_shift_index(dir, 1, &ip1, &jp1, &kp1);
        prj_z4c_shift_index(dir, 2, &ip2, &jp2, &kp2);
        prj_z4c_shift_index(dir, 3, &ip3, &jp3, &kp3);
        prj_z4c_shift_index(dir, 4, &ip4, &jp4, &kp4);
        dl = prj_z4c_get(z, q_var, im4, jm4, km4) / 60.0 -
            2.0 * prj_z4c_get(z, q_var, im3, jm3, km3) / 15.0 +
            0.5 * prj_z4c_get(z, q_var, im2, jm2, km2) -
            4.0 * prj_z4c_get(z, q_var, im1, jm1, km1) / 3.0 +
            7.0 * q0 / 12.0 + 2.0 * prj_z4c_get(z, q_var, ip1, jp1, kp1) / 5.0 -
            prj_z4c_get(z, q_var, ip2, jp2, kp2) / 30.0;
        dr = -prj_z4c_get(z, q_var, ip4, jp4, kp4) / 60.0 +
            2.0 * prj_z4c_get(z, q_var, ip3, jp3, kp3) / 15.0 -
            0.5 * prj_z4c_get(z, q_var, ip2, jp2, kp2) +
            4.0 * prj_z4c_get(z, q_var, ip1, jp1, kp1) / 3.0 -
            7.0 * q0 / 12.0 - 2.0 * prj_z4c_get(z, q_var, im1, jm1, km1) / 5.0 +
            prj_z4c_get(z, q_var, im2, jm2, km2) / 30.0;
    }
#endif
    return (beta < 0.0 ? beta * dl : beta * dr) * idx[dir];
}

static double prj_z4c_Diss(const double *z, int v, int dir, const double idx[3],
    int i, int j, int k)
{
#if PRJ_NGHOST_Z4C == 2
    const int n = 5;
    int offs[5] = {-2, -1, 0, 1, 2};
    double coef[5] = {1.0, -4.0, 6.0, -4.0, 1.0};
#elif PRJ_NGHOST_Z4C == 3
    const int n = 7;
    int offs[7] = {-3, -2, -1, 0, 1, 2, 3};
    double coef[7] = {1.0, -6.0, 15.0, -20.0, 15.0, -6.0, 1.0};
#else
    const int n = 9;
    int offs[9] = {-4, -3, -2, -1, 0, 1, 2, 3, 4};
    double coef[9] = {1.0, -8.0, 28.0, -56.0, 70.0, -56.0, 28.0, -8.0, 1.0};
#endif
    double out = 0.0;
    int p;

    for (p = 0; p < n; ++p) {
        int ii = i, jj = j, kk = k;
        prj_z4c_shift_index(dir, offs[p], &ii, &jj, &kk);
        out += coef[p] * prj_z4c_get(z, v, ii, jj, kk);
    }
    return out * idx[dir];
}

static double prj_z4c_diss_coeff(const prj_z4c_params *params)
{
    double coeff;

    if (params == 0 || params->diss == 0.0) {
        return 0.0;
    }
    coeff = params->diss * pow(2.0, -2.0 * (double)PRJ_NGHOST_Z4C);
#if (PRJ_NGHOST_Z4C % 2) == 0
    coeff = -coeff;
#endif
    return coeff;
}

static double prj_z4c_det3(const double g[3][3])
{
    return g[0][0] * (g[1][1] * g[2][2] - g[1][2] * g[2][1])
        - g[0][1] * (g[1][0] * g[2][2] - g[1][2] * g[2][0])
        + g[0][2] * (g[1][0] * g[2][1] - g[1][1] * g[2][0]);
}

static void prj_z4c_inv3(const double g[3][3], double inv[3][3], double *det_out)
{
    double det = prj_z4c_det3(g);
    double odet;

    if (!isfinite(det) || fabs(det) < 1.0e-300) {
        det = det < 0.0 ? -1.0e-300 : 1.0e-300;
    }
    odet = 1.0 / det;
    inv[0][0] = (g[1][1] * g[2][2] - g[1][2] * g[2][1]) * odet;
    inv[0][1] = (g[0][2] * g[2][1] - g[0][1] * g[2][2]) * odet;
    inv[0][2] = (g[0][1] * g[1][2] - g[0][2] * g[1][1]) * odet;
    inv[1][0] = inv[0][1];
    inv[1][1] = (g[0][0] * g[2][2] - g[0][2] * g[2][0]) * odet;
    inv[1][2] = (g[0][2] * g[1][0] - g[0][0] * g[1][2]) * odet;
    inv[2][0] = inv[0][2];
    inv[2][1] = inv[1][2];
    inv[2][2] = (g[0][0] * g[1][1] - g[0][1] * g[1][0]) * odet;
    if (det_out != 0) {
        *det_out = det;
    }
}

static void prj_z4c_load_metric(const double *z, int i, int j, int k, double g[3][3])
{
    int a, b;

    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            g[a][b] = prj_z4c_get(z, prj_z4c_g_var(a, b), i, j, k);
        }
    }
}

static double prj_z4c_guarded_chi(const prj_z4c_params *opt, double chi)
{
    double floor = opt != 0 ? opt->chi_min_floor : 1.0e-12;
    double div_floor = opt != 0 ? opt->chi_div_floor : -1000.0;

    if (!isfinite(chi)) {
        return floor;
    }
    if (chi <= div_floor) {
        chi = div_floor;
    }
    if (chi <= floor) {
        chi = floor;
    }
    return chi;
}

static void prj_z4c_load_A(const double *z, int i, int j, int k, double A[3][3]);

static int prj_z4c_fd1_stencil_in_storage(int dir, int i, int j, int k)
{
    int off[6];
    double c[6];
    int n;
    int p;

    prj_z4c_fd1_coeff(&n, off, c);
    (void)c;
    for (p = 0; p < n; ++p) {
        int ii = i;
        int jj = j;
        int kk = k;

        prj_z4c_shift_index(dir, off[p], &ii, &jj, &kk);
        if (!prj_z4c_cell_in_storage(ii, jj, kk)) {
            return 0;
        }
    }
    return 1;
}

int prj_z4c_load_hydro_geom(const prj_mesh *mesh, const prj_block *block,
    int stage, int i, int j, int k, prj_z4c_hydro_geom *geom)
{
    const prj_z4c_params *opt;
    const double *z;
    double g[3][3];
    double gu[3][3];
    double det_g;
    double idx[3];
    double chi;
    double chi_guarded;
    double chi_power;
    double inv_factor;
    double factor;
    double dchi[3] = {0.0, 0.0, 0.0};
    int deriv_ok[3];
    int a;
    int b;
    int dir;

    if (!prj_z4c_runtime_enabled(mesh) || block == 0 || geom == 0 ||
        !prj_z4c_cell_in_storage(i, j, k)) {
        return 0;
    }
    z = prj_block_z4c_stage_const(block, prj_stage_slot_from_bf_arg(stage));
    if (z == 0) {
        return 0;
    }
    opt = &mesh->z4c_params;
    for (dir = 0; dir < 3; ++dir) {
        if (block->dx[dir] == 0.0) {
            return 0;
        }
        idx[dir] = 1.0 / block->dx[dir];
        deriv_ok[dir] = prj_z4c_fd1_stencil_in_storage(dir, i, j, k);
    }

    prj_z4c_load_metric(z, i, j, k, g);
    prj_z4c_inv3(g, gu, &det_g);
    if (!isfinite(det_g) || det_g <= 0.0) {
        return 0;
    }

    chi = prj_z4c_get(z, PRJ_Z4C_CHI, i, j, k);
    chi_guarded = prj_z4c_guarded_chi(opt, chi);
    chi_power = opt->chi_psi_power;
    if (!isfinite(chi_power) || fabs(chi_power) < 1.0e-300) {
        chi_power = -4.0;
    }
    inv_factor = pow(chi_guarded, -4.0 / chi_power);
    if (!isfinite(inv_factor) || inv_factor <= 0.0) {
        return 0;
    }
    factor = 1.0 / inv_factor;

    geom->sqrt_gamma = sqrt(det_g * factor * factor * factor);
    if (!isfinite(geom->sqrt_gamma) || geom->sqrt_gamma <= 0.0) {
        return 0;
    }
    geom->alpha = prj_z4c_get(z, PRJ_Z4C_ALPHA, i, j, k);
    if (!isfinite(geom->alpha) || geom->alpha <= 0.0) {
        geom->alpha = 1.0e-12;
    }
    for (a = 0; a < 3; ++a) {
        geom->beta[a] = prj_z4c_get(z, prj_z4c_beta_var(a), i, j, k);
        if (!isfinite(geom->beta[a])) {
            geom->beta[a] = 0.0;
        }
    }
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            geom->gamma[a][b] = factor * g[a][b];
            geom->gamma_inv[a][b] = inv_factor * gu[a][b];
        }
    }

    {
        double A[3][3];
        double Ktrace = prj_z4c_get(z, PRJ_Z4C_KHAT, i, j, k) +
            2.0 * prj_z4c_get(z, PRJ_Z4C_THETA, i, j, k);

        prj_z4c_load_A(z, i, j, k, A);
        for (a = 0; a < 3; ++a) {
            for (b = 0; b < 3; ++b) {
                geom->K_dd[a][b] = factor * A[a][b] +
                    (1.0 / 3.0) * geom->gamma[a][b] * Ktrace;
            }
        }
    }

    for (dir = 0; dir < 3; ++dir) {
        geom->dalpha[dir] = 0.0;
        for (a = 0; a < 3; ++a) {
            geom->dbeta[dir][a] = 0.0;
        }
        for (a = 0; a < 3; ++a) {
            for (b = 0; b < 3; ++b) {
                geom->dgamma[dir][a][b] = 0.0;
            }
        }
        if (!deriv_ok[dir]) {
            continue;
        }
        dchi[dir] = prj_z4c_Dx(z, PRJ_Z4C_CHI, dir, idx, i, j, k);
        geom->dalpha[dir] = prj_z4c_Dx(z, PRJ_Z4C_ALPHA, dir, idx, i, j, k);
        for (a = 0; a < 3; ++a) {
            geom->dbeta[dir][a] =
                prj_z4c_Dx(z, prj_z4c_beta_var(a), dir, idx, i, j, k);
        }
        for (a = 0; a < 3; ++a) {
            for (b = 0; b < 3; ++b) {
                double dg = prj_z4c_Dx(z, prj_z4c_g_var(a, b), dir, idx, i, j, k);
                double dfactor = factor * (4.0 / chi_power) *
                    dchi[dir] / chi_guarded;

                geom->dgamma[dir][a][b] = factor * dg + dfactor * g[a][b];
            }
        }
    }
    return 1;
}

int prj_z4c_load_hydro_metric_geom(const prj_mesh *mesh, const prj_block *block,
    int stage, int i, int j, int k, prj_z4c_hydro_geom *geom)
{
    const prj_z4c_params *opt;
    const double *z;
    double g[3][3];
    double gu[3][3];
    double det_g;
    double chi;
    double chi_guarded;
    double chi_power;
    double inv_factor;
    double factor;
    int a;
    int b;

    if (!prj_z4c_runtime_enabled(mesh) || block == 0 || geom == 0 ||
        !prj_z4c_cell_in_storage(i, j, k)) {
        return 0;
    }
    z = prj_block_z4c_stage_const(block, prj_stage_slot_from_bf_arg(stage));
    if (z == 0) {
        return 0;
    }
    opt = &mesh->z4c_params;

    prj_z4c_load_metric(z, i, j, k, g);
    prj_z4c_inv3(g, gu, &det_g);
    if (!isfinite(det_g) || det_g <= 0.0) {
        return 0;
    }

    chi = prj_z4c_get(z, PRJ_Z4C_CHI, i, j, k);
    chi_guarded = prj_z4c_guarded_chi(opt, chi);
    chi_power = opt->chi_psi_power;
    if (!isfinite(chi_power) || fabs(chi_power) < 1.0e-300) {
        chi_power = -4.0;
    }
    inv_factor = pow(chi_guarded, -4.0 / chi_power);
    if (!isfinite(inv_factor) || inv_factor <= 0.0) {
        return 0;
    }
    factor = 1.0 / inv_factor;

    geom->sqrt_gamma = sqrt(det_g * factor * factor * factor);
    if (!isfinite(geom->sqrt_gamma) || geom->sqrt_gamma <= 0.0) {
        return 0;
    }
    geom->alpha = prj_z4c_get(z, PRJ_Z4C_ALPHA, i, j, k);
    if (!isfinite(geom->alpha) || geom->alpha <= 0.0) {
        geom->alpha = 1.0e-12;
    }
    for (a = 0; a < 3; ++a) {
        geom->beta[a] = prj_z4c_get(z, prj_z4c_beta_var(a), i, j, k);
        if (!isfinite(geom->beta[a])) {
            geom->beta[a] = 0.0;
        }
    }
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            geom->gamma[a][b] = factor * g[a][b];
            geom->gamma_inv[a][b] = inv_factor * gu[a][b];
        }
    }
    return 1;
}

static void prj_z4c_load_minkowski_hydro_geom(prj_z4c_hydro_geom *geom)
{
    int a;
    int b;

    memset(geom, 0, sizeof(*geom));
    geom->alpha = 1.0;
    geom->sqrt_gamma = 1.0;
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            geom->gamma[a][b] = a == b ? 1.0 : 0.0;
            geom->gamma_inv[a][b] = a == b ? 1.0 : 0.0;
        }
    }
}

static void prj_z4c_load_A(const double *z, int i, int j, int k, double A[3][3])
{
    int a, b;

    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            A[a][b] = prj_z4c_get(z, prj_z4c_A_var(a, b), i, j, k);
        }
    }
}

static void prj_z4c_enforce_cell(double *z, int i, int j, int k)
{
    double g[3][3];
    double gu[3][3];
    double A[3][3];
    double det;
    double scale;
    double trA = 0.0;
    int a, b;

    prj_z4c_load_metric(z, i, j, k, g);
    det = prj_z4c_det3(g);
    if (!isfinite(det) || det <= 0.0) {
        for (a = 0; a < 3; ++a) {
            for (b = 0; b < 3; ++b) {
                g[a][b] = a == b ? 1.0 : 0.0;
            }
        }
        det = 1.0;
    }
    scale = pow(det, -1.0 / 3.0);
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            g[a][b] *= scale;
        }
    }
    prj_z4c_set(z, PRJ_Z4C_GXX, i, j, k, g[0][0]);
    prj_z4c_set(z, PRJ_Z4C_GXY, i, j, k, g[0][1]);
    prj_z4c_set(z, PRJ_Z4C_GXZ, i, j, k, g[0][2]);
    prj_z4c_set(z, PRJ_Z4C_GYY, i, j, k, g[1][1]);
    prj_z4c_set(z, PRJ_Z4C_GYZ, i, j, k, g[1][2]);
    prj_z4c_set(z, PRJ_Z4C_GZZ, i, j, k, g[2][2]);

    prj_z4c_inv3(g, gu, 0);
    prj_z4c_load_A(z, i, j, k, A);
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            trA += gu[a][b] * A[a][b];
        }
    }
    for (a = 0; a < 3; ++a) {
        for (b = a; b < 3; ++b) {
            prj_z4c_set(z, prj_z4c_A_var(a, b), i, j, k,
                A[a][b] - (1.0 / 3.0) * trA * g[a][b]);
        }
    }
}

static void prj_z4c_enforce_range(prj_mesh *mesh, const prj_mpi *mpi, int stage,
    int ilo, int ihi, int jlo, int jhi, int klo, int khi)
{
    int bidx;

    if (!prj_z4c_runtime_enabled(mesh)) {
        return;
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        double *z;
        int i, j, k;

        if (!prj_z4c_local_block(mpi, block)) {
            continue;
        }
        z = prj_block_z4c_stage(block, stage);
        for (i = ilo; i < ihi; ++i) {
            for (j = jlo; j < jhi; ++j) {
                for (k = klo; k < khi; ++k) {
                    if (mesh->z4c_params.floor_chi &&
                        z[Z4CIDX(PRJ_Z4C_CHI, i, j, k)] < mesh->z4c_params.chi_min_floor) {
                        z[Z4CIDX(PRJ_Z4C_CHI, i, j, k)] = mesh->z4c_params.chi_min_floor;
                    }
                    prj_z4c_enforce_cell(z, i, j, k);
                }
            }
        }
    }
}

double prj_z4c_calc_dt_seconds(const prj_mesh *mesh, const prj_mpi *mpi, double cfl)
{
    double dt_min = HUGE_VAL;
    int bidx;

    if (!prj_z4c_runtime_enabled(mesh)) {
        return HUGE_VAL;
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        const prj_block *block = &mesh->blocks[bidx];
        const prj_z4c_params *opt = &mesh->z4c_params;
        const double *z;
        int i, j, k;

        if (!prj_z4c_local_block(mpi, block)) {
            continue;
        }
        z = prj_block_z4c_stage_const(block, 0);
        if (z == 0) {
            continue;
        }
        for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
            for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    double beta[3];
                    double alpha;
                    double lapse_f;
                    double gauge_speed;
                    double denom;
                    double dt_cell;
                    int a;

                    alpha = prj_z4c_get(z, PRJ_Z4C_ALPHA, i, j, k);
                    lapse_f = opt->lapse_oplog * opt->lapse_harmonicf +
                        opt->lapse_harmonic * alpha;
                    if (!isfinite(alpha) || alpha <= 0.0 ||
                        !isfinite(lapse_f) || lapse_f < 0.0) {
                        denom = 1.0e300;
                    } else {
                        gauge_speed = sqrt(alpha * lapse_f);
                        denom = 0.0;
                        for (a = 0; a < 3; ++a) {
                            beta[a] = prj_z4c_get(z, prj_z4c_beta_var(a), i, j, k);
                            if (!isfinite(beta[a])) {
                                denom = 1.0e300;
                                break;
                            }
                            denom += (gauge_speed + fabs(beta[a])) / block->dx[a];
                        }
                    }
                    if (!isfinite(denom) || denom <= 0.0) {
                        dt_cell = HUGE_VAL;
                    } else {
                        dt_cell = cfl / (PRJ_CLIGHT * denom);
                    }
                    if (dt_cell < dt_min) {
                        dt_min = dt_cell;
                    }
                }
            }
        }
    }
    return prj_mpi_min_dt(mpi, dt_min);
}

void prj_z4c_init_mesh_flat(prj_mesh *mesh, const prj_mpi *mpi)
{
    int bidx;

    if (!prj_z4c_runtime_enabled(mesh)) {
        return;
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        int stage;

        if (!prj_z4c_local_block(mpi, block)) {
            continue;
        }
        prj_fill(block->z4c, (size_t)PRJ_BLOCK_NSTAGES * (size_t)PRJ_NZ4C *
            (size_t)PRJ_BLOCK_NCELLS_Z4C, 0.0);
        prj_fill(block->z4c_rhs, (size_t)PRJ_BLOCK_NSTAGES * (size_t)PRJ_NZ4C *
            (size_t)PRJ_BLOCK_NCELLS_Z4C, 0.0);
        for (stage = 0; stage < PRJ_BLOCK_NSTAGES; ++stage) {
            double *z = prj_block_z4c_stage(block, stage);
            int i, j, k;

            for (i = -PRJ_NGHOST_Z4C; i < PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C; ++i) {
                for (j = -PRJ_NGHOST_Z4C; j < PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C; ++j) {
                    for (k = -PRJ_NGHOST_Z4C; k < PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C; ++k) {
                        z[Z4CIDX(PRJ_Z4C_CHI, i, j, k)] = 1.0;
                        z[Z4CIDX(PRJ_Z4C_GXX, i, j, k)] = 1.0;
                        z[Z4CIDX(PRJ_Z4C_GYY, i, j, k)] = 1.0;
                        z[Z4CIDX(PRJ_Z4C_GZZ, i, j, k)] = 1.0;
                        z[Z4CIDX(PRJ_Z4C_ALPHA, i, j, k)] = 1.0;
                    }
                }
            }
        }
    }
    mesh->z4c_initialized = 1;
}

/* Linear interpolation of a monotonically-increasing radial table. Clamps to
   the endpoints outside [r[0], r[npts-1]]. */
static double prj_z4c_interp_radial(const double *r_grid, const double *val,
    int npts, double r)
{
    int lo;
    int hi;

    if (npts <= 0) {
        return 0.0;
    }
    if (r <= r_grid[0]) {
        return val[0];
    }
    if (r >= r_grid[npts - 1]) {
        return val[npts - 1];
    }
    lo = 0;
    hi = npts - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;

        if (r_grid[mid] > r) {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    {
        double denom = r_grid[hi] - r_grid[lo];
        double w = denom > 0.0 ? (r - r_grid[lo]) / denom : 0.0;

        return (1.0 - w) * val[lo] + w * val[hi];
    }
}

/* Conformally-flat spherical initial data: fill chi and alpha from a 1D radial
   table (indexed by grid radius measured from x_com), keeping the conformal
   metric flat and the extrinsic curvature, shift, and Theta zero. Mirrors the
   structure of prj_z4c_init_mesh_flat. */
void prj_z4c_init_mesh_spherical(prj_mesh *mesh, const prj_mpi *mpi,
    const double x_com[3], const double *r_grid, const double *chi_grid,
    const double *alpha_grid, int npts)
{
    int bidx;

    if (!prj_z4c_runtime_enabled(mesh)) {
        return;
    }
    if (r_grid == 0 || chi_grid == 0 || alpha_grid == 0 || npts <= 0) {
        prj_z4c_fail("prj_z4c_init_mesh_spherical: missing radial table");
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        int stage;

        if (!prj_z4c_local_block(mpi, block)) {
            continue;
        }
        prj_fill(block->z4c, (size_t)PRJ_BLOCK_NSTAGES * (size_t)PRJ_NZ4C *
            (size_t)PRJ_BLOCK_NCELLS_Z4C, 0.0);
        prj_fill(block->z4c_rhs, (size_t)PRJ_BLOCK_NSTAGES * (size_t)PRJ_NZ4C *
            (size_t)PRJ_BLOCK_NCELLS_Z4C, 0.0);
        for (stage = 0; stage < PRJ_BLOCK_NSTAGES; ++stage) {
            double *z = prj_block_z4c_stage(block, stage);
            int i, j, k;

            for (i = -PRJ_NGHOST_Z4C; i < PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C; ++i) {
                for (j = -PRJ_NGHOST_Z4C; j < PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C; ++j) {
                    for (k = -PRJ_NGHOST_Z4C; k < PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C; ++k) {
                        double dx1 = block->xmin[0] + ((double)i + 0.5) * block->dx[0] - x_com[0];
                        double dx2 = block->xmin[1] + ((double)j + 0.5) * block->dx[1] - x_com[1];
                        double dx3 = block->xmin[2] + ((double)k + 0.5) * block->dx[2] - x_com[2];
                        double r = sqrt(dx1 * dx1 + dx2 * dx2 + dx3 * dx3);
                        double chi = prj_z4c_interp_radial(r_grid, chi_grid, npts, r);
                        double alpha = prj_z4c_interp_radial(r_grid, alpha_grid, npts, r);

                        z[Z4CIDX(PRJ_Z4C_CHI, i, j, k)] = chi;
                        z[Z4CIDX(PRJ_Z4C_GXX, i, j, k)] = 1.0;
                        z[Z4CIDX(PRJ_Z4C_GXY, i, j, k)] = 0.0;
                        z[Z4CIDX(PRJ_Z4C_GXZ, i, j, k)] = 0.0;
                        z[Z4CIDX(PRJ_Z4C_GYY, i, j, k)] = 1.0;
                        z[Z4CIDX(PRJ_Z4C_GYZ, i, j, k)] = 0.0;
                        z[Z4CIDX(PRJ_Z4C_GZZ, i, j, k)] = 1.0;
                        z[Z4CIDX(PRJ_Z4C_KHAT, i, j, k)] = 0.0;
                        z[Z4CIDX(PRJ_Z4C_AXX, i, j, k)] = 0.0;
                        z[Z4CIDX(PRJ_Z4C_AXY, i, j, k)] = 0.0;
                        z[Z4CIDX(PRJ_Z4C_AXZ, i, j, k)] = 0.0;
                        z[Z4CIDX(PRJ_Z4C_AYY, i, j, k)] = 0.0;
                        z[Z4CIDX(PRJ_Z4C_AYZ, i, j, k)] = 0.0;
                        z[Z4CIDX(PRJ_Z4C_AZZ, i, j, k)] = 0.0;
                        z[Z4CIDX(PRJ_Z4C_GAMX, i, j, k)] = 0.0;
                        z[Z4CIDX(PRJ_Z4C_GAMY, i, j, k)] = 0.0;
                        z[Z4CIDX(PRJ_Z4C_GAMZ, i, j, k)] = 0.0;
                        z[Z4CIDX(PRJ_Z4C_THETA, i, j, k)] = 0.0;
                        z[Z4CIDX(PRJ_Z4C_ALPHA, i, j, k)] = alpha;
                        z[Z4CIDX(PRJ_Z4C_BETAX, i, j, k)] = 0.0;
                        z[Z4CIDX(PRJ_Z4C_BETAY, i, j, k)] = 0.0;
                        z[Z4CIDX(PRJ_Z4C_BETAZ, i, j, k)] = 0.0;
                    }
                }
            }
        }
    }
    mesh->z4c_initialized = 1;
}

void prj_z4c_init_punctures(prj_mesh *mesh, const prj_mpi *mpi, int npunctures,
    const double centers_cm[][3], const double masses_cm[],
    const double momenta_cm[][3], double floor_radius_cm)
{
    int bidx;
    double eps2;

    if (!prj_z4c_runtime_enabled(mesh)) {
        return;
    }
    if (npunctures <= 0 || centers_cm == 0 || masses_cm == 0) {
        prj_z4c_fail("prj_z4c_init_punctures: missing puncture data");
    }
    eps2 = floor_radius_cm > 0.0 ? floor_radius_cm * floor_radius_cm : 0.0;

    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        int stage;

        if (!prj_z4c_local_block(mpi, block)) {
            continue;
        }
        prj_fill(block->z4c, (size_t)PRJ_BLOCK_NSTAGES * (size_t)PRJ_NZ4C *
            (size_t)PRJ_BLOCK_NCELLS_Z4C, 0.0);
        prj_fill(block->z4c_rhs, (size_t)PRJ_BLOCK_NSTAGES * (size_t)PRJ_NZ4C *
            (size_t)PRJ_BLOCK_NCELLS_Z4C, 0.0);
        for (stage = 0; stage < PRJ_BLOCK_NSTAGES; ++stage) {
            double *z = prj_block_z4c_stage(block, stage);
            int i, j, k;

            for (i = -PRJ_NGHOST_Z4C; i < PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C; ++i) {
                for (j = -PRJ_NGHOST_Z4C; j < PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C; ++j) {
                    for (k = -PRJ_NGHOST_Z4C; k < PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C; ++k) {
                        double x[3];
                        double psi = 1.0;
                        double Aconf[3][3] = {{0.0, 0.0, 0.0},
                                               {0.0, 0.0, 0.0},
                                               {0.0, 0.0, 0.0}};
                        double psi_to_chi;
                        double psi_m6;
                        int p;
                        int a;
                        int b;

                        x[0] = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
                        x[1] = block->xmin[1] + ((double)j + 0.5) * block->dx[1];
                        x[2] = block->xmin[2] + ((double)k + 0.5) * block->dx[2];

                        for (p = 0; p < npunctures; ++p) {
                            double rvec[3];
                            double r2;
                            double reff;

                            rvec[0] = x[0] - centers_cm[p][0];
                            rvec[1] = x[1] - centers_cm[p][1];
                            rvec[2] = x[2] - centers_cm[p][2];
                            r2 = rvec[0] * rvec[0] + rvec[1] * rvec[1] + rvec[2] * rvec[2];
                            reff = sqrt(r2 + eps2);
                            if (reff < 1.0e-300) {
                                reff = 1.0e-300;
                            }
                            psi += 0.5 * masses_cm[p] / reff;

                            if (momenta_cm != 0) {
                                double n[3] = {0.0, 0.0, 0.0};
                                double Pdotn = 0.0;
                                double r_for_n = sqrt(r2);
                                double coeff;

                                if (r_for_n > 1.0e-300) {
                                    n[0] = rvec[0] / r_for_n;
                                    n[1] = rvec[1] / r_for_n;
                                    n[2] = rvec[2] / r_for_n;
                                }
                                Pdotn = momenta_cm[p][0] * n[0] +
                                    momenta_cm[p][1] * n[1] +
                                    momenta_cm[p][2] * n[2];
                                coeff = 1.5 / (reff * reff);
                                for (a = 0; a < 3; ++a) {
                                    for (b = a; b < 3; ++b) {
                                        double delta_ab = a == b ? 1.0 : 0.0;

                                        Aconf[a][b] += coeff *
                                            (momenta_cm[p][a] * n[b] +
                                             momenta_cm[p][b] * n[a] -
                                             (delta_ab - n[a] * n[b]) * Pdotn);
                                    }
                                }
                            }
                        }

                        psi_to_chi = pow(psi, mesh->z4c_params.chi_psi_power);
                        psi_m6 = pow(psi, -6.0);
                        z[Z4CIDX(PRJ_Z4C_CHI, i, j, k)] = psi_to_chi;
                        z[Z4CIDX(PRJ_Z4C_GXX, i, j, k)] = 1.0;
                        z[Z4CIDX(PRJ_Z4C_GXY, i, j, k)] = 0.0;
                        z[Z4CIDX(PRJ_Z4C_GXZ, i, j, k)] = 0.0;
                        z[Z4CIDX(PRJ_Z4C_GYY, i, j, k)] = 1.0;
                        z[Z4CIDX(PRJ_Z4C_GYZ, i, j, k)] = 0.0;
                        z[Z4CIDX(PRJ_Z4C_GZZ, i, j, k)] = 1.0;
                        z[Z4CIDX(PRJ_Z4C_KHAT, i, j, k)] = 0.0;
                        z[Z4CIDX(PRJ_Z4C_AXX, i, j, k)] = psi_m6 * Aconf[0][0];
                        z[Z4CIDX(PRJ_Z4C_AXY, i, j, k)] = psi_m6 * Aconf[0][1];
                        z[Z4CIDX(PRJ_Z4C_AXZ, i, j, k)] = psi_m6 * Aconf[0][2];
                        z[Z4CIDX(PRJ_Z4C_AYY, i, j, k)] = psi_m6 * Aconf[1][1];
                        z[Z4CIDX(PRJ_Z4C_AYZ, i, j, k)] = psi_m6 * Aconf[1][2];
                        z[Z4CIDX(PRJ_Z4C_AZZ, i, j, k)] = psi_m6 * Aconf[2][2];
                        z[Z4CIDX(PRJ_Z4C_GAMX, i, j, k)] = 0.0;
                        z[Z4CIDX(PRJ_Z4C_GAMY, i, j, k)] = 0.0;
                        z[Z4CIDX(PRJ_Z4C_GAMZ, i, j, k)] = 0.0;
                        z[Z4CIDX(PRJ_Z4C_THETA, i, j, k)] = 0.0;
                        z[Z4CIDX(PRJ_Z4C_ALPHA, i, j, k)] = pow(psi, -2.0);
                        z[Z4CIDX(PRJ_Z4C_BETAX, i, j, k)] = 0.0;
                        z[Z4CIDX(PRJ_Z4C_BETAY, i, j, k)] = 0.0;
                        z[Z4CIDX(PRJ_Z4C_BETAZ, i, j, k)] = 0.0;
                    }
                }
            }
        }
    }
    mesh->z4c_initialized = 1;
}

static int prj_z4c_parity(int var, int normal_dir)
{
    if (normal_dir == X1DIR) {
        return (var == PRJ_Z4C_GXY || var == PRJ_Z4C_GXZ ||
            var == PRJ_Z4C_AXY || var == PRJ_Z4C_AXZ ||
            var == PRJ_Z4C_GAMX || var == PRJ_Z4C_BETAX) ? -1 : 1;
    }
    if (normal_dir == X2DIR) {
        return (var == PRJ_Z4C_GXY || var == PRJ_Z4C_GYZ ||
            var == PRJ_Z4C_AXY || var == PRJ_Z4C_AYZ ||
            var == PRJ_Z4C_GAMY || var == PRJ_Z4C_BETAY) ? -1 : 1;
    }
    return (var == PRJ_Z4C_GXZ || var == PRJ_Z4C_GYZ ||
        var == PRJ_Z4C_AXZ || var == PRJ_Z4C_AYZ ||
        var == PRJ_Z4C_GAMZ || var == PRJ_Z4C_BETAZ) ? -1 : 1;
}

static double prj_z4c_extrapolate(const double *z, int var, int i, int j, int k,
    int offi, int offj, int offk, int delta, int order)
{
    double f0 = prj_z4c_get(z, var, i, j, k);
    double f1 = prj_z4c_get(z, var, i + offi, j + offj, k + offk);

    if (order <= 2) {
        return f0 + (double)delta * (f0 - f1);
    } else {
        double f2 = prj_z4c_get(z, var, i + 2 * offi, j + 2 * offj, k + 2 * offk);

        if (order == 3) {
            return 0.5 * (f0 * (1.0 + delta) * (2.0 + delta) +
                (double)delta * (f2 + (double)delta * f2 - 2.0 * f1 * (2.0 + delta)));
        } else {
            double f3 = prj_z4c_get(z, var, i + 3 * offi, j + 3 * offj, k + 3 * offk);

            return (-3.0 * f1 * (double)delta * (2.0 + delta) * (3.0 + delta) +
                f0 * (1.0 + delta) * (2.0 + delta) * (3.0 + delta) +
                (double)delta * (1.0 + delta) *
                (-f3 * (2.0 + delta) + 3.0 * f2 * (3.0 + delta))) / 6.0;
        }
    }
}

static int prj_z4c_axis_bc_type(const prj_bc *bc, int axis, int outer)
{
    if (axis == X1DIR) {
        return outer ? bc->bc_x1_outer : bc->bc_x1_inner;
    }
    if (axis == X2DIR) {
        return outer ? bc->bc_x2_outer : bc->bc_x2_inner;
    }
    return outer ? bc->bc_x3_outer : bc->bc_x3_inner;
}

static double prj_z4c_axis_domain_face(const prj_mesh *mesh, int axis, int outer)
{
    if (axis == X1DIR) {
        return outer ? mesh->coord.x1max : mesh->coord.x1min;
    }
    if (axis == X2DIR) {
        return outer ? mesh->coord.x2max : mesh->coord.x2min;
    }
    return outer ? mesh->coord.x3max : mesh->coord.x3min;
}

static void prj_z4c_apply_physical_axis(const prj_mesh *mesh, const prj_bc *bc,
    prj_block *block, double *z, int axis, int outer)
{
    double block_face = outer ? block->xmax[axis] : block->xmin[axis];
    double domain_face = prj_z4c_axis_domain_face(mesh, axis, outer);
    double tol = 1.0e-8 * PRJ_MAX(fabs(block->dx[axis]), 1.0);
    int bc_type;
    int order;
    int p;

    if (fabs(block_face - domain_face) > tol) {
        return;
    }
    bc_type = prj_z4c_axis_bc_type(bc, axis, outer);
    order = mesh->z4c_extrap_order;
    for (p = 0; p < PRJ_NGHOST_Z4C; ++p) {
        int ghost = outer ? PRJ_BLOCK_SIZE + p : -p - 1;
        int mirror = outer ? PRJ_BLOCK_SIZE - p - 1 : p;
        int base = outer ? PRJ_BLOCK_SIZE - 1 : 0;
        int off = outer ? -1 : 1;
        int i, j, k;

        for (i = -PRJ_NGHOST_Z4C; i < PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C; ++i) {
            for (j = -PRJ_NGHOST_Z4C; j < PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C; ++j) {
                for (k = -PRJ_NGHOST_Z4C; k < PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C; ++k) {
                    int dst_idx[3] = {i, j, k};
                    int src_idx[3] = {i, j, k};
                    int var;

                    dst_idx[axis] = ghost;
                    src_idx[axis] = mirror;
                    for (var = 0; var < PRJ_NZ4C; ++var) {
                        double value;

                        if (bc_type == PRJ_BC_REFLECT) {
                            value = (double)prj_z4c_parity(var, axis) *
                                prj_z4c_get(z, var, src_idx[0], src_idx[1], src_idx[2]);
                        } else {
                            int base_idx[3] = {i, j, k};
                            int off_idx[3] = {0, 0, 0};

                            base_idx[axis] = base;
                            off_idx[axis] = off;
                            value = prj_z4c_extrapolate(z, var,
                                base_idx[0], base_idx[1], base_idx[2],
                                off_idx[0], off_idx[1], off_idx[2], p + 1, order);
                        }
                        prj_z4c_set(z, var, dst_idx[0], dst_idx[1], dst_idx[2], value);
                    }
                }
            }
        }
    }
}

static void prj_z4c_apply_physical_bcs(prj_mesh *mesh, const prj_mpi *mpi,
    const prj_bc *bc, int stage)
{
    int bidx;

    if (mesh == 0 || bc == 0) {
        return;
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        double *z;
        int axis;

        if (!prj_z4c_local_block(mpi, block)) {
            continue;
        }
        z = prj_block_z4c_stage(block, stage);
        for (axis = 0; axis < 3; ++axis) {
            prj_z4c_apply_physical_axis(mesh, bc, block, z, axis, 0);
            prj_z4c_apply_physical_axis(mesh, bc, block, z, axis, 1);
        }
    }
}

static double prj_z4c_restrict_value(const double *src, int var, int i, int j, int k)
{
    double sum = 0.0;
    int di, dj, dk;

    for (di = 0; di < 2; ++di) {
        for (dj = 0; dj < 2; ++dj) {
            for (dk = 0; dk < 2; ++dk) {
                sum += prj_z4c_get(src, var, 2 * i + di, 2 * j + dj, 2 * k + dk);
            }
        }
    }
    return 0.125 * sum;
}

static double prj_z4c_prolong_value(const double *src, int var, int i, int j, int k,
    const double target[3])
{
    double c = prj_z4c_get(src, var, i, j, k);
    double sx = 0.5 * (prj_z4c_get(src, var, i + 1, j, k) -
        prj_z4c_get(src, var, i - 1, j, k));
    double sy = 0.5 * (prj_z4c_get(src, var, i, j + 1, k) -
        prj_z4c_get(src, var, i, j - 1, k));
    double sz = 0.5 * (prj_z4c_get(src, var, i, j, k + 1) -
        prj_z4c_get(src, var, i, j, k - 1));

    return c + target[0] * sx + target[1] * sy + target[2] * sz;
}

static void prj_z4c_clear_block_aux(prj_block *block)
{
    if (block == 0) {
        return;
    }
    if (block->z4c_rhs != 0) {
        prj_fill(block->z4c_rhs, (size_t)PRJ_BLOCK_NSTAGES * (size_t)PRJ_NZ4C *
            (size_t)PRJ_BLOCK_NCELLS_Z4C, 0.0);
    }
}

static void prj_z4c_amr_prolong_axis_weights(int parent_index, int fine_odd,
    int *base, int *n, double w[5])
{
#if PRJ_NGHOST_Z4C == 2
    static const double pro[3] = {0.15625, 0.9375, -0.09375};
    int q;

    *n = 3;
    *base = parent_index - 1;
    for (q = 0; q < *n; ++q) {
        w[q] = fine_odd ? pro[*n - 1 - q] : pro[q];
    }
#elif PRJ_NGHOST_Z4C == 4
    static const double pro[5] = {
        -0.02197265625, 0.205078125, 0.9228515625, -0.123046875, 0.01708984375
    };
    int q;

    *n = 5;
    *base = parent_index - 2;
    for (q = 0; q < *n; ++q) {
        w[q] = fine_odd ? pro[*n - 1 - q] : pro[q];
    }
#else
    (void)parent_index;
    (void)fine_odd;
    *base = 0;
    *n = 0;
    memset(w, 0, 5U * sizeof(*w));
#endif
}

static double prj_z4c_amr_prolong_value(const double *src, int var,
    int i, int j, int k, int child_oct)
{
#if PRJ_NGHOST_Z4C == 2 || PRJ_NGHOST_Z4C == 4
    int gi = ((child_oct & 1) ? PRJ_BLOCK_SIZE : 0) + i;
    int gj = ((child_oct & 2) ? PRJ_BLOCK_SIZE : 0) + j;
    int gk = ((child_oct & 4) ? PRJ_BLOCK_SIZE : 0) + k;
    int pi = gi / 2;
    int pj = gj / 2;
    int pk = gk / 2;
    int ibase, jbase, kbase;
    int ni, nj, nk;
    double wi[5], wj[5], wk[5];
    double out = 0.0;
    int ii, jj, kk;

    prj_z4c_amr_prolong_axis_weights(pi, gi & 1, &ibase, &ni, wi);
    prj_z4c_amr_prolong_axis_weights(pj, gj & 1, &jbase, &nj, wj);
    prj_z4c_amr_prolong_axis_weights(pk, gk & 1, &kbase, &nk, wk);
    for (ii = 0; ii < ni; ++ii) {
        for (jj = 0; jj < nj; ++jj) {
            for (kk = 0; kk < nk; ++kk) {
                out += wi[ii] * wj[jj] * wk[kk] *
                    prj_z4c_get(src, var, ibase + ii, jbase + jj, kbase + kk);
            }
        }
    }
    return out;
#else
    double target[3] = {
        (i % 2 == 0) ? -0.25 : 0.25,
        (j % 2 == 0) ? -0.25 : 0.25,
        (k % 2 == 0) ? -0.25 : 0.25
    };
    int pi = i / 2 + ((child_oct & 1) ? PRJ_BLOCK_SIZE / 2 : 0);
    int pj = j / 2 + ((child_oct & 2) ? PRJ_BLOCK_SIZE / 2 : 0);
    int pk = k / 2 + ((child_oct & 4) ? PRJ_BLOCK_SIZE / 2 : 0);

    return prj_z4c_prolong_value(src, var, pi, pj, pk, target);
#endif
}

static void prj_z4c_amr_restrict_axis_weights(int fine_even, int *base,
    int *n, double w[5])
{
#if PRJ_NGHOST_Z4C == 2
    static const double res[3] = {0.375, 0.75, -0.125};
    int left = fine_even < PRJ_BLOCK_SIZE / 2;
    int q;

    *n = 3;
    *base = left ? fine_even : fine_even - 1;
    for (q = 0; q < *n; ++q) {
        w[q] = left ? res[q] : res[*n - 1 - q];
    }
#elif PRJ_NGHOST_Z4C == 4
    static const double res[5] = {
        -0.0390625, 0.46875, 0.703125, -0.15625, 0.0234375
    };
    static const double res_edge[5] = {
        0.2734375, 1.09375, -0.546875, 0.21875, -0.0390625
    };
    const double *rw = (fine_even == 0 || fine_even == PRJ_BLOCK_SIZE - 2) ?
        res_edge : res;
    int left = fine_even < PRJ_BLOCK_SIZE / 2;
    int q;

    *n = 5;
    *base = left ? fine_even - 1 : fine_even - 2;
    if (fine_even == 0) {
        *base += 1;
    } else if (fine_even == PRJ_BLOCK_SIZE - 2) {
        *base -= 1;
    }
    for (q = 0; q < *n; ++q) {
        w[q] = left ? rw[q] : rw[*n - 1 - q];
    }
#else
    (void)fine_even;
    *base = 0;
    *n = 0;
    memset(w, 0, 5U * sizeof(*w));
#endif
}

static double prj_z4c_amr_restrict_value(const double *src, int var, int i, int j, int k)
{
#if PRJ_NGHOST_Z4C == 2 || PRJ_NGHOST_Z4C == 4
    int fi = 2 * i;
    int fj = 2 * j;
    int fk = 2 * k;
    int ibase, jbase, kbase;
    int ni, nj, nk;
    double wi[5], wj[5], wk[5];
    double out = 0.0;
    int ii, jj, kk;

    prj_z4c_amr_restrict_axis_weights(fi, &ibase, &ni, wi);
    prj_z4c_amr_restrict_axis_weights(fj, &jbase, &nj, wj);
    prj_z4c_amr_restrict_axis_weights(fk, &kbase, &nk, wk);
    for (ii = 0; ii < ni; ++ii) {
        for (jj = 0; jj < nj; ++jj) {
            for (kk = 0; kk < nk; ++kk) {
                out += wi[ii] * wj[jj] * wk[kk] *
                    prj_z4c_get(src, var, ibase + ii, jbase + jj, kbase + kk);
            }
        }
    }
    return out;
#else
    return prj_z4c_restrict_value(src, var, i, j, k);
#endif
}

static int prj_z4c_fill_kind_from_rel_level(int rel_level)
{
    if (rel_level == 0) {
        return PRJ_BOUNDARY_FILL_SAME_LEVEL;
    }
    return rel_level < 0 ? PRJ_BOUNDARY_FILL_RESTRICTION :
        PRJ_BOUNDARY_FILL_PROLONGATION;
}

static double prj_z4c_sample_slot_value(const double *src, const prj_neighbor *slot,
    int var, int i, int j, int k)
{
    if (slot->rel_level == 0) {
        return prj_z4c_get(src, var,
            i + slot->send_loc_start_z4c[0],
            j + slot->send_loc_start_z4c[1],
            k + slot->send_loc_start_z4c[2]);
    }
    if (slot->rel_level < 0) {
        return prj_z4c_restrict_value(src, var,
            i + slot->send_loc_start_z4c[0] / 2,
            j + slot->send_loc_start_z4c[1] / 2,
            k + slot->send_loc_start_z4c[2] / 2);
    } else {
        double target[3];
        int ai = i + slot->recv_loc_start_z4c[0];
        int aj = j + slot->recv_loc_start_z4c[1];
        int ak = k + slot->recv_loc_start_z4c[2];

        target[0] = (ai % 2 == 0) ? -0.25 : 0.25;
        target[1] = (aj % 2 == 0) ? -0.25 : 0.25;
        target[2] = (ak % 2 == 0) ? -0.25 : 0.25;
        return prj_z4c_prolong_value(src, var,
            i / 2 + slot->send_loc_start_z4c[0],
            j / 2 + slot->send_loc_start_z4c[1],
            k / 2 + slot->send_loc_start_z4c[2], target);
    }
}

static void prj_z4c_local_send(prj_mesh *mesh, const prj_mpi *mpi, int stage, int fill_kind)
{
    int bidx;

    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        double *src;
        int n;

        if (!prj_z4c_local_block(mpi, block)) {
            continue;
        }
        src = prj_block_z4c_stage(block, stage);
        for (n = 0; n < 56; ++n) {
            const prj_neighbor *slot = &block->slot[n];
            int nid = slot->id;
            prj_block *neighbor;
            double *dst;
            int i, j, k;

            if (nid < 0 || nid >= mesh->nblocks) {
                continue;
            }
            neighbor = &mesh->blocks[nid];
            if (!prj_z4c_local_block(mpi, neighbor) || neighbor->rank != block->rank ||
                prj_z4c_fill_kind_from_rel_level(slot->rel_level) != fill_kind) {
                continue;
            }
            dst = prj_block_z4c_stage(neighbor, stage);
            for (i = 0; i < slot->recv_loc_end_z4c[0] - slot->recv_loc_start_z4c[0]; ++i) {
                for (j = 0; j < slot->recv_loc_end_z4c[1] - slot->recv_loc_start_z4c[1]; ++j) {
                    for (k = 0; k < slot->recv_loc_end_z4c[2] - slot->recv_loc_start_z4c[2]; ++k) {
                        int di = i + slot->recv_loc_start_z4c[0];
                        int dj = j + slot->recv_loc_start_z4c[1];
                        int dk = k + slot->recv_loc_start_z4c[2];
                        int var;

                        if (!prj_z4c_cell_in_storage(di, dj, dk)) {
                            continue;
                        }
                        for (var = 0; var < PRJ_NZ4C; ++var) {
                            prj_z4c_set(dst, var, di, dj, dk,
                                prj_z4c_sample_slot_value(src, slot, var, i, j, k));
                        }
                    }
                }
            }
        }
    }
}

static int prj_z4c_decode_cell_index(int code, int *i, int *j, int *k)
{
    *k = code % PRJ_BS_Z4C - PRJ_NGHOST_Z4C;
    code /= PRJ_BS_Z4C;
    *j = code % PRJ_BS_Z4C - PRJ_NGHOST_Z4C;
    code /= PRJ_BS_Z4C;
    *i = code - PRJ_NGHOST_Z4C;
    return 0;
}

static int prj_z4c_buffer_record_count(const int *sizes)
{
    int count = 0;

    if (sizes == 0) {
        return 0;
    }
    while (sizes[count] >= 0) {
        count += 1;
    }
    return count;
}

static int prj_z4c_buffer_record_total(const int *sizes)
{
    int total = 0;
    int n = prj_z4c_buffer_record_count(sizes);
    int i;

    for (i = 0; i < n; ++i) {
        total += sizes[i];
    }
    return total;
}

static void prj_z4c_mpi_exchange(prj_mesh *mesh, prj_mpi *mpi, int stage, int fill_kind)
{
#if defined(PRJ_ENABLE_MPI)
    MPI_Request *requests;
    double **send_values;
    double **recv_values;
    int request_count = 0;
    int nb;

    if (mesh == 0 || mpi == 0 || mpi->totrank <= 1 || mpi->neighbor_number == 0) {
        return;
    }
    requests = (MPI_Request *)mpi->request_buffer;
    if (requests == 0 || mpi->request_capacity < 2 * mpi->neighbor_number) {
        prj_z4c_fail("prj_z4c_mpi_exchange: missing MPI request buffer");
    }
    send_values = (double **)prj_calloc((size_t)mpi->neighbor_number, sizeof(*send_values));
    recv_values = (double **)prj_calloc((size_t)mpi->neighbor_number, sizeof(*recv_values));
    if (send_values == 0 || recv_values == 0) {
        prj_z4c_fail("prj_z4c_mpi_exchange: allocation failed");
    }

    for (nb = 0; nb < mpi->neighbor_number; ++nb) {
        prj_mpi_buffer *buffer = &mpi->neighbor_buffer[nb];
        int send_records = buffer->cell_send_count_z4c_by_kind[fill_kind];
        int recv_records = buffer->cell_recv_count_z4c_by_kind[fill_kind];
        int pos = 0;
        int bidx;

        if (send_records > 0) {
            send_values[nb] = (double *)prj_malloc((size_t)send_records *
                (size_t)PRJ_NZ4C * sizeof(**send_values));
        }
        if (recv_records > 0) {
            recv_values[nb] = (double *)prj_malloc((size_t)recv_records *
                (size_t)PRJ_NZ4C * sizeof(**recv_values));
        }
        for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
            prj_block *block = &mesh->blocks[bidx];
            double *src;
            int n;

            if (!prj_z4c_local_block(mpi, block)) {
                continue;
            }
            src = prj_block_z4c_stage(block, stage);
            for (n = 0; n < 56; ++n) {
                const prj_neighbor *slot = &block->slot[n];
                int nid = slot->id;
                int i, j, k;

                if (nid < 0 || nid >= mesh->nblocks ||
                    mesh->blocks[nid].rank != buffer->receiver_rank ||
                    prj_z4c_fill_kind_from_rel_level(slot->rel_level) != fill_kind) {
                    continue;
                }
                for (i = 0; i < slot->recv_loc_end_z4c[0] - slot->recv_loc_start_z4c[0]; ++i) {
                    for (j = 0; j < slot->recv_loc_end_z4c[1] - slot->recv_loc_start_z4c[1]; ++j) {
                        for (k = 0; k < slot->recv_loc_end_z4c[2] - slot->recv_loc_start_z4c[2]; ++k) {
                            int var;

                            for (var = 0; var < PRJ_NZ4C; ++var) {
                                if (pos >= send_records * PRJ_NZ4C) {
                                    prj_z4c_fail("prj_z4c_mpi_exchange: send buffer overflow");
                                }
                                send_values[nb][pos++] =
                                    prj_z4c_sample_slot_value(src, slot, var, i, j, k);
                            }
                        }
                    }
                }
            }
        }
        if (pos != send_records * PRJ_NZ4C) {
            prj_z4c_fail("prj_z4c_mpi_exchange: send count mismatch");
        }
        MPI_Irecv(recv_values[nb], recv_records * PRJ_NZ4C, MPI_DOUBLE,
            buffer->receiver_rank, PRJ_Z4C_MPI_TAG, MPI_COMM_WORLD,
            &requests[request_count++]);
        MPI_Isend(send_values[nb], send_records * PRJ_NZ4C, MPI_DOUBLE,
            buffer->receiver_rank, PRJ_Z4C_MPI_TAG, MPI_COMM_WORLD,
            &requests[request_count++]);
    }
    if (request_count > 0) {
        MPI_Waitall(request_count, requests, MPI_STATUSES_IGNORE);
    }
    for (nb = 0; nb < mpi->neighbor_number; ++nb) {
        prj_mpi_buffer *buffer = &mpi->neighbor_buffer[nb];
        int total = prj_z4c_buffer_record_total(buffer->cell_data_size_recv_z4c);
        int pos = 0;
        int r;

        for (r = 0; r < total; ++r) {
            int block_id = buffer->cell_data_idx_recv_z4c[0][r];
            int code = buffer->cell_data_idx_recv_z4c[1][r];
            int sample_kind = buffer->cell_data_idx_recv_z4c[2][r];
            int i, j, k, var;
            prj_block *block;
            double *dst;

            if (sample_kind != fill_kind) {
                continue;
            }
            prj_z4c_decode_cell_index(code, &i, &j, &k);
            if (!prj_z4c_cell_in_storage(i, j, k)) {
                pos += PRJ_NZ4C;
                continue;
            }
            if (block_id < 0 || block_id >= mesh->nblocks) {
                pos += PRJ_NZ4C;
                continue;
            }
            block = &mesh->blocks[block_id];
            if (!prj_z4c_local_block(mpi, block)) {
                pos += PRJ_NZ4C;
                continue;
            }
            dst = prj_block_z4c_stage(block, stage);
            for (var = 0; var < PRJ_NZ4C; ++var) {
                prj_z4c_set(dst, var, i, j, k, recv_values[nb][pos++]);
            }
        }
        free(send_values[nb]);
        free(recv_values[nb]);
    }
    free(send_values);
    free(recv_values);
#else
    (void)mesh;
    (void)mpi;
    (void)stage;
    (void)fill_kind;
#endif
}

void prj_z4c_fill_ghosts(prj_mesh *mesh, prj_mpi *mpi, const prj_bc *bc, int stage)
{
    int pass;
    int fill_kinds[3] = {
        PRJ_BOUNDARY_FILL_SAME_LEVEL,
        PRJ_BOUNDARY_FILL_RESTRICTION,
        PRJ_BOUNDARY_FILL_PROLONGATION
    };

    if (!prj_z4c_runtime_enabled(mesh)) {
        return;
    }
    if (stage < 0 || stage >= PRJ_BLOCK_NSTAGES) {
        prj_z4c_fail("prj_z4c_fill_ghosts: invalid stage");
    }
    if (bc != 0) {
        prj_z4c_apply_physical_bcs(mesh, mpi, bc, stage);
    }
    for (pass = 0; pass < 3; ++pass) {
        prj_z4c_local_send(mesh, mpi, stage, fill_kinds[pass]);
        prj_z4c_mpi_exchange(mesh, mpi, stage, fill_kinds[pass]);
    }
    if (bc != 0) {
        prj_z4c_apply_physical_bcs(mesh, mpi, bc, stage);
    }
}

#if PRJ_USE_RADIATION_M1
static void prj_z4c_lower_vec(const double g[3][3], const double v_con[3],
    double v_cov[3])
{
    int a, b;

    for (a = 0; a < 3; ++a) {
        v_cov[a] = 0.0;
        for (b = 0; b < 3; ++b) {
            v_cov[a] += g[a][b] * v_con[b];
        }
    }
}
#endif

static void prj_z4c_build_tmunu_cell_cgs(const prj_rad *rad, const prj_block *block,
    const prj_z4c_hydro_geom *geom, const double *W_mhd, const double *W_rad,
    int i, int j, int k, double t_cgs[PRJ_NTMUNU])
{
    double stress[3][3];
    double mom[3] = {0.0, 0.0, 0.0};
    int a, b;

    (void)rad;
    if (geom == 0) {
        prj_z4c_fail("prj_z4c_build_tmunu_cell_cgs: missing geometry");
    }
    memset(t_cgs, 0, (size_t)PRJ_NTMUNU * sizeof(*t_cgs));
    memset(stress, 0, sizeof(stress));

    if (W_mhd != 0) {
        prj_eos_gr_geom egeom;
        prj_eos_grmhd_state state;
        double Wc[PRJ_NVAR_PRIM];
        double pressure;
        int status;
        int v;

        pressure = block != 0 && block->eosvar != 0 ?
            block->eosvar[EIDX(PRJ_EOSVAR_PRESSURE, i, j, k)] : 0.0;
        for (a = 0; a < 3; ++a) {
            for (b = 0; b < 3; ++b) {
                egeom.gamma[a][b] = geom->gamma[a][b];
            }
        }
        for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
            Wc[v] = W_mhd[WIDX(v, i, j, k)];
        }
        if (Wc[PRJ_PRIM_RHO] != 0.0 || pressure != 0.0) {
            status = prj_eos_grmhd_state_from_prim(0, &egeom, Wc, pressure,
                &state, PRJ_EOS_CTX_MAIN);
            if (status != PRJ_EOS_GR_OK) {
                prj_z4c_fail("prj_z4c_build_tmunu_cell_cgs: invalid GRMHD state");
            }
            t_cgs[PRJ_TMUNU_E] += state.E;
            for (a = 0; a < 3; ++a) {
                mom[a] += state.S_cov[a];
                for (b = 0; b < 3; ++b) {
                    stress[a][b] += state.stress_cov[a][b];
                }
            }
        }
    }

#if PRJ_USE_RADIATION_M1
    if (rad != 0 && W_rad != 0) {
        int field, group;

        for (field = 0; field < PRJ_NRAD; ++field) {
            for (group = 0; group < PRJ_NEGROUP; ++group) {
                double Erad = W_rad[WIDX(PRJ_RAD_PRIM_E(field, group), i, j, k)] *
                    RAD_SCALE;
                double Frad[3] = {
                    W_rad[WIDX(PRJ_RAD_PRIM_F1(field, group), i, j, k)] * RAD_SCALE,
                    W_rad[WIDX(PRJ_RAD_PRIM_F2(field, group), i, j, k)] * RAD_SCALE,
                    W_rad[WIDX(PRJ_RAD_PRIM_F3(field, group), i, j, k)] * RAD_SCALE
                };
                double Frad_cov[3];
                double P[3][3];
                double Pcov[3][3];

                prj_rad_m1_pressure(rad, Erad, Frad[0], Frad[1], Frad[2], P);
                prj_z4c_lower_vec(geom->gamma, Frad, Frad_cov);
                memset(Pcov, 0, sizeof(Pcov));
                for (a = 0; a < 3; ++a) {
                    for (b = 0; b < 3; ++b) {
                        int c, d;

                        for (c = 0; c < 3; ++c) {
                            for (d = 0; d < 3; ++d) {
                                Pcov[a][b] += geom->gamma[a][c] *
                                    geom->gamma[b][d] * P[c][d];
                            }
                        }
                    }
                }
                t_cgs[PRJ_TMUNU_E] += Erad;
                for (a = 0; a < 3; ++a) {
                    mom[a] += Frad_cov[a] / PRJ_CLIGHT;
                    for (b = 0; b < 3; ++b) {
                        stress[a][b] += Pcov[a][b];
                    }
                }
            }
        }
    }
#else
    (void)W_rad;
#endif

    for (a = 0; a < 3; ++a) {
        t_cgs[prj_z4c_tmunu_mom_var(a)] = mom[a];
        for (b = a; b < 3; ++b) {
            t_cgs[prj_z4c_tmunu_stress_var(a, b)] = stress[a][b];
        }
    }
}

static void prj_z4c_build_tmunu_cell_geo(const prj_rad *rad, const prj_block *block,
    const prj_z4c_hydro_geom *geom, const double *W_mhd, const double *W_rad,
    int i, int j, int k, double t_geo[PRJ_NTMUNU])
{
    double t_cgs[PRJ_NTMUNU];
    int v;

    prj_z4c_build_tmunu_cell_cgs(rad, block, geom, W_mhd, W_rad,
        i, j, k, t_cgs);
    for (v = 0; v < PRJ_NTMUNU; ++v) {
        t_geo[v] = prj_z4c_tmunu_cgs_to_geo(t_cgs[v]);
    }
}

static void prj_z4c_compute_rhs_cell(const prj_mesh *mesh, const prj_block *block,
    const double *z, const double *tmunu, double *rhs, int i, int j, int k, double tau_cm)
{
    const prj_z4c_params *opt = &mesh->z4c_params;
    double idx[3] = {1.0 / block->dx[0], 1.0 / block->dx[1], 1.0 / block->dx[2]};
    double g[3][3], gu[3][3], A[3][3], Auu[3][3], AA_dd[3][3];
    double R_dd[3][3], Rphi_dd[3][3], Ddalpha_dd[3][3], Ddphi_dd[3][3];
    double Gamma_ddd[3][3][3], Gamma_udd[3][3][3];
    double dg[3][3][3], ddg[3][3][3][3];
    double dalpha[3], dchi[3], dphi[3], dKhat[3], dTheta[3];
    double ddalpha[3][3], ddchi[3][3], dbeta[3][3], dGam[3][3];
    double ddbeta[3][3][3];
    double Gamma_u[3], DA_u[3], LGam[3], Lbeta[3], Lg[3][3], LA[3][3];
    double T_S_dd[3][3], T_S_d[3], T_E = 0.0, T_S = 0.0;
    double Lalpha = 0.0, Lchi = 0.0, LKhat = 0.0, LTheta = 0.0;
    double detg, chi, chi_guarded, oopsi4, K, AA = 0.0, R = 0.0, Ht, Ddalpha = 0.0;
    double div_beta = 0.0;
    double alpha;
    int a, b, c, d, e, var;

    memset(g, 0, sizeof(g));
    memset(gu, 0, sizeof(gu));
    memset(A, 0, sizeof(A));
    memset(Auu, 0, sizeof(Auu));
    memset(AA_dd, 0, sizeof(AA_dd));
    memset(R_dd, 0, sizeof(R_dd));
    memset(Rphi_dd, 0, sizeof(Rphi_dd));
    memset(Ddalpha_dd, 0, sizeof(Ddalpha_dd));
    memset(Ddphi_dd, 0, sizeof(Ddphi_dd));
    memset(Gamma_ddd, 0, sizeof(Gamma_ddd));
    memset(Gamma_udd, 0, sizeof(Gamma_udd));
    memset(dg, 0, sizeof(dg));
    memset(ddg, 0, sizeof(ddg));
    memset(dalpha, 0, sizeof(dalpha));
    memset(dchi, 0, sizeof(dchi));
    memset(dphi, 0, sizeof(dphi));
    memset(dKhat, 0, sizeof(dKhat));
    memset(dTheta, 0, sizeof(dTheta));
    memset(ddalpha, 0, sizeof(ddalpha));
    memset(ddchi, 0, sizeof(ddchi));
    memset(dbeta, 0, sizeof(dbeta));
    memset(dGam, 0, sizeof(dGam));
    memset(ddbeta, 0, sizeof(ddbeta));
    memset(Gamma_u, 0, sizeof(Gamma_u));
    memset(DA_u, 0, sizeof(DA_u));
    memset(LGam, 0, sizeof(LGam));
    memset(Lbeta, 0, sizeof(Lbeta));
    memset(Lg, 0, sizeof(Lg));
    memset(LA, 0, sizeof(LA));
    memset(T_S_dd, 0, sizeof(T_S_dd));
    memset(T_S_d, 0, sizeof(T_S_d));

    prj_z4c_load_metric(z, i, j, k, g);
    prj_z4c_load_A(z, i, j, k, A);
    prj_z4c_inv3(g, gu, &detg);
    (void)detg;
    chi = prj_z4c_get(z, PRJ_Z4C_CHI, i, j, k);
    chi_guarded = chi > opt->chi_div_floor ? chi : opt->chi_div_floor;
    if (chi_guarded <= opt->chi_min_floor || !isfinite(chi_guarded)) {
        chi_guarded = opt->chi_min_floor;
    }
    oopsi4 = pow(chi_guarded, -4.0 / opt->chi_psi_power);
    K = prj_z4c_get(z, PRJ_Z4C_KHAT, i, j, k) +
        2.0 * prj_z4c_get(z, PRJ_Z4C_THETA, i, j, k);
    alpha = prj_z4c_get(z, PRJ_Z4C_ALPHA, i, j, k);
    if (tmunu != 0) {
        T_E = tmunu[PRJ_TMUNU_E];
        for (a = 0; a < 3; ++a) {
            T_S_d[a] = tmunu[prj_z4c_tmunu_mom_var(a)];
            for (b = 0; b < 3; ++b) {
                T_S_dd[a][b] = tmunu[prj_z4c_tmunu_stress_var(a, b)];
                T_S += oopsi4 * gu[a][b] * T_S_dd[a][b];
            }
        }
    }

    for (a = 0; a < 3; ++a) {
        dalpha[a] = prj_z4c_Dx(z, PRJ_Z4C_ALPHA, a, idx, i, j, k);
        dchi[a] = prj_z4c_Dx(z, PRJ_Z4C_CHI, a, idx, i, j, k);
        dKhat[a] = prj_z4c_Dx(z, PRJ_Z4C_KHAT, a, idx, i, j, k);
        dTheta[a] = prj_z4c_Dx(z, PRJ_Z4C_THETA, a, idx, i, j, k);
        Lalpha += prj_z4c_Lx(z, prj_z4c_beta_var(a), PRJ_Z4C_ALPHA, a, idx, i, j, k);
        Lchi += prj_z4c_Lx(z, prj_z4c_beta_var(a), PRJ_Z4C_CHI, a, idx, i, j, k);
        LKhat += prj_z4c_Lx(z, prj_z4c_beta_var(a), PRJ_Z4C_KHAT, a, idx, i, j, k);
        LTheta += prj_z4c_Lx(z, prj_z4c_beta_var(a), PRJ_Z4C_THETA, a, idx, i, j, k);
    }
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            dbeta[b][a] = prj_z4c_Dx(z, prj_z4c_beta_var(a), b, idx, i, j, k);
            dGam[b][a] = prj_z4c_Dx(z, prj_z4c_Gam_var(a), b, idx, i, j, k);
            Lbeta[b] += prj_z4c_Lx(z, prj_z4c_beta_var(a), prj_z4c_beta_var(b),
                a, idx, i, j, k);
            LGam[b] += prj_z4c_Lx(z, prj_z4c_beta_var(a), prj_z4c_Gam_var(b),
                a, idx, i, j, k);
        }
    }
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            for (c = 0; c < 3; ++c) {
                dg[c][a][b] = prj_z4c_Dx(z, prj_z4c_g_var(a, b), c, idx, i, j, k);
            }
        }
    }
    for (a = 0; a < 3; ++a) {
        ddalpha[a][a] = prj_z4c_Dxx(z, PRJ_Z4C_ALPHA, a, idx, i, j, k);
        ddchi[a][a] = prj_z4c_Dxx(z, PRJ_Z4C_CHI, a, idx, i, j, k);
        for (b = a + 1; b < 3; ++b) {
            ddalpha[a][b] = ddalpha[b][a] =
                prj_z4c_Dxy(z, PRJ_Z4C_ALPHA, a, b, idx, i, j, k);
            ddchi[a][b] = ddchi[b][a] =
                prj_z4c_Dxy(z, PRJ_Z4C_CHI, a, b, idx, i, j, k);
        }
    }
    for (c = 0; c < 3; ++c) {
        for (a = 0; a < 3; ++a) {
            ddbeta[a][a][c] = prj_z4c_Dxx(z, prj_z4c_beta_var(c), a, idx, i, j, k);
            for (b = a + 1; b < 3; ++b) {
                ddbeta[a][b][c] = ddbeta[b][a][c] =
                    prj_z4c_Dxy(z, prj_z4c_beta_var(c), a, b, idx, i, j, k);
            }
        }
    }
    for (c = 0; c < 3; ++c) {
        for (d = 0; d < 3; ++d) {
            for (a = 0; a < 3; ++a) {
                ddg[a][a][c][d] = prj_z4c_Dxx(z, prj_z4c_g_var(c, d), a, idx, i, j, k);
                for (b = a + 1; b < 3; ++b) {
                    ddg[a][b][c][d] = ddg[b][a][c][d] =
                        prj_z4c_Dxy(z, prj_z4c_g_var(c, d), a, b, idx, i, j, k);
                }
            }
        }
    }
    for (a = 0; a < 3; ++a) {
        for (b = a; b < 3; ++b) {
            for (c = 0; c < 3; ++c) {
                Lg[a][b] += prj_z4c_Lx(z, prj_z4c_beta_var(c),
                    prj_z4c_g_var(a, b), c, idx, i, j, k);
                LA[a][b] += prj_z4c_Lx(z, prj_z4c_beta_var(c),
                    prj_z4c_A_var(a, b), c, idx, i, j, k);
            }
        }
    }

    for (c = 0; c < 3; ++c) {
        for (a = 0; a < 3; ++a) {
            for (b = 0; b < 3; ++b) {
                Gamma_ddd[c][a][b] = 0.5 * (dg[a][b][c] + dg[b][a][c] - dg[c][a][b]);
            }
        }
    }
    for (c = 0; c < 3; ++c) {
        for (a = 0; a < 3; ++a) {
            for (b = 0; b < 3; ++b) {
                for (d = 0; d < 3; ++d) {
                    Gamma_udd[c][a][b] += gu[c][d] * Gamma_ddd[d][a][b];
                }
            }
        }
    }
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            for (c = 0; c < 3; ++c) {
                Gamma_u[a] += gu[b][c] * Gamma_udd[a][b][c];
            }
        }
    }
    for (a = 0; a < 3; ++a) {
        for (b = a; b < 3; ++b) {
            for (c = 0; c < 3; ++c) {
                R_dd[a][b] += 0.5 * (g[c][a] * dGam[b][c] +
                    g[c][b] * dGam[a][c] +
                    Gamma_u[c] * (Gamma_ddd[a][b][c] + Gamma_ddd[b][a][c]));
            }
            for (c = 0; c < 3; ++c) {
                for (d = 0; d < 3; ++d) {
                    R_dd[a][b] -= 0.5 * gu[c][d] * ddg[c][d][a][b];
                }
            }
            for (c = 0; c < 3; ++c) {
                for (d = 0; d < 3; ++d) {
                    for (e = 0; e < 3; ++e) {
                        R_dd[a][b] += gu[c][d] * (
                            Gamma_udd[e][c][a] * Gamma_ddd[b][e][d] +
                            Gamma_udd[e][c][b] * Gamma_ddd[a][e][d] +
                            Gamma_udd[e][a][d] * Gamma_ddd[e][c][b]);
                    }
                }
            }
            R_dd[b][a] = R_dd[a][b];
        }
    }

    for (a = 0; a < 3; ++a) {
        dphi[a] = dchi[a] / (chi_guarded * opt->chi_psi_power);
    }
    for (a = 0; a < 3; ++a) {
        for (b = a; b < 3; ++b) {
            Ddphi_dd[a][b] = ddchi[a][b] / (chi_guarded * opt->chi_psi_power) -
                opt->chi_psi_power * dphi[a] * dphi[b];
            for (c = 0; c < 3; ++c) {
                Ddphi_dd[a][b] -= Gamma_udd[c][a][b] * dphi[c];
            }
            Ddphi_dd[b][a] = Ddphi_dd[a][b];
        }
    }
    for (a = 0; a < 3; ++a) {
        for (b = a; b < 3; ++b) {
            Rphi_dd[a][b] = 4.0 * dphi[a] * dphi[b] - 2.0 * Ddphi_dd[a][b];
            for (c = 0; c < 3; ++c) {
                for (d = 0; d < 3; ++d) {
                    Rphi_dd[a][b] -= 2.0 * g[a][b] * gu[c][d] *
                        (Ddphi_dd[c][d] + 2.0 * dphi[c] * dphi[d]);
                }
            }
            Rphi_dd[b][a] = Rphi_dd[a][b];
        }
    }
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            Ddalpha_dd[a][b] = ddalpha[a][b] -
                2.0 * (dphi[a] * dalpha[b] + dphi[b] * dalpha[a]);
            for (c = 0; c < 3; ++c) {
                Ddalpha_dd[a][b] -= Gamma_udd[c][a][b] * dalpha[c];
                for (d = 0; d < 3; ++d) {
                    Ddalpha_dd[a][b] += 2.0 * g[a][b] * gu[c][d] *
                        dphi[c] * dalpha[d];
                }
            }
            Ddalpha += oopsi4 * gu[a][b] * Ddalpha_dd[a][b];
        }
    }
    for (a = 0; a < 3; ++a) {
        for (b = a; b < 3; ++b) {
            for (c = 0; c < 3; ++c) {
                for (d = 0; d < 3; ++d) {
                    AA_dd[a][b] += gu[c][d] * A[a][c] * A[d][b];
                }
            }
            AA_dd[b][a] = AA_dd[a][b];
        }
    }
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            AA += gu[a][b] * AA_dd[a][b];
            for (c = 0; c < 3; ++c) {
                for (d = 0; d < 3; ++d) {
                    Auu[a][b] += gu[a][c] * gu[b][d] * A[c][d];
                }
            }
        }
    }
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            DA_u[a] -= 1.5 * Auu[a][b] * dchi[b] / chi_guarded;
            DA_u[a] -= (1.0 / 3.0) * gu[a][b] * (2.0 * dKhat[b] + dTheta[b]);
        }
        for (b = 0; b < 3; ++b) {
            for (c = 0; c < 3; ++c) {
                DA_u[a] += Gamma_udd[a][b][c] * Auu[b][c];
            }
        }
    }
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            R += oopsi4 * gu[a][b] * (R_dd[a][b] + Rphi_dd[a][b]);
        }
    }
    Ht = R + (2.0 / 3.0) * K * K - AA;

    for (a = 0; a < 3; ++a) {
        div_beta += dbeta[a][a];
    }
    Lchi += (1.0 / 6.0) * opt->chi_psi_power * chi_guarded * div_beta;
    for (a = 0; a < 3; ++a) {
        LGam[a] += (2.0 / 3.0) * Gamma_u[a] * div_beta;
        for (b = 0; b < 3; ++b) {
            double ddbeta_d = 0.0;

            for (c = 0; c < 3; ++c) {
                ddbeta_d += (1.0 / 3.0) * ddbeta[b][c][c];
            }
            LGam[a] += gu[a][b] * ddbeta_d - Gamma_u[b] * dbeta[b][a];
            for (c = 0; c < 3; ++c) {
                LGam[a] += gu[b][c] * ddbeta[b][c][a];
            }
        }
    }
    for (a = 0; a < 3; ++a) {
        for (b = a; b < 3; ++b) {
            Lg[a][b] -= (2.0 / 3.0) * g[a][b] * div_beta;
            LA[a][b] -= (2.0 / 3.0) * A[a][b] * div_beta;
            for (c = 0; c < 3; ++c) {
                Lg[a][b] += dbeta[a][c] * g[b][c] + dbeta[b][c] * g[a][c];
                LA[a][b] += dbeta[b][c] * A[a][c] + dbeta[a][c] * A[b][c];
            }
        }
    }

    rhs[Z4CIDX(PRJ_Z4C_KHAT, i, j, k)] =
        -Ddalpha + alpha * (AA + (1.0 / 3.0) * K * K) + LKhat +
        opt->damp_kappa1_inv_cm * (1.0 - opt->damp_kappa2) * alpha *
        prj_z4c_get(z, PRJ_Z4C_THETA, i, j, k) +
        4.0 * M_PI * alpha * (T_S + T_E);
    rhs[Z4CIDX(PRJ_Z4C_CHI, i, j, k)] =
        Lchi - (1.0 / 6.0) * opt->chi_psi_power * chi_guarded * alpha * K;
    rhs[Z4CIDX(PRJ_Z4C_THETA, i, j, k)] =
        (LTheta + alpha * (0.5 * Ht - (2.0 + opt->damp_kappa2) *
        opt->damp_kappa1_inv_cm * prj_z4c_get(z, PRJ_Z4C_THETA, i, j, k)) -
        8.0 * M_PI * alpha * T_E) * (double)opt->use_z4c;
    for (a = 0; a < 3; ++a) {
        double value = 2.0 * alpha * DA_u[a] + LGam[a] -
            2.0 * alpha * opt->damp_kappa1_inv_cm *
            (prj_z4c_get(z, prj_z4c_Gam_var(a), i, j, k) - Gamma_u[a]);

        for (b = 0; b < 3; ++b) {
            value -= 2.0 * Auu[a][b] * dalpha[b];
            value -= 16.0 * M_PI * alpha * gu[a][b] * T_S_d[b];
        }
        rhs[Z4CIDX(prj_z4c_Gam_var(a), i, j, k)] = value;
    }
    for (a = 0; a < 3; ++a) {
        for (b = a; b < 3; ++b) {
            double rhs_A;

            rhs[Z4CIDX(prj_z4c_g_var(a, b), i, j, k)] =
                -2.0 * alpha * A[a][b] + Lg[a][b];
            rhs_A = oopsi4 * (-Ddalpha_dd[a][b] +
                alpha * (R_dd[a][b] + Rphi_dd[a][b]));
            rhs_A -= (1.0 / 3.0) * g[a][b] * (-Ddalpha + alpha * R);
            rhs_A += alpha * (K * A[a][b] - 2.0 * AA_dd[a][b]);
            rhs_A += LA[a][b];
            rhs_A -= 8.0 * M_PI * alpha *
                (oopsi4 * T_S_dd[a][b] - (1.0 / 3.0) * T_S * g[a][b]);
            rhs[Z4CIDX(prj_z4c_A_var(a, b), i, j, k)] = rhs_A;
        }
    }
    {
        double f = opt->lapse_oplog * opt->lapse_harmonicf +
            opt->lapse_harmonic * alpha;
        double rhs_alpha = opt->lapse_advect * Lalpha - f * alpha *
            prj_z4c_get(z, PRJ_Z4C_KHAT, i, j, k);

        if (opt->slow_start_lapse) {
            double W2 = chi > opt->chi_min_floor ? chi : opt->chi_min_floor;
            double W = sqrt(W2);
            double tscale = tau_cm / opt->ssl_damping_time_cm;

            rhs_alpha += opt->ssl_damping_amp_inv_cm * (W - alpha) *
                pow(W, (double)opt->ssl_damping_index) * exp(-0.5 * tscale * tscale);
        }
        rhs[Z4CIDX(PRJ_Z4C_ALPHA, i, j, k)] = rhs_alpha;
    }
    for (a = 0; a < 3; ++a) {
        double value = opt->shift_Gamma * prj_z4c_get(z, prj_z4c_Gam_var(a), i, j, k) +
            opt->shift_advect * Lbeta[a] -
            opt->shift_eta_inv_cm * prj_z4c_get(z, prj_z4c_beta_var(a), i, j, k);

        value += opt->shift_alpha2Gamma * alpha * alpha *
            prj_z4c_get(z, prj_z4c_Gam_var(a), i, j, k);
        for (b = 0; b < 3; ++b) {
            value += opt->shift_H * alpha * chi_guarded *
                (0.5 * alpha * dchi[b] - dalpha[b]) * gu[a][b];
        }
        rhs[Z4CIDX(prj_z4c_beta_var(a), i, j, k)] = value;
    }
    {
        double diss_coeff = prj_z4c_diss_coeff(opt);

        if (diss_coeff == 0.0) {
            return;
        }
        for (var = 0; var < PRJ_NZ4C; ++var) {
            for (a = 0; a < 3; ++a) {
                rhs[Z4CIDX(var, i, j, k)] += diss_coeff *
                    prj_z4c_Diss(z, var, a, idx, i, j, k);
            }
        }
    }
}

void prj_z4c_compute_rhs(prj_mesh *mesh, const prj_mpi *mpi,
    const prj_rad *rad, int state_stage, int rhs_stage, double tau_cm)
{
    int bidx;
    int use_full_dynamic_gr;

    if (!prj_z4c_runtime_enabled(mesh)) {
        return;
    }
    if (state_stage < 0 || state_stage >= PRJ_BLOCK_NSTAGES ||
        rhs_stage < 0 || rhs_stage >= PRJ_BLOCK_NSTAGES) {
        prj_z4c_fail("prj_z4c_compute_rhs: invalid stage");
    }
    use_full_dynamic_gr = prj_eos_full_dynamic_gr_enabled(mesh);
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        const double *z;
        const double *W_mhd;
        const double *W_rad;
        double *rhs;
        int i, j, k;

        if (!prj_z4c_local_block(mpi, block)) {
            continue;
        }
        z = prj_block_z4c_stage_const(block, state_stage);
        W_mhd = prj_block_mhd_stage_const(block, state_stage);
        W_rad = prj_block_rad_stage_const(block, state_stage);
        rhs = prj_block_z4c_rhs_stage(block, rhs_stage);
        if (z == 0 || rhs == 0) {
            prj_z4c_fail("prj_z4c_compute_rhs: missing stage storage");
        }
        prj_fill(rhs, (size_t)PRJ_NZ4C * (size_t)PRJ_BLOCK_NCELLS_Z4C, 0.0);
        for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
            for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    prj_z4c_hydro_geom geom;
                    double tmunu[PRJ_NTMUNU];

                    if (use_full_dynamic_gr) {
                        if (!prj_z4c_load_hydro_geom(mesh, block, state_stage,
                                i, j, k, &geom)) {
                            prj_z4c_fail("prj_z4c_compute_rhs: failed to load matter geometry");
                        }
                    } else {
                        prj_z4c_load_minkowski_hydro_geom(&geom);
                    }
                    prj_z4c_build_tmunu_cell_geo(rad, block, &geom,
                        W_mhd, W_rad, i, j, k, tmunu);
                    prj_z4c_compute_rhs_cell(mesh, block, z, tmunu, rhs,
                        i, j, k, tau_cm);
                }
            }
        }
    }
}

static int prj_z4c_sommerfeld_face_enabled(const prj_mesh *mesh, const prj_bc *bc,
    const prj_block *block, int axis, int outer)
{
    double block_face;
    double domain_face;
    double tol;
    int bc_type;

    if (mesh == 0 || bc == 0 || block == 0) {
        return 0;
    }
    block_face = outer ? block->xmax[axis] : block->xmin[axis];
    domain_face = prj_z4c_axis_domain_face(mesh, axis, outer);
    tol = 1.0e-8 * PRJ_MAX(fabs(block->dx[axis]), 1.0);
    if (fabs(block_face - domain_face) > tol) {
        return 0;
    }
    bc_type = prj_z4c_axis_bc_type(bc, axis, outer);
    if (bc_type == PRJ_BC_OUTFLOW) {
        return 1;
    }
    if (bc_type == PRJ_BC_USER && mesh->z4c_params.user_Sbc != 0) {
        return 1;
    }
    return 0;
}

static double prj_z4c_sommerfeld_rhs_var(const prj_block *block, const double *z,
    int var, int i, int j, int k, const double s[3], double r, double speed)
{
    double idx[3] = {1.0 / block->dx[0], 1.0 / block->dx[1], 1.0 / block->dx[2]};
    double adv = 0.0;
    int axis;

    for (axis = 0; axis < 3; ++axis) {
        adv += s[axis] * prj_z4c_Dx2nd(z, var, axis, idx, i, j, k);
    }
    return -speed * prj_z4c_get(z, var, i, j, k) / r - speed * adv;
}

static void prj_z4c_apply_sommerfeld_cell(const prj_block *block, const double *z,
    double *rhs, int i, int j, int k)
{
    double x[3];
    double r;
    double s[3];
    double sqrt2 = sqrt(2.0);
    int a;

    x[0] = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
    x[1] = block->xmin[1] + ((double)j + 0.5) * block->dx[1];
    x[2] = block->xmin[2] + ((double)k + 0.5) * block->dx[2];
    r = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
    if (r < 1.0e-300) {
        r = 1.0e-300;
    }
    for (a = 0; a < 3; ++a) {
        s[a] = x[a] / r;
    }

    rhs[Z4CIDX(PRJ_Z4C_THETA, i, j, k)] =
        prj_z4c_sommerfeld_rhs_var(block, z, PRJ_Z4C_THETA, i, j, k, s, r, 1.0);
    rhs[Z4CIDX(PRJ_Z4C_KHAT, i, j, k)] =
        prj_z4c_sommerfeld_rhs_var(block, z, PRJ_Z4C_KHAT, i, j, k, s, r, sqrt2);
    for (a = 0; a < 3; ++a) {
        rhs[Z4CIDX(prj_z4c_Gam_var(a), i, j, k)] =
            prj_z4c_sommerfeld_rhs_var(block, z, prj_z4c_Gam_var(a),
                i, j, k, s, r, 1.0);
    }
    rhs[Z4CIDX(PRJ_Z4C_AXX, i, j, k)] =
        prj_z4c_sommerfeld_rhs_var(block, z, PRJ_Z4C_AXX, i, j, k, s, r, 1.0);
    rhs[Z4CIDX(PRJ_Z4C_AXY, i, j, k)] =
        prj_z4c_sommerfeld_rhs_var(block, z, PRJ_Z4C_AXY, i, j, k, s, r, 1.0);
    rhs[Z4CIDX(PRJ_Z4C_AXZ, i, j, k)] =
        prj_z4c_sommerfeld_rhs_var(block, z, PRJ_Z4C_AXZ, i, j, k, s, r, 1.0);
    rhs[Z4CIDX(PRJ_Z4C_AYY, i, j, k)] =
        prj_z4c_sommerfeld_rhs_var(block, z, PRJ_Z4C_AYY, i, j, k, s, r, 1.0);
    rhs[Z4CIDX(PRJ_Z4C_AYZ, i, j, k)] =
        prj_z4c_sommerfeld_rhs_var(block, z, PRJ_Z4C_AYZ, i, j, k, s, r, 1.0);
    rhs[Z4CIDX(PRJ_Z4C_AZZ, i, j, k)] =
        prj_z4c_sommerfeld_rhs_var(block, z, PRJ_Z4C_AZZ, i, j, k, s, r, 1.0);
}

static void prj_z4c_apply_sommerfeld_face(const prj_mesh *mesh, const prj_bc *bc,
    const prj_block *block, const double *z, double *rhs, int axis, int outer)
{
    int ilo = 0, ihi = PRJ_BLOCK_SIZE;
    int jlo = 0, jhi = PRJ_BLOCK_SIZE;
    int klo = 0, khi = PRJ_BLOCK_SIZE;
    int i, j, k;

    if (!prj_z4c_sommerfeld_face_enabled(mesh, bc, block, axis, outer)) {
        return;
    }
    if (axis == X1DIR) {
        ilo = outer ? PRJ_BLOCK_SIZE - 1 : 0;
        ihi = ilo + 1;
    } else if (axis == X2DIR) {
        jlo = outer ? PRJ_BLOCK_SIZE - 1 : 0;
        jhi = jlo + 1;
    } else {
        klo = outer ? PRJ_BLOCK_SIZE - 1 : 0;
        khi = klo + 1;
    }
    for (i = ilo; i < ihi; ++i) {
        for (j = jlo; j < jhi; ++j) {
            for (k = klo; k < khi; ++k) {
                prj_z4c_apply_sommerfeld_cell(block, z, rhs, i, j, k);
            }
        }
    }
}

void prj_z4c_apply_sommerfeld_rhs(prj_mesh *mesh, const prj_mpi *mpi,
    const prj_bc *bc, int state_stage, int rhs_stage)
{
    int bidx;

    if (!prj_z4c_runtime_enabled(mesh) || bc == 0) {
        return;
    }
    if (state_stage < 0 || state_stage >= PRJ_BLOCK_NSTAGES ||
        rhs_stage < 0 || rhs_stage >= PRJ_BLOCK_NSTAGES) {
        prj_z4c_fail("prj_z4c_apply_sommerfeld_rhs: invalid stage");
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        const double *z;
        double *rhs;
        int axis;

        if (!prj_z4c_local_block(mpi, block)) {
            continue;
        }
        z = prj_block_z4c_stage_const(block, state_stage);
        rhs = prj_block_z4c_rhs_stage(block, rhs_stage);
        for (axis = 0; axis < 3; ++axis) {
            prj_z4c_apply_sommerfeld_face(mesh, bc, block, z, rhs, axis, 0);
            prj_z4c_apply_sommerfeld_face(mesh, bc, block, z, rhs, axis, 1);
        }
    }
}

void prj_z4c_update_linear(prj_mesh *mesh, const prj_mpi *mpi,
    int dst_stage, int a_stage, double a_w, int b_stage, double b_w,
    int rhs_stage, double dtau_cm)
{
    int bidx;

    if (!prj_z4c_runtime_enabled(mesh)) {
        return;
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        double *dst;
        const double *a_src;
        const double *b_src;
        const double *rhs;
        int i, j, k, var;

        if (!prj_z4c_local_block(mpi, block)) {
            continue;
        }
        dst = prj_block_z4c_stage(block, dst_stage);
        a_src = prj_block_z4c_stage_const(block, a_stage);
        b_src = prj_block_z4c_stage_const(block, b_stage);
        rhs = prj_block_z4c_rhs_stage_const(block, rhs_stage);
        for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
            for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    for (var = 0; var < PRJ_NZ4C; ++var) {
                        dst[Z4CIDX(var, i, j, k)] =
                            a_w * a_src[Z4CIDX(var, i, j, k)] +
                            b_w * b_src[Z4CIDX(var, i, j, k)] +
                            dtau_cm * rhs[Z4CIDX(var, i, j, k)];
                    }
                }
            }
        }
    }
}

void prj_z4c_update_linear_cell(prj_block *block,
    int dst_stage, int a_stage, double a_w, int b_stage, double b_w,
    int rhs_stage, double dtau_cm, int i, int j, int k)
{
    double *dst = prj_block_z4c_stage(block, dst_stage);
    const double *a_src = prj_block_z4c_stage_const(block, a_stage);
    const double *b_src = prj_block_z4c_stage_const(block, b_stage);
    const double *rhs = prj_block_z4c_rhs_stage_const(block, rhs_stage);
    int var;

    for (var = 0; var < PRJ_NZ4C; ++var) {
        dst[Z4CIDX(var, i, j, k)] =
            a_w * a_src[Z4CIDX(var, i, j, k)] +
            b_w * b_src[Z4CIDX(var, i, j, k)] +
            dtau_cm * rhs[Z4CIDX(var, i, j, k)];
    }
}

void prj_z4c_finalize_stage(prj_mesh *mesh, prj_mpi *mpi, const prj_bc *bc, int stage)
{
    if (!prj_z4c_runtime_enabled(mesh)) {
        return;
    }
    prj_z4c_enforce_range(mesh, mpi, stage, 0, PRJ_BLOCK_SIZE, 0, PRJ_BLOCK_SIZE,
        0, PRJ_BLOCK_SIZE);
    prj_z4c_fill_ghosts(mesh, mpi, bc, stage);
    prj_z4c_enforce_range(mesh, mpi, stage, -PRJ_NGHOST_Z4C, PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C,
        -PRJ_NGHOST_Z4C, PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C,
        -PRJ_NGHOST_Z4C, PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C);
}

void prj_z4c_save_stage(prj_mesh *mesh, const prj_mpi *mpi, int dst_stage, int src_stage)
{
    int bidx;

    if (!prj_z4c_runtime_enabled(mesh)) {
        return;
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        double *dst;
        const double *src;

        if (!prj_z4c_local_block(mpi, block)) {
            continue;
        }
        dst = prj_block_z4c_stage(block, dst_stage);
        src = prj_block_z4c_stage_const(block, src_stage);
        memcpy(dst, src, (size_t)PRJ_NZ4C * (size_t)PRJ_BLOCK_NCELLS_Z4C * sizeof(*dst));
    }
}

void prj_z4c_blend_with_saved(prj_mesh *mesh, prj_mpi *mpi, const prj_bc *bc,
    int saved_stage, double saved_weight)
{
    int bidx;
    double current_weight = 1.0 - saved_weight;

    if (!prj_z4c_runtime_enabled(mesh)) {
        return;
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        double *z;
        const double *saved;
        int i, j, k, var;

        if (!prj_z4c_local_block(mpi, block)) {
            continue;
        }
        z = prj_block_z4c_stage(block, 0);
        saved = prj_block_z4c_stage_const(block, saved_stage);
        for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
            for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    for (var = 0; var < PRJ_NZ4C; ++var) {
                        z[Z4CIDX(var, i, j, k)] =
                            saved_weight * saved[Z4CIDX(var, i, j, k)] +
                            current_weight * z[Z4CIDX(var, i, j, k)];
                    }
                }
            }
        }
    }
    prj_z4c_finalize_stage(mesh, mpi, bc, 0);
}

void prj_z4c_amr_prolongate_child(const prj_block *parent, prj_block *child, int child_oct)
{
    int i, j, k, var, stage;

    if (parent == 0 || child == 0 || parent->z4c == 0 || child->z4c == 0) {
        return;
    }
    for (stage = 0; stage < PRJ_BLOCK_NSTAGES; ++stage) {
        const double *src = prj_block_z4c_stage_const(parent, stage);
        double *dst = prj_block_z4c_stage(child, stage);

        for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
            for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    for (var = 0; var < PRJ_NZ4C; ++var) {
                        dst[Z4CIDX(var, i, j, k)] =
                            prj_z4c_amr_prolong_value(src, var, i, j, k, child_oct);
                    }
                }
            }
        }
    }
    prj_z4c_clear_block_aux(child);
}

void prj_z4c_amr_restrict_parent(const prj_block *children[8], prj_block *parent)
{
    int oct, stage;

    if (children == 0 || parent == 0 || parent->z4c == 0) {
        return;
    }
    for (stage = 0; stage < PRJ_BLOCK_NSTAGES; ++stage) {
        double *dst = prj_block_z4c_stage(parent, stage);

        for (oct = 0; oct < 8; ++oct) {
            const prj_block *child = children[oct];
            const double *src;
            int ioff = (oct & 1) ? PRJ_BLOCK_SIZE / 2 : 0;
            int joff = (oct & 2) ? PRJ_BLOCK_SIZE / 2 : 0;
            int koff = (oct & 4) ? PRJ_BLOCK_SIZE / 2 : 0;
            int i, j, k, var;

            if (child == 0 || child->z4c == 0) {
                continue;
            }
            src = prj_block_z4c_stage_const(child, stage);
            for (i = 0; i < PRJ_BLOCK_SIZE / 2; ++i) {
                for (j = 0; j < PRJ_BLOCK_SIZE / 2; ++j) {
                    for (k = 0; k < PRJ_BLOCK_SIZE / 2; ++k) {
                        for (var = 0; var < PRJ_NZ4C; ++var) {
                            dst[Z4CIDX(var, i + ioff, j + joff, k + koff)] =
                                prj_z4c_amr_restrict_value(src, var, i, j, k);
                        }
                    }
                }
            }
        }
    }
    prj_z4c_clear_block_aux(parent);
}

#else

int prj_z4c_runtime_enabled(const prj_mesh *mesh)
{
    (void)mesh;
    return 0;
}

int prj_z4c_load_hydro_geom(const prj_mesh *mesh, const prj_block *block,
    int stage, int i, int j, int k, prj_z4c_hydro_geom *geom)
{
    (void)mesh;
    (void)block;
    (void)stage;
    (void)i;
    (void)j;
    (void)k;
    (void)geom;
    return 0;
}

int prj_z4c_load_hydro_metric_geom(const prj_mesh *mesh, const prj_block *block,
    int stage, int i, int j, int k, prj_z4c_hydro_geom *geom)
{
    (void)mesh;
    (void)block;
    (void)stage;
    (void)i;
    (void)j;
    (void)k;
    (void)geom;
    return 0;
}

double prj_z4c_calc_dt_seconds(const prj_mesh *mesh, const prj_mpi *mpi, double cfl)
{
    (void)mesh;
    (void)mpi;
    (void)cfl;
    return HUGE_VAL;
}

void prj_z4c_init_mesh_flat(prj_mesh *mesh, const prj_mpi *mpi)
{
    (void)mesh;
    (void)mpi;
}

void prj_z4c_init_mesh_spherical(prj_mesh *mesh, const prj_mpi *mpi,
    const double x_com[3], const double *r_grid, const double *chi_grid,
    const double *alpha_grid, int npts)
{
    (void)mesh;
    (void)mpi;
    (void)x_com;
    (void)r_grid;
    (void)chi_grid;
    (void)alpha_grid;
    (void)npts;
}

void prj_z4c_init_punctures(prj_mesh *mesh, const prj_mpi *mpi, int npunctures,
    const double centers_cm[][3], const double masses_cm[],
    const double momenta_cm[][3], double floor_radius_cm)
{
    (void)mesh;
    (void)mpi;
    (void)npunctures;
    (void)centers_cm;
    (void)masses_cm;
    (void)momenta_cm;
    (void)floor_radius_cm;
}

void prj_z4c_fill_ghosts(prj_mesh *mesh, prj_mpi *mpi, const prj_bc *bc, int stage)
{
    (void)mesh;
    (void)mpi;
    (void)bc;
    (void)stage;
}

void prj_z4c_compute_rhs(prj_mesh *mesh, const prj_mpi *mpi,
    const prj_rad *rad, int state_stage, int rhs_stage, double tau_cm)
{
    (void)mesh;
    (void)mpi;
    (void)rad;
    (void)state_stage;
    (void)rhs_stage;
    (void)tau_cm;
}

void prj_z4c_apply_sommerfeld_rhs(prj_mesh *mesh, const prj_mpi *mpi,
    const prj_bc *bc, int state_stage, int rhs_stage)
{
    (void)mesh;
    (void)mpi;
    (void)bc;
    (void)state_stage;
    (void)rhs_stage;
}

void prj_z4c_update_linear(prj_mesh *mesh, const prj_mpi *mpi,
    int dst_stage, int a_stage, double a_w, int b_stage, double b_w,
    int rhs_stage, double dtau_cm)
{
    (void)mesh;
    (void)mpi;
    (void)dst_stage;
    (void)a_stage;
    (void)a_w;
    (void)b_stage;
    (void)b_w;
    (void)rhs_stage;
    (void)dtau_cm;
}

void prj_z4c_update_linear_cell(prj_block *block,
    int dst_stage, int a_stage, double a_w, int b_stage, double b_w,
    int rhs_stage, double dtau_cm, int i, int j, int k)
{
    (void)block;
    (void)dst_stage;
    (void)a_stage;
    (void)a_w;
    (void)b_stage;
    (void)b_w;
    (void)rhs_stage;
    (void)dtau_cm;
    (void)i;
    (void)j;
    (void)k;
}

void prj_z4c_finalize_stage(prj_mesh *mesh, prj_mpi *mpi, const prj_bc *bc, int stage)
{
    (void)mesh;
    (void)mpi;
    (void)bc;
    (void)stage;
}

void prj_z4c_save_stage(prj_mesh *mesh, const prj_mpi *mpi, int dst_stage, int src_stage)
{
    (void)mesh;
    (void)mpi;
    (void)dst_stage;
    (void)src_stage;
}

void prj_z4c_blend_with_saved(prj_mesh *mesh, prj_mpi *mpi, const prj_bc *bc,
    int saved_stage, double saved_weight)
{
    (void)mesh;
    (void)mpi;
    (void)bc;
    (void)saved_stage;
    (void)saved_weight;
}

void prj_z4c_amr_prolongate_child(const prj_block *parent, prj_block *child, int child_oct)
{
    (void)parent;
    (void)child;
    (void)child_oct;
}

void prj_z4c_amr_restrict_parent(const prj_block *children[8], prj_block *parent)
{
    (void)children;
    (void)parent;
}

#endif
