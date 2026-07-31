#ifndef PRJ_RADIATION_H
#define PRJ_RADIATION_H

#include <math.h>

#define PRJ_CLIGHT 2.99792458e10
#define PRJ_HPLANCK 6.62607015e-27

#if PRJ_NRAD > 0 && PRJ_MIXED_PRECISION_FLUX
/* Radiation E/F are stored internally in RAD_SCALE*erg units, so they are
   already O(1) and fit in single precision without any further rescaling. */
static inline float prj_rad_mixed_pack(double x)
{
    return (float)x;
}

static inline double prj_rad_mixed_unpack(float x)
{
    return (double)x;
}

static inline double prj_rad_mixed_round(double x)
{
    return prj_rad_mixed_unpack(prj_rad_mixed_pack(x));
}
#endif

void prj_rad_init(prj_rad *rad);
#if PRJ_USE_RADIATION_FSA
void prj_rad_fsa_calculate_directions(prj_rad *rad);
void prj_rad_fsa_free_geometry(prj_rad *rad);
#if PRJ_USE_RADIAL_FRAME_FSA
void prj_rad_fsa_refresh_block_geometry(const prj_rad *rad, prj_block *block);
void prj_rad_fsa_refresh_mesh_geometry(const prj_rad *rad, prj_mesh *mesh,
    const prj_mpi *mpi);
#endif
void prj_rad_fsa_rotated_dir(const prj_block *block, int i, int j, int k,
    const double n0[3], double n[3]);
void prj_rad_fsa_rotated_angle_dir(const prj_rad *rad, const prj_block *block,
    int angle, int i, int j, int k, double n[3]);
#endif
void prj_rad_prim2cons(const double *W, double *U);
void prj_rad_cons2prim(const double *U, double *W);
void prj_rad_flux(const prj_rad *rad, const double *WL, const double *WR,
    double lapse, const double *chi_face,
    double dx_dir, double v_face, double *flux);
#if PRJ_USE_RADIATION_FSA
void prj_rad_flux_fsa(const prj_rad *rad, const prj_block *block,
    const double *WL, const double *WR, double lapse, int dir, double v_face,
    int il, int jl, int kl, int ir, int jr, int kr, double *flux);
#endif
void prj_rad_energy_update(prj_rad *rad, prj_eos *eos, double *u, double dt, double lapse, double *final_temperature, double *kappa_out);
void prj_rad_momentum_update(prj_rad *rad, prj_eos *eos, double *u, double dt, double lapse, double temperature, const double *kappa_in);
#if PRJ_USE_RADIATION_FSA
void prj_rad_fsa_clamp_intensities(double *u);
void prj_rad_energy_update_fsa(prj_rad *rad, prj_eos *eos, double *u, double dt, double lapse, double *final_temperature, double *final_ye, double *kappa_out, double *eta_out);
void prj_rad_energy_momentum_update_fsa(prj_rad *rad, const prj_block *block,
    int i, int j, int k, prj_eos *eos, double *u, double dt, double lapse);
#if DO_FFC
void prj_rad_ffc_fsa(prj_rad *rad, double *u, double dt);
#endif
#endif
void prj_rad_m1_pressure(const prj_rad *rad, double E, double F1, double F2, double F3, double P[3][3]);
#if PRJ_DYNAMIC_GR && PRJ_USE_RADIATION_M1
typedef struct prj_rad_gr_m1_closure_ctx {
    double gamma[3][3];
    double gamma_inv[3][3];
    double dgamma[3][3][3];
    double K_dd[3][3];
    double vcon[3];
    double dvdx[3][3];
} prj_rad_gr_m1_closure_ctx;

/* Per-side fluid kinematics for GR M1 moment projections. These are invariant
 * across all energy groups of a face side, so callers that need the same frame
 * repeatedly build this once with prj_rad_gr_m1_prepare_side and pass it to the
 * _cached helpers. */
typedef struct prj_rad_gr_m1_side_data {
    double vcon[3];
    double vcov[3];
    double u_cov[3];
    double wlor;
    double beta2;
} prj_rad_gr_m1_side_data;

typedef struct prj_rad_grm1_closure_coeffs {
    double a_coef;
    double b_coef;
    double s;
    double f2;
} prj_rad_grm1_closure_coeffs;

static inline void prj_rad_grm1_closure_coeffs_from_F2(double E, double F2,
    prj_rad_grm1_closure_coeffs *coeffs)
{
    double cE;
    double f2;
    double s;

    if (coeffs == 0) {
        return;
    }
    if (!isfinite(F2) || F2 < 0.0) {
        F2 = 0.0;
    }
    cE = PRJ_CLIGHT * E;
    f2 = cE > 0.0 ? F2 / (cE * cE) : 0.0;
    if (!isfinite(f2) || f2 < 0.0) {
        f2 = 0.0;
    } else if (f2 > 1.0) {
        f2 = 1.0;
    }
    s = sqrt(4.0 - 3.0 * f2);
    coeffs->a_coef = E * (s - 1.0) / 3.0;
    coeffs->b_coef = E > 0.0 ?
        3.0 / ((2.0 + s) * PRJ_CLIGHT * PRJ_CLIGHT * E) : 0.0;
    coeffs->s = s;
    coeffs->f2 = f2;
}

void prj_rad_gr_m1_prepare_side(const prj_rad_gr_m1_closure_ctx *ctx,
    prj_rad_gr_m1_side_data *side);

/* Spatial M1 pressure tensor in the Eulerian (lab) frame, algebraic Levermore
 * closure: P^{ij} = E a(f) gamma^{ij} + b(f) F^i F^j/(c^2 E), f = |F|/(cE).
 * gamma_inv is the spatial inverse 3-metric; Fcon the contravariant spatial
 * flux; Fmag = sqrt(F_i F^i).  Callers guarantee E > 0 and f <= 1. */
void prj_rad_grm1_pressure_lab(double E, const double Fcon[3],
    double Fmag, const double gamma_inv[3][3], double P[3][3]);

/* Build R^{alpha beta} = E n^a n^b + (n^a F^b + n^b F^a)/c + P^{ab} from the
 * Eulerian moments with the algebraic lab-frame Levermore closure (spatial
 * block from prj_rad_grm1_pressure_lab). `Fcon` is the physical contravariant
 * lab-frame flux; the stress-tensor time-space leg is F^i/c internally. */
int prj_rad_grm1_build_R(const double g_cov[4][4], const double g_con[4][4],
    double alpha, double E, const double Fcon[3], double Rcon[4][4]);

/* Contract the third radiation moment with the velocity gradient,
 * drift^alpha = M^{alpha beta gamma} u_{beta;gamma}.  The third moment is not
 * materialized; the contraction is formed from J, H^alpha, L^{alpha beta}, and
 * the analytic Levermore closure.  `ucon`, `Rcon`, and `ducov` are in the same
 * basis as g_ab/g^ab. */
int prj_rad_grm1_freq_drift(const double g_cov[4][4],
    const double g_con[4][4], const double ucon[4],
    const double Rcon[4][4], const double ducov[4][4], double drift[4]);

void prj_rad_gr_m1_pressure(const prj_rad *rad,
    const prj_rad_gr_m1_closure_ctx *ctx, double E, const double Fcov[3],
    double P[3][3]);

void prj_rad_gr_m1_pressure_cached(const prj_rad *rad,
    const prj_rad_gr_m1_closure_ctx *ctx,
    const prj_rad_gr_m1_side_data *side, double E, const double Fcov[3],
    double P[3][3]);

void prj_rad_gr_m1_pressure_fbar_cached(const prj_rad *rad,
    const prj_rad_gr_m1_closure_ctx *ctx,
    const prj_rad_gr_m1_side_data *side, double E, const double Fcov[3],
    double P[3][3], double *fbar_out, double *J0_out, double H0_out[3]);
#endif
#if PRJ_NRAD > 0
void prj_rad_m1_wavespeeds(double E, double F1, double F2, double F3,
    double *lam_min, double *lam_max);
void prj_rad_m1_wavespeeds_with_fluxmag(double E, double F1, double Fmag, double inv_Fmag,
    double f, double *lam_min, double *lam_max);
#endif
void prj_rad_freq_flux_apply(const prj_rad *rad, const prj_block *block,
    const double *W_state, double *u, int ic, int jc, int kc, double lapse, double dt);
#if PRJ_DYNAMIC_GR && PRJ_USE_RADIATION_M1
struct prj_z4c_hydro_geom;
void prj_rad_freq_flux_apply_gr_m1_geom(const prj_rad *rad,
    const prj_mesh *mesh, const prj_block *block, int z4c_stage,
    const struct prj_z4c_hydro_geom *geom, const double *W_state, double *u,
    int ic, int jc, int kc, double dt,
    const double observer_time_derivative[4]);
void prj_rad_freq_flux_apply_gr_m1(const prj_rad *rad, const prj_mesh *mesh,
    const prj_block *block, int z4c_stage, const double *W_state, double *u,
    int ic, int jc, int kc, double dt, const double observer_time_derivative[4]);
void prj_rad_gr_m1_matter_update(prj_rad *rad, prj_eos *eos,
    const prj_mesh *mesh, const prj_block *block, int z4c_stage, double *u,
    double *prim, int i, int j, int k, double dt,
    double *final_temperature);
void prj_rad_gr_m1_matter_update_geom(prj_rad *rad, prj_eos *eos,
    const prj_mesh *mesh, const prj_block *block, int z4c_stage,
    const struct prj_z4c_hydro_geom *geom, double *u, double *prim,
    int i, int j, int k, double dt, double *final_temperature);
#endif
void prj_rad_ang_flux_apply(const prj_rad *rad, const prj_block *block,
    const double *W_state, double *u, int ic, int jc, int kc, double lapse, double dt);

#endif
