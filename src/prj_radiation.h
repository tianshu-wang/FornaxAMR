#ifndef PRJ_RADIATION_H
#define PRJ_RADIATION_H

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
    double opacity;
    int have_shear;
} prj_rad_gr_m1_closure_ctx;

void prj_rad_gr_m1_pressure(const prj_rad *rad,
    const prj_rad_gr_m1_closure_ctx *ctx, double E, const double Fcov[3],
    double P[3][3]);
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
void prj_rad_freq_flux_apply_gr_m1(const prj_rad *rad, const prj_mesh *mesh,
    const prj_block *block, int z4c_stage, const double *W_state, double *u,
    int ic, int jc, int kc, double dt);
void prj_rad_gr_m1_matter_update(prj_rad *rad, prj_eos *eos,
    const prj_mesh *mesh, const prj_block *block, int z4c_stage, double *u,
    int i, int j, int k, double dt, double *final_temperature);
#endif
void prj_rad_ang_flux_apply(const prj_rad *rad, const prj_block *block,
    const double *W_state, double *u, int ic, int jc, int kc, double lapse, double dt);

#endif
