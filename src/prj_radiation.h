#ifndef PRJ_RADIATION_H
#define PRJ_RADIATION_H

#define PRJ_CLIGHT 2.99792458e10
#define PRJ_HPLANCK 6.62607015e-27

#if PRJ_NRAD > 0 && PRJ_MIXED_PRECISION
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
void prj_rad_prim2cons(const double *W, double *U);
void prj_rad_cons2prim(const double *U, double *W);
void prj_rad_flux(const prj_rad *rad, const double *WL, const double *WR,
    double lapse, const double *chi_face,
    double dx_dir, double v_face, double *flux);
void prj_rad_energy_update(prj_rad *rad, prj_eos *eos, double *u, double dt, double lapse, double *final_temperature, double *kappa_out);
void prj_rad_momentum_update(prj_rad *rad, prj_eos *eos, double *u, double dt, double lapse, double temperature, const double *kappa_in);
void prj_rad_m1_pressure(const prj_rad *rad, double E, double F1, double F2, double F3, double P[3][3]);
#if PRJ_NRAD > 0
void prj_rad_m1_wavespeeds(double E, double F1, double F2, double F3,
    double *lam_min, double *lam_max);
void prj_rad_m1_wavespeeds_with_fluxmag(double E, double F1, double Fmag, double inv_Fmag,
    double f, double *lam_min, double *lam_max);
#endif
void prj_rad_freq_flux_apply(const prj_rad *rad, const prj_block *block,
    const double *W_state, double *u, int ic, int jc, int kc, double lapse, double dt);

#endif
