#ifndef PRJ_GRAVITY_H
#define PRJ_GRAVITY_H

#define PRJ_GNEWT 6.67430e-8
#define PRJ_GRAVITY_USE_GR 1
#ifndef LMAX
#define LMAX 4
#endif
#if LMAX < 1
#error "LMAX must be at least 1"
#endif
#define PRJ_YLM_INDEX(l, m) ((l) * (l) + (l) + (m))

void prj_gravity_init(prj_sim *sim);
void prj_gravity_rebuild_grid(prj_sim *sim);
void prj_gravity_free(prj_grav *grav);
void prj_gravity_monopole_reduce(prj_mesh *mesh, int stage);
void prj_gravity_monopole_integrate(prj_mesh *mesh);
int prj_gravity_update_center_of_mass(prj_mesh *mesh, double x_com_err_tol);
void prj_gravity_cache_block(prj_block *block);
void prj_gravity_cache_mesh(prj_mesh *mesh);
double prj_gravity_block_accel_at(const prj_block *block, int i, int j, int k);
double prj_gravity_block_lapse_at(const prj_block *block, int i, int j, int k);
const prj_grav *prj_gravity_active_monopole(void);
void prj_gravity_real_spherical_harmonics_all(double x1, double x2, double x3, double *out);
double prj_gravity_real_spherical_harmonic(int l, int m, double x1, double x2, double x3);
double prj_gravity_interp_accel(const prj_grav *grav, double r);
double prj_gravity_interp_lapse(const prj_grav *grav, double r);
int prj_gravity_apply(void);

#endif
