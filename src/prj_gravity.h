#ifndef PRJ_GRAVITY_H
#define PRJ_GRAVITY_H

#define PRJ_GNEWT 6.67430e-8
#define PRJ_GRAVITY_USE_GR 1

void prj_gravity_init(prj_sim *sim);
void prj_gravity_rebuild_grid(prj_sim *sim);
void prj_gravity_free(prj_grav_mono *grav_mono);
void prj_gravity_monopole_reduce(prj_mesh *mesh, int stage);
void prj_gravity_monopole_integrate(prj_mesh *mesh);
void prj_gravity_cache_block(prj_block *block);
void prj_gravity_cache_mesh(prj_mesh *mesh);
double prj_gravity_block_accel_at(const prj_block *block, int i, int j, int k);
double prj_gravity_block_lapse_at(const prj_block *block, int i, int j, int k);
const prj_grav_mono *prj_gravity_active_monopole(void);
double prj_gravity_interp_accel(const prj_grav_mono *grav_mono, double r);
double prj_gravity_interp_lapse(const prj_grav_mono *grav_mono, double r);
int prj_gravity_apply(void);

#endif
