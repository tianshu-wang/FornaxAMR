#ifndef PRJ_GRAVITY_H
#define PRJ_GRAVITY_H

#define PRJ_GNEWT 6.67430e-8
#define PRJ_GRAVITY_USE_GR 1

void prj_gravity_init(prj_sim *sim);
void prj_gravity_free(prj_grav_mono *grav_mono);
void prj_gravity_monopole_reduce(prj_mesh *mesh, prj_eos *eos);
void prj_gravity_monopole_integrate(prj_mesh *mesh);
const prj_grav_mono *prj_gravity_active_monopole(void);
int prj_gravity_apply(void);

#endif
