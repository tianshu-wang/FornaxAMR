#ifndef PRJ_FLUX_H
#define PRJ_FLUX_H

void prj_flux_update(prj_eos *eos, prj_rad *rad, const prj_mesh *mesh,
    prj_block *block, double *W, double *eosvar, double *flux[3], int use_bf1);
void prj_flux_fill_transport_opacity_active(prj_mesh *mesh, prj_rad *rad,
    const prj_mpi *mpi, int stage);
void prj_flux_fill_transport_opacity_halo(prj_mesh *mesh, prj_rad *rad,
    const prj_mpi *mpi, int stage);
void prj_flux_div(const prj_block *block, int i, int j, int k, double *fluxdiv);

#endif
