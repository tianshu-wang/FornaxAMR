#ifndef PRJ_EOS_H
#define PRJ_EOS_H

enum prj_eos_quantity {
    PRJ_EOS_EINT = 0,
    PRJ_EOS_PRESSURE = 1,
    PRJ_EOS_GAMMA = 2,
    PRJ_EOS_TEMPERATURE = 3,
    PRJ_EOS_NQUANT = 4
};

enum prj_eosvar_quantity {
    PRJ_EOSVAR_PRESSURE = 0,
    PRJ_EOSVAR_TEMPERATURE = 1,
    PRJ_EOSVAR_GAMMA = 2
};

void prj_eos_init(prj_eos *eos);
void prj_eos_rty(prj_eos *eos, double rho, double T, double ye, double *eos_quantities);
void prj_eos_rey(prj_eos *eos, double rho, double eint, double ye, double *eos_quantities);
void prj_eos_fill_block(prj_eos *eos, prj_block *block, double *W);
void prj_eos_fill_active_cells(prj_mesh *mesh, prj_eos *eos, int stage);
void prj_eos_fill_mesh(prj_mesh *mesh, prj_eos *eos, int stage);
void prj_eos_prim2cons(prj_eos *eos, double *W, double *U);
void prj_eos_cons2prim(prj_eos *eos, double *U, double *W);

#endif
