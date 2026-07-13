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

/* Where an EOS lookup was invoked from, so an out-of-range failure can report
 * whether it happened during AMR adaptation or the regular main loop. Passed
 * explicitly through the call chain (no static/global state). */
enum prj_eos_call_ctx {
    PRJ_EOS_CTX_MAIN = 0,
    PRJ_EOS_CTX_AMR = 1
};

typedef struct prj_eos_gr_geom {
    double gamma[3][3];
} prj_eos_gr_geom;

enum prj_eos_gr_status {
    PRJ_EOS_GR_OK = 0,
    PRJ_EOS_GR_NULL_ARG = 1,
    PRJ_EOS_GR_BAD_METRIC = 2,
    PRJ_EOS_GR_BAD_STATE = 3,
    PRJ_EOS_GR_NO_CONVERGE = 4
};

void prj_eos_init(prj_eos *eos, const prj_mpi *mpi);
void prj_eos_rty(prj_eos *eos, double rho, double T, double ye, double *eos_quantities,
    enum prj_eos_call_ctx ctx);
double prj_eos_rty_eint(prj_eos *eos, double rho, double T, double ye,
    double *deint_dlnT, double *deint_dYe, enum prj_eos_call_ctx ctx);
double prj_eos_rty_geteta(prj_eos *eos, double rho, double T, double ye,
    enum prj_eos_call_ctx ctx);
void prj_eos_rey(prj_eos *eos, double rho, double eint, double ye, double *eos_quantities,
    enum prj_eos_call_ctx ctx);
double prj_eos_low_temp_eint(prj_eos *eos, double rho, double ye, enum prj_eos_call_ctx ctx);
void prj_eos_fill_block(prj_eos *eos, prj_block *block, double *W, enum prj_eos_call_ctx ctx);
void prj_eos_fill_active_cells(prj_mesh *mesh, prj_eos *eos, const prj_mpi *mpi, int stage,
    enum prj_eos_call_ctx ctx);
void prj_eos_fill_mesh(prj_mesh *mesh, prj_eos *eos, const prj_mpi *mpi, int stage,
    enum prj_eos_call_ctx ctx);
void prj_eos_fill_ghost_cons(prj_mesh *mesh, prj_eos *eos, const prj_mpi *mpi, int stage,
    enum prj_eos_call_ctx ctx);
void prj_eos_prim2cons(prj_eos *eos, double *W, double *U);
void prj_eos_cons2prim(prj_eos *eos, double *U, double *W);
int prj_eos_gr_prim2cons(prj_eos *eos, const prj_eos_gr_geom *geom,
    const double *W, double *U, enum prj_eos_call_ctx ctx);
int prj_eos_gr_cons2prim(prj_eos *eos, const prj_eos_gr_geom *geom,
    const double *U, double *W, enum prj_eos_call_ctx ctx);
int prj_eos_full_dynamic_gr_enabled(const prj_mesh *mesh);
void prj_eos_cell_prim2cons(prj_eos *eos, const prj_mesh *mesh,
    const prj_block *block, int z4c_stage, int i, int j, int k,
    const double *W, double *U, enum prj_eos_call_ctx ctx);
void prj_eos_cell_cons2prim(prj_eos *eos, const prj_mesh *mesh,
    const prj_block *block, int z4c_stage, int i, int j, int k,
    const double *U, double *W, enum prj_eos_call_ctx ctx);

#endif
