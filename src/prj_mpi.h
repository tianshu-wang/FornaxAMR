#ifndef PRJ_MPI_H
#define PRJ_MPI_H

void prj_mpi_barrier(const prj_mpi *mpi);
void prj_mpi_init(int *argc, char ***argv, prj_mpi *mpi);
void prj_mpi_decompose(prj_mesh *mesh, const prj_mpi *mpi);
void prj_mpi_prepare(prj_mesh *mesh, prj_mpi *mpi);
void prj_mpi_exchange_ghosts(prj_mesh *mesh, prj_mpi *mpi, int stage, int fill_kind);
void prj_mpi_exchange_ghosts_and_bf(prj_mesh *mesh, prj_mpi *mpi, int stage, int fill_kind, int use_bf1);
void prj_mpi_post_ghosts_and_bf(prj_mesh *mesh, prj_mpi *mpi, int stage, int fill_kind, int use_bf1);
void prj_mpi_wait_ghosts_and_bf(prj_mesh *mesh, prj_mpi *mpi, int stage, int fill_kind, int use_bf1);
void prj_mpi_exchange_fluxes(prj_mesh *mesh, prj_mpi *mpi);
void prj_mpi_exchange_fluxes_and_emf(prj_mesh *mesh, prj_mpi *mpi);
#if PRJ_MHD
void prj_mpi_exchange_bf(prj_mesh *mesh, prj_mpi *mpi, int use_bf1, int fill_kind);
void prj_mpi_exchange_emf(prj_mesh *mesh, prj_mpi *mpi);
void prj_mpi_exchange_amr_mhd_prolongate_bf(const prj_mesh *mesh, prj_mpi *mpi);
int prj_mpi_amr_mhd_prolongate_bf_one(const prj_mpi *mpi, const prj_block *parent,
    const prj_neighbor *slot, prj_block *child, int child_oct,
    int use_bf1, int dir);
#endif
double prj_mpi_min_dt(const prj_mpi *mpi, double local_dt);
double prj_mpi_global_sum(const prj_mpi *mpi, double local_val);
void prj_mpi_rebalance(prj_mesh *mesh, prj_mpi *mpi);
int prj_mpi_get_neighbor_rank(const prj_mesh *mesh, int neighbor_block_id);
void prj_mpi_finalize(prj_mpi *mpi);

#endif
