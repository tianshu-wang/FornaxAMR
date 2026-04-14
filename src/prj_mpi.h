#ifndef PRJ_MPI_H
#define PRJ_MPI_H

void prj_mpi_init(int *argc, char ***argv, prj_mpi *mpi);
void prj_mpi_decompose(prj_mesh *mesh);
void prj_mpi_prepare(prj_mesh *mesh, prj_mpi *mpi);
void prj_mpi_exchange_ghosts(prj_mesh *mesh, prj_mpi *mpi, int stage, int fill_kind);
void prj_mpi_exchange_fluxes(prj_mesh *mesh, prj_mpi *mpi);
double prj_mpi_min_dt(double local_dt);
double prj_mpi_global_sum(double local_val);
void prj_mpi_rebalance(prj_mesh *mesh);
int prj_mpi_get_neighbor_rank(const prj_mesh *mesh, int neighbor_block_id);
prj_mpi *prj_mpi_current(void);
void prj_mpi_finalize(void);

#endif
