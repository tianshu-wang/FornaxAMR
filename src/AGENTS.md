# PRJ Source Files

## Core Data Structures

- **prj_defs.h** - Constants: `PRJ_BLOCK_SIZE=16`, `PRJ_NGHOST=2`, `PRJ_NDIM=3`. Enums: cons/prim vars, BC types, AMR estimators, EOS kinds.
- **prj_types.h** - Structures: `prj_block`, `prj_mesh`, `prj_sim`, `prj_eos`, `prj_rad`, `prj_grav_mono`, `prj_mpi`, `prj_coord`, `prj_neighbor`, `prj_bc`.

## Modules

- **prj_mesh** - Block allocation, geometry setup, neighbor computation, AMR base marking.
- **prj_amr** - Tag cells, refine/coarsen blocks, enforce 2:1 balance, prolong/restrict.
- **prj_boundary** - Physical BCs (outflow, reflect), fill ghost zones, MPI receive.
- **prj_eos** - Ideal gas or table lookup (rty, rey), prim↔cons conversion, fill mesh.
- **prj_reconstruct** - Slope limiting (minmod-like) for primitive variables.
- **prj_riemann** - HLLE/HLLC solvers, shock detection, flux communication.
- **prj_flux** - Compute interface fluxes, divergence to update.
- **prj_src** - Source terms: geometry, monopole gravity, radiation velocity gradient.
- **prj_timeint** - SSP RK2: calc_dt, stage1, stage2, step.
- **prj_gravity** - Monopole approximation, grid rebuilding, lapse/interp, cache blocks.
- **prj_mpi** - Init, domain decompose, ghost/flux exchange, rebalance, min_dt reduction.
- **prj_io** - Parameter parsing, HDF5 dump/restart writing.
- **prj_radiation** - M1 closure, prim↔cons, energy/momentum update, wavespeeds.
- **prj_rad3_opac** - Opacity lookup tables (kappa, sigma, delta, eta) for radiation.
- **prj_utils** - `prj_fill()` utility.