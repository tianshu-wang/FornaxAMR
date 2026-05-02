# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**FornaxAMR** — an adaptive mesh refinement (AMR) hydrodynamics code written in C99 with MPI parallelization, designed for astrophysical simulations (core-collapse supernovae, blast waves). Optional features: monopole gravity, M1 radiation transport, MHD.

## Build

```bash
make                        # default build (auto-detects mpicc or cc)
make setup=cc               # build with a setup config (e.g., setups/cc.mk)
make DEBUG=1                # debug build (-g -O0 -DPRJ_DEBUG)
make clean
```

Feature flags (override defaults from `setup.mk`):
```bash
make GRAVITY=1 RADIATION=1 MHD=1
```

Machine-specific configs live in `machines/*.mk` (use `macos.mk` for local dev). Symlink or copy to `machine.mk` before building. Setup configs (`setups/*.mk`) set default feature flags per problem type.

## Run

```bash
./prj --param params/sedov.txt
mpirun -n 4 ./prj --param params/cc.txt
```

The `--param` flag is required. Output goes to `output/`. Visualization: `python analysis/visualize.py`.

## Problems

| Name | Description | Notes |
|------|-------------|-------|
| `sedov` | Sedov-Taylor blast wave | |
| `magnetized_sedov` | Sedov with uniform B-field | needs `MHD=1` |
| `cc` | Core collapse | needs tabular EOS at `../eos_tmp/SFHoEOS__ye__0.035_0.56_50__logT_-4.793_2.176_500__logrho_-8.699_15.5_500_extend.dat` |
| `magnetized_cc` | Core collapse with dipole B-field | |
| `ccsn` | Core-collapse supernova | |
| `magnetized_ccsn` | CCSN with dipole B-field | |
| `shock1d` | 1D shock tube | |
| `general` | Generic/custom problem | |

## Test

```bash
make test    # compiles and runs all tests in tests/
```

## Code Architecture

### Core Data Structures (`src/prj_defs.h`, `src/prj_types.h`)

- `prj_block` — one AMR block: 16³ active cells + 2-layer ghost zones (20³ storage), holds primitive (`W`, `W1`), conservative (`U`), fluxes, and MHD face fields.
- `prj_mesh` — collection of blocks plus AMR parameters, EOS, boundary conditions.
- `prj_sim` — top-level simulation state.
- `prj_eos` — equation of state (ideal gas or tabular).
- Constants: `PRJ_BLOCK_SIZE=16`, `PRJ_NGHOST=2`, `PRJ_NDIM=3`.

### Module Responsibilities (`src/`)

| Module | Role |
|--------|------|
| `prj_mesh` | Block allocation, geometry, neighbor finding (AMR-aware face/edge/corner), base AMR marking |
| `prj_amr` | Cell tagging, refine/coarsen, 2:1 balance, prolongation (fine→coarse) and restriction (coarse→fine) |
| `prj_boundary` | Physical BCs (outflow, reflective), ghost zone filling, MPI receive coordination |
| `prj_eos` | Ideal gas and table-based lookups; prim↔cons conversion; fill EOS vars over mesh |
| `prj_reconstruct` | Slope limiting (minmod-like) for primitive variable reconstruction at interfaces |
| `prj_riemann` | HLLE/HLLC Riemann solvers; shock detection (Minoshima & Miyoshi method) |
| `prj_flux` | Interface flux computation; flux divergence for conservative update |
| `prj_src` | Source terms: geometric (cylindrical/spherical), monopole gravity, radiation velocity gradient |
| `prj_timeint` | SSP RK2 time integration: `calc_dt`, `stage1`, `stage2`, `step` |
| `prj_gravity` | Monopole gravity: grid rebuilding, lapse interpolation, block caching |
| `prj_mpi` | Domain decomposition, ghost/flux exchange, load rebalancing, global min-dt reduction |
| `prj_io` | Parameter file parsing, HDF5 dump and restart I/O |
| `prj_radiation` | M1 closure radiation: prim↔cons, energy/momentum update, wavespeeds |
| `prj_rad3_opac` | Opacity lookup tables (κ, σ, δ, η) |
| `prj_mhd` | Magnetic field evolution, face/edge tracking for AMR divergence-free reconstruction |

### Execution Flow (`src/main.c`)

1. Parse `--param` file → select problem → initialize mesh, EOS, MPI decomposition
2. Set initial conditions via `problems/prj_problem_<name>.c`
3. Main loop per timestep:
   - Fill ghost zones → fill EOS vars → compute CFL dt
   - SSP RK2: reconstruct → Riemann solve → flux divergence → update U
   - Every `amr_interval` steps: tag → refine/coarsen → rebalance MPI
   - Write HDF5 snapshots and restart files at configured intervals

### Adding a New Problem

Create `problems/prj_problem_<name>.c` implementing the problem init function, add it to `SRCS` in the Makefile, and add a matching `params/<name>.txt` parameter file.
