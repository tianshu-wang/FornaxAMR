# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build Commands

```bash
# Build the main executable
make all

# Run unit tests
make test

# Clean build artifacts
make clean
```

**Build configuration** is split across three files:
- `machine.mk` — symlink or copy from `machines/<hostname>.mk` to select the target machine
- `setup.mk` — symlink or copy from `setups/<problem>.mk` to select the problem
- `Makefile` — master build logic

**Feature flags** (set in `setups/*.mk` or on the command line):
| Flag | Default | Effect |
|------|---------|--------|
| `GRAVITY=1` | off | Enable monopole gravity |
| `RADIATION=1` | off | Enable multigroup radiation transport |
| `MHD=1` | off | Enable magnetohydrodynamics |
| `NRAD=<n>` | — | Number of radiation energy groups |
| `DEBUG=1` | off | `-O0 -g` instead of `-O3` |
| `BLOCK_SIZE=<n>` | 16 | Cells per block edge |

EOS tables live outside the repo at `../eos_tmp/` (read permissions are configured in `.claude/settings.local.json`).

## Architecture

This is a 3D block-structured AMR (Adaptive Mesh Refinement) astrophysics simulation code written in C99. It targets core-collapse supernovae and related phenomena on HPC clusters using MPI.

### Key data structures (`src/prj_types.h`)

- **`prj_block`** — a single 16³ AMR block with ghost layers; holds primitive variables (ρ, v, P, B), conserved vars, and face fluxes.
- **`prj_mesh`** — octree collection of blocks with parent-child refinement hierarchy and neighbor topology (56 neighbors per block).
- **`prj_sim`** — top-level simulation object; owns the mesh, EOS, radiation state, gravity, BCs, and time parameters.

### Physics modules (`src/`)

| Module | Purpose |
|--------|---------|
| `prj_timeint` | 2-stage Runge-Kutta driver |
| `prj_flux` | Flux calculation using PPM reconstruction |
| `prj_riemann` | HLLE / HLLC / HLLD Riemann solvers |
| `prj_mhd` | MHD: magnetic face values, EMF, constrained transport |
| `prj_radiation` | M1 closure multigroup radiation transport |
| `prj_gravity` | Monopole gravity on 1D radial grid |
| `prj_eos` | Ideal-gas or tabulated EOS (T-ρ-Ye lookups) |
| `prj_amr` | Refinement tagging, prolongation, restriction, 2-to-1 enforcement |
| `prj_boundary` | Ghost-cell fill: physical BCs + MPI exchange |
| `prj_mpi` | Domain decomposition, structured neighbor communication |
| `prj_io` | HDF5 dumps and restarts |

### Time-step loop (in `prj_timeint.c`)

```
For each RK stage:
  1. Fill ghosts   (MPI exchange + physical BCs)
  2. Prims → cons  (EOS lookup)
  3. Reconstruction + Riemann solve → face fluxes
  4. Source terms  (gravity, radiation, geometric)
  5. Advance conserved vars
  6. Update primitive vars
After RK: AMR adapt → I/O (both interval-triggered)
```

### Problem setups (`problems/`)

Each file provides `prj_problem_init()` and optional per-step hooks. Select via `setup.mk`. Available problems: `general`, `sedov`, `sedov_offcenter`, `shock1d`, `cc` (core-collapse), `ccsn`.

Runtime parameters (mesh size, physics toggles, output cadence, etc.) are read from plain-text files in `params/` at startup.

### Compile-time configuration (`src/prj_defs.h`)

Physics features and array dimensions are controlled by `#define` macros. When adding a new physics module or changing array shapes, update `prj_defs.h` and the corresponding `prj_types.h` struct fields together.

### Post-processing

`analysis/visualize.py` — Python/matplotlib script that reads HDF5 output dumps using h5py.
