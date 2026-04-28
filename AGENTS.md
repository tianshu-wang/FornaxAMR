# PRJ - AMR Hydrodynamics Code

## Build

```bash
# Default build (uses mpicc or cc)
make

# Build with setup (e.g., cc setup for core collapse problems)
make setup=cc

# Debug build
make DEBUG=1
```

## Run

```bash
# Required: --param <file>
./prj --param params/sedov.txt
```

## Problems

- `sedov` - Sedov-Taylor explosion
- `cc` - Core collapse (requires EOS table at `../eos_tmp/SFHoEOS__ye__0.035_0.56_50__logT_-4.793_2.176_500__logrho_-8.699_15.5_500_extend.dat`)
- `ccsn` - Core-collapse supernova
- `shock1d` - 1D shock tube
- `general` - Generic problem

## Test

```bash
make test
```

## Notes

- MPI required 
- HDF5 required for output
- Machine configs in `machines/*.mk` for HPC systems. Use macos.mk for Mac OS.
- Setup configs in `setups/*.mk`.
- Output written to `output/` directory
- Visualization: `python analysis/visualize.py`
