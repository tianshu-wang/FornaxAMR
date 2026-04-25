# AGENTS.md

## Build and Run

```bash
# Build with a preset (e.g., sedov, cc, ccsn)
make -f Makefile -include setups/<name>.mk

# Run (requires --param)
mpirun -n <N> ./prj --param params/<name>.txt

# Resolution / max_level overrides
./prj --param params/<name>.txt --resolution 32 --max_level 2
```

Presets map to both `.mk` files (build config) and `.txt` files (runtime params):
- `sedov`, `magnetized_sedov`
- `cc`, `magnetized_cc` (gravity)
- `ccsn`, `magnetized_ccsn` (gravity + radiation)
- `shock1d`
- `sedov_offcenter`

ccsn params reference sibling dirs: `../eos_tmp/SFHoEOS_*.dat`, `../s9.0.swbj15.fornax`, `../opacbin.extendT.param`.

## Testing

```bash
make test   # runs each binary in tests/ in sequence
```

## Build Config

- Compiler: `CC` env var (defaults to `mpicc` if found, else `cc`)
- Feature flags via `-include setups/*.mk`: `GRAVITY`, `RADIATION`, `MHD`, `NRAD`, `NEGROUP`, `DUMP_SINGLE_PRECISION`
- `MHD_DEBUG=0` disables debug print statements (enabled by default when MHD=1)
- machine.mk auto-detects HDF5 via brew prefix or pkg-config; override with `HDF5_CFLAGS` / `HDF5_LIBS`
- `DEBUG=1` enables `-O0 -g -DPRJ_DEBUG`

## Code Layout

- `src/` — core solver modules (mesh, mhd, radiation, amr, eos, etc.)
- `problems/` — problem-specific initial conditions (sedov, ccsn, cc, shock1d)
- `params/` — runtime parameter files
- `setups/` — build preset .mk files
- `analysis/visualize.py` — visualization script