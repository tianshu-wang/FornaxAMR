# Radiation verification problems

These reproduce Sections 4.1--4.7 of Anninos & Fragile (2020). Each build
selects analytic microphysics at compile time and requires no opacity or EOS
table. The benchmark setup files also set `INELASTIC_SCATTERING := 0`; normal
builds retain the default value of `1`.

```sh
make clean
make machine=macos SETUP_MK=tests/rad/free_streaming.mk
mpirun -n 1 ./prj --param tests/rad/free_streaming.txt
python tests/rad/analyze_free_streaming.py output/dump_00001.h5
```

Replace the setup, parameter, and analysis names with `diffusive_source` or
`picket_fence` for the other cases. All parameter files use a uniform mesh and
one root block along the two symmetric dimensions.

Sections 4.4--4.7 are represented by twelve compile-time configurations:

- `sphere_thin`, `sphere_thick`
- `doppler_15`, `doppler_25`
- `grav_redshift_15`, `grav_redshift_25`, `grav_doppler_15`, `grav_doppler_25`
- `shock_case1` through `shock_case4`

Build any configuration with, for example:

```sh
make clean
make machine=macos SETUP_MK=tests/rad/doppler_25.mk
mpirun -n 1 ./prj --param tests/rad/doppler_25.txt
python tests/rad/analyze_doppler.py output/dump_00010.h5 --groups 25 \
  --previous output/dump_00009.h5 --coarse path/to/doppler_15_final.h5
```

The sphere and redshift calculations opt into exact one-dimensional spherical
shell divergence. The sphere and redshift material profiles are prescribed
reservoirs: their density, temperature, and velocity do not evolve, while their
analytic LTE absorption and emission continue to act on radiation. The gravity
setups use the analytic uniform-sphere metric callback. Shock Case 3 selects
the analytic ideal EOS with gamma 2; all other ideal-EOS builds retain gamma
5/3. Shock tubes additionally use full relativistic hydro/radiation transport
on compile-time frozen flat spacetime (`DYNAMIC_GR=1`,
`FIXED_SPACETIME=1`). These options are off by default.

Analysis entry points are `analyze_sphere.py`, `analyze_doppler.py`,
`analyze_grav_redshift.py`, and `analyze_rad_shock.py`. Each accepts an explicit
dump and exits nonzero when its benchmark criteria fail. Pass the preceding
sphere, Doppler, or redshift dump with `--previous` for stationarity; pass the
matching 15-group final dump to a 25-group redshift analysis with `--coarse` to
enforce convergence. Pass the initial shock dump with `--initial` for the
special-relativistic conservation check.

The analysis scripts require NumPy, h5py, and Matplotlib. The committed Picket
Fence CSV records its DOI and normalization, but its current rows are
provisional Figure 3 marker reconstructions pending access to the original
Su--Olson Case B transport table.
