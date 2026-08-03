# Radiation verification problems

These files document the planned reproductions of Sections 4.1--4.7 of Anninos
& Fragile (2020). The parameter files, Python analyses, and reference data are
retained while the C problems are migrated to the general `EOS=USER` and
`OPAC=USER` interfaces.

```sh
make machine=macos SETUP_MK=tests/rad/free_streaming.mk
```

Every radiation setup currently stops immediately with `Radiation benchmark
pending USER problem migration`; this prevents an incomplete benchmark from
silently compiling with production physics.

Sections 4.4--4.7 are represented by twelve compile-time configurations:

- `sphere_thin`, `sphere_thick`
- `doppler_15`, `doppler_25`
- `grav_redshift_15`, `grav_redshift_25`, `grav_doppler_15`, `grav_doppler_25`
- `shock_case1` through `shock_case4`

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
