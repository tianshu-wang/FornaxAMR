# Radiation verification problems

These reproduce Sections 4.1--4.3 of Anninos & Fragile (2020). Each build
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

The analysis scripts require NumPy, h5py, and Matplotlib. The committed Picket
Fence CSV records its DOI and normalization, but its current rows are
provisional Figure 3 marker reconstructions pending access to the original
Su--Olson Case B transport table.
