# GR primitive and conserved variable reference

This note summarizes the meaning, index convention, and units of the hydro,
MHD, and M1 radiation variables in the full dynamic-GR path. The main code
locations are:

- `src/prj_defs.h`: primitive/conserved slot macros.
- `src/prj_eos.c`: GR primitive-to-conserved and conserved-to-primitive maps.
- `src/prj_radiation.c`: M1 radiation primitive/conserved copies and GR M1
  tensor helpers.

## Storage convention

The primitive state `W` is local and undensitized. In full dynamic GR, the
conserved state `U` or cell-local `u` is densitized by `sqrt(gamma)`:

```text
U[v] = sqrt(gamma) * Uloc[v]
```

Here `gamma` is the determinant of the spatial metric `gamma_ij`.

Radiation values are additionally stored in internal scaled units. For physical
CGS radiation units:

```text
E_or_F_physical = RAD_SCALE * E_or_F_code
```

`RAD_SCALE` defaults to `1e25`.

## Hydro primitives

| slot | meaning | units | index type |
|---|---|---:|---|
| `PRJ_PRIM_RHO` | rest-mass density `rho` | `g cm^-3` | scalar |
| `PRJ_PRIM_V1..V3` | Eulerian normal-frame physical 3-velocity `v^i` | `cm s^-1` | contravariant spatial vector |
| `PRJ_PRIM_EINT` | specific internal energy `epsilon` | `erg g^-1` | scalar |
| `PRJ_PRIM_YE` | electron fraction `Ye` | dimensionless | scalar |

The code forms the dimensionless velocity

```text
beta^i = v^i / c
```

and lowers it with the spatial metric:

```text
beta_i = gamma_ij beta^j
```

## MHD primitives

These slots exist only when `PRJ_MHD=1`.

| slot | meaning | units | index type |
|---|---|---:|---|
| `PRJ_PRIM_B1..B3` | Eulerian magnetic field `B^i` | code magnetic-field units | contravariant spatial vector |

The code uses a magnetic-field normalization where magnetic pressure/energy
density is built from `B^2/2`. Thus `B^2` has pressure or energy-density units
(`erg cm^-3`) in the MHD equations. Dump output writes `B` unscaled, with no
extra `sqrt(4*pi)` conversion.

## Hydro/MHD conserved variables

| slot | local GR meaning before densitization | stored full-GR value | units | index type |
|---|---|---|---:|---|
| `PRJ_CONS_RHO` | `D = rho Wlor` | `sqrt(gamma) D` | `g cm^-3` | scalar density |
| `PRJ_CONS_MOM1..MOM3` | `S_i / c` | `sqrt(gamma) S_i / c` | `g cm^-2 s^-1` | covariant spatial vector density |
| `PRJ_CONS_ETOT` | `tau` | `sqrt(gamma) tau` | `erg cm^-3` | scalar energy density |
| `PRJ_CONS_YE` | `D Ye` | `sqrt(gamma) D Ye` | `g cm^-3` | scalar density |
| `PRJ_CONS_B1..B3` | `B^i` | `sqrt(gamma) B^i` | code magnetic-field units | contravariant spatial vector density |

`Wlor` is the Lorentz factor. The conserved momentum is covariant: the code
computes `S_i` using `beta_i` and stores `S_i / c`. Raising it, when needed,
uses `gamma^ij`.

The GRMHD state helper uses:

```text
h = rho c^2 + rho epsilon + P
S_i = (h Wlor^2 + B^2) beta_i - (B.beta) B_i
tau = (rho epsilon + P) Wlor^2
    + rho c^2 Wlor (Wlor - 1)
    - P + B^2 - 0.5 * ((B.beta)^2 + B^2 / Wlor^2)
```

where `B_i = gamma_ij B^j`.

## M1 radiation primitives

For M1 radiation, primitive and conserved macros intentionally use the same
slot numbers. The physical meaning depends on whether the array is primitive
or conserved.

| slot | primitive meaning | physical CGS units | stored code units | index type |
|---|---|---:|---:|---|
| `PRJ_PRIM_RAD_E(field,g)` | group-integrated Eulerian radiation energy density `E_g` | `erg cm^-3` | `(erg cm^-3) / RAD_SCALE` | scalar |
| `PRJ_PRIM_RAD_F1..F3(field,g)` | group-integrated Eulerian radiation flux `F_{g,i}` | `erg cm^-2 s^-1` | `(erg cm^-2 s^-1) / RAD_SCALE` | covariant spatial vector |

`field` is the radiation species/component index and `g` is the energy-group
index.

The M1 flux is stored covariantly as `F_i`. When the code needs the
contravariant flux, it raises with the inverse spatial metric:

```text
F^i = gamma^ij F_j
```

The GR M1 stress-tensor algebra uses `F^i / c` for the time-space leg, so
`F/c` has the same units as `E`.

## M1 radiation conserved variables

In full dynamic GR, the cell-local conserved radiation state stores the
metric-densitized versions:

| slot | stored full-GR value | physical CGS units before scaling | stored code units |
|---|---|---:|---:|
| `PRJ_CONS_RAD_E(field,g)` | `sqrt(gamma) E_g` | `erg cm^-3` | `sqrt(gamma) * (erg cm^-3) / RAD_SCALE` |
| `PRJ_CONS_RAD_F1..F3(field,g)` | `sqrt(gamma) F_{g,i}` | `erg cm^-2 s^-1` | `sqrt(gamma) * (erg cm^-2 s^-1) / RAD_SCALE` |

These are group-integrated quantities. If a spectral per-particle-energy value
is needed, divide by the energy-bin width:

```text
E_spectral_g   = E_g / degroup_erg[field][g]
F_spectral_g,i = F_g,i / degroup_erg[field][g]
```

`degroup_erg[field][g]` is the group width in erg.

## FSA radiation note

When `PRJ_USE_RADIATION_FSA=1`, radiation uses angular intensity slots instead
of M1 `E/F_i` slots:

```text
PRJ_PRIM_RAD_I(field,g,angle)
PRJ_CONS_RAD_I(field,g,angle)
```

These are angular-cell radiation energy/intensity-like scalar quantities in
the same internal `RAD_SCALE`-scaled convention. In full dynamic GR, the
conserved value is densitized by `sqrt(gamma)`.

## Quick answers

- Primitive velocity `v^i`: contravariant.
- Primitive magnetic field `B^i`: contravariant.
- Conserved momentum `S_i / c`: covariant.
- M1 radiation flux `F_{g,i}`: covariant.
- M1 radiation energy `E_g`: group-integrated energy density,
  `erg cm^-3` physically, divided by `RAD_SCALE` in storage.
- M1 radiation flux `F_{g,i}`: group-integrated energy flux,
  `erg cm^-2 s^-1` physically, divided by `RAD_SCALE` in storage.
