# GR M1 notation mapping: Anninos & Fragile paper to PRJ code

Source paper: `/Users/tianshuw/Downloads/arXiv-2007.12195v1/radm1.tex`,
especially Sections 2, 3.1, 3.2, and 3.3.

Code scope: the dynamic-GR M1 implementation in `src/prj_radiation.c`,
`src/prj_radiation.h`, `src/prj_flux.c`, `src/prj_src.c`,
`src/prj_defs.h`, `src/prj_types.h`, and related Z4c stress-energy assembly.

## Main convention differences

The paper mostly uses units with `c = 1`. The code keeps explicit
`PRJ_CLIGHT` factors. In stress-tensor algebra the code uses `F^i / c` as the
`M^{0i}` leg, while the stored/evolved radiation flux variables are the
physical fluxes `F_i`.

The paper describes radiation-frame primitives `(E_R, u_R^i)`. This code's GR
M1 path does not store those as primitives. It stores Eulerian-normal M1
moments: lab-frame energy `E` and covariant lab-frame flux `F_i`. Fluid-frame
`J`, `H^alpha`, `L^{alpha beta}`, and the pressure tensor are derived when
needed.

Radiation values in `W_rad`, `U_rad`, and dumps are in internal scaled units.
Multiply by `RAD_SCALE` to recover physical CGS radiation energy/flux units.

## Indexing and storage

| Paper notation | Meaning | Code representation |
|---|---|---|
| Species index `i` | Radiation component/species | `field`, `0 <= field < PRJ_NRAD`; flattened with `field * PRJ_NEGROUP + group` |
| Frequency/bin index `(\nu)` or `n` | Frequency group/bin | `group`, `0 <= group < PRJ_NEGROUP`; code uses particle energy grids, despite variable names like `nu_face` |
| `\nu_l`, `\nu_u` | Bin domain limits | `rad->emin[field]`, `rad->emax[field]` |
| `\nu_{g+1/2}` | Bin face/edge | `rad->eedge[field][gf]` |
| `\nu_g` | Bin center | `rad->egroup[field][group]`; also `rad->egroup_erg[field][group]` |
| `\delta\nu_g` | Bin width | `rad->eedge[field][g+1] - rad->eedge[field][g]`; erg version `rad->degroup_erg[field][g]` |
| group-integrated moment | Paper often writes spectral moments; tests also discuss binned integrals | Code stores group-integrated moments. Frequency-flux code temporarily divides by `dnu` with `inv_dnu` to reconstruct spectral face values. |

Radiation variable indices:

| Concept | Full prim/cons index | Radiation-only index | Storage array |
|---|---|---|---|
| Lab energy `E_g` | `PRJ_PRIM_RAD_E(field, group)` / `PRJ_CONS_RAD_E(field, group)` | `PRJ_RAD_PRIM_E(field, group)` / `PRJ_RAD_CONS_E(field, group)` | `block->W_rad`, `block->U_rad`, `u[...]` |
| Covariant flux `F_{1,g}` | `PRJ_PRIM_RAD_F1` / `PRJ_CONS_RAD_F1` | `PRJ_RAD_PRIM_F1` / `PRJ_RAD_CONS_F1` | same |
| Covariant flux `F_{2,g}` | `PRJ_PRIM_RAD_F2` / `PRJ_CONS_RAD_F2` | `PRJ_RAD_PRIM_F2` / `PRJ_RAD_CONS_F2` | same |
| Covariant flux `F_{3,g}` | `PRJ_PRIM_RAD_F3` / `PRJ_CONS_RAD_F3` | `PRJ_RAD_PRIM_F3` / `PRJ_RAD_CONS_F3` | same |

`W_rad` and `U_rad` use radiation-only indices with `WIDX(...)` or
`RADVIDX(...)`. A cell-local conservative vector `u[PRJ_NVAR_CONS]` uses the
full `PRJ_CONS_RAD_*` indices. The macros intentionally alias prim/cons
variable numbers; physical conservative meaning comes from the storage context.

## Fluid and geometry variables

| Paper notation | Meaning | Code representation |
|---|---|---|
| `\rho` | Rest-mass density | `W[PRJ_PRIM_RHO]`; conservative `u[PRJ_CONS_RHO]` |
| `\epsilon` | Specific internal energy | `W[PRJ_PRIM_EINT]`; `eint` |
| `P_\mathrm{gas}` | Gas pressure | EOS output, usually `block->eosvar[EIDX(PRJ_EOSVAR_PRESSURE,...)]` or `eos_q[...]` |
| `Y_e` | Electron fraction | `W[PRJ_PRIM_YE]`; conservative `u[PRJ_CONS_YE] = rho * Ye` |
| `u^\alpha` | Fluid 4-velocity | In GR M1 closure, normal-frame components: `ucon[0] = side->wlor`, `ucon[a+1] = side->wlor * side->vcon[a]` |
| `u_\alpha` | Fluid covariant 4-velocity | `ucov[0] = -side->wlor`, `ucov[a+1] = side->u_cov[a]` |
| `w = \alpha u^0` | Lorentz factor | `side->wlor` |
| `V^i = u^i/u^0` | Transport velocity | Not stored directly in GR M1 closure; code stores physical 3-velocity in `W[PRJ_PRIM_V1+a]`, then `vcon[a] = W[...] / PRJ_CLIGHT` |
| `\alpha` | Lapse | `geom.alpha` |
| `\beta^i` | Shift | `geom.beta[i]` |
| `\gamma_{ij}` | Spatial metric | `geom.gamma[i][j]`, `ctx->gamma[i][j]` |
| `\gamma^{ij}` | Inverse spatial metric | `geom.gamma_inv[i][j]`, `ctx->gamma_inv[i][j]` |
| `g_{\alpha\beta}` | 4-metric | Full coordinate metric in `prj_flux_gr_m1_metric4`; Eulerian normal-frame metric in `prj_rad_grm1_normal_metric_from_ctx` and `prj_src_gr_m1_normal_metric` |
| `\Gamma^\beta_{\alpha\gamma}` | Connection | Expanded into 3+1 terms using `geom.dalpha`, `geom.dbeta`, `geom.dgamma`, `geom.K_dd` in source and frequency-drift paths |
| `h_{\alpha\beta}=g_{\alpha\beta}+u_\alpha u_\beta` | Fluid projection | Built in `prj_rad_grm1_prepare_m3_data` as `hcov[4][4]`; inverse/projection `hcon[4][4]` |

Do not confuse paper `W`/`w` with code array `W`. In the paper excerpt,
`w = alpha u^0` is the Lorentz factor; in the code, `W` usually means the
primitive array.

## Radiation moments and stress tensor

| Paper notation | Meaning | Code representation |
|---|---|---|
| `R^{\alpha\beta}_{(\nu)}` | Spectral radiation stress-energy tensor | Local `Rcon[4][4]`, built by `prj_rad_grm1_build_R(g_cov,g_con,alpha,E,Fcon,Rcon)` |
| `E_{(\nu)}` | Eulerian/lab radiation energy density | Stored `W_rad[WIDX(PRJ_RAD_PRIM_E(field,g),...)]`; cell-local `E`, `Eg[g]`; conservative GR value is `sqrt_gamma * E` in `u[PRJ_CONS_RAD_E]` |
| `F^i_{(\nu)}` | Eulerian/lab contravariant flux | `Fcon[i]`, raised from stored `Fcov[i]`; stress-tensor leg is `Fhat_con[i] = Fcon[i] / PRJ_CLIGHT` |
| `F_{j(\nu)}` | Eulerian/lab covariant flux | Stored/evolved M1 flux variables: `PRJ_RAD_PRIM_F1/F2/F3`, local `Fcov[3]`; conservative GR value is `sqrt_gamma * Fcov[j]` |
| `P^{ij}_{(\nu)}` | Eulerian pressure/stress tensor | `Pcon[3][3]` or `P[3][3]`, taken from the spatial block `Rcon[i+1][j+1]` in the Eulerian normal-frame build |
| `P_{ij}` | Covariant pressure/stress | `Pcov[3][3]` when lowering for Z4c stress-energy |
| `J_{(\nu)}` | Fluid-frame radiation energy density | `m3.J` from `prj_rad_grm1_prepare_m3_data`, computed as `R^{ab} u_a u_b` |
| `H^\alpha_{(\nu)}` | Fluid-frame radiation flux/momentum | `m3.H[4]` from `prj_rad_grm1_prepare_m3_data`; units match energy density because stress-tensor flux legs use `F/c` |
| `L^{\alpha\beta}_{(\nu)}` | Fluid-frame pressure tensor | `m3.L[4][4]` from `prj_rad_grm1_prepare_m3_data` |
| `E_{R(\nu)}`, `u_R^\alpha` | Radiation-rest-frame energy and velocity | Algebraically eliminated inside `prj_rad_grm1_build_R`; not returned or stored |

## M1 closure

| Paper notation | Meaning | Code representation |
|---|---|---|
| `\xi` | Fluid-frame flux factor, `sqrt(H_alpha H^alpha) / J` | `fbar`; computed from `m3.Hnorm2` and `m3.J` after `R^{ab}` is built |
| `\chi(\xi)` | Eddington factor | Analytic Levermore formula `prj_rad_m1_chi_exact`; GR transport wave speeds call it through `prj_flux_gr_m1_chi` |
| Eulerian flux factor `|F|/(cE)` | Not the same as `xi` when fluid moves | Used only to clamp admissible input states before building `R^{ab}` |
| `[P^{ij}]_\mathrm{thin}` | Free-streaming pressure | The finite null limit of `prj_rad_grm1_build_R`, `R^{ab}=R^{0a}R^{0b}/R^{00}` |
| `[P^{ij}]_\mathrm{thick}` | Optically thick pressure | Isotropic radiation-fluid limit of Eq. (4), recovered by the same `build_R` formula |
| `(3 chi - 1)/2` | Thin interpolation weight | `thin_w` |
| `3(1 - chi)/2` | Thick interpolation weight | `thick_w` |
| `P^{ij} = thin_w P_thin + thick_w P_thick` | M1 pressure closure | Represented implicitly by `prj_rad_grm1_build_R`; no implicit `fbar` root solve is active in GR M1 |

## Conservative equations and face fluxes

| Paper notation | Meaning | Code representation |
|---|---|---|
| `mathbf U = [D, cal E, cal S_j, cal R_(nu), cal R_j(nu)]` | Conservative state | Hydro: `u[PRJ_CONS_RHO]`, `u[PRJ_CONS_ETOT]`, `u[PRJ_CONS_MOM1..3]`; radiation: `u[PRJ_CONS_RAD_E/F1/F2/F3(field,g)]` |
| `cal R_(nu) = -sqrt(-g) R^0_0` | Conserved radiation energy | `U[0] = sqrt(gamma) alpha^2 Rcon[0][0]` in `prj_flux_gr_m1_state_flux` |
| `cal R_{j(nu)} = sqrt(-g) R^0_j` | Conserved radiation momentum | `U[1+j] = c alpha sqrt(gamma) R^0_j` via `prj_flux_gr_m1_R_mixed` |
| `-sqrt(-g) R^i_0` | Radiation energy flux through face | Code uses the opposite overall sign from paper Eq. (19): for x1-local face, `F[0] = c alpha^2 sqrt(gamma) Rcon[0][1]` |
| `sqrt(-g) R^i_j` | Radiation momentum flux through face | For x1-local face, `F[1+j] = c^2 alpha sqrt(gamma) R^1_j` |
| HLL face solver | HRSC flux from left/right states | `prj_flux_gr_m1`, using `UL/UR`, `FphysL/FphysR`, `sL/sR`, and an optical-depth diffusion limiter |
| `lambda_\pm` | Characteristic speeds | `prj_flux_gr_m1_wavespeeds`; thin/thick blended using analytic Levermore `chi`; multiplied by `PRJ_CLIGHT` |

The face-flux implementation is written in an x1-local frame inside
`prj_flux_gr_m1_state_flux`; surrounding flux code handles face directions.

## Curvature and frequency-shift sources

| Paper notation | Meaning | Code representation |
|---|---|---|
| `-sqrt(-g) R^alpha_beta Gamma^beta_{0 alpha}` | Radiation energy curvature source | `prj_src_gr_m1_z4c`: `energy_src` built from `Fcon . dalpha / alpha` and `Pcon : K_dd` |
| `+sqrt(-g) R^alpha_beta Gamma^beta_{j alpha}` | Radiation momentum curvature source | `prj_src_gr_m1_z4c`: `mom_src` built from `E*dalpha`, `Fcov*dbeta`, and `Pcon*dgamma` |
| `M^{alpha beta gamma}` | Third moment for frequency shifts | Not materialized; `prj_rad_grm1_freq_drift()` and `prj_rad_grm1_freq_drift_3p1()` contract/project it from local `J`, `H`, `L` |
| `N^{alpha beta gamma}` | Third angular moment closure | Not materialized; thin/thick pieces are evaluated from analytic Levermore `chi(xi)` when moment components or contractions are needed |
| `partial_nu[nu ...]` | Frequency-space advection | `energy_face[gf]`, `momentum_face[gf][a]` in `prj_rad_freq_flux_apply_gr_m1` |
| Energy drift integrand | Paper's energy-space source projected into lab energy equation | `energy_drift` / local `Acon`, then `Acon_spec[g] = Acon / dnu` |
| Momentum drift integrand | Paper's momentum-space source projected into lab momentum equation | `momentum_drift_cov[3]` / local `Mq_cov[g][a]`, then `Mq_spec[g][a] = Mq_cov[g][a] / dnu` |
| Boundary zero frequency flux | Conservative frequency-domain closure | Code sets outer frequency faces to zero by only filling internal `gf = 1..PRJ_NEGROUP-1` |

The GR frequency-shift code uses the Cardall/Endeve/Mezzacappa 3+1 projected
form, not the paper's Shibata contraction written literally. The physical role
is the same: Doppler and gravitational redistribution among energy groups.

## Matter coupling and opacities

| Paper notation | Meaning | Code representation |
|---|---|---|
| `G^\mu_(nu)` / `G_alpha` | Radiation-matter four-force | Not stored as a single explicit vector in the GR M1 path; effects are applied in `prj_rad_gr_m1_matter_update` |
| `kappa^a_(nu)` | Absorption/emission opacity per mass in paper formula | Code `kappa[idx]` is the density-multiplied absorption coefficient used in `dt * c * kappa`; from `rad->absopac` or `block->kappa_cell` |
| `kappa^s_(nu)` | Scattering opacity per mass in paper formula | Code `sigma[idx]` is the density-multiplied scattering coefficient; from `rad->scaopac` or `block->sigma_cell` |
| Scattering transport correction | Not explicit in paper Eq. `Galpha` | Code `delta[idx]`, with effective momentum opacity `kappa + sigma * (1 - delta / 3)` |
| `B_(nu)` | Equilibrium emissivity/source function | Code uses emissivity `eta[idx]` from `rad->emis`; group integrated and scaled for the implicit energy update |
| `rho kappa` in paper source | Interaction coefficient with units `1/length` | Already folded into code `kappa`, `sigma`, `eta` interpolation outputs |
| Matter temperature `T` | Gas temperature | `T`, from EOS `PRJ_EOS_TEMPERATURE` |
| Lepton coupling | Not central in paper's photon/neutrino formula excerpt | Code uses `rad->x_e[field][group]` in energy/lepton residuals |

GR M1 matter coupling flow in code:

1. Clamp lab `E,F_i`.
2. Build `R^{ab}` with `prj_rad_grm1_build_R`.
3. Decompose that tensor to fluid-frame `J,H` with `prj_rad_grm1_prepare_m3_data`.
4. Put `J` and `c H^i` into a temporary non-GR M1 source vector.
5. Apply inelastic, energy, and momentum updates.
6. Convert `dJ,dH` back to lab conserved changes and back-react on hydro
   energy/momentum.

This is not the full paper Newton solve over `(5 + 4 N_B)` primitives. The code
freezes the fluid frame velocity during the GR source update; the source file
also notes that a future fully implicit GR-M1 solve should include velocity in
strongly momentum-coupled cells.

## Quick lookup by common code names

| Code name | Paper symbol/meaning |
|---|---|
| `E`, `Eg[g]` | `E_(nu)` lab-frame group-integrated energy |
| `Fcov[3]` | `F_j(nu)` stored covariant lab flux |
| `Fcon[3]` | `F^i_(nu)` raised lab flux |
| `Rcon[4][4]` | `R^{alpha beta}_{(nu)}` stress tensor |
| `Rcon[0][i+1]` | `R^{0i}_{(nu)}`, equal to the 3+1 Eq. (2) combination and carrying `F^i/c` units |
| `P`, `Pcon`, `Pcon[3][3]` | `P^{ij}_{(nu)}` |
| `Pcov` | `P_{ij}` |
| `m3.J` | `J_(nu)` |
| `m3.H[4]` | `H^alpha_(nu)` |
| `m3.L[4][4]` | `L^{alpha beta}_{(nu)}` |
| `prj_rad_grm1_freq_drift()` | `M^{alpha beta gamma}_{(nu)} u_{beta;gamma}` / third-moment frequency-drift contraction |
| `fbar` | Fluid-frame flux factor `xi` used in `chi(xi)` |
| `chi` in closure functions | Eddington factor `chi` |
| `chi_ext` | Total opacity-like coefficient for the HLL diffusion limiter; not Eddington `chi` |
| `kappa` | Density-multiplied absorption coefficient |
| `sigma` | Density-multiplied scattering coefficient |
| `delta` | Scattering anisotropy correction; appears as `1 - delta/3` |
| `eta` | Group emissivity/source |
| `field` | Radiation species/component index |
| `group`, `g`, `gf` | Energy group index / energy face index |
| `nu_face` | Energy-group face array, `rad->eedge[field]` |
| `dt_geom` | `dt * alpha * sqrt_gamma` multiplier for GR frequency update |
| `rad_rhs` | Explicit radiation RHS, including geometry/curvature source terms |
