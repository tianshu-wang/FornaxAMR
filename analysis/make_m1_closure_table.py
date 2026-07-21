#!/usr/bin/env python3
"""Precompute the GR M1 fluid-frame closure table for FornaxAMR.

For the no-radiation-viscosity GR M1 closure the self-consistent fluid-frame
flux factor -- and hence the Eddington factor chi -- is a *universal* function of
three dimensionless scalars, independent of the spatial metric, the radiation
energy scale E, and the opacity:

    f    = |F| / (c E)                       Eulerian flux factor        [0, 1]
    mu   = (F . v) / (|F| |v|)               cos(angle) between F and v  [-1, 1]
    beta = |v| / c                           fluid speed                 [0, 1)

This holds because the spatial metric can always be orthonormalized (covariance)
and the closure is then rotationally invariant, so only |F|, |v| and their
relative angle survive.  Working in that flat orthonormal frame
(gamma_ij = delta_ij, c = 1, E = 1), the pressure is built with NO per-cell
linear solve:

    P^{ij}(fbar) = w_thin(chi) Pthin^{ij}(E,F) + w_thick(chi) Pthick^{ij}(E,F,v)

where the optically-thick tensor Pthick is the closed form of eqs (29)-(31)
(fluid-frame isotropy, valid only without radiation viscosity).  fbar is then
root-found so derived_fbar( P(fbar) ) == fbar.

It tabulates the residual fbar - f, where f is the Eulerian flux factor, and
the resulting Eddington factor chi = chi_closure(fbar) over a
(NF x NMU x NBETA) grid and writes them in a plain-text format that C99 can
parse with fscanf.  Storing fbar - f makes the beta=0 slab exactly small and
interpolates the smooth correction rather than the leading Eulerian flux.

The closure function chi(xi) is pluggable; the equation being solved uses it too.
"Levermore" is provided as the default and matches prj_rad_m1_chi_exact() in the
C code:  chi(f) = (3 + 4 f^2) / (5 + 2 sqrt(4 - 3 f^2)).
"""

import argparse
import sys
import numpy as np


# --------------------------------------------------------------------------
# Closure functions chi(xi).  Add new closures here; the solver is agnostic.
# --------------------------------------------------------------------------
def chi_levermore(f):
    """Levermore (1984) closure; matches prj_rad_m1_chi_exact() in the C code."""
    f = np.clip(np.asarray(f, dtype=float), 0.0, 1.0)
    return (3.0 + 4.0 * f * f) / (5.0 + 2.0 * np.sqrt(4.0 - 3.0 * f * f))


CLOSURES = {
    "Levermore": chi_levermore,
}


# --------------------------------------------------------------------------
# Vectorized closure.  Every array below carries a leading "cell" axis of
# length N so a whole slab of (f, mu) at fixed beta is handled at once.
# --------------------------------------------------------------------------
def thick_pressure_closed(E, Fcon, vcon, wlor):
    """Closed-form optically-thick pressure tensor P_TH^{ij} -- no linear solve.

    Implements eqs (29)-(31) of the Foucart et al. M1 formulation: in the
    optically thick limit the pressure is isotropic in the fluid rest frame,
    L^{ab} = (1/3) J h^{ab}, which lets J and the fluid-frame flux H be solved
    for E, F, v in closed form.  This holds ONLY without radiation viscosity
    (the -(4 lbar/15) sigma term breaks fluid-frame isotropy).

    Flat metric, c = 1: raise/lower are the identity, so F^i v_i, F_j, v_j are
    plain contractions and H^i = H_j.  Fcon is the stress-tensor flux (F/c).
    Shapes: E (N,), Fcon/vcon (N,3), wlor (N,) -> P (N,3,3).
    """
    W2 = wlor * wlor
    Fv = np.einsum("na,na->n", Fcon, vcon)                        # (28): F^i v_i
    # (29): fluid-frame energy density J
    J = 3.0 / (2.0 * W2 + 1.0) * ((2.0 * W2 - 1.0) * E - 2.0 * W2 * Fv)
    # (30): fluid-frame flux H^i  (flat metric: H^i = gamma^{ij} H_j = H_j)
    Hcon = (Fcon / wlor[:, None]
            + (wlor[:, None] * vcon) / (2.0 * W2 + 1.0)[:, None]
            * ((4.0 * W2 + 1.0) * Fv - 4.0 * W2 * E)[:, None])
    # (31): P_TH^{ij} = (4/3) W^2 J v^i v^j + W (H^i v^j + v^i H^j) + (1/3) J g^{ij}
    P = (4.0 / 3.0) * (W2 * J)[:, None, None] * vcon[:, :, None] * vcon[:, None, :]
    P += wlor[:, None, None] * (Hcon[:, :, None] * vcon[:, None, :]
                                + vcon[:, :, None] * Hcon[:, None, :])
    P += (1.0 / 3.0) * J[:, None, None] * np.eye(3)[None, :, :]
    return P


def _derived_fbar(P, E, wlor, u_cov, Fhat_con, J0):
    """Fluid-frame flux factor from the pressure tensor.

    Mirrors prj_rad_gr_m1_derived_fbar (flat metric, gamma_ij = delta_ij).
    Shapes: P (N,3,3); E (N,); wlor (N,); u_cov,Fhat_con (N,3); J0 (N,) -> (N,).
    """
    R0 = -E * wlor + np.einsum("na,na->n", Fhat_con, u_cov)
    Rcon = -wlor[:, None] * Fhat_con + np.einsum("nab,nb->na", P, u_cov)
    J = J0 + np.einsum("nab,na,nb->n", P, u_cov, u_cov)

    num = (wlor * wlor - 1.0) * R0 * R0
    num += -2.0 * wlor * R0 * np.einsum("na,na->n", u_cov, Rcon)
    # (delta_ab + u_a u_b) Rcon_a Rcon_b = |Rcon|^2 + (u . Rcon)^2
    uR = np.einsum("na,na->n", u_cov, Rcon)
    num += np.einsum("na,na->n", Rcon, Rcon) + uR * uR

    num = np.where(np.isfinite(num) & (num >= 0.0), num, 0.0)
    denom = np.abs(J)
    out = np.where(denom > 0.0, np.sqrt(num) / np.maximum(denom, 1e-300), 0.0)
    return np.clip(out, 0.0, 1.0)


def solve_slab(f, mu, beta, chi_fn, n_bisect=80):
    """Solve the closure for a slab of cells sharing beta.

    f, mu are 1-D arrays of equal length N (already broadcast for this beta).
    Returns (fbar, chi) at the self-consistent fluid-frame flux factor.  The
    pressure is the explicit M1 blend of the thin closure and the closed-form
    optically-thick tensor (eqs 29-31) -- no per-cell linear solve.
    """
    N = f.shape[0]
    E = np.ones(N)
    c = 1.0

    # --- geometry / velocity: v along z, F in the x-z plane at angle mu ------
    beta = min(beta, 1.0)
    beta2 = beta * beta
    if beta2 >= 1.0:  # mirror the C clip in prepare_side
        beta2 = 1.0 - 1.0e-12
        beta = np.sqrt(beta2)
    vcon = np.zeros((N, 3))
    vcon[:, 2] = beta
    wlor = np.full(N, 1.0 / np.sqrt(1.0 - beta2))
    u_cov = wlor[:, None] * vcon  # flat metric: v_cov = v_con, u_cov = W v

    # --- flux vector: |F| = f (since c = E = 1), angle mu to v (=z) ----------
    Fcon = np.zeros((N, 3))
    Fcon[:, 0] = f * np.sqrt(np.clip(1.0 - mu * mu, 0.0, 1.0))
    Fcon[:, 2] = f * mu
    F2 = np.einsum("na,na->n", Fcon, Fcon)
    Fhat_con = Fcon / c

    # --- thin pressure Pthin (explicit in E, F) ------------------------------
    # Pthin[a][b] = E F^a F^b / F2   (0 where F2 == 0)
    safeF2 = np.where(F2 > 0.0, F2, 1.0)
    Pthin = (E / safeF2)[:, None, None] * Fcon[:, :, None] * Fcon[:, None, :]
    Pthin = np.where((F2 > 0.0)[:, None, None], Pthin, 0.0)

    # --- closed-form optically-thick pressure (eqs 29-31, no linear solve) ---
    Pthick = thick_pressure_closed(E, Fhat_con, vcon, wlor)

    # --- E/F-only piece needed by derived_fbar: J0 = E W^2 - 2 W^2 (F.v)/c ---
    FdotV = np.einsum("na,na->n", Fhat_con, vcon)  # v_cov = v_con
    J0 = E * wlor * wlor - 2.0 * wlor * wlor * FdotV

    # --- outer root-find: derived_fbar(P(fbar)) - fbar = 0 on [0, 1] ---------
    def residual(fbar):
        chi = chi_fn(fbar)
        w_thin = 0.5 * (3.0 * chi - 1.0)
        w_thick = 1.5 * (1.0 - chi)
        P = w_thin[:, None, None] * Pthin + w_thick[:, None, None] * Pthick
        return _derived_fbar(P, E, wlor, u_cov, Fhat_con, J0) - fbar

    lo = np.zeros(N)
    hi = np.ones(N)
    # The M1 residual is monotone decreasing on [0, 1] (g(0) >= 0, g(1) <= 0):
    # g(mid) > 0 => root to the right (lo = mid), else root to the left (hi = mid).
    # Fully vectorizable across the slab and correct at the g(0)=0 boundary.
    for _ in range(n_bisect):
        mid = 0.5 * (lo + hi)
        gm = residual(mid)
        right = gm > 0.0
        lo = np.where(right, mid, lo)
        hi = np.where(right, hi, mid)
    fbar = 0.5 * (lo + hi)
    return fbar, chi_fn(fbar)


# --------------------------------------------------------------------------
def build_grids(nf, nmu, nbeta_log, beta_min, beta_max):
    f_grid = np.linspace(0.0, 1.0, nf)
    mu_grid = np.linspace(-1.0, 1.0, nmu)
    beta_grid = np.concatenate([[0.0],
                                np.logspace(np.log10(beta_min),
                                            np.log10(beta_max), nbeta_log)])
    return f_grid, mu_grid, beta_grid


def build_table(f_grid, mu_grid, beta_grid, chi_fn, verbose=True):
    nf, nmu, nbeta = len(f_grid), len(mu_grid), len(beta_grid)
    fbar_table = np.empty((nf, nmu, nbeta))
    chi_table = np.empty((nf, nmu, nbeta))
    ff, mm = np.meshgrid(f_grid, mu_grid, indexing="ij")
    ff, mm = ff.ravel(), mm.ravel()
    for ib, beta in enumerate(beta_grid):
        fbar, chi = solve_slab(ff, mm, beta, chi_fn)
        fbar_table[:, :, ib] = fbar.reshape(nf, nmu)
        chi_table[:, :, ib] = chi.reshape(nf, nmu)
        if verbose:
            sys.stderr.write("\r  beta %3d/%3d (=%.4g)" % (ib + 1, nbeta, beta))
            sys.stderr.flush()
    if verbose:
        sys.stderr.write("\n")
    return fbar_table, chi_table


def write_table(path, closure_name, f_grid, mu_grid, beta_grid,
                fbar_table, chi_table):
    """Plain-text, C99-fscanf-friendly.  Layout (whitespace separated tokens):

        magic/version line
        keyword value lines for metadata
        'f_grid'    then nf doubles
        'mu_grid'   then nmu doubles
        'beta_grid' then nbeta doubles
        'fbar_minus_f_data' then nf*nmu*nbeta doubles
        'chi_data'          then nf*nmu*nbeta doubles

        fbar is reconstructed as:
            fbar = f_grid[if] + fbar_minus_f_data[idx]

        Both data arrays are C row-major with beta fastest:
            idx = (if*nmu + imu)*nbeta + ibeta
    """
    nf, nmu, nbeta = chi_table.shape
    if fbar_table.shape != chi_table.shape:
        raise ValueError("fbar_table and chi_table shapes differ")

    def dump(fh, arr):
        flat = np.asarray(arr).ravel()
        for i in range(0, flat.size, 8):
            fh.write(" ".join("%.17g" % v for v in flat[i:i + 8]))
            fh.write("\n")

    fbar_minus_f_table = fbar_table - f_grid[:, None, None]

    with open(path, "w") as fh:
        fh.write("FORNAX_M1_CLOSURE_TABLE 3\n")
        fh.write("closure %s\n" % closure_name)
        fh.write("ndim 3\n")
        fh.write("nf %d\n" % nf)
        fh.write("nmu %d\n" % nmu)
        fh.write("nbeta %d\n" % nbeta)
        fh.write("fbar_storage fbar_minus_f\n")
        fh.write("layout data[(if*nmu+imu)*nbeta+ibeta]  # ibeta fastest\n")
        fh.write("f_grid\n"); dump(fh, f_grid)
        fh.write("mu_grid\n"); dump(fh, mu_grid)
        fh.write("beta_grid\n"); dump(fh, beta_grid)
        fh.write("fbar_minus_f_data\n"); dump(fh, fbar_minus_f_table)
        fh.write("chi_data\n"); dump(fh, chi_table)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--closure", default="Levermore", choices=sorted(CLOSURES))
    ap.add_argument("--nf", type=int, default=100, help="flux-factor points in [0,1]")
    ap.add_argument("--nmu", type=int, default=100, help="mu points in [-1,1]")
    ap.add_argument("--nbeta", type=int, default=100,
                    help="log-spaced beta points (a beta=0 row is added -> nbeta+1)")
    ap.add_argument("--beta-min", type=float, default=1.0e-3)
    ap.add_argument("--beta-max", type=float, default=1.0)
    ap.add_argument("-o", "--output", default="m1_closure_table.txt")
    args = ap.parse_args()

    chi_fn = CLOSURES[args.closure]
    f_grid, mu_grid, beta_grid = build_grids(args.nf, args.nmu, args.nbeta,
                                             args.beta_min, args.beta_max)
    sys.stderr.write("Building %s closure table: nf=%d nmu=%d nbeta=%d "
                     "(1 beta=0 row + %d log-spaced)\n"
                     % (args.closure, len(f_grid), len(mu_grid),
                        len(beta_grid), args.nbeta))
    fbar_table, chi_table = build_table(f_grid, mu_grid, beta_grid, chi_fn)

    # sanity: at beta=0 the fluid frame is the Eulerian frame, so fbar == f and
    # chi == chi(f).
    err_fbar = np.max(np.abs(fbar_table[:, :, 0] - f_grid[:, None]))
    err_chi = np.max(np.abs(chi_table[:, :, 0] - chi_fn(f_grid)[:, None]))
    sys.stderr.write("beta=0 consistency  max|fbar - f| = %.3e\n" % err_fbar)
    sys.stderr.write("beta=0 consistency  max|chi - chi(f)| = %.3e\n" % err_chi)
    sys.stderr.write("fbar range [%.6f, %.6f]\n" %
                     (fbar_table.min(), fbar_table.max()))
    fbar_minus_f_table = fbar_table - f_grid[:, None, None]
    sys.stderr.write("(fbar - f) range [%.6f, %.6f]\n" %
                     (fbar_minus_f_table.min(), fbar_minus_f_table.max()))
    sys.stderr.write("chi range [%.6f, %.6f]\n" %
                     (chi_table.min(), chi_table.max()))

    write_table(args.output, args.closure, f_grid, mu_grid, beta_grid,
                fbar_table, chi_table)
    sys.stderr.write("wrote %s\n" % args.output)


if __name__ == "__main__":
    main()
