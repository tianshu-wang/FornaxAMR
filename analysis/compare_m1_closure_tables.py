#!/usr/bin/env python3
"""Compare fast approximations for the GR M1 closure fbar table.

The script builds candidate approximations for

    fbar(f, mu, beta)

where f = |F|/(cE), mu = cos(F, v), and beta = |v|/c.  It validates each
candidate against direct root solves and reports whether the max absolute fbar
error reaches the requested target.

The default candidates are speed-first double tables:

  raw      : trilinear interpolation of fbar
  delta_f  : trilinear interpolation of fbar - f
  edge     : trilinear interpolation of fbar - [(1-f) fbar(f=0,beta) + f]
  plusblend: delta_f plus a blended high-resolution mu=+1 face correction

over several velocity coordinates:

  beta, eta=atanh(beta), wbeta=W*beta, logw=log(W)

and over a few axis placements:

  linear : uniform in the physical variable
  zeroN  : clustered toward zero
  edgeN  : f = 1 - (1 - s)^N, clustered toward f=1
  bothN  : clustered toward both f=0 and f=1
  plusN  : mu = 1 - 2 (1 - s)^N, clustered toward mu=1
  edgeN  : symmetric mu endpoint clustering

Axis candidates can be polar triples:

  FAXIS:MUAXIS:VAXIS

or flux-component quadruples:

  cart:XAXIS:YAXIS:VAXIS      with x=f*mu, y=f*sqrt(1-mu^2)
  cart2:XAXIS:YAXIS:VAXIS     with x=f*mu, y=f^2*(1-mu^2)
  cartnorm:XAXIS:YAXIS:VAXIS  with x=f*mu, y=fperp/sqrt(1-x^2)
  thickcart*:...               same y choices, but x is centered at the
                                boosted optically-thick point
  thickpolar:FAXIS:MUAXIS:VAXIS
                              polar table with f centered at
                              f_iso(beta)/mu where reachable
  isopolar:FAXIS:MUAXIS:VAXIS
                              polar table with f centered at f_iso(beta);
                              center1 is normalized f - f_iso(beta)

Chebyshev fits can be included with --include-cheb, but those are arithmetic
heavy at runtime; they are mainly useful as a size/accuracy reference.
"""

import argparse
import os
import sys
from dataclasses import dataclass

import numpy as np
from numpy.polynomial import chebyshev as cheb


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

sys.dont_write_bytecode = True

import make_m1_closure_table as closure_gen  # noqa: E402


@dataclass(frozen=True)
class GridSpec:
    nf: int
    nmu: int
    nq: int


@dataclass(frozen=True)
class AxisCandidate:
    mode: str
    xcoord: str
    ycoord: str
    qcoord: str


@dataclass
class CandidateResult:
    method: str
    axes: str
    spec: str
    bytes_used: int
    corner_loads: str
    max_err: float
    p99_err: float
    rms_err: float
    passed: bool
    worst_f: float
    worst_mu: float
    worst_beta: float
    worst_exact: float
    worst_pred: float


def parse_csv(text):
    return [item.strip() for item in text.split(",") if item.strip()]


def parse_dims(text):
    specs = []
    for item in parse_csv(text):
        parts = item.lower().split("x")
        if len(parts) != 3:
            raise ValueError("dimension spec must look like NFxNMUxNQ")
        nf, nmu, nq = (int(p) for p in parts)
        if nf < 2 or nmu < 2 or nq < 2:
            raise ValueError("all dimensions must be >= 2")
        specs.append(GridSpec(nf, nmu, nq))
    return specs


def parse_triplet(text):
    parts = parse_csv(text)
    if len(parts) != 3:
        raise ValueError("expected three comma-separated integers")
    vals = tuple(int(p) for p in parts)
    if any(v < 1 for v in vals):
        raise ValueError("all entries must be positive")
    return vals


def parse_axis_candidates(text):
    candidates = []
    for item in parse_csv(text):
        parts = item.split(":")
        if len(parts) == 3:
            candidates.append(AxisCandidate("polar", *(p.strip() for p in parts)))
        elif len(parts) == 4:
            mode, xcoord, ycoord, qcoord = (p.strip() for p in parts)
            if mode not in ("polar", "isopolar", "thickpolar",
                            "cart", "cart2", "cartnorm",
                            "thickcart", "thickcart2", "thickcartnorm"):
                raise ValueError("unknown axis candidate mode %r" % mode)
            candidates.append(AxisCandidate(mode, xcoord, ycoord, qcoord))
        else:
            raise ValueError("axis candidates must look like FAXIS:MUAXIS:VAXIS "
                             "or MODE:XAXIS:YAXIS:VAXIS")
    return candidates


def axis_power(name, prefix):
    if not name.startswith(prefix):
        return None
    power = int(name[len(prefix):])
    if power < 1:
        raise ValueError("%s power must be positive" % prefix)
    return power


def beta_clip(beta):
    return np.clip(beta, 0.0, 1.0 - 1.0e-12)


def beta_max_clip(beta_max):
    if beta_max <= 0.0:
        raise ValueError("beta_max must be positive")
    return min(beta_max, 1.0 - 1.0e-12)


def beta_to_q(beta, coord):
    beta = beta_clip(np.asarray(beta, dtype=float))
    if coord == "beta":
        return beta
    if coord == "eta":
        return np.arctanh(beta)
    if coord == "wbeta":
        return beta / np.sqrt(np.maximum(1.0 - beta * beta, 1.0e-300))
    if coord == "logw":
        return -0.5 * np.log1p(-beta * beta)
    raise ValueError("unknown velocity coordinate %r" % coord)


def q_to_beta(q, coord):
    q = np.maximum(np.asarray(q, dtype=float), 0.0)
    if coord == "beta":
        return beta_clip(q)
    if coord == "eta":
        return beta_clip(np.tanh(q))
    if coord == "wbeta":
        return beta_clip(q / np.sqrt(1.0 + q * q))
    if coord == "logw":
        return beta_clip(np.sqrt(np.maximum(0.0, 1.0 - np.exp(-2.0 * q))))
    raise ValueError("unknown velocity coordinate %r" % coord)


def thick_lab_flux(beta):
    beta = beta_clip(np.asarray(beta, dtype=float))
    return 4.0 * beta / (3.0 + beta * beta)


def center_to_s(x, center, coord):
    x = np.clip(np.asarray(x, dtype=float), -1.0, 1.0)
    center = np.clip(np.asarray(center, dtype=float), -1.0 + 1.0e-14,
                     1.0 - 1.0e-14)
    power = axis_power(coord, "center")
    if power is None:
        return mu_to_s(x, coord)

    left = x < center
    s = np.empty(np.broadcast_shapes(np.shape(x), np.shape(center)))
    xb = np.broadcast_to(x, s.shape)
    cb = np.broadcast_to(center, s.shape)
    left_dist = np.clip((cb - xb) / np.maximum(cb + 1.0, 1.0e-300),
                        0.0, 1.0)
    right_dist = np.clip((xb - cb) / np.maximum(1.0 - cb, 1.0e-300),
                         0.0, 1.0)
    s[left] = 0.5 * (1.0 - np.power(left_dist[left], 1.0 / power))
    s[~left] = 0.5 * (1.0 + np.power(right_dist[~left], 1.0 / power))
    return s


def s_to_center(s, center, coord):
    s = np.clip(np.asarray(s, dtype=float), 0.0, 1.0)
    center = np.clip(np.asarray(center, dtype=float), -1.0 + 1.0e-14,
                     1.0 - 1.0e-14)
    power = axis_power(coord, "center")
    if power is None:
        return s_to_mu(s, coord)

    left = s < 0.5
    out = np.empty(np.broadcast_shapes(np.shape(s), np.shape(center)))
    sb = np.broadcast_to(s, out.shape)
    cb = np.broadcast_to(center, out.shape)
    out[left] = cb[left] - (cb[left] + 1.0) * np.power(1.0 - 2.0 * sb[left],
                                                       power)
    out[~left] = cb[~left] + (1.0 - cb[~left]) * np.power(2.0 * sb[~left] - 1.0,
                                                          power)
    return np.clip(out, -1.0, 1.0)


def f_to_s(f, coord):
    f = np.clip(np.asarray(f, dtype=float), 0.0, 1.0)
    if coord == "linear":
        return f
    power = axis_power(coord, "zero")
    if power is not None:
        return np.power(f, 1.0 / power)
    power = axis_power(coord, "edge")
    if power is not None:
        return 1.0 - np.power(1.0 - f, 1.0 / power)
    power = axis_power(coord, "both")
    if power is not None:
        x = 2.0 * f - 1.0
        signed_s = np.sign(x) * (1.0 - np.power(1.0 - np.abs(x),
                                                 1.0 / power))
        return 0.5 * (signed_s + 1.0)
    raise ValueError("unknown f-axis coordinate %r" % coord)


def s_to_f(s, coord):
    s = np.clip(np.asarray(s, dtype=float), 0.0, 1.0)
    if coord == "linear":
        return s
    power = axis_power(coord, "zero")
    if power is not None:
        return np.power(s, power)
    power = axis_power(coord, "edge")
    if power is not None:
        return 1.0 - np.power(1.0 - s, power)
    power = axis_power(coord, "both")
    if power is not None:
        x = 2.0 * s - 1.0
        return 0.5 * (np.sign(x) *
                      (1.0 - np.power(1.0 - np.abs(x), power)) + 1.0)
    raise ValueError("unknown f-axis coordinate %r" % coord)


def mu_to_s(mu, coord):
    mu = np.clip(np.asarray(mu, dtype=float), -1.0, 1.0)
    if coord == "linear":
        return 0.5 * (mu + 1.0)
    power = axis_power(coord, "plus")
    if power is not None:
        return 1.0 - np.power(0.5 * (1.0 - mu), 1.0 / power)
    power = axis_power(coord, "edge")
    if power is not None:
        x = np.abs(mu)
        signed_s = np.sign(mu) * (1.0 - np.power(1.0 - x, 1.0 / power))
        return 0.5 * (signed_s + 1.0)
    raise ValueError("unknown mu-axis coordinate %r" % coord)


def s_to_mu(s, coord):
    s = np.clip(np.asarray(s, dtype=float), 0.0, 1.0)
    if coord == "linear":
        return 2.0 * s - 1.0
    power = axis_power(coord, "plus")
    if power is not None:
        return 1.0 - 2.0 * np.power(1.0 - s, power)
    power = axis_power(coord, "edge")
    if power is not None:
        x = 2.0 * s - 1.0
        return np.sign(x) * (1.0 - np.power(1.0 - np.abs(x), power))
    raise ValueError("unknown mu-axis coordinate %r" % coord)


def thickpolar_center_f(mu, beta):
    mu = np.asarray(mu, dtype=float)
    beta = np.asarray(beta, dtype=float)
    f0 = np.where(mu > 0.0, thick_lab_flux(beta) / np.maximum(mu, 1.0e-300),
                  1.0)
    return np.clip(f0, 0.0, 1.0)


def f_to_s_thickpolar(f, mu, beta, coord):
    z = 2.0 * np.clip(np.asarray(f, dtype=float), 0.0, 1.0) - 1.0
    center = 2.0 * thickpolar_center_f(mu, beta) - 1.0
    return center_to_s(z, center, coord)


def s_to_f_thickpolar(s, mu, beta, coord):
    center = 2.0 * thickpolar_center_f(mu, beta) - 1.0
    z = s_to_center(s, center, coord)
    return 0.5 * (z + 1.0)


def f_to_s_isopolar(f, beta, coord):
    z = 2.0 * np.clip(np.asarray(f, dtype=float), 0.0, 1.0) - 1.0
    center = 2.0 * thick_lab_flux(beta) - 1.0
    return center_to_s(z, center, coord)


def s_to_f_isopolar(s, beta, coord):
    center = 2.0 * thick_lab_flux(beta) - 1.0
    z = s_to_center(s, center, coord)
    return 0.5 * (z + 1.0)


def exact_fbar_points(f, mu, beta, chi_fn, n_bisect):
    """Vectorized direct solve for arbitrary beta at each point."""
    f = np.clip(np.asarray(f, dtype=float), 0.0, 1.0)
    mu = np.clip(np.asarray(mu, dtype=float), -1.0, 1.0)
    beta = beta_clip(np.asarray(beta, dtype=float))
    if f.shape != mu.shape or f.shape != beta.shape:
        raise ValueError("f, mu, beta arrays must have the same shape")

    n = f.size
    e = np.ones(n)
    beta2 = beta * beta
    wlor = 1.0 / np.sqrt(1.0 - beta2)

    vcon = np.zeros((n, 3))
    vcon[:, 2] = beta
    u_cov = wlor[:, None] * vcon

    fcon = np.zeros((n, 3))
    fcon[:, 0] = f * np.sqrt(np.clip(1.0 - mu * mu, 0.0, 1.0))
    fcon[:, 2] = f * mu
    f2 = np.einsum("na,na->n", fcon, fcon)
    fhat_con = fcon

    safef2 = np.where(f2 > 0.0, f2, 1.0)
    pthin = (e / safef2)[:, None, None] * fcon[:, :, None] * fcon[:, None, :]
    pthin = np.where((f2 > 0.0)[:, None, None], pthin, 0.0)
    pthick = closure_gen.thick_pressure_closed(e, fhat_con, vcon, wlor)

    fdotv = np.einsum("na,na->n", fhat_con, vcon)
    j0 = e * wlor * wlor - 2.0 * wlor * wlor * fdotv

    def residual(fbar):
        chi = chi_fn(fbar)
        w_thin = 0.5 * (3.0 * chi - 1.0)
        w_thick = 1.5 * (1.0 - chi)
        p = w_thin[:, None, None] * pthin + w_thick[:, None, None] * pthick
        return closure_gen._derived_fbar(p, e, wlor, u_cov, fhat_con, j0) - fbar

    lo = np.zeros(n)
    hi = np.ones(n)
    for _ in range(n_bisect):
        mid = 0.5 * (lo + hi)
        gm = residual(mid)
        right = gm > 0.0
        lo = np.where(right, mid, lo)
        hi = np.where(right, hi, mid)
    return 0.5 * (lo + hi)


def exact_slab(f_grid, mu_grid, beta, chi_fn, n_bisect):
    ff, mm = np.meshgrid(f_grid, mu_grid, indexing="ij")
    bb = np.full(ff.size, beta)
    return exact_fbar_points(ff.ravel(), mm.ravel(), bb, chi_fn,
                             n_bisect).reshape(len(f_grid), len(mu_grid))


def flux_components(f, mu):
    f = np.clip(np.asarray(f, dtype=float), 0.0, 1.0)
    mu = np.clip(np.asarray(mu, dtype=float), -1.0, 1.0)
    x = f * mu
    y = f * np.sqrt(np.clip(1.0 - mu * mu, 0.0, 1.0))
    return x, y


def component_kind(mode):
    if mode.startswith("thick"):
        return mode[len("thick"):]
    return mode


def component_nodes_to_f_mu(mode, x, ycoord_value):
    x = np.clip(np.asarray(x, dtype=float), -1.0, 1.0)
    ycoord_value = np.clip(np.asarray(ycoord_value, dtype=float), 0.0, 1.0)
    kind = component_kind(mode)
    if kind == "cart":
        y = ycoord_value
    elif kind == "cart2":
        y = np.sqrt(ycoord_value)
    elif kind == "cartnorm":
        y = ycoord_value * np.sqrt(np.clip(1.0 - x * x, 0.0, 1.0))
    else:
        raise ValueError("unknown component mode %r" % mode)

    radius = np.sqrt(x * x + y * y)
    f = np.minimum(radius, 1.0)
    mu = np.where(radius > 0.0, x / np.maximum(radius, 1.0e-300), 0.0)
    return f, np.clip(mu, -1.0, 1.0)


def component_interp_scaled_coords(mode, xcoord, ycoord, f, mu, beta):
    x, y = flux_components(f, mu)
    kind = component_kind(mode)
    if kind == "cart":
        ytab = y
    elif kind == "cart2":
        ytab = y * y
    elif kind == "cartnorm":
        denom = np.sqrt(np.clip(1.0 - x * x, 0.0, 1.0))
        ytab = np.where(denom > 0.0, y / np.maximum(denom, 1.0e-300), 0.0)
    else:
        raise ValueError("unknown component mode %r" % mode)

    if mode.startswith("thick"):
        sx = center_to_s(x, thick_lab_flux(beta), xcoord)
    else:
        sx = mu_to_s(x, xcoord)
    return sx, f_to_s(np.clip(ytab, 0.0, 1.0), ycoord)


def make_validation_points(n_random, beta_max, seed):
    rng = np.random.default_rng(seed)
    n_random = max(0, int(n_random))
    beta_max = beta_max_clip(beta_max)

    f = rng.random(n_random)
    mu = rng.uniform(-1.0, 1.0, n_random)
    beta = np.empty(n_random)

    n1 = n_random // 3
    n2 = n_random // 3
    n3 = n_random - n1 - n2
    beta[:n1] = rng.uniform(0.0, beta_max, n1)
    umax = beta_to_q(beta_max, "wbeta")
    beta[n1:n1 + n2] = q_to_beta(rng.uniform(0.0, umax, n2), "wbeta")
    etamax = beta_to_q(beta_max, "eta")
    beta[n1 + n2:] = q_to_beta(rng.uniform(0.0, etamax, n3), "eta")

    f_edge = np.array([0.0, 1.0e-10, 1.0e-8, 1.0e-4, 0.05, 0.25,
                       0.5, 0.75, 0.9, 0.95, 0.975, 0.99,
                       0.995, 0.999, 1.0 - 1.0e-8, 1.0])
    mu_edge = np.array([-1.0, -0.75, -0.25, 0.0, 0.25, 0.75,
                        0.9, 0.99, 0.999, 1.0])
    fixed_beta_edge = np.array([0.0, 1.0e-8, 1.0e-4, 1.0e-3, 1.0e-2,
                                0.05, 0.1, 0.25, 0.5])
    scaled_beta_edge = np.array([0.75 * beta_max, 0.9 * beta_max,
                                 0.95 * beta_max, 0.99 * beta_max,
                                 beta_max])
    beta_edge = np.unique(np.concatenate([
        fixed_beta_edge[fixed_beta_edge <= beta_max],
        scaled_beta_edge,
    ]))
    ff, mm, bb = np.meshgrid(f_edge, mu_edge, beta_edge, indexing="ij")

    f = np.concatenate([f, ff.ravel()])
    mu = np.concatenate([mu, mm.ravel()])
    beta = np.concatenate([beta, bb.ravel()])
    return f, mu, beta


def build_trilinear_grid(spec, fcoord, mucoord, qcoord, beta_max, chi_fn,
                         n_bisect):
    f_grid = s_to_f(np.linspace(0.0, 1.0, spec.nf), fcoord)
    mu_grid = s_to_mu(np.linspace(0.0, 1.0, spec.nmu), mucoord)
    qmax = float(beta_to_q(beta_max, qcoord))
    q_grid = np.linspace(0.0, qmax, spec.nq)
    beta_grid = q_to_beta(q_grid, qcoord)

    table = np.empty((spec.nf, spec.nmu, spec.nq))
    for iq, beta in enumerate(beta_grid):
        table[:, :, iq] = exact_slab(f_grid, mu_grid, beta, chi_fn, n_bisect)
    return f_grid, mu_grid, q_grid, beta_grid, table


def build_isopolar_grid(spec, fcoord, mucoord, qcoord, beta_max, chi_fn,
                        n_bisect):
    s_grid = np.linspace(0.0, 1.0, spec.nf)
    mu_grid = s_to_mu(np.linspace(0.0, 1.0, spec.nmu), mucoord)
    qmax = float(beta_to_q(beta_max, qcoord))
    q_grid = np.linspace(0.0, qmax, spec.nq)
    beta_grid = q_to_beta(q_grid, qcoord)

    table = np.empty((spec.nf, spec.nmu, spec.nq))
    f_node_table = np.empty_like(table)
    for iq, beta in enumerate(beta_grid):
        f_grid = s_to_f_isopolar(s_grid, beta, fcoord)
        table[:, :, iq] = exact_slab(f_grid, mu_grid, beta, chi_fn, n_bisect)
        f_node_table[:, :, iq] = f_grid[:, None]
    return None, mu_grid, q_grid, beta_grid, f_node_table, table


def build_thickpolar_grid(spec, fcoord, mucoord, qcoord, beta_max, chi_fn,
                          n_bisect):
    s_grid = np.linspace(0.0, 1.0, spec.nf)
    mu_grid = s_to_mu(np.linspace(0.0, 1.0, spec.nmu), mucoord)
    qmax = float(beta_to_q(beta_max, qcoord))
    q_grid = np.linspace(0.0, qmax, spec.nq)
    beta_grid = q_to_beta(q_grid, qcoord)

    table = np.empty((spec.nf, spec.nmu, spec.nq))
    f_node_table = np.empty_like(table)
    ss, mm = np.meshgrid(s_grid, mu_grid, indexing="ij")
    for iq, beta in enumerate(beta_grid):
        f_nodes = s_to_f_thickpolar(ss, mm, beta, fcoord)
        bb = np.full(f_nodes.size, beta)
        table[:, :, iq] = exact_fbar_points(f_nodes.ravel(), mm.ravel(),
                                            bb, chi_fn, n_bisect).reshape(
                                                spec.nf, spec.nmu)
        f_node_table[:, :, iq] = f_nodes
    return None, mu_grid, q_grid, beta_grid, f_node_table, table


def build_component_grid(spec, mode, xcoord, ycoord, qcoord, beta_max, chi_fn,
                         n_bisect):
    y_grid = s_to_f(np.linspace(0.0, 1.0, spec.nmu), ycoord)
    qmax = float(beta_to_q(beta_max, qcoord))
    q_grid = np.linspace(0.0, qmax, spec.nq)
    beta_grid = q_to_beta(q_grid, qcoord)

    table = np.empty((spec.nf, spec.nmu, spec.nq))
    f_node_table = np.empty_like(table)
    s_grid = np.linspace(0.0, 1.0, spec.nf)
    for iq, beta in enumerate(beta_grid):
        if mode.startswith("thick"):
            x_grid = s_to_center(s_grid, thick_lab_flux(beta), xcoord)
        else:
            x_grid = s_to_mu(s_grid, xcoord)
        xx, yy = np.meshgrid(x_grid, y_grid, indexing="ij")
        f_nodes, mu_nodes = component_nodes_to_f_mu(mode, xx, yy)
        bb = np.full(f_nodes.size, beta)
        table[:, :, iq] = exact_fbar_points(f_nodes.ravel(), mu_nodes.ravel(),
                                            bb, chi_fn, n_bisect).reshape(
                                                spec.nf, spec.nmu)
        f_node_table[:, :, iq] = f_nodes
    return None, y_grid, q_grid, beta_grid, f_node_table, table


def interp_axis_uniform(x, xmin, xmax, n):
    if xmax == xmin:
        scaled = np.zeros_like(x, dtype=float)
    else:
        scaled = (np.asarray(x, dtype=float) - xmin) * ((n - 1) / (xmax - xmin))
    scaled = np.clip(scaled, 0.0, n - 1)
    i0 = np.floor(scaled).astype(np.int64)
    i0 = np.minimum(i0, n - 2)
    t = scaled - i0
    return i0, i0 + 1, t


def interp3_scaled(table, qmax, sx, sy, q):
    nf, nmu, nq = table.shape
    i0, i1, tx = interp_axis_uniform(sx, 0.0, 1.0, nf)
    j0, j1, ty = interp_axis_uniform(sy, 0.0, 1.0, nmu)
    k0, k1, tz = interp_axis_uniform(q, 0.0, qmax, nq)

    c000 = table[i0, j0, k0]
    c100 = table[i1, j0, k0]
    c010 = table[i0, j1, k0]
    c110 = table[i1, j1, k0]
    c001 = table[i0, j0, k1]
    c101 = table[i1, j0, k1]
    c011 = table[i0, j1, k1]
    c111 = table[i1, j1, k1]

    c00 = c000 * (1.0 - tx) + c100 * tx
    c10 = c010 * (1.0 - tx) + c110 * tx
    c01 = c001 * (1.0 - tx) + c101 * tx
    c11 = c011 * (1.0 - tx) + c111 * tx
    c0 = c00 * (1.0 - ty) + c10 * ty
    c1 = c01 * (1.0 - ty) + c11 * ty
    return c0 * (1.0 - tz) + c1 * tz


def interp3_uniform(table, qmax, fcoord, mucoord, qcoord, f, mu, beta):
    q = beta_to_q(beta, qcoord)
    sf = f_to_s(f, fcoord)
    smu = mu_to_s(mu, mucoord)
    return interp3_scaled(table, qmax, sf, smu, q)


def interp1_uniform(values, qmax, qcoord, beta):
    q = beta_to_q(beta, qcoord)
    n = len(values)
    k0, k1, t = interp_axis_uniform(q, 0.0, qmax, n)
    return values[k0] * (1.0 - t) + values[k1] * t


def interp2_scaled(table, sx, q, qmax):
    nf, nq = table.shape
    i0, i1, tx = interp_axis_uniform(sx, 0.0, 1.0, nf)
    k0, k1, tz = interp_axis_uniform(q, 0.0, qmax, nq)

    c00 = table[i0, k0]
    c10 = table[i1, k0]
    c01 = table[i0, k1]
    c11 = table[i1, k1]
    c0 = c00 * (1.0 - tx) + c10 * tx
    c1 = c01 * (1.0 - tx) + c11 * tx
    return c0 * (1.0 - tz) + c1 * tz


def face_s_to_f(s, beta, coord):
    center = 2.0 * thick_lab_flux(beta) - 1.0
    return 0.5 * (s_to_center(s, center, coord) + 1.0)


def face_f_to_s(f, beta, coord):
    center = 2.0 * thick_lab_flux(beta) - 1.0
    z = 2.0 * np.clip(np.asarray(f, dtype=float), 0.0, 1.0) - 1.0
    return center_to_s(z, center, coord)


def build_plus_face(nf, nq, qcoord, beta_max, chi_fn, n_bisect, fcoord):
    s_grid = np.linspace(0.0, 1.0, nf)
    qmax = float(beta_to_q(beta_max, qcoord))
    q_grid = np.linspace(0.0, qmax, nq)
    beta_grid = q_to_beta(q_grid, qcoord)
    table = np.empty((nf, nq))
    for iq, beta in enumerate(beta_grid):
        f_nodes = face_s_to_f(s_grid, beta, fcoord)
        table[:, iq] = exact_fbar_points(
            f_nodes, np.ones_like(f_nodes), np.full_like(f_nodes, beta),
            chi_fn, n_bisect)
    return q_grid, table


def plus_face_prediction(face_table, q_grid, qcoord, fcoord, f, beta):
    qmax = q_grid[-1]
    q = beta_to_q(beta, qcoord)
    sx = face_f_to_s(f, beta, fcoord)
    return interp2_scaled(face_table, sx, q, qmax)


def polar_plus_row_prediction(fbar_table, q_grid, fcoord, qcoord, f, beta):
    qmax = q_grid[-1]
    nf, _, _ = fbar_table.shape
    sf = f_to_s(f, fcoord)
    q = beta_to_q(beta, qcoord)
    plus_row = fbar_table[:, -1, :]
    return interp2_scaled(plus_row, sf, q, qmax)


def plus_blend_prediction(fbar_table, q_grid, face_table, face_q_grid,
                          fcoord, mucoord, qcoord, face_fcoord, f, mu, beta,
                          width, power, xwidth, xpower):
    base = candidate_prediction("delta_f", fbar_table, q_grid, fcoord, mucoord,
                                qcoord, f, mu, beta)
    face = plus_face_prediction(face_table, face_q_grid, qcoord, face_fcoord,
                                f, beta)
    coarse_face = polar_plus_row_prediction(fbar_table, q_grid, fcoord, qcoord,
                                            f, beta)
    x, fperp = flux_components(f, mu)
    w_perp = 1.0 / (1.0 + np.power(fperp / max(width, 1.0e-300), power))
    xdist = np.abs(x - thick_lab_flux(beta))
    w_para = 1.0 / (1.0 + np.power(xdist / max(xwidth, 1.0e-300), xpower))
    weight = w_perp * w_para
    weight = np.where(mu > 0.0, weight, 0.0)
    return np.clip(base + weight * (face - coarse_face), 0.0, 1.0)


def candidate_prediction(method, fbar_table, q_grid, fcoord, mucoord, qcoord,
                         f, mu, beta):
    qmax = q_grid[-1]
    nf = fbar_table.shape[0]
    f_grid = s_to_f(np.linspace(0.0, 1.0, nf), fcoord)
    f_node = f_grid[:, None, None]

    if method == "raw":
        pred = interp3_uniform(fbar_table, qmax, fcoord, mucoord, qcoord,
                               f, mu, beta)
    elif method == "delta_f":
        pred = f + interp3_uniform(fbar_table - f_node, qmax, fcoord, mucoord,
                                   qcoord, f, mu, beta)
    elif method == "edge":
        fbar0 = fbar_table[0, 0, :].copy()
        base_nodes = (1.0 - f_node) * fbar0[None, None, :] + f_node
        residual = fbar_table - base_nodes
        base = (1.0 - f) * interp1_uniform(fbar0, qmax, qcoord, beta) + f
        pred = base + interp3_uniform(residual, qmax, fcoord, mucoord, qcoord,
                                      f, mu, beta)
    else:
        raise ValueError("unknown table method %r" % method)
    return np.clip(pred, 0.0, 1.0)


def isopolar_candidate_prediction(method, fbar_table, f_node, q_grid, fcoord,
                                  mucoord, qcoord, f, mu, beta):
    qmax = q_grid[-1]
    sf = f_to_s_isopolar(f, beta, fcoord)
    smu = mu_to_s(mu, mucoord)
    q = beta_to_q(beta, qcoord)

    if method == "raw":
        pred = interp3_scaled(fbar_table, qmax, sf, smu, q)
    elif method == "delta_f":
        pred = f + interp3_scaled(fbar_table - f_node, qmax, sf, smu, q)
    else:
        raise ValueError("isopolar supports raw and delta_f, not %r" %
                         method)
    return np.clip(pred, 0.0, 1.0)


def thickpolar_candidate_prediction(method, fbar_table, f_node, q_grid, fcoord,
                                    mucoord, qcoord, f, mu, beta):
    qmax = q_grid[-1]
    sf = f_to_s_thickpolar(f, mu, beta, fcoord)
    smu = mu_to_s(mu, mucoord)
    q = beta_to_q(beta, qcoord)

    if method == "raw":
        pred = interp3_scaled(fbar_table, qmax, sf, smu, q)
    elif method == "delta_f":
        pred = f + interp3_scaled(fbar_table - f_node, qmax, sf, smu, q)
    else:
        raise ValueError("thickpolar supports raw and delta_f, not %r" %
                         method)
    return np.clip(pred, 0.0, 1.0)


def component_candidate_prediction(method, fbar_table, f_node, q_grid, mode,
                                   xcoord, ycoord, qcoord, f, mu, beta):
    qmax = q_grid[-1]
    sx, sy = component_interp_scaled_coords(mode, xcoord, ycoord, f, mu, beta)
    q = beta_to_q(beta, qcoord)

    if method == "raw":
        pred = interp3_scaled(fbar_table, qmax, sx, sy, q)
    elif method == "delta_f":
        pred = f + interp3_scaled(fbar_table - f_node, qmax, sx, sy, q)
    else:
        raise ValueError("component tables support raw and delta_f, not %r" %
                         method)
    return np.clip(pred, 0.0, 1.0)


def summarize(method, axes, spec, bytes_used, corner_loads, pred, exact,
              target, f, mu, beta):
    err = np.abs(pred - exact)
    worst = int(np.argmax(err))
    return CandidateResult(
        method=method,
        axes=axes,
        spec=spec,
        bytes_used=bytes_used,
        corner_loads=corner_loads,
        max_err=float(np.max(err)),
        p99_err=float(np.percentile(err, 99.0)),
        rms_err=float(np.sqrt(np.mean(err * err))),
        passed=bool(np.max(err) <= target),
        worst_f=float(f[worst]),
        worst_mu=float(mu[worst]),
        worst_beta=float(beta[worst]),
        worst_exact=float(exact[worst]),
        worst_pred=float(pred[worst]),
    )


def fit_cheb_candidate(coord, beta_max, deg, train, chi_fn, n_bisect):
    nf, nmu, nq = train
    df, dmu, dq = deg
    xf, _ = cheb.chebgauss(nf)
    xmu, _ = cheb.chebgauss(nmu)
    xq, _ = cheb.chebgauss(nq)
    f_grid = 0.5 * (xf + 1.0)
    mu_grid = xmu
    qmax = float(beta_to_q(beta_max, coord))
    q_grid = 0.5 * (xq + 1.0) * qmax
    beta_grid = q_to_beta(q_grid, coord)

    values = np.empty((nf, nmu, nq))
    for iq, beta in enumerate(beta_grid):
        values[:, :, iq] = exact_slab(f_grid, mu_grid, beta, chi_fn, n_bisect)

    xx, yy, zz = np.meshgrid(xf, xmu, xq, indexing="ij")
    vand = cheb.chebvander3d(xx.ravel(), yy.ravel(), zz.ravel(), deg)
    coef_flat, _, _, _ = np.linalg.lstsq(vand, values.ravel(), rcond=None)
    return coef_flat.reshape((df + 1, dmu + 1, dq + 1)), qmax


def eval_cheb(coef, qmax, coord, f, mu, beta):
    xf = 2.0 * np.clip(f, 0.0, 1.0) - 1.0
    xmu = np.clip(mu, -1.0, 1.0)
    q = np.clip(beta_to_q(beta, coord), 0.0, qmax)
    xq = 2.0 * q / qmax - 1.0 if qmax > 0.0 else np.zeros_like(q)
    return np.clip(cheb.chebval3d(xf, xmu, xq, coef), 0.0, 1.0)


def format_bytes(nbytes):
    if nbytes < 1024:
        return "%d B" % nbytes
    if nbytes < 1024 * 1024:
        return "%.1f KiB" % (nbytes / 1024.0)
    return "%.2f MiB" % (nbytes / (1024.0 * 1024.0))


def print_results(results):
    print("")
    print("%-10s %-20s %-11s %10s %8s %12s %12s %12s %s" %
          ("method", "axes", "grid/deg", "bytes", "loads",
           "max_err", "p99_err", "rms_err", "pass"))
    print("-" * 114)
    for r in sorted(results, key=lambda x: (not x.passed, x.max_err)):
        print("%-10s %-20s %-11s %10s %8s %12.4e %12.4e %12.4e %s" %
              (r.method, r.axes, r.spec, format_bytes(r.bytes_used),
               r.corner_loads, r.max_err, r.p99_err, r.rms_err,
               "yes" if r.passed else "no"))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--target", type=float, default=1.0e-4,
                    help="required max absolute fbar error")
    ap.add_argument("--dims", default="24x24x32,32x32x48,40x40x56",
                    help="comma-separated NFxNMUxNQ trilinear table sizes")
    ap.add_argument("--axis-candidates",
                    default=("linear:linear:beta,linear:linear:eta,"
                             "linear:linear:wbeta,edge2:plus2:wbeta,"
                             "edge4:plus4:wbeta,edge4:plus4:eta,"
                             "both8:plus8:eta,"
                             "both4:edge4:eta,cart:edge4:zero4:eta,"
                             "cart2:edge4:both4:eta,"
                             "cartnorm:edge4:both4:eta,"
                             "thickcart:center4:zero4:eta,"
                             "thickcartnorm:center4:zero4:eta,"
                             "isopolar:center1:plus8:eta,"
                             "thickpolar:center4:plus4:eta"),
                    help=("comma-separated FAXIS:MUAXIS:VAXIS candidates; "
                          "or MODE:XAXIS:YAXIS:VAXIS for component tables. "
                          "f/y axes: linear, zeroN, edgeN, bothN; "
                          "mu/x axes: linear, plusN, edgeN, centerN; "
                          "v axes: beta, eta, wbeta, logw"))
    ap.add_argument("--methods", default="raw,delta_f,edge",
                    help=("comma-separated methods: raw, delta_f, edge, "
                          "plusblend"))
    ap.add_argument("--closure", default="Levermore",
                    choices=sorted(closure_gen.CLOSURES))
    ap.add_argument("--beta-max", type=float, default=0.999,
                    help="maximum beta covered by the table and validation set")
    ap.add_argument("--validation-points", type=int, default=8192,
                    help="number of random validation points before edge cases")
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--bisect", type=int, default=80,
                    help="bisection iterations for exact root solves")
    ap.add_argument("--face-fcoord", default="center2",
                    help="f coordinate for the supplemental mu=+1 face")
    ap.add_argument("--face-nf", type=int, default=0,
                    help="face f points for plusblend; 0 reuses table NF")
    ap.add_argument("--face-nq", type=int, default=0,
                    help="face velocity points for plusblend; 0 reuses table NQ")
    ap.add_argument("--face-width", type=float, default=0.02,
                    help="perpendicular flux width for plusblend correction")
    ap.add_argument("--face-power", type=float, default=4.0,
                    help="perpendicular flux falloff power for plusblend")
    ap.add_argument("--face-xwidth", type=float, default=0.06,
                    help="parallel distance width around f_iso for plusblend")
    ap.add_argument("--face-xpower", type=float, default=4.0,
                    help="parallel distance falloff power for plusblend")
    ap.add_argument("--include-cheb", action="store_true",
                    help="also fit/evaluate a tensor Chebyshev approximation")
    ap.add_argument("--cheb-coords", default="eta,wbeta")
    ap.add_argument("--cheb-deg", default="10,10,12",
                    help="Chebyshev degrees df,dmu,dq")
    ap.add_argument("--cheb-train", default="14,14,18",
                    help="Chebyshev training nodes nf,nmu,nq")
    args = ap.parse_args()

    beta_max = beta_max_clip(args.beta_max)
    specs = parse_dims(args.dims)
    axis_candidates = parse_axis_candidates(args.axis_candidates)
    methods = parse_csv(args.methods)
    chi_fn = closure_gen.CLOSURES[args.closure]

    f_val, mu_val, beta_val = make_validation_points(args.validation_points,
                                                     beta_max, args.seed)
    exact = exact_fbar_points(f_val, mu_val, beta_val, chi_fn, args.bisect)
    print("validation points: %d" % exact.size)
    print("target max |dfbar|: %.3e" % args.target)
    print("beta_max: %.12g" % beta_max)

    results = []
    for cand in axis_candidates:
        qmax = float(beta_to_q(beta_max, cand.qcoord))
        if not np.isfinite(qmax) or qmax <= 0.0:
            raise ValueError("invalid qmax for coord %s" % cand.qcoord)
        axes = "%s/%s/%s" % (cand.xcoord, cand.ycoord, cand.qcoord)
        if cand.mode != "polar":
            axes = "%s:%s" % (cand.mode, axes)
        for spec in specs:
            label = "%dx%dx%d" % (spec.nf, spec.nmu, spec.nq)
            sys.stderr.write("building %s table %s...\n" % (axes, label))
            if cand.mode == "polar":
                _, _, q_grid, _, fbar_table = build_trilinear_grid(
                    spec, cand.xcoord, cand.ycoord, cand.qcoord, beta_max,
                    chi_fn, args.bisect)
                f_node = None
                face_q_grid = None
                face_table = None
                if "plusblend" in methods:
                    face_nf = args.face_nf if args.face_nf > 0 else spec.nf
                    face_nq = args.face_nq if args.face_nq > 0 else spec.nq
                    face_q_grid, face_table = build_plus_face(
                        face_nf, face_nq, cand.qcoord, beta_max, chi_fn,
                        args.bisect, args.face_fcoord)
            elif cand.mode == "isopolar":
                _, _, q_grid, _, f_node, fbar_table = build_isopolar_grid(
                    spec, cand.xcoord, cand.ycoord, cand.qcoord, beta_max,
                    chi_fn, args.bisect)
                face_q_grid = None
                face_table = None
            elif cand.mode == "thickpolar":
                _, _, q_grid, _, f_node, fbar_table = build_thickpolar_grid(
                    spec, cand.xcoord, cand.ycoord, cand.qcoord, beta_max,
                    chi_fn, args.bisect)
                face_q_grid = None
                face_table = None
            else:
                _, _, q_grid, _, f_node, fbar_table = build_component_grid(
                    spec, cand.mode, cand.xcoord, cand.ycoord, cand.qcoord,
                    beta_max, chi_fn, args.bisect)
                face_q_grid = None
                face_table = None
            for method in methods:
                if method == "plusblend" and cand.mode != "polar":
                    continue
                if cand.mode != "polar" and method == "edge":
                    continue
                if cand.mode == "polar":
                    if method == "plusblend":
                        pred = plus_blend_prediction(
                            fbar_table, q_grid, face_table, face_q_grid,
                            cand.xcoord, cand.ycoord, cand.qcoord,
                            args.face_fcoord, f_val, mu_val, beta_val,
                            args.face_width, args.face_power,
                            args.face_xwidth, args.face_xpower)
                    else:
                        pred = candidate_prediction(method, fbar_table, q_grid,
                                                    cand.xcoord, cand.ycoord,
                                                    cand.qcoord, f_val, mu_val,
                                                    beta_val)
                elif cand.mode == "isopolar":
                    pred = isopolar_candidate_prediction(
                        method, fbar_table, f_node, q_grid, cand.xcoord,
                        cand.ycoord, cand.qcoord, f_val, mu_val, beta_val)
                elif cand.mode == "thickpolar":
                    pred = thickpolar_candidate_prediction(
                        method, fbar_table, f_node, q_grid, cand.xcoord,
                        cand.ycoord, cand.qcoord, f_val, mu_val, beta_val)
                else:
                    pred = component_candidate_prediction(
                        method, fbar_table, f_node, q_grid, cand.mode,
                        cand.xcoord, cand.ycoord, cand.qcoord, f_val, mu_val,
                        beta_val)
                extra = spec.nq if method == "edge" else 0
                if method == "plusblend":
                    extra += face_table.size
                bytes_used = (fbar_table.size + extra) * 8
                if method == "edge":
                    loads = "10"
                elif method == "plusblend":
                    loads = "16"
                else:
                    loads = "8"
                results.append(summarize(method, axes, label, bytes_used,
                                         loads, pred, exact, args.target,
                                         f_val, mu_val, beta_val))

    if args.include_cheb:
        deg = parse_triplet(args.cheb_deg)
        train = parse_triplet(args.cheb_train)
        for coord in parse_csv(args.cheb_coords):
            label = "%d,%d,%d" % deg
            sys.stderr.write("fitting Chebyshev %s deg=%s train=%s...\n" %
                             (coord, label, train))
            coef, qmax = fit_cheb_candidate(coord, beta_max, deg, train,
                                            chi_fn, args.bisect)
            pred = eval_cheb(coef, qmax, coord, f_val, mu_val, beta_val)
            bytes_used = coef.size * 8
            results.append(summarize("cheb", coord, label, bytes_used,
                                     str(coef.size), pred, exact, args.target,
                                     f_val, mu_val, beta_val))

    print_results(results)
    passed = [r for r in results if r.passed]
    if passed:
        best = min(passed, key=lambda r: (r.bytes_used, r.max_err))
        print("\nsmallest passing candidate: %s %s %s, %s, max_err=%.4e" %
              (best.method, best.axes, best.spec, format_bytes(best.bytes_used),
               best.max_err))
    else:
        best = min(results, key=lambda r: r.max_err)
        print("\nno candidate reached target; best max_err was %.4e for %s %s %s" %
              (best.max_err, best.method, best.axes, best.spec))
        print("worst point for that candidate: f=%.8g mu=%.8g beta=%.8g "
              "exact=%.8g pred=%.8g" %
              (best.worst_f, best.worst_mu, best.worst_beta,
               best.worst_exact, best.worst_pred))


if __name__ == "__main__":
    main()
