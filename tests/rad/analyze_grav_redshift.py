#!/usr/bin/env python3
"""Section 4.6 gravitational and combined redshift oracle."""
import argparse
import numpy as np
import matplotlib.pyplot as plt
from common import load_line
from redshift_common import metrics


def errors(d, doppler):
    r, mean, expected, lum = metrics(d, True, doppler)
    rel = np.abs(mean / expected - 1.0)
    err = np.max(rel[r > 1.1e6])
    peak = np.max(rel[(r > 7.0e6) & (r < 9.0e6)]) if doppler else 0.0
    return r, mean, expected, lum, err, peak


p = argparse.ArgumentParser()
p.add_argument("dump")
p.add_argument("--groups", type=int, choices=(15, 25), required=True)
p.add_argument("--doppler", action="store_true")
p.add_argument("--previous", help="penultimate dump for the 0.5% stationarity check")
p.add_argument("--coarse", help="matching 15-group final dump for convergence")
p.add_argument("--output", default="grav_redshift.png")
a = p.parse_args()
d = load_line(a.dump, a.groups)
r, mean, expected, lum, err, peak = errors(d, a.doppler)
finite = np.all(np.isfinite(mean)) and np.all(d["E"] >= 0)
stationary = 0.0
if a.previous:
    old = load_line(a.previous, a.groups)
    e0, e1 = old["E"].sum(0), d["E"].sum(0)
    mask = e1 > 1.0e-12 * np.max(e1)
    stationary = np.max(np.abs(e1[mask] - e0[mask]) / np.maximum(e1[mask], 1.0e-300))
improved = True
if a.coarse:
    if a.groups != 25:
        p.error("--coarse is only meaningful for a 25-group result")
    _, _, _, _, coarse_err, coarse_peak = errors(load_line(a.coarse, 15), a.doppler)
    improved = err < coarse_err and (not a.doppler or peak < coarse_peak)
print(f"groups={a.groups} doppler={a.doppler} max_error={err:.6e} "
      f"peak={peak:.6e} stationarity={stationary:.6e} improved={improved}")
plt.plot(r / 1.0e5, mean / expected, label="PRJ/analytic")
plt.axhline(1, color="k"); plt.xlabel("r [km]"); plt.ylabel("mean-energy ratio")
plt.tight_layout(); plt.savefig(a.output, dpi=160)
if not (finite and stationary <= 0.005 and improved and
        err <= (0.033 if a.groups == 15 else 0.017) and peak <= 0.012):
    raise SystemExit(1)
