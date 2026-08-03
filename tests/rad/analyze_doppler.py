#!/usr/bin/env python3
"""Section 4.5 spectral and stationarity oracle."""
import argparse
import numpy as np
import matplotlib.pyplot as plt
from common import load_line
from redshift_common import metrics


def errors(d):
    r, mean, expected, lum = metrics(d, False, True)
    rel = np.abs(mean / expected - 1.0)
    near = (r > 1.0e6) & (r < 7.0e6)
    peak = (r >= 7.0e6) & (r <= 9.0e6)
    far = r > 1.0e7
    return r, mean, expected, lum, np.max(rel[near]), np.max(rel[peak]), \
        np.std(lum[far]) / max(abs(np.mean(lum[far])), 1.0e-300)


p = argparse.ArgumentParser()
p.add_argument("dump")
p.add_argument("--groups", type=int, choices=(15, 25), required=True)
p.add_argument("--previous", help="penultimate dump for the 0.5% stationarity check")
p.add_argument("--coarse", help="15-group final dump; required to test convergence when supplied")
p.add_argument("--output", default="doppler.png")
a = p.parse_args()
d = load_line(a.dump, a.groups)
r, mean, expected, lum, nearerr, peakerr, lscatter = errors(d)
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
    _, _, _, _, coarse_near, coarse_peak, _ = errors(load_line(a.coarse, 15))
    improved = nearerr < coarse_near and peakerr < coarse_peak
print(f"groups={a.groups} near={nearerr:.6e} peak={peakerr:.6e} "
      f"luminosity_scatter={lscatter:.6e} stationarity={stationary:.6e} improved={improved}")
far = r > 1.0e7
fig, ax = plt.subplots(1, 2, figsize=(10, 4))
ax[0].plot(r / 1.0e5, lum / np.mean(lum[far])); ax[0].axhline(1, color="k")
ax[1].plot(r / 1.0e5, mean / expected); ax[1].axhline(1, color="k")
ax[0].set(xlabel="r [km]", ylabel="L/Lout")
ax[1].set(xlabel="r [km]", ylabel="mean-energy ratio")
fig.tight_layout(); fig.savefig(a.output, dpi=160)
if not (finite and stationary <= 0.005 and improved and lscatter <= 0.03 and
        nearerr <= (0.055 if a.groups == 15 else 0.032) and
        peakerr <= (0.02 if a.groups == 15 else 0.01)):
    raise SystemExit(1)
