#!/usr/bin/env python3
"""Invariant/end-state oracle for the four Section 4.7 shock tubes."""
import argparse
import numpy as np
import matplotlib.pyplot as plt
from common import load_line

C = 2.99792458e10
TABLE = {
    1: ((1, 3e-5, .0015, 1e-8), (2.4, 1.61e-4, 6.25e-3, 2.51e-7)),
    2: ((1, 4e-3, .25, 2e-5), (3.11, .04512, .0804, 3.46e-3)),
    3: ((1, 60, 10, 2), (8, 2.34e3, 1.25, 1.14e3)),
    4: ((1, 6e-3, .69, .18), (3.65, 3.59e-2, .189, 1.30)),
}


def radiation_beta(E, F):
    f = np.clip(F / (C * np.maximum(E, 1.0e-300)), -1.0, 1.0)
    return np.where(np.abs(f) < 1.0e-14, 0.0,
                    (2.0 - np.sqrt(np.maximum(4.0 - 3.0*f*f, 0.0))) / f)


def conserved(d, gamma_ad):
    rho, pressure, v = d["rho"], d["pressure"], d["vel"][0]
    beta = v / C; lor = 1.0 / np.sqrt(np.maximum(1.0 - beta*beta, 1.0e-30))
    h = 1.0 + gamma_ad / (gamma_ad - 1.0) * pressure / (rho*C*C)
    E, F = d["E"].sum(0), d["F"].sum(0)
    return np.array([np.sum(rho*lor),
                     np.sum(rho*h*lor*lor*v + F/(C*C)),
                     np.sum(rho*h*lor*lor*C*C - pressure + E)])


def total_variation_ratio(y, left, right):
    return np.sum(np.abs(np.diff(y))) / max(abs(right-left), 1.0e-300)


p = argparse.ArgumentParser(); p.add_argument("dump"); p.add_argument("--initial")
p.add_argument("--case", type=int, choices=(1, 2, 3, 4), required=True)
p.add_argument("--output", default="rad_shock.png"); a = p.parse_args()
d = load_line(a.dump, 18); x = d["x"]; E = d["E"].sum(0); F = d["F"].sum(0)
left, right = TABLE[a.case]; lm = x < -16; rm = x > 16
rho_unit = np.mean(d["rho"][lm]) / left[0]; punit = rho_unit*C*C
beta = d["vel"][0] / C; ux = beta / np.sqrt(np.maximum(1-beta*beta, 1e-30))
rad_beta = radiation_beta(E, F)
rad_ux = rad_beta / np.sqrt(np.maximum(1-rad_beta*rad_beta, 1e-30))
vals = ((np.mean(d["rho"][lm])/rho_unit, np.mean(d["pressure"][lm])/punit,
         np.mean(ux[lm]), np.mean(E[lm])/(rho_unit*C*C)),
        (np.mean(d["rho"][rm])/rho_unit, np.mean(d["pressure"][rm])/punit,
         np.mean(ux[rm]), np.mean(E[rm])/(rho_unit*C*C)))
expected = []
for state in TABLE[a.case]:
    lor = np.sqrt(1 + state[2]*state[2])
    expected.append((state[0], state[1], state[2], state[3]*(4*lor*lor-1)/3))
enderr = max(abs(vals[s][q]/expected[s][q]-1) for s in range(2) for q in range(4))
frameerr = max(abs(np.mean(rad_ux[lm])/left[2]-1), abs(np.mean(rad_ux[rm])/right[2]-1))
finite = (np.all(np.isfinite(E)) and np.all(E >= 0) and np.all(d["rho"] > 0) and
          np.all(d["pressure"] > 0) and np.all(np.isfinite(d["temp"])) and
          np.all(d["temp"] > 0) and np.all(np.abs(beta) < 1) and
          np.all(np.abs(F) <= C*E*(1+1e-12)))
drift = 0.0
if a.initial:
    old = load_line(a.initial, 18); gamma_ad = 2.0 if a.case == 3 else 5.0/3.0
    now, before = conserved(d, gamma_ad), conserved(old, gamma_ad)
    drift = np.max(np.abs(now-before) / np.maximum(np.abs(before), 1e-300))
fields = (d["rho"]/rho_unit, d["pressure"]/punit, ux,
          E/(rho_unit*C*C), rad_ux)
tv = max(total_variation_ratio(y, np.mean(y[lm]), np.mean(y[rm])) for y in fields)
rho_left, rho_right = np.mean(d["rho"][lm]), np.mean(d["rho"][rm])
rho_mid = 0.5 * (rho_left + rho_right)
crossings = np.flatnonzero((d["rho"][:-1]-rho_mid) * (d["rho"][1:]-rho_mid) <= 0)
if crossings.size:
    q = crossings[np.argmin(np.abs(0.5*(x[crossings]+x[crossings+1])))]
    frac = (rho_mid-d["rho"][q]) / (d["rho"][q+1]-d["rho"][q])
    mid = x[q] + frac * (x[q+1]-x[q])
else:
    mid = x[np.argmin(abs(d["rho"]-rho_mid))]
midpoint = abs(mid)/np.mean(np.diff(x))
print(f"case={a.case} endpoint={enderr:.6e} conservation={drift:.6e} "
      f"frame={frameerr:.6e} topology={tv:.6e} midpoint_cells={midpoint:.3f}")
fig, ax = plt.subplots(5, 1, sharex=True, figsize=(6, 10))
labels = ("rho", "P", "u^x", "R", "u_rad^x")
for z, y, label in zip(ax, fields, labels): z.plot(x, y); z.set_ylabel(label)
ax[-1].set_xlabel("x"); fig.tight_layout(); fig.savefig(a.output, dpi=160)
if not (finite and enderr <= .03 and drift <= .01 and frameerr <= .03 and
        tv <= 1.25 and midpoint <= 2):
    raise SystemExit(1)
