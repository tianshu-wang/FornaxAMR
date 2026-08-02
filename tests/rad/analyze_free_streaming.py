#!/usr/bin/env python3
import argparse
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import load_line

C = 2.99792458e10
AR = 7.5657e-15

p = argparse.ArgumentParser()
p.add_argument("dump")
p.add_argument("--output", default="free_streaming.png")
a = p.parse_args()
d = load_line(a.dump, 3)
edges = np.logspace(-2, 4, 4)
width = np.diff(edges)
hot, cold = AR*1e24, AR*1e16
front = C*d["time"]
exact = np.where(d["x"] < front, hot, cold)
total = d["E"].sum(axis=0)
dx = np.median(np.diff(d["x"]))
mask = np.abs(d["x"]-front) > 2*dx
l1 = np.mean(np.abs(total[mask]-exact[mask])/exact[mask])
mid = 0.5*(hot+cold)
front_num = d["x"][np.argmin(np.abs(total-mid))]
front_cells = abs(front_num-front)/dx
causal = np.all(np.abs(d["F"]) <= C*d["E"]*(1+1e-8))
finite = np.all(np.isfinite(d["E"])) and np.all(d["E"] >= 0)
fig, ax = plt.subplots(figsize=(6,4))
for g in range(3): ax.semilogy(d["x"], d["E"][g]/width[g], label=f"group {g+1}")
ax.semilogy(d["x"], exact, "k--", label="analytic total")
ax.set(xlabel="x (mean free paths)", ylabel="radiation energy density")
ax.legend(); fig.tight_layout(); fig.savefig(a.output, dpi=150)
ok = l1 <= .02 and front_cells <= 2 and causal and finite
print(f"mean_relative_error={l1:.6e} front_error_cells={front_cells:.3f} causal={causal} finite={finite}")
raise SystemExit(0 if ok else 1)
