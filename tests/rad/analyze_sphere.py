#!/usr/bin/env python3
"""Analyze either homogeneous radiating-sphere configuration (Section 4.4)."""
import argparse
import numpy as np
import matplotlib.pyplot as plt
from common import load_line

C = 2.99792458e10

def analytic(r, rstar, kappa, source):
    mu, weight = np.polynomial.legendre.leggauss(512)
    energy = np.empty_like(r); flux = np.empty_like(r)
    for n, radius in enumerate(r):
        g = np.sqrt(max(0.0, 1.0 - (rstar / max(radius, 1e-300)) ** 2)) if radius else 0.0
        if radius <= rstar:
            path = radius * mu + np.sqrt(np.maximum(0.0, rstar*rstar-radius*radius*(1-mu*mu)))
        else:
            path = np.where(mu >= g, 2*np.sqrt(np.maximum(0.0,
                rstar*rstar-radius*radius*(1-mu*mu))), 0.0)
        intensity = source * (1.0 - np.exp(-kappa * path))
        energy[n] = 0.5 * np.sum(weight * intensity)
        flux[n] = 0.5 * C * np.sum(weight * mu * intensity)
    return energy, flux

p = argparse.ArgumentParser()
p.add_argument("dump"); p.add_argument("--previous")
p.add_argument("--regime", choices=("thin", "thick"), required=True)
p.add_argument("--output", default="sphere.png")
a = p.parse_args()
d = load_line(a.dump, 1); x=d["x"]; E=d["E"][0]; F=d["F"][0]
if a.regime == "thin": R, B, kappa, tol = 3e6, .8, 5.04e-6, .05
else: R, B, kappa, tol = 5e6, 10., 2.8e-4, .03
Ea, Fa = analytic(x, 1e6, kappa, B)
l1 = np.mean(np.abs(E-Ea))/max(np.mean(np.abs(Ea)), 1e-300)
outer = x > 1.2e6
flux_ratio = np.divide(F, C*E, out=np.zeros_like(F), where=E>0)
asym = np.mean(np.abs(E[outer]-Ea[outer]))/max(np.mean(Ea[outer]),1e-300)
stationary = 0.0
if a.previous:
    old=load_line(a.previous,1); stationary=np.max(np.abs(E-old["E"][0])/np.maximum(E,1e-300))
finite=np.all(np.isfinite(E)) and np.all(E>=0) and np.all(np.abs(F)<=C*E*(1+1e-12))
print(f"regime={a.regime} l1={l1:.6e} outer={asym:.6e} stationary={stationary:.6e}")
fig,ax=plt.subplots(2,1,sharex=True,figsize=(6,7)); ax[0].plot(x/R,E,"r.",ms=2,label="PRJ"); ax[0].plot(x/R,Ea,"k-",label="analytic"); ax[1].plot(x/R,flux_ratio,".",ms=2); ax[1].plot(x/R,np.divide(Fa,C*Ea,out=np.zeros_like(Fa),where=Ea>0),"k-"); ax[0].set_ylabel("E"); ax[1].set_ylabel("F/(cE)"); ax[1].set_xlabel("r/R"); ax[0].legend(); fig.tight_layout(); fig.savefig(a.output,dpi=160)
if not (finite and l1 <= tol and asym <= .02 and (not a.previous or stationary <= .005)): raise SystemExit(1)
