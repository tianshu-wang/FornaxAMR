#!/usr/bin/env python3
import argparse
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import load_line

C=2.99792458e10; R=2e9; KS=2.5e-6; RHO=9e14
p=argparse.ArgumentParser(); p.add_argument("dump"); p.add_argument("--output",default="diffusive_source.png"); a=p.parse_args()
d=load_line(a.dump,3); t=200*R/C+d["time"]; that=C*t
profile=(KS/that)**1.5*np.exp(-3*KS*d["x"]**2/(4*that))
exact=RHO*C*C*profile; exact_f=d["x"]/(2*t)*exact
num=d["E"].sum(axis=0); flux=d["F"].sum(axis=0)
mask=exact > exact.max()*1e-6
rel=np.abs(num[mask]-exact[mask])/exact[mask]
l1=float(np.mean(rel)); linf=float(np.max(rel))
sign=np.all(flux[d["x"]>0] >= -1e-12*np.max(np.abs(flux)))
finite=np.all(np.isfinite(num)) and np.all(num>=0)
fig,ax=plt.subplots(2,1,figsize=(6,6),sharex=True)
ax[0].plot(d["x"]/R,num/(RHO*C*C),label="numerical"); ax[0].plot(d["x"]/R,profile,"k--",label="analytic")
ax[1].plot(d["x"]/R,flux/(RHO*C**3),label="numerical"); ax[1].plot(d["x"]/R,exact_f/(RHO*C**3),"k--",label="analytic")
ax[0].set_ylabel("E/(rho c^2)"); ax[1].set_ylabel("F/(rho c^3)"); ax[1].set_xlabel("r/R"); ax[0].legend(); fig.tight_layout(); fig.savefig(a.output,dpi=150)
ok=l1<=.01 and linf<=.03 and sign and finite
print(f"physical_time_R_over_c={t*C/R:.6f} mean_relative_error={l1:.6e} max_relative_error={linf:.6e} outward_flux={sign} finite={finite}")
raise SystemExit(0 if ok else 1)
