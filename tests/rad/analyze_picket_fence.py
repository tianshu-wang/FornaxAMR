#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import load_line

p=argparse.ArgumentParser(); p.add_argument("dump"); p.add_argument("--reference",default=str(Path(__file__).with_name("picket_fence_reference.csv"))); p.add_argument("--output",default="picket_fence.png"); a=p.parse_args()
d=load_line(a.dump,20); ref=np.genfromtxt(a.reference,delimiter=",",names=True,comments="#")
x=11*d["x"]; norm=7.5657e-15*1e12
u1=d["E"][0::2].sum(axis=0)/norm; u2=d["E"][1::2].sum(axis=0)/norm
theta=d["temp"]/1.0e3
vals=[theta,u1,u2]; names=["T/T0","U1","U2"]; refs=[ref["theta"],ref["U1"],ref["U2"]]
errors=[]
for val,rv in zip(vals,refs):
    interp=np.interp(ref["x"],x,val); mask=rv>1e-4; errors.append(float(np.max(np.abs(interp[mask]-rv[mask])/rv[mask])))
tau=d["time"]*2.99792458e10*11
finite=all(np.all(np.isfinite(v)) for v in vals) and np.all(d["E"]>=0)
static=np.max(np.abs(d["rho"]-1))<1e-10 and np.max(np.abs(d["vel"]))<1e-10
fig,ax=plt.subplots(2,1,figsize=(6,6),sharex=True)
ax[0].plot(x,theta,label="numerical"); ax[0].plot(ref["x"],ref["theta"],"o",mfc="none",label="transport")
ax[1].plot(x,u1,label="U1"); ax[1].plot(x,u2,label="U2"); ax[1].plot(ref["x"],ref["U1"],"o",mfc="none"); ax[1].plot(ref["x"],ref["U2"],"o",mfc="none")
ax[0].set_ylabel("T/T0"); ax[1].set_ylabel("U"); ax[1].set_xlabel("optical coordinate"); ax[0].legend(); ax[1].legend(); fig.tight_layout(); fig.savefig(a.output,dpi=150)
ok=max(errors)<=.10 and finite and static and tau<=10.0+1e-12
print(f"tau={tau:.6f} max_errors="+",".join(f"{n}:{e:.6e}" for n,e in zip(names,errors))+f" finite={finite} static={static}")
raise SystemExit(0 if ok else 1)
