"""Shared Section 4.5/4.6 spectral diagnostics."""
import numpy as np
C=2.99792458e10; G=6.67430e-8
def velocity(r):
    return np.where(r<=7e6,0,np.where(r<=8e6,-.2*C*(r-7e6)/1e6,-.2*C*(8e6/r)**2))
def lapse(r):
    rs=1e6; rho=9e14; mass=4*np.pi*rho*rs**3/3
    phi=np.where(r<=rs,-G*mass*(3*rs*rs-r*r)/(2*rs**3),-G*mass/r)
    return np.sqrt(1+2*phi/C**2)
def metrics(d, gravity=False, doppler=True):
    n=d["E"].shape[0]; edge=np.geomspace(1,50,n+1); ec=np.sqrt(edge[:-1]*edge[1:])
    E=d["E"].sum(0); F=d["F"].sum(0); number=(d["E"]/ec[:,None]).sum(0)
    mean=E/np.maximum(number,1e-300); r=d["x"]
    v=velocity(r) if doppler else np.zeros_like(r); alpha=lapse(r) if gravity else np.ones_like(r)
    w=1/np.sqrt(1-(v/C)**2); alpha_star=lapse(np.array([1e6]))[0] if gravity else 1
    expected=15.76*alpha_star/(alpha*w*(1+v/C)); lum=4*np.pi*r*r*F
    return r,mean,expected,lum
