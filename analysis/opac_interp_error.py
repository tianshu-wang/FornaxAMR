#!/usr/bin/env python3
"""Adjacent-node ratio + linear-vs-log interpolation error for the opacity table.

Reads the raw binary opacity table and, for each quantity, computes the
distribution of adjacent-node log-differences Delta = ln(y_{n+1}/y_n) along the
rho, T and Ye axes, then the midpoint relative error incurred by replacing the
current log (geometric) interpolation with plain linear-in-value interpolation:

    err_mid = (1+r)/(2*sqrt(r)) - 1 = cosh(Delta/2) - 1 ,  r = exp(Delta)

Stored values for abs/sca/emis are natural-log (lookup does exp(val+rl));
sdelta is stored linear, so its "log" interpretation does not apply and it is
reported separately on the linear values for reference.
"""
import numpy as np

PATH = "/Users/tianshuw/vibe-coding/test-codex/opacity.SFHo.juo.horo.brem1.extendedT.bin"
NG, NRO, NT, NYE = 40, 50, 50, 30
FIELD = NG * NRO * NT * NYE          # 3,000,000 floats per record
IRECL = 4 * FIELD + 10               # bytes per Fortran-style record
QUANT = ["absopac", "scaopac", "emis", "sdelta"]   # record_base/3 = 0,1,2,3
AXES = {"rho": 3, "T": 2, "Ye": 1}   # C-order shape (ng, Ye, T, rho)

def read_record(f, rec):
    f.seek(rec * IRECL)
    a = np.fromfile(f, dtype="<f4", count=FIELD).astype(np.float64)
    return a.reshape(NG, NYE, NT, NRO)   # (ng, Ye, T, rho)

def pct(x, ps):
    return np.percentile(x, ps)

def main():
    ps = [50, 90, 99, 99.9, 100]
    err_thresh = [0.01, 0.05, 0.20, 0.50]
    with open(PATH, "rb") as f:
        for q, name in enumerate(QUANT):
            # stack the 3 species for this quantity
            data = np.stack([read_record(f, q * 3 + s) for s in range(3)], axis=0)
            islog = name != "sdelta"
            data = np.log(np.abs(data)) if not islog else data  # work in log(value)
            floor = data.min()
            is_floor = data <= floor + 0.5          # sentinel "zero" clamp
            print("=" * 78)
            print(f"{name}   (stored {'log' if islog else 'linear'})   "
                  f"log-value range [{floor:.3f}, {data.max():.3f}]  "
                  f"floored cells={is_floor.mean()*100:.1f}%")
            for axname, ax in AXES.items():
                lo_f = np.take(is_floor, range(0, data.shape[ax]-1), axis=ax)
                hi_f = np.take(is_floor, range(1, data.shape[ax]),   axis=ax)
                d_all = np.diff(data, axis=ax)          # Delta = ln(r)
                active = ~(lo_f | hi_f)                 # both endpoints above floor
                frac_floor_edge = 1.0 - active.mean()
                d = d_all[active]
                absd = np.abs(d)
                err = np.cosh(d/2.0) - 1.0              # midpoint lin-vs-log error
                ar = pct(absd, ps)
                er = pct(err, ps)
                frac = [np.mean(err > t) for t in err_thresh]
                print(f"  {axname:3s}  active N={d.size:>8d}  "
                      f"(floor-touching edges {frac_floor_edge*100:.1f}%)")
                print(f"        |Δln|    p50={ar[0]:.3f} p90={ar[1]:.3f} p99={ar[2]:.3f} "
                      f"p99.9={ar[3]:.3f} max={ar[4]:.3f}")
                print(f"        ratio r  p50={np.exp(ar[0]):.2f} p90={np.exp(ar[1]):.2f} "
                      f"p99={np.exp(ar[2]):.2f} max={np.exp(ar[4]):.1f}")
                print(f"        err_mid  p50={er[0]*100:.3f}% p90={er[1]*100:.3f}% "
                      f"p99={er[2]*100:.2f}% p99.9={er[3]*100:.1f}% max={er[4]*100:.1f}%")
                print(f"        frac active edges err> "
                      + "  ".join(f"{int(t*100)}%:{frac[i]*100:.2f}%"
                                  for i, t in enumerate(err_thresh)))
    print("=" * 78)

if __name__ == "__main__":
    main()
