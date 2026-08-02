#!/usr/bin/env python3
"""Plot 1D profiles of a Sod shock tube dump.

Reads one HDF5 dump written by prj (output/dump_XXXXX.h5), extracts the line
along x at the mid (y, z) cell of every active block, and plots density,
velocity, and pressure vs x.

Usage:
    python plot_shocktube.py [dump.h5]   # default: latest output/dump_*.h5
Simplified from analysis/visualize.py.
"""
import glob
import sys

import h5py
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# MetaData column layout (see prj_io_fill_metadata).
M_ID, M_ACTIVE = 0, 3
M_XMIN, M_DX = slice(4, 7), slice(10, 13)
# Data / Eos variable indices.
DENSITY, V1 = 0, 1
PRESSURE = 0


def latest_dump():
    files = sorted(glob.glob("output/dump_*.h5"))
    if not files:
        sys.exit("no output/dump_*.h5 found; run the simulation first")
    return files[-1]


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else latest_dump()
    with h5py.File(path, "r") as h5:
        time = float(h5.attrs["time"])
        bs = int(h5.attrs["block_size"])
        data = h5["Data"]
        eos = h5["Eos"]
        meta = h5["MetaData"][:]

        jc = kc = bs // 2
        xs, rho, vel, pres = [], [], [], []
        for bid in range(meta.shape[0]):
            row = meta[bid]
            if row[M_ACTIVE] < 0.5 or row[M_ID] < 0:
                continue
            xmin = row[M_XMIN]
            dx = row[M_DX]
            d = np.asarray(data[bid, DENSITY, :]).reshape(bs, bs, bs)
            v = np.asarray(data[bid, V1, :]).reshape(bs, bs, bs)
            p = np.asarray(eos[bid, PRESSURE, :]).reshape(bs, bs, bs)
            xc = xmin[0] + (np.arange(bs) + 0.5) * dx[0]
            xs.append(xc)
            rho.append(d[:, jc, kc])
            vel.append(v[:, jc, kc])
            pres.append(p[:, jc, kc])

    if not xs:
        sys.exit("no active blocks in dump")

    xs = np.concatenate(xs)
    order = np.argsort(xs)
    xs = xs[order]
    rho = np.concatenate(rho)[order]
    vel = np.concatenate(vel)[order]
    pres = np.concatenate(pres)[order]

    fig, axes = plt.subplots(3, 1, figsize=(6, 7), sharex=True)
    for ax, y, name in zip(axes, (rho, vel, pres),
                           ("density", "velocity", "pressure")):
        ax.plot(xs, y, "-o", ms=3, lw=1)
        ax.set_ylabel(name)
        ax.grid(True, alpha=0.3)
    axes[0].set_title("Sod shock tube  t = %.3f" % time)
    axes[-1].set_xlabel("x")

    out = path.rsplit(".", 1)[0] + "_shocktube.png"
    fig.savefig(out, dpi=130, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
