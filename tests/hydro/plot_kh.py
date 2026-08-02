#!/usr/bin/env python3
"""Plot the density field of a 2D Kelvin-Helmholtz dump.

Reads one HDF5 dump written by prj (output/dump_XXXXX.h5), slices the xy plane
at mid-z, and draws each active AMR block at its own resolution with
pcolormesh so mixed refinement levels stitch seamlessly.

Usage:
    python plot_kh.py [dump.h5]        # default: latest output/dump_*.h5
Simplified from analysis/visualize.py (drops MHD/radiation/km-scaling).
"""
import glob
import sys

import h5py
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# MetaData column layout (see prj_io_fill_metadata).
M_ID, M_LEVEL, M_ACTIVE = 0, 2, 3
M_XMIN, M_XMAX, M_DX = slice(4, 7), slice(7, 10), slice(10, 13)
DENSITY = 0  # Data variable index


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
        data = h5["Data"]           # [nblk][nvar][bs^3]
        meta = h5["MetaData"][:]     # [nblk][ncol]

        blocks = []
        vmin, vmax = np.inf, -np.inf
        for bid in range(meta.shape[0]):
            row = meta[bid]
            if row[M_ACTIVE] < 0.5 or row[M_ID] < 0:
                continue
            xmin = row[M_XMIN]
            dx = row[M_DX]
            vol = np.asarray(data[bid, DENSITY, :]).reshape(bs, bs, bs)
            plane = vol[:, :, bs // 2]          # (x, y) at mid-z
            blocks.append((xmin, dx, plane))
            vmin = min(vmin, plane.min())
            vmax = max(vmax, plane.max())

    if not blocks:
        sys.exit("no active blocks in dump")

    fig, ax = plt.subplots(figsize=(6, 6))
    for xmin, dx, plane in blocks:
        xedges = xmin[0] + np.arange(bs + 1) * dx[0]
        yedges = xmin[1] + np.arange(bs + 1) * dx[1]
        ax.pcolormesh(xedges, yedges, plane.T, shading="flat",
                      vmin=vmin, vmax=vmax, cmap="viridis")

    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("KH density  t = %.3f  (%d blocks)" % (time, len(blocks)))
    sm = plt.cm.ScalarMappable(cmap="viridis",
                               norm=plt.Normalize(vmin=vmin, vmax=vmax))
    fig.colorbar(sm, ax=ax, label="density")

    out = path.rsplit(".", 1)[0] + "_kh.png"
    fig.savefig(out, dpi=130, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
