#!/usr/bin/env python3
"""Print compact per-variable summaries for a PRJ Z4c dump/restart file."""

import sys

import numpy as np

try:
    import h5py
except ModuleNotFoundError:
    h5py = None


Z4C_NAMES = [
    "chi",
    "gxx", "gxy", "gxz", "gyy", "gyz", "gzz",
    "Khat",
    "Axx", "Axy", "Axz", "Ayy", "Ayz", "Azz",
    "Gamx", "Gamy", "Gamz",
    "Theta",
    "alpha",
    "betax", "betay", "betaz",
]


def main(argv):
    if len(argv) != 2:
        print(f"usage: {argv[0]} output/dump_00000.h5", file=sys.stderr)
        return 2
    if h5py is None:
        print("h5py is not installed; build and use analysis/z4c_dump_summary.c instead", file=sys.stderr)
        return 2

    with h5py.File(argv[1], "r") as h5:
        if "Z4c" not in h5:
            print(f"{argv[1]}: no Z4c dataset", file=sys.stderr)
            return 1
        z4c = h5["Z4c"][...]
        step = h5.attrs.get("step", -1)
        time = h5.attrs.get("time", np.nan)

    print(f"# file={argv[1]} step={step} time={time:.17e} shape={z4c.shape}")
    for var, name in enumerate(Z4C_NAMES):
        data = z4c[:, var, :]
        finite = np.isfinite(data)
        finite_count = int(finite.sum())
        if finite_count == 0:
            print(f"{name:>7s} finite=0")
            continue
        data_finite = data[finite]
        print(
            f"{name:>7s} finite={finite_count:d} "
            f"min={data_finite.min():.9e} "
            f"max={data_finite.max():.9e} "
            f"maxabs={np.abs(data_finite).max():.9e}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
