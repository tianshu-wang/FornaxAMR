#!/usr/bin/env python3
"""Slice-plane visualizer that reads restart_*.h5 files instead of dump_*.h5.

Restart files store the full simulation state needed to resume a run, so the
on-disk layout differs from the dump format that ``visualize.py`` consumes:

  * ``Data``  -> double[nblocks][PRJ_NVAR_PRIM][BS**3]
                 The three spatial indices are flattened into one contiguous
                 BS**3 axis (active cells only, C-order i*BS*BS + j*BS + k).
                 Stored in DOUBLE precision (dumps are single precision).
  * ``MetaData`` -> double[nblocks][PRJ_IO_METADATA_SIZE]
                 col 0 = id, col 2 = level, col 3 = active, cols 4..6 = xmin.
  * ``Bf`` (MHD) -> face fields, not used here.

There are no ``variable_names`` attributes and no ``eos`` dataset, so:
  * variables are addressed by their fixed primitive index (set below);
  * ``pressure`` / ``temperature`` / ``gamma`` CANNOT be plotted from a restart
    file (they are EOS-derived and never written to restarts).

See the efficiency note at the bottom of this file for how this compares to the
dump-based visualize.py.
"""

from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm, Normalize, SymLogNorm


# User settings. Choose a primitive variable name, or derived variables "B"/"vr".
VARIABLE = "vr"
OUTPUT_DIR = Path("output")
PLANES = ("xy", "yz", "xz")

# Set these to None to show the full domain. Ranges are in km.
X_RANGE = [-1e4, 1e4]
Y_RANGE = [-1e4, 1e4]

# Choose from "normalize", "lognorm", or "symlognorm".
COLOR_SCALE = "normalize"
COLOR_VMIN = None
COLOR_VMAX = None
SYMLOG_LINTHRESH = 1.0e-6
SYMLOG_LINSCALE = 1.0
CM_PER_KM = 1.0e5
AXIS_NAMES = ("x", "y", "z")

# --- Build-configuration constants (NOT stored in the restart file) ----------
# The restart format does not record whether MHD/radiation were compiled in, so
# the variable->index mapping must be supplied here to match the build that
# produced the files. Defaults mirror src/prj_defs.h.
HAS_MHD = False
NRAD = 0          # PRJ_NRAD
NEGROUP = 12      # PRJ_NEGROUP
RAD_GROUP_STRIDE = 4  # 1 + PRJ_NDIM
# B is written to dumps scaled by sqrt(4*pi); restarts store it unscaled. Apply
# the same factor so colour scales match visualize.py.
MHD_DUMP_SCALE = np.sqrt(4.0 * np.pi)

# Metadata column layout (prj_io_fill_metadata in src/prj_io.c).
META_COL_ID = 0
META_COL_LEVEL = 2
META_COL_ACTIVE = 3
META_COL_XMIN = 4  # cols 4,5,6


def build_prim_index_map():
    """variable name -> primitive index in the restart Data dataset."""
    index_map = {
        "density": 0,
        "v1": 1,
        "v2": 2,
        "v3": 3,
        "eint": 4,
        "ye": 5,
    }
    nhydro = 6
    if HAS_MHD:
        index_map.update({"B1": 6, "B2": 7, "B3": 8})
        nhydro = 9
    # Radiation primitives: PRJ_PRIM_RAD_E(field, group) = nhydro +
    # ((field*NEGROUP + group) * RAD_GROUP_STRIDE), with components E,F1,F2,F3.
    component_names = ("E", "F1", "F2", "F3")
    for field in range(NRAD):
        for group in range(NEGROUP):
            base = nhydro + ((field * NEGROUP + group) * RAD_GROUP_STRIDE)
            for comp in range(RAD_GROUP_STRIDE):
                name = f"rad{field}_g{group}_{component_names[comp]}"
                index_map[name] = base + comp
    return index_map


# Variables that come from the EOS dataset in dumps and are absent in restarts.
EOS_ONLY_VARIABLES = ("pressure", "temperature", "gamma")


def plot_path_for(output_dir: Path, variable: str, plane: str, dump_id: str) -> Path:
    return output_dir / "plots_restart" / variable / f"{plane}-{dump_id}.png"


def load_restart_files(output_dir: Path):
    files = sorted(output_dir.glob("restart_*.h5"))
    if not files:
        raise FileNotFoundError("no restart files found in output/")
    return files


def resolve_variable(index_map, variable: str):
    if variable not in index_map:
        available = ", ".join(sorted(index_map))
        raise KeyError(f"variable {variable!r} not found; available variables: {available}")
    return index_map[variable]


def magnetic_field_strength(B1, B2, B3):
    return np.sqrt(B1 * B1 + B2 * B2 + B3 * B3)


def block_extent(xmin, level, coord, root_nx):
    root_dx = np.array(
        [
            (coord[1] - coord[0]) / root_nx[0],
            (coord[3] - coord[2]) / root_nx[1],
            (coord[5] - coord[4]) / root_nx[2],
        ]
    )
    return root_dx / (2 ** level)


def plane_axes(plane: str):
    if plane == "yz":
        return 0, 1, 2
    if plane == "xz":
        return 1, 0, 2
    return 2, 0, 1


def reference_variable_name(variable: str) -> str:
    if variable == "B":
        return "B1"
    if variable == "vr":
        return "v1"
    return variable


def read_block_volume(data, index_map, variable: str, bid: int, block_size: int) -> np.ndarray:
    """Read one (block, variable) row from Data and reshape to [i, j, k].

    The restart layout flattens the spatial indices, so the smallest read that
    yields any plane is the full block volume (BS**3 values).
    """
    var_idx = resolve_variable(index_map, variable)
    flat = data[bid, var_idx, :]
    vol = flat.reshape(block_size, block_size, block_size)
    if variable in ("B1", "B2", "B3"):
        vol = vol * MHD_DUMP_SCALE
    return vol


def slice_plane(vol: np.ndarray, plane: str, cell: int) -> np.ndarray:
    if plane == "yz":
        return vol[cell, :, :]
    if plane == "xz":
        return vol[:, cell, :]
    return vol[:, :, cell]


def orient_plane_coordinate(values, axis, normal_axis, axis_a):
    if axis == normal_axis:
        return values
    if axis == axis_a:
        return values[:, None]
    return values[None, :]


def radial_velocity_from_plane(xmin, level, plane, cell, coord, root_nx, v1, v2, v3):
    normal_axis, axis_a, axis_b = plane_axes(plane)
    extent = block_extent(xmin, level, coord, root_nx)
    n = v1.shape[0]
    dx = extent / n
    axes = [None, None, None]

    axes[normal_axis] = xmin[normal_axis] + (cell + 0.5) * dx[normal_axis]
    axes[axis_a] = xmin[axis_a] + (np.arange(v1.shape[0]) + 0.5) * dx[axis_a]
    axes[axis_b] = xmin[axis_b] + (np.arange(v1.shape[1]) + 0.5) * dx[axis_b]

    x = orient_plane_coordinate(axes[0], 0, normal_axis, axis_a)
    y = orient_plane_coordinate(axes[1], 1, normal_axis, axis_a)
    z = orient_plane_coordinate(axes[2], 2, normal_axis, axis_a)
    radius = np.sqrt(x * x + y * y + z * z)
    numerator = v1 * x + v2 * y + v3 * z
    vr = np.empty_like(v1, dtype=np.result_type(v1, v2, v3, float))

    np.divide(numerator, radius, out=vr, where=radius > 0.0)
    vr[radius <= 0.0] = 0.0
    return vr


def read_plane_values(data, index_map, variable, bid, plane, cell, xmin, level,
                      coord, root_nx, block_size):
    if variable == "B":
        return magnetic_field_strength(
            slice_plane(read_block_volume(data, index_map, "B1", bid, block_size), plane, cell),
            slice_plane(read_block_volume(data, index_map, "B2", bid, block_size), plane, cell),
            slice_plane(read_block_volume(data, index_map, "B3", bid, block_size), plane, cell),
        )
    if variable == "vr":
        return radial_velocity_from_plane(
            xmin, level, plane, cell, coord, root_nx,
            slice_plane(read_block_volume(data, index_map, "v1", bid, block_size), plane, cell),
            slice_plane(read_block_volume(data, index_map, "v2", bid, block_size), plane, cell),
            slice_plane(read_block_volume(data, index_map, "v3", bid, block_size), plane, cell),
        )
    return slice_plane(read_block_volume(data, index_map, variable, bid, block_size), plane, cell)


def collect_plane_blocks(data, index_map, variable, coords, levels, plane, coord,
                         root_nx, block_size):
    normal_axis, axis_a, axis_b = plane_axes(plane)
    slices = []
    plane_value = 0.0
    tol = 1.0e-12
    domain_max = coord[2 * normal_axis + 1]
    n = block_size
    for bid, xmin in enumerate(coords):
        level = int(levels[bid])
        extent = block_extent(xmin, level, coord, root_nx)
        xmax = xmin + extent
        intersects = (
            xmin[normal_axis] - tol <= plane_value < xmax[normal_axis] - tol or
            (abs(plane_value - domain_max) <= tol and
             xmin[normal_axis] - tol <= plane_value <= xmax[normal_axis] + tol)
        )
        if not intersects:
            continue
        if X_RANGE is not None and (xmax[axis_a] <= X_RANGE[0] or xmin[axis_a] >= X_RANGE[1]):
            continue
        if Y_RANGE is not None and (xmax[axis_b] <= Y_RANGE[0] or xmin[axis_b] >= Y_RANGE[1]):
            continue
        dx = extent / n
        cell = int(np.clip(np.floor((plane_value - xmin[normal_axis]) / dx[normal_axis]), 0, n - 1))
        plane_values = read_plane_values(
            data, index_map, variable, bid, plane, cell, xmin, level, coord, root_nx, block_size
        )
        xedges = xmin[axis_a] + np.arange(n + 1) * dx[axis_a]
        yedges = xmin[axis_b] + np.arange(n + 1) * dx[axis_b]
        slices.append((xmin, xmax, xedges, yedges, plane_values))
    return slices


def draw_block_grid(ax, xmin, xmax, plane: str) -> None:
    _, axis_a, axis_b = plane_axes(plane)
    ax.plot([xmin[axis_a], xmax[axis_a]], [xmin[axis_b], xmin[axis_b]], color="grey", alpha=0.35, linewidth=0.8)
    ax.plot([xmin[axis_a], xmax[axis_a]], [xmax[axis_b], xmax[axis_b]], color="grey", alpha=0.35, linewidth=0.8)
    ax.plot([xmin[axis_a], xmin[axis_a]], [xmin[axis_b], xmax[axis_b]], color="grey", alpha=0.35, linewidth=0.8)
    ax.plot([xmax[axis_a], xmax[axis_a]], [xmin[axis_b], xmax[axis_b]], color="grey", alpha=0.35, linewidth=0.8)


def color_limits(plane_blocks):
    vmin = min(np.min(plane_values) for _, _, _, _, plane_values in plane_blocks)
    vmax = max(np.max(plane_values) for _, _, _, _, plane_values in plane_blocks)
    if np.isclose(vmin, vmax):
        pad = max(1.0e-12, 1.0e-6 * max(1.0, abs(vmin)))
        vmin -= pad
        vmax += pad
    if COLOR_VMIN is not None:
        vmin = COLOR_VMIN
    if COLOR_VMAX is not None:
        vmax = COLOR_VMAX
    return vmin, vmax


def build_norm(vmin: float, vmax: float):
    scale = COLOR_SCALE.lower()
    if scale == "normalize":
        return Normalize(vmin=vmin, vmax=vmax)
    if scale == "lognorm":
        if vmin <= 0.0:
            raise ValueError("LogNorm requires positive color limits; set COLOR_VMIN > 0 or choose another scale")
        return LogNorm(vmin=vmin, vmax=vmax)
    if scale == "symlognorm":
        return SymLogNorm(linthresh=SYMLOG_LINTHRESH, linscale=SYMLOG_LINSCALE, vmin=vmin, vmax=vmax)
    raise ValueError(f"unknown COLOR_SCALE '{COLOR_SCALE}'")


def plot_plane(output_dir, dump_id, dump_time, plane_blocks, variable, plane) -> None:
    plot_path = plot_path_for(output_dir, variable, plane, dump_id)
    plot_path.parent.mkdir(parents=True, exist_ok=True)
    if plot_path.exists():
        return

    fig, ax = plt.subplots(figsize=(6, 6))
    pcm = None
    if not plane_blocks:
        raise RuntimeError(f"no blocks intersect the {plane} plane")
    vmin, vmax = color_limits(plane_blocks)
    norm = build_norm(vmin, vmax)
    for xmin, xmax, xedges, yedges, plane_values in plane_blocks:
        pcm = ax.pcolormesh(xedges, yedges, plane_values.T, shading="flat", cmap="viridis", norm=norm)
        draw_block_grid(ax, xmin, xmax, plane)
    fig.colorbar(pcm, ax=ax, label=variable)
    _, axis_a, axis_b = plane_axes(plane)
    ax.set_title(f"{variable} at time = {dump_time:g} s")
    ax.set_xlabel(f"{AXIS_NAMES[axis_a]} [km]")
    ax.set_ylabel(f"{AXIS_NAMES[axis_b]} [km]")
    ax.set_aspect("equal")
    if X_RANGE is not None:
        ax.set_xlim(X_RANGE)
    if Y_RANGE is not None:
        ax.set_ylim(Y_RANGE)
    fig.savefig(plot_path, dpi=150)
    plt.close(fig)


def read_active_metadata(h5):
    """Return (coords[km], levels) for active leaf blocks only."""
    meta = h5["MetaData"]
    ids = meta[:, META_COL_ID]
    active = meta[:, META_COL_ACTIVE]
    keep = np.flatnonzero((active >= 0.5) & (ids >= 0.0))
    levels = meta[:, META_COL_LEVEL][keep].astype(int)
    coords = meta[:, META_COL_XMIN:META_COL_XMIN + 3][keep] / CM_PER_KM
    return keep, coords, levels


def main() -> None:
    variable = VARIABLE
    output_dir = OUTPUT_DIR
    if variable in EOS_ONLY_VARIABLES:
        raise SystemExit(
            f"{variable!r} is EOS-derived and is not stored in restart files; "
            "use visualize.py on dump_*.h5 for that variable."
        )
    index_map = build_prim_index_map()
    restart_files = load_restart_files(output_dir)
    checked_missing_plots = False

    for restart_path in restart_files:
        dump_id = restart_path.stem.split("_")[-1]
        planes = [p for p in PLANES if not plot_path_for(output_dir, variable, p, dump_id).exists()]
        if not planes:
            continue
        checked_missing_plots = True
        with h5py.File(restart_path, "r") as h5:
            if "coord" not in h5.attrs or "root_nx" not in h5.attrs:
                continue
            dump_time = float(h5.attrs["time"])
            coord = np.asarray(h5.attrs["coord"], dtype=float) / CM_PER_KM
            root_nx = h5.attrs["root_nx"]
            block_size = int(h5.attrs["block_size"])
            keep, coords, levels = read_active_metadata(h5)
            data = h5["Data"]
            for plane in planes:
                # collect_plane_blocks indexes Data by the kept (active) block ids
                plane_blocks = collect_plane_blocks(
                    KeptData(data, keep), index_map, variable, coords, levels,
                    plane, coord, root_nx, block_size
                )
                plot_plane(output_dir, dump_id, dump_time, plane_blocks, variable, plane)
        print(dump_id)
    if not checked_missing_plots:
        print("nothing to do (all plots already exist)")


class KeptData:
    """Maps compacted active-block indices back to original Data rows."""

    def __init__(self, data, keep):
        self._data = data
        self._keep = keep

    def __getitem__(self, key):
        bid, var_idx, cell_slice = key
        return self._data[int(self._keep[bid]), var_idx, cell_slice]


if __name__ == "__main__":
    main()


# -----------------------------------------------------------------------------
# Efficiency note vs visualize.py (dump-based) -- measured on a matched state
# (step 500, 29359 active blocks, MHD on, BS=8; dump 722 MB, restart 2.9 GB)
# -----------------------------------------------------------------------------
# The restart reader moves a lot MORE DATA but is NOT slower in wall time:
#
#   * Bytes touched for a centered 3-plane set (804 intersecting blocks x 3 vel
#     components): dump ~0.6 MB vs restart ~9.9 MB == 16x more for restart.
#     That 16x = 8x (full BS**3 block volume vs one BS*BS plane, forced by the
#     flattened cell axis) x 2x (float64 vs float32).
#
#   * Yet warm-cache read time was comparable, restart even marginally faster
#     (~530 ms vs ~690 ms for the 3-plane set). When the file is in page cache,
#     bytes moved are cheap and HDF5 hyperslab-SELECTION overhead dominates:
#     the dump issues many small strided plane selections, while the restart
#     does one contiguous BS**3 slab read + a fast in-memory numpy reshape/slice.
#
#   * End-to-end the two scripts are indistinguishable (~3.3 s each) because
#     matplotlib import + rendering dominates both.
#
# Where the restart format actually costs more:
#   - Cold cache / data that does not fit in RAM: the 16x bytes-touched (and 4x
#     larger file) translate directly into disk-read time.
#   - Storage / I/O at scale: 2.9 GB vs 722 MB for the same state.
#   - Restarts also carry inactive/non-leaf blocks, a 639-wide MetaData table,
#     and (MHD) face fields; we skip inactive blocks on read, so per-plane cost
#     is unaffected, but the file and the MetaData read are heavier.
#
# Capability difference (not efficiency): pressure/temperature/gamma are absent
# from restarts, and MHD/radiation presence must be configured manually above
# because the restart file does not record the build's variable layout.
