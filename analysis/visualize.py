#!/usr/bin/env python3

from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm, Normalize, SymLogNorm


# User settings. Choose a dump variable name, or derived variables "B" or "vr".
VARIABLE = "temperature"
OUTPUT_DIR = Path("output")
PLANES = ("xy", "yz", "xz")

# Set these to None to show the full domain.
# Ranges are in km.
X_RANGE = [-6e2,6e2]
Y_RANGE = [-6e2,6e2]

# Choose from "normalize", "lognorm", or "symlognorm".
COLOR_SCALE = "lognorm"
COLOR_VMIN = None #1e7
COLOR_VMAX = None #1e15
SYMLOG_LINTHRESH = 1.0e-6
SYMLOG_LINSCALE = 1.0
CM_PER_KM = 1.0e5
AXIS_NAMES = ("x", "y", "z")


def plot_path_for(output_dir: Path, variable: str, plane: str, dump_id: str) -> Path:
    return output_dir / "plots" / variable / f"{plane}-{dump_id}.png"


def load_dump_files(output_dir: Path):
    dumps = sorted(output_dir.glob("dump_*.h5"))
    if not dumps:
        raise FileNotFoundError("no dump files found in output/")
    return dumps


def decode_variable_name(value) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8").rstrip("\x00")
    return str(value).rstrip("\x00")


def decode_variable_names(attr):
    values = np.asarray(attr)
    if values.shape == ():
        text = decode_variable_name(values.item())
        return [part.strip() for part in text.split(",") if part.strip()]
    return [decode_variable_name(value) for value in values.tolist()]


def build_variable_map(h5):
    variable_map = {}
    for dataset_name, obj in h5.items():
        if not isinstance(obj, h5py.Dataset) or "variable_names" not in obj.attrs:
            continue
        names = decode_variable_names(obj.attrs["variable_names"])
        if obj.ndim < 2 or obj.shape[1] != len(names):
            raise ValueError(f"dataset {dataset_name!r} has inconsistent variable_names")
        for var_idx, name in enumerate(names):
            if name in variable_map:
                raise ValueError(f"duplicate dump variable name {name!r}")
            variable_map[name] = (dataset_name, var_idx)
    if not variable_map:
        raise ValueError("dump file has no grouped datasets with variable_names attributes")
    return variable_map


def resolve_variable(variable_map, variable: str):
    if variable not in variable_map:
        available = ", ".join(sorted(variable_map))
        raise KeyError(f"variable {variable!r} not found in dump; available variables: {available}")
    return variable_map[variable]


def magnetic_field_strength(B1: np.ndarray, B2: np.ndarray, B3: np.ndarray) -> np.ndarray:
    return np.sqrt(B1 * B1 + B2 * B2 + B3 * B3)


def block_extent(xmin: np.ndarray, level: int, coord: np.ndarray, root_nx: np.ndarray) -> np.ndarray:
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


def read_dataset_plane(h5, variable_map, variable: str, bid: int, plane: str, cell: int) -> np.ndarray:
    dataset_name, var_idx = resolve_variable(variable_map, variable)
    dataset = h5[dataset_name]
    if dataset.ndim != 5:
        raise ValueError(f"variable {variable!r} is not cell-centered")
    if plane == "yz":
        return dataset[bid, var_idx, cell, :, :]
    if plane == "xz":
        return dataset[bid, var_idx, :, cell, :]
    return dataset[bid, var_idx, :, :, cell]


def orient_plane_coordinate(values, axis: int, normal_axis: int, axis_a: int):
    if axis == normal_axis:
        return values
    if axis == axis_a:
        return values[:, None]
    return values[None, :]


def radial_velocity_from_plane(
    xmin: np.ndarray,
    level: int,
    plane: str,
    cell: int,
    coord: np.ndarray,
    root_nx: np.ndarray,
    v1: np.ndarray,
    v2: np.ndarray,
    v3: np.ndarray,
) -> np.ndarray:
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


def read_plane_values(
    h5,
    variable_map,
    variable: str,
    bid: int,
    plane: str,
    cell: int,
    xmin: np.ndarray,
    level: int,
    coord: np.ndarray,
    root_nx: np.ndarray,
) -> np.ndarray:
    if variable == "pressure":
        return read_dataset_plane(h5, variable_map, "pressure", bid, plane, cell)
    if variable == "B":
        return magnetic_field_strength(
            read_dataset_plane(h5, variable_map, "B1", bid, plane, cell),
            read_dataset_plane(h5, variable_map, "B2", bid, plane, cell),
            read_dataset_plane(h5, variable_map, "B3", bid, plane, cell),
        )
    if variable == "vr":
        return radial_velocity_from_plane(
            xmin,
            level,
            plane,
            cell,
            coord,
            root_nx,
            read_dataset_plane(h5, variable_map, "v1", bid, plane, cell),
            read_dataset_plane(h5, variable_map, "v2", bid, plane, cell),
            read_dataset_plane(h5, variable_map, "v3", bid, plane, cell),
        )
    return read_dataset_plane(h5, variable_map, variable, bid, plane, cell)


def collect_plane_blocks_from_h5(
    h5,
    variable_map,
    variable: str,
    coords: np.ndarray,
    levels: np.ndarray,
    plane: str,
    coord: np.ndarray,
    root_nx: np.ndarray,
):
    normal_axis, axis_a, axis_b = plane_axes(plane)
    slices = []
    plane_value = 0.0
    tol = 1.0e-12
    domain_max = coord[2 * normal_axis + 1]
    ref_dataset_name, _ = resolve_variable(variable_map, reference_variable_name(variable))
    ref_dataset = h5[ref_dataset_name]
    if ref_dataset.ndim != 5:
        raise ValueError(f"variable {variable!r} is not cell-centered")
    n = ref_dataset.shape[2]
    for bid, xmin in enumerate(coords):
        level = int(levels[bid])
        extent = block_extent(xmin, level, coord, root_nx)
        xmax = xmin + extent
        intersects = (
            xmin[normal_axis] - tol <= plane_value < xmax[normal_axis] - tol or
            (abs(plane_value - domain_max) <= tol and xmin[normal_axis] - tol <= plane_value <= xmax[normal_axis] + tol)
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
            h5, variable_map, variable, bid, plane, cell, xmin, level, coord, root_nx
        )
        xedges = xmin[axis_a] + np.arange(n + 1) * dx[axis_a]
        yedges = xmin[axis_b] + np.arange(n + 1) * dx[axis_b]
        slices.append((xmin, xmax, xedges, yedges, plane_values))
    return slices


def draw_block_grid(ax, xmin: np.ndarray, xmax: np.ndarray, plane: str) -> None:
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
        return SymLogNorm(
            linthresh=SYMLOG_LINTHRESH,
            linscale=SYMLOG_LINSCALE,
            vmin=vmin,
            vmax=vmax,
        )
    raise ValueError(f"unknown COLOR_SCALE '{COLOR_SCALE}'")


def plot_plane(output_dir: Path, dump_id: str, dump_time: float, plane_blocks, variable: str, plane: str) -> None:
    plot_path = plot_path_for(output_dir, variable, plane, dump_id)
    plot_dir = plot_path.parent
    plot_dir.mkdir(parents=True, exist_ok=True)

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


def main() -> None:
    variable = VARIABLE
    output_dir = OUTPUT_DIR
    dump_files = load_dump_files(output_dir)
    checked_missing_plots = False
    found_compatible_dump = False

    for dump_path in dump_files:
        dump_id = dump_path.stem.split("_")[-1]
        planes = [plane for plane in PLANES if not plot_path_for(output_dir, variable, plane, dump_id).exists()]
        if not planes:
            continue
        checked_missing_plots = True
        with h5py.File(dump_path, "r") as h5:
            if "coord" not in h5.attrs or "root_nx" not in h5.attrs:
                continue
            found_compatible_dump = True
            variable_map = build_variable_map(h5)
            dump_time = float(h5.attrs["time"])
            coords = h5["coordinate"][...] / CM_PER_KM
            levels = h5["level"][...]
            coord = np.asarray(h5.attrs["coord"], dtype=float) / CM_PER_KM
            root_nx = h5.attrs["root_nx"]
            for plane in planes:
                plane_blocks = collect_plane_blocks_from_h5(
                    h5, variable_map, variable, coords, levels, plane, coord, root_nx
                )
                plot_plane(output_dir, dump_id, dump_time, plane_blocks, variable, plane)
        print(dump_id)
    if checked_missing_plots and not found_compatible_dump:
        raise FileNotFoundError("no compatible dump files with geometry metadata found in output/")


if __name__ == "__main__":
    main()
