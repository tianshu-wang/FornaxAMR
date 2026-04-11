#!/usr/bin/env python3

from pathlib import Path
import sys

import h5py
import matplotlib.pyplot as plt
import numpy as np


def load_dump_files(output_dir: Path):
    compatible = []
    dumps = sorted(output_dir.glob("dump_*.h5"))
    if not dumps:
        raise FileNotFoundError("no dump files found in output/")
    for dump in dumps:
        with h5py.File(dump, "r") as h5:
            if "coord" in h5.attrs and "root_nx" in h5.attrs:
                compatible.append(dump)
    if not compatible:
        raise FileNotFoundError("no compatible dump files with geometry metadata found in output/")
    return compatible


def pressure_from_state(density: np.ndarray, eint: np.ndarray) -> np.ndarray:
    gamma = 5.0 / 3.0
    return (gamma - 1.0) * density * eint


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


def collect_plane_blocks(coords: np.ndarray, levels: np.ndarray, values: np.ndarray, plane: str, coord: np.ndarray, root_nx: np.ndarray):
    normal_axis, axis_a, axis_b = plane_axes(plane)
    slices = []
    plane_value = 0.0
    tol = 1.0e-12
    domain_max = coord[2 * normal_axis + 1]
    for bid, xmin in enumerate(coords):
        extent = block_extent(xmin, int(levels[bid]), coord, root_nx)
        xmax = xmin + extent
        intersects = (
            xmin[normal_axis] - tol <= plane_value < xmax[normal_axis] - tol or
            (abs(plane_value - domain_max) <= tol and xmin[normal_axis] - tol <= plane_value <= xmax[normal_axis] + tol)
        )
        if not intersects:
            continue
        block = values[bid]
        n = block.shape[0]
        dx = extent / n
        cell = int(np.clip(np.floor((plane_value - xmin[normal_axis]) / dx[normal_axis]), 0, n - 1))
        if plane == "yz":
            plane_values = block[cell, :, :]
        elif plane == "xz":
            plane_values = block[:, cell, :]
        else:
            plane_values = block[:, :, cell]
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


def plot_plane(output_dir: Path, step: int, coords: np.ndarray, levels: np.ndarray, values: np.ndarray, variable: str, plane: str, coord: np.ndarray, root_nx: np.ndarray) -> None:
    plot_dir = output_dir / "plots" / variable
    plot_dir.mkdir(parents=True, exist_ok=True)
    plot_path = plot_dir / f"{plane}-{step:05d}.png"

    if plot_path.exists():
        return

    fig, ax = plt.subplots(figsize=(6, 6))
    pcm = None
    plane_blocks = collect_plane_blocks(coords, levels, values, plane, coord, root_nx)
    if not plane_blocks:
        raise RuntimeError(f"no blocks intersect the {plane} plane")
    vmin = min(np.min(plane_values) for _, _, _, _, plane_values in plane_blocks)
    vmax = max(np.max(plane_values) for _, _, _, _, plane_values in plane_blocks)
    if np.isclose(vmin, vmax):
        pad = max(1.0e-12, 1.0e-6 * max(1.0, abs(vmin)))
        vmin -= pad
        vmax += pad
    for xmin, xmax, xedges, yedges, plane_values in plane_blocks:
        pcm = ax.pcolormesh(xedges, yedges, plane_values.T, shading="flat", cmap="viridis", vmin=vmin, vmax=vmax)
        draw_block_grid(ax, xmin, xmax, plane)
    fig.colorbar(pcm, ax=ax, label=variable)
    ax.set_title(f"{variable} {plane} step {step:05d}")
    ax.set_aspect("equal")
    fig.savefig(plot_path, dpi=150)
    plt.close(fig)


def main() -> None:
    variable = sys.argv[1] if len(sys.argv) > 1 else "pressure"
    output_dir = Path("output")
    dump_files = load_dump_files(output_dir)

    for dump_path in dump_files:
        step = int(dump_path.stem.split("_")[-1])
        with h5py.File(dump_path, "r") as h5:
            coords = h5["coordinate"][...]
            levels = h5["level"][...]
            coord = h5.attrs["coord"]
            root_nx = h5.attrs["root_nx"]
            if variable == "pressure":
                values = pressure_from_state(h5["density"][...], h5["eint"][...])
            else:
                values = h5[variable][...]
        for plane in ("xy", "yz", "xz"):
            plot_plane(output_dir, step, coords, levels, values, variable, plane, coord, root_nx)


if __name__ == "__main__":
    main()
