#!/usr/bin/env python3
"""Diagnose pressure-scale-height AMR changes in PRJ dump files.

This script does not require any PRJ internals at runtime, but it mirrors the
dump layout and the pressure_scale_height estimator:

    indicator = dx * rho * |g| / pressure

It is intended to be run in an environment with numpy and h5py.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import h5py
import numpy as np


BLOCK_SIZE = 8
NBINS_GRAVITY = 1024
EPS = 1.0e-300


BlockKey = Tuple[int, int, int, int]


@dataclass
class DumpMeta:
    path: Path
    dump_index: int
    step: int
    time: float
    coord: np.ndarray
    root_nx: np.ndarray
    levels: np.ndarray
    origins: np.ndarray
    keys: List[BlockKey]


@dataclass
class BlockStats:
    max_i: np.ndarray
    max_boundary_i: np.ndarray
    mean_i: np.ndarray
    max_rho: np.ndarray
    max_pressure: np.ndarray
    max_g: np.ndarray
    max_radius: np.ndarray
    argmax_flat: np.ndarray


@dataclass
class PostParentEquivRecord:
    value: float
    rho: float
    pressure: float
    g: float
    radius: float
    child_level: int
    child_key: BlockKey


def dump_path(dump_dir: Path, index: int) -> Path:
    return dump_dir / f"dump_{index:05d}.h5"


def parse_param_file(path: Optional[Path]) -> Dict[str, str]:
    if path is None:
        return {}
    params: Dict[str, str] = {}
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.split("#", 1)[0].split("//", 1)[0].strip()
            if not line or "=" not in line:
                continue
            key, value = line.split("=", 1)
            params[key.strip()] = value.strip()
    return params


def param_float(params: Dict[str, str], key: str) -> Optional[float]:
    if key not in params:
        return None
    try:
        return float(params[key])
    except ValueError:
        return None


def param_int(params: Dict[str, str], key: str) -> Optional[int]:
    value = param_float(params, key)
    if value is None:
        return None
    return int(value)


def read_scalar_attr(handle: h5py.File, name: str, default=None):
    if name not in handle.attrs:
        return default
    value = handle.attrs[name]
    if isinstance(value, np.ndarray):
        return value.item() if value.shape == () else value
    return value


def load_meta(path: Path) -> DumpMeta:
    with h5py.File(path, "r") as handle:
        levels = np.asarray(handle["level"][:], dtype=np.int64)
        origins = np.asarray(handle["coordinate"][:], dtype=np.float64)
        coord = np.asarray(handle.attrs["coord"], dtype=np.float64)
        root_nx = np.asarray(handle.attrs["root_nx"], dtype=np.int64)
        dump_index = int(read_scalar_attr(handle, "dump_index", -1))
        step = int(read_scalar_attr(handle, "step", -1))
        time = float(read_scalar_attr(handle, "time", math.nan))
    keys = make_block_keys(levels, origins, coord, root_nx)
    return DumpMeta(path, dump_index, step, time, coord, root_nx, levels, origins, keys)


def axis_minima(coord: np.ndarray) -> np.ndarray:
    return np.array([coord[0], coord[2], coord[4]], dtype=np.float64)


def axis_maxima(coord: np.ndarray) -> np.ndarray:
    return np.array([coord[1], coord[3], coord[5]], dtype=np.float64)


def axis_spans(coord: np.ndarray) -> np.ndarray:
    return axis_maxima(coord) - axis_minima(coord)


def block_widths(levels: np.ndarray, coord: np.ndarray, root_nx: np.ndarray) -> np.ndarray:
    root_width = axis_spans(coord) / root_nx.astype(np.float64)
    scale = np.ldexp(np.ones_like(levels, dtype=np.float64), -levels.astype(np.int64))
    return scale[:, None] * root_width[None, :]


def cell_widths(levels: np.ndarray, coord: np.ndarray, root_nx: np.ndarray) -> np.ndarray:
    return block_widths(levels, coord, root_nx) / float(BLOCK_SIZE)


def make_block_keys(
    levels: np.ndarray,
    origins: np.ndarray,
    coord: np.ndarray,
    root_nx: np.ndarray,
) -> List[BlockKey]:
    xmin = axis_minima(coord)
    widths = block_widths(levels, coord, root_nx)
    grid = np.rint((origins - xmin[None, :]) / widths).astype(np.int64)
    return [
        (int(levels[i]), int(grid[i, 0]), int(grid[i, 1]), int(grid[i, 2]))
        for i in range(levels.shape[0])
    ]


def level_hist(levels: np.ndarray) -> Dict[int, int]:
    return {int(k): int(v) for k, v in sorted(Counter(levels.tolist()).items())}


def effective_gravity_level(
    meta: DumpMeta,
    requested_max_level: Optional[int],
    min_dx: Optional[float],
) -> int:
    root_cell = float(np.max(axis_spans(meta.coord) / meta.root_nx.astype(np.float64) / BLOCK_SIZE))
    if requested_max_level is None:
        level = int(np.max(meta.levels))
    else:
        level = int(requested_max_level)

    if min_dx is not None and min_dx > 0.0:
        min_dx_level = 0
        cell_size = root_cell
        while cell_size > min_dx:
            cell_size *= 0.5
            min_dx_level += 1
        if requested_max_level is None:
            level = min_dx_level
        else:
            level = min(level, min_dx_level)
    return max(level, 0)


def gravity_rf(meta: DumpMeta, gravity_level: int) -> np.ndarray:
    rmax = 0.5 * float(np.min(axis_spans(meta.coord)))
    root_cell = float(np.max(axis_spans(meta.coord) / meta.root_nx.astype(np.float64) / BLOCK_SIZE))
    min_cell = math.ldexp(root_cell, -gravity_level)
    dr_min = 0.5 * min_cell
    if dr_min <= 0.0 or rmax <= dr_min:
        dr_min = rmax / float(NBINS_GRAVITY)

    rf = np.empty(NBINS_GRAVITY + 1, dtype=np.float64)
    rf[0] = 0.0
    log_span = math.log(rmax / dr_min)
    for i in range(1, NBINS_GRAVITY + 1):
        rf[i] = dr_min * math.exp(log_span * float(i - 1) / float(NBINS_GRAVITY - 1))
    return rf


def interp_accel_like_prj(accel: np.ndarray, rf: np.ndarray, radius: np.ndarray) -> np.ndarray:
    centers = 0.5 * (rf[:-1] + rf[1:])
    out = np.interp(radius, centers, accel, left=accel[0], right=accel[-1])
    first = radius <= centers[0]
    if np.any(first):
        if centers[0] > 0.0:
            out[first] = accel[0] * (radius[first] / centers[0])
        else:
            out[first] = accel[0]
    return out


def cell_radii_for_chunk(meta: DumpMeta, start: int, stop: int, x_com: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    levels = meta.levels[start:stop]
    origins = meta.origins[start:stop]
    dcell = cell_widths(levels, meta.coord, meta.root_nx)
    q = np.arange(BLOCK_SIZE, dtype=np.float64) + 0.5
    x = origins[:, 0, None, None, None] + dcell[:, 0, None, None, None] * q[None, :, None, None]
    y = origins[:, 1, None, None, None] + dcell[:, 1, None, None, None] * q[None, None, :, None]
    z = origins[:, 2, None, None, None] + dcell[:, 2, None, None, None] * q[None, None, None, :]
    radius = np.sqrt((x - x_com[0]) ** 2 + (y - x_com[1]) ** 2 + (z - x_com[2]) ** 2)
    dx = np.max(dcell, axis=1)
    return radius, dx


def indicator_chunk(
    meta: DumpMeta,
    handle: h5py.File,
    rf: np.ndarray,
    start: int,
    stop: int,
    x_com: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    rho = np.asarray(handle["hydro"][start:stop, 0, :, :, :], dtype=np.float64)
    pressure = np.asarray(handle["eos"][start:stop, 0, :, :, :], dtype=np.float64)
    accel = np.asarray(handle["accel"][:], dtype=np.float64)
    radius, dx = cell_radii_for_chunk(meta, start, stop, x_com)
    g = np.abs(interp_accel_like_prj(accel, rf, radius))
    indicator = dx[:, None, None, None] * rho * g / np.maximum(pressure, EPS)
    return indicator, rho, pressure, g, radius


def compute_block_stats(
    meta: DumpMeta,
    rf: np.ndarray,
    x_com: np.ndarray,
    chunk_size: int,
) -> BlockStats:
    nblocks = meta.levels.shape[0]
    max_i = np.empty(nblocks, dtype=np.float64)
    max_boundary_i = np.empty(nblocks, dtype=np.float64)
    mean_i = np.empty(nblocks, dtype=np.float64)
    max_rho = np.empty(nblocks, dtype=np.float64)
    max_pressure = np.empty(nblocks, dtype=np.float64)
    max_g = np.empty(nblocks, dtype=np.float64)
    max_radius = np.empty(nblocks, dtype=np.float64)
    argmax_flat = np.empty(nblocks, dtype=np.int64)
    boundary_mask = np.zeros((BLOCK_SIZE, BLOCK_SIZE, BLOCK_SIZE), dtype=bool)
    boundary_mask[0, :, :] = True
    boundary_mask[-1, :, :] = True
    boundary_mask[:, 0, :] = True
    boundary_mask[:, -1, :] = True
    boundary_mask[:, :, 0] = True
    boundary_mask[:, :, -1] = True
    boundary_flat = boundary_mask.ravel()

    with h5py.File(meta.path, "r") as handle:
        for start in range(0, nblocks, chunk_size):
            stop = min(start + chunk_size, nblocks)
            indicator, rho, pressure, g, radius = indicator_chunk(meta, handle, rf, start, stop, x_com)
            flat_i = indicator.reshape(stop - start, -1)
            flat_rho = rho.reshape(stop - start, -1)
            flat_p = pressure.reshape(stop - start, -1)
            flat_g = g.reshape(stop - start, -1)
            flat_r = radius.reshape(stop - start, -1)
            arg = np.argmax(flat_i, axis=1)
            rows = np.arange(stop - start)
            max_i[start:stop] = flat_i[rows, arg]
            max_boundary_i[start:stop] = np.max(flat_i[:, boundary_flat], axis=1)
            mean_i[start:stop] = np.mean(flat_i, axis=1)
            max_rho[start:stop] = flat_rho[rows, arg]
            max_pressure[start:stop] = flat_p[rows, arg]
            max_g[start:stop] = flat_g[rows, arg]
            max_radius[start:stop] = flat_r[rows, arg]
            argmax_flat[start:stop] = arg

    return BlockStats(max_i, max_boundary_i, mean_i, max_rho, max_pressure, max_g, max_radius, argmax_flat)


def quantiles(values: np.ndarray) -> Dict[str, float]:
    values = values[np.isfinite(values)]
    if values.size == 0:
        return {name: math.nan for name in ("min", "p01", "p05", "p50", "p95", "p99", "max")}
    qs = np.quantile(values, [0.0, 0.01, 0.05, 0.5, 0.95, 0.99, 1.0])
    return {
        "min": float(qs[0]),
        "p01": float(qs[1]),
        "p05": float(qs[2]),
        "p50": float(qs[3]),
        "p95": float(qs[4]),
        "p99": float(qs[5]),
        "max": float(qs[6]),
    }


def same_geometry(a: DumpMeta, b: DumpMeta) -> bool:
    return (
        a.levels.shape == b.levels.shape
        and np.array_equal(a.levels, b.levels)
        and np.allclose(a.origins, b.origins, rtol=0.0, atol=1.0e-8)
        and a.keys == b.keys
    )


def compare_same_mesh(
    before: DumpMeta,
    after: DumpMeta,
    rf_before: np.ndarray,
    rf_after: np.ndarray,
    x_com: np.ndarray,
    chunk_size: int,
    block_mask: Optional[np.ndarray] = None,
) -> Dict[str, object]:
    if not same_geometry(before, after):
        raise RuntimeError(f"{before.path.name} and {after.path.name} do not have identical block geometry")

    ratios: Dict[str, List[np.ndarray]] = {
        "rho": [],
        "pressure": [],
        "g": [],
        "indicator": [],
    }
    log_terms: Dict[str, List[np.ndarray]] = {
        "log_rho": [],
        "log_g": [],
        "minus_log_pressure": [],
        "log_indicator": [],
    }

    nblocks = before.levels.shape[0]
    with h5py.File(before.path, "r") as h0, h5py.File(after.path, "r") as h1:
        for start in range(0, nblocks, chunk_size):
            stop = min(start + chunk_size, nblocks)
            if block_mask is not None:
                local_mask = block_mask[start:stop]
                if not np.any(local_mask):
                    continue
            else:
                local_mask = slice(None)

            i0, rho0, p0, g0, _ = indicator_chunk(before, h0, rf_before, start, stop, x_com)
            i1, rho1, p1, g1, _ = indicator_chunk(after, h1, rf_after, start, stop, x_com)

            rho0 = rho0[local_mask]
            p0 = p0[local_mask]
            g0 = g0[local_mask]
            i0 = i0[local_mask]
            rho1 = rho1[local_mask]
            p1 = p1[local_mask]
            g1 = g1[local_mask]
            i1 = i1[local_mask]

            valid = (rho0 > 0.0) & (p0 > 0.0) & (g0 > 0.0) & (i0 > 0.0) & (p1 > 0.0)
            if not np.any(valid):
                continue
            rrho = rho1[valid] / rho0[valid]
            rp = p1[valid] / p0[valid]
            rg = g1[valid] / g0[valid]
            ri = i1[valid] / i0[valid]
            ratios["rho"].append(rrho.ravel())
            ratios["pressure"].append(rp.ravel())
            ratios["g"].append(rg.ravel())
            ratios["indicator"].append(ri.ravel())
            log_terms["log_rho"].append(np.log(rrho).ravel())
            log_terms["log_g"].append(np.log(rg).ravel())
            log_terms["minus_log_pressure"].append((-np.log(rp)).ravel())
            log_terms["log_indicator"].append(np.log(ri).ravel())

    ratio_summary = {}
    log_summary = {}
    for name, chunks in ratios.items():
        values = np.concatenate(chunks) if chunks else np.array([], dtype=np.float64)
        ratio_summary[name] = quantiles(values)
    for name, chunks in log_terms.items():
        values = np.concatenate(chunks) if chunks else np.array([], dtype=np.float64)
        log_summary[name] = quantiles(values)

    med_abs = {
        "rho": abs(log_summary["log_rho"]["p50"]),
        "g": abs(log_summary["log_g"]["p50"]),
        "pressure": abs(log_summary["minus_log_pressure"]["p50"]),
    }
    dominant = max(med_abs, key=med_abs.get)
    return {
        "ratio_quantiles": ratio_summary,
        "log_contribution_quantiles": log_summary,
        "dominant_by_median_log_abs": dominant,
    }


def descendant_counts(post: DumpMeta) -> Dict[BlockKey, int]:
    counts: Dict[BlockKey, int] = defaultdict(int)
    for level, ix, iy, iz in post.keys:
        for ancestor_level in range(level):
            shift = level - ancestor_level
            counts[(ancestor_level, ix >> shift, iy >> shift, iz >> shift)] += 1
    return counts


def refined_parent_mask(pre: DumpMeta, post: DumpMeta) -> Tuple[np.ndarray, np.ndarray]:
    post_same = set(post.keys)
    desc = descendant_counts(post)
    mask = np.zeros(len(pre.keys), dtype=bool)
    counts = np.zeros(len(pre.keys), dtype=np.int64)
    for i, key in enumerate(pre.keys):
        ndesc = desc.get(key, 0)
        counts[i] = ndesc
        mask[i] = key not in post_same and ndesc > 0
    return mask, counts


def aggregate_post_descendant_max(pre: DumpMeta, post: DumpMeta, post_stats: BlockStats) -> Dict[BlockKey, float]:
    out: Dict[BlockKey, float] = {}
    for post_idx, key in enumerate(post.keys):
        level, ix, iy, iz = key
        value = float(post_stats.max_i[post_idx])
        for ancestor_level in range(level):
            shift = level - ancestor_level
            ancestor = (ancestor_level, ix >> shift, iy >> shift, iz >> shift)
            if ancestor not in out or value > out[ancestor]:
                out[ancestor] = value
    return out


def aggregate_post_descendant_parent_equiv_max(
    pre: DumpMeta,
    post: DumpMeta,
    post_stats: BlockStats,
) -> Dict[BlockKey, float]:
    out: Dict[BlockKey, float] = {}
    for post_idx, key in enumerate(post.keys):
        level, ix, iy, iz = key
        child_value = float(post_stats.max_i[post_idx])
        for ancestor_level in range(level):
            shift = level - ancestor_level
            ancestor = (ancestor_level, ix >> shift, iy >> shift, iz >> shift)
            parent_equiv = math.ldexp(child_value, shift)
            if ancestor not in out or parent_equiv > out[ancestor]:
                out[ancestor] = parent_equiv
    return out


def aggregate_post_descendant_parent_equiv_records(
    post: DumpMeta,
    post_stats: BlockStats,
) -> Dict[BlockKey, PostParentEquivRecord]:
    out: Dict[BlockKey, PostParentEquivRecord] = {}
    for post_idx, key in enumerate(post.keys):
        level, ix, iy, iz = key
        child_value = float(post_stats.max_i[post_idx])
        for ancestor_level in range(level):
            shift = level - ancestor_level
            ancestor = (ancestor_level, ix >> shift, iy >> shift, iz >> shift)
            parent_equiv = math.ldexp(child_value, shift)
            if ancestor not in out or parent_equiv > out[ancestor].value:
                out[ancestor] = PostParentEquivRecord(
                    value=parent_equiv,
                    rho=float(post_stats.max_rho[post_idx]),
                    pressure=float(post_stats.max_pressure[post_idx]),
                    g=float(post_stats.max_g[post_idx]),
                    radius=float(post_stats.max_radius[post_idx]),
                    child_level=level,
                    child_key=key,
                )
    return out


def print_absolute_indicator_summary(label: str, meta: DumpMeta, stats: BlockStats, threshold: float) -> None:
    print(f"\n{label} absolute pressure-scale-height indicator")
    for name, values in (("block max(I)", stats.max_i), ("boundary-cell max(I)", stats.max_boundary_i)):
        q = quantiles(values)
        print(
            f"  {name:20s}: p50={q['p50']:.6e} p90={np.quantile(values, 0.90):.6e} "
            f"p99={q['p99']:.6e} max={q['max']:.6e}  count>thr={int(np.count_nonzero(values > threshold))}"
        )
        by_level = {}
        for level in sorted(set(int(x) for x in meta.levels.tolist())):
            mask = meta.levels == level
            by_level[level] = int(np.count_nonzero(values[mask] > threshold))
        print(f"    count>thr by level: {by_level}")


def write_refined_parent_csv(
    path: Path,
    pre: DumpMeta,
    base_stats: BlockStats,
    pre_stats: BlockStats,
    post_desc_max: Dict[BlockKey, float],
    post_parent_equiv_max: Dict[BlockKey, float],
    post_parent_equiv_records: Dict[BlockKey, PostParentEquivRecord],
    refined_mask: np.ndarray,
    descendant_count: np.ndarray,
    threshold: float,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    widths = block_widths(pre.levels, pre.coord, pre.root_nx)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "pre_block_index",
            "level",
            "ix",
            "iy",
            "iz",
            "xmin",
            "ymin",
            "zmin",
            "block_width",
            "descendant_blocks_in_post",
            "maxI_base",
            "maxI_pre",
            "maxI_pre_over_threshold",
            "maxI_pre_over_base",
            "maxI_post_descendant",
            "maxI_post_parent_equiv",
            "maxI_post_parent_equiv_over_threshold",
            "maxI_post_parent_equiv_over_pre",
            "post_rho_at_parent_equiv_maxI",
            "post_pressure_at_parent_equiv_maxI",
            "post_g_at_parent_equiv_maxI",
            "post_radius_at_parent_equiv_maxI",
            "post_child_level_at_parent_equiv_maxI",
            "rho_post_over_pre_at_maxI",
            "pressure_post_over_pre_at_maxI",
            "g_post_over_pre_at_maxI",
            "rho_at_pre_maxI",
            "pressure_at_pre_maxI",
            "g_at_pre_maxI",
            "radius_at_pre_maxI",
        ])
        order = np.argsort(pre_stats.max_i[refined_mask])[::-1]
        refined_indices = np.flatnonzero(refined_mask)[order]
        for idx in refined_indices:
            key = pre.keys[idx]
            base = float(base_stats.max_i[idx])
            pre_i = float(pre_stats.max_i[idx])
            post_parent_equiv = post_parent_equiv_max.get(key, math.nan)
            post_record = post_parent_equiv_records.get(key)
            if post_record is None:
                post_rho = math.nan
                post_pressure = math.nan
                post_g = math.nan
                post_radius = math.nan
                post_child_level = math.nan
            else:
                post_rho = post_record.rho
                post_pressure = post_record.pressure
                post_g = post_record.g
                post_radius = post_record.radius
                post_child_level = post_record.child_level
            writer.writerow([
                int(idx),
                key[0],
                key[1],
                key[2],
                key[3],
                float(pre.origins[idx, 0]),
                float(pre.origins[idx, 1]),
                float(pre.origins[idx, 2]),
                float(np.max(widths[idx])),
                int(descendant_count[idx]),
                base,
                pre_i,
                pre_i / threshold if threshold > 0.0 else math.nan,
                pre_i / base if base > 0.0 else math.nan,
                post_desc_max.get(key, math.nan),
                post_parent_equiv,
                post_parent_equiv / threshold if threshold > 0.0 else math.nan,
                post_parent_equiv / pre_i if pre_i > 0.0 else math.nan,
                post_rho,
                post_pressure,
                post_g,
                post_radius,
                post_child_level,
                post_rho / float(pre_stats.max_rho[idx]) if pre_stats.max_rho[idx] > 0.0 else math.nan,
                post_pressure / float(pre_stats.max_pressure[idx]) if pre_stats.max_pressure[idx] > 0.0 else math.nan,
                post_g / float(pre_stats.max_g[idx]) if pre_stats.max_g[idx] > 0.0 else math.nan,
                float(pre_stats.max_rho[idx]),
                float(pre_stats.max_pressure[idx]),
                float(pre_stats.max_g[idx]),
                float(pre_stats.max_radius[idx]),
            ])


def print_ratio_table(title: str, summary: Dict[str, object]) -> None:
    print(f"\n{title}")
    print("  ratio quantiles after/before:")
    for name, q in summary["ratio_quantiles"].items():
        print(
            f"    {name:9s}: p50={q['p50']:.6e} p95={q['p95']:.6e} "
            f"p99={q['p99']:.6e} max={q['max']:.6e}"
        )
    logs = summary["log_contribution_quantiles"]
    print("  median log contributions to log(I_after/I_before):")
    print(
        "    log(rho)={:.6e}  log(g)={:.6e}  -log(P)={:.6e}  total={:.6e}".format(
            logs["log_rho"]["p50"],
            logs["log_g"]["p50"],
            logs["minus_log_pressure"]["p50"],
            logs["log_indicator"]["p50"],
        )
    )
    print(f"  dominant by median absolute log term: {summary['dominant_by_median_log_abs']}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dump_dir", type=Path, help="Directory containing dump_XXXXX.h5 files")
    parser.add_argument("--param", type=Path, default=None, help="optional PRJ parfile for max_level/min_dx/threshold")
    parser.add_argument("--base", type=int, default=0, help="baseline dump index")
    parser.add_argument("--pre", type=int, default=3, help="last dump before refinement jump")
    parser.add_argument("--post", type=int, default=4, help="first dump after refinement jump")
    parser.add_argument("--threshold", type=float, default=None, help="pressure_scale_height refine threshold")
    parser.add_argument("--max-level", type=int, default=None, help="configured mesh max_level; default infers active max")
    parser.add_argument("--min-dx", type=float, default=None, help="configured min_dx, if any")
    parser.add_argument("--x-com", type=float, nargs=3, default=(0.0, 0.0, 0.0), help="gravity center of mass")
    parser.add_argument("--chunk-size", type=int, default=256, help="number of blocks to process per HDF5 chunk")
    parser.add_argument("--top", type=int, default=20, help="number of top refined parent blocks to print")
    parser.add_argument("--out-dir", type=Path, default=None, help="directory for CSV/JSON diagnostics")
    args = parser.parse_args()

    dump_dir = args.dump_dir
    out_dir = args.out_dir if args.out_dir is not None else dump_dir / "psh_diagnostics"
    x_com = np.asarray(args.x_com, dtype=np.float64)
    params = parse_param_file(args.param)
    par_max_level = param_int(params, "max_level")
    par_min_dx = param_float(params, "min_dx")
    par_threshold = param_float(params, "amr_refine_thresh")
    par_use_multipole = param_int(params, "use_multipole_gravity")
    requested_max_level = args.max_level
    if requested_max_level is None and par_max_level is not None and par_max_level >= 0:
        requested_max_level = par_max_level
    if requested_max_level is not None and requested_max_level < 0:
        requested_max_level = None
    min_dx = args.min_dx if args.min_dx is not None else par_min_dx
    threshold = args.threshold if args.threshold is not None else (par_threshold if par_threshold is not None else 0.25)

    base = load_meta(dump_path(dump_dir, args.base))
    pre = load_meta(dump_path(dump_dir, args.pre))
    post = load_meta(dump_path(dump_dir, args.post))

    print("Dump overview")
    for meta in (base, pre, post):
        print(
            f"  dump_{meta.dump_index:05d}: time={meta.time:.9e} step={meta.step} "
            f"nblocks={meta.levels.size} levels={level_hist(meta.levels)}"
        )

    gravity_level = effective_gravity_level(pre, requested_max_level, min_dx)
    if requested_max_level is None and min_dx is None:
        print(
            "\nWARNING: gravity grid level was inferred from max active level. "
            "For exact |g| reconstruction, pass the parfile max_level/min_dx."
        )
    print(
        f"Using threshold={threshold:.6e}, gravity radial grid level={gravity_level}, "
        f"x_com={tuple(float(x) for x in x_com)}"
    )
    if par_use_multipole is None or par_use_multipole != 0:
        print(
            "WARNING: dump files store only the radial monopole accel[] table. "
            "If multipole gravity was enabled, the full per-cell |g| used by AMR "
            "cannot be reconstructed exactly from these dumps."
        )

    rf_base = gravity_rf(base, gravity_level)
    rf_pre = gravity_rf(pre, gravity_level)
    rf_post = gravity_rf(post, gravity_level)

    print("\nComputing block pressure-scale-height statistics...")
    base_stats = compute_block_stats(base, rf_base, x_com, args.chunk_size)
    pre_stats = compute_block_stats(pre, rf_pre, x_com, args.chunk_size)
    post_stats = compute_block_stats(post, rf_post, x_com, args.chunk_size)
    print_absolute_indicator_summary(f"dump_{args.base:05d}", base, base_stats, threshold)
    print_absolute_indicator_summary(f"dump_{args.pre:05d}", pre, pre_stats, threshold)
    print_absolute_indicator_summary(f"dump_{args.post:05d}", post, post_stats, threshold)

    refined_mask, desc_count = refined_parent_mask(pre, post)
    nref = int(np.count_nonzero(refined_mask))
    print(f"Refined parent regions from dump_{args.pre:05d} to dump_{args.post:05d}: {nref}")
    if nref > 0:
        vals = pre_stats.max_i[refined_mask]
        print(
            "  pre-dump parent max(I) quantiles: "
            f"p50={np.quantile(vals, 0.50):.6e} "
            f"p90={np.quantile(vals, 0.90):.6e} "
            f"p99={np.quantile(vals, 0.99):.6e} "
            f"max={np.max(vals):.6e}"
        )
        for frac in (1.0, 0.9, 0.75, 0.5):
            count = int(np.count_nonzero(vals >= frac * threshold))
            print(f"  parents with pre max(I) >= {frac:.2f} threshold: {count}/{nref}")

    if same_geometry(base, pre):
        all_summary = compare_same_mesh(base, pre, rf_base, rf_pre, x_com, args.chunk_size)
        print_ratio_table(f"dump_{args.base:05d} -> dump_{args.pre:05d}, all unchanged-geometry cells", all_summary)
        if nref > 0:
            ref_summary = compare_same_mesh(
                base, pre, rf_base, rf_pre, x_com, args.chunk_size, block_mask=refined_mask
            )
            print_ratio_table(
                f"dump_{args.base:05d} -> dump_{args.pre:05d}, only regions refined by dump_{args.post:05d}",
                ref_summary,
            )
    else:
        all_summary = None
        ref_summary = None
        print(f"\nSkipping same-mesh ratio analysis: dump_{args.base:05d} and dump_{args.pre:05d} differ.")

    post_desc_max = aggregate_post_descendant_max(pre, post, post_stats)
    post_parent_equiv_max = aggregate_post_descendant_parent_equiv_max(pre, post, post_stats)
    post_parent_equiv_records = aggregate_post_descendant_parent_equiv_records(post, post_stats)
    if nref > 0:
        post_parent_equiv_values = np.array(
            [post_parent_equiv_max.get(pre.keys[idx], math.nan) for idx in np.flatnonzero(refined_mask)],
            dtype=np.float64,
        )
        finite = post_parent_equiv_values[np.isfinite(post_parent_equiv_values)]
        if finite.size > 0:
            print(
                "\nDump-post descendant values rescaled to their pre-dump parent dx: "
                f"p50={np.quantile(finite, 0.50):.6e} "
                f"p90={np.quantile(finite, 0.90):.6e} "
                f"p99={np.quantile(finite, 0.99):.6e} "
                f"max={np.max(finite):.6e} "
                f"count>thr={int(np.count_nonzero(finite > threshold))}/{finite.size}"
            )
    csv_path = out_dir / f"refined_parents_{args.pre:05d}_to_{args.post:05d}.csv"
    write_refined_parent_csv(
        csv_path, pre, base_stats, pre_stats, post_desc_max, post_parent_equiv_max,
        post_parent_equiv_records,
        refined_mask, desc_count, threshold
    )
    print(f"\nWrote refined parent table: {csv_path}")

    if nref > 0:
        print(f"\nTop {min(args.top, nref)} refined parent candidates by pre-dump max(I):")
        order = np.argsort(pre_stats.max_i[refined_mask])[::-1]
        refined_indices = np.flatnonzero(refined_mask)[order[: args.top]]
        for rank, idx in enumerate(refined_indices, start=1):
            key = pre.keys[idx]
            print(
                f"  {rank:2d}. block={idx:6d} key={key} desc={desc_count[idx]:4d} "
                f"Ipre={pre_stats.max_i[idx]:.6e} ({pre_stats.max_i[idx] / threshold:.3f} thr) "
                f"Ipost_parent_equiv={post_parent_equiv_max.get(key, math.nan):.6e} "
                f"I0={base_stats.max_i[idx]:.6e} r={pre_stats.max_radius[idx]:.6e} "
                f"rho={pre_stats.max_rho[idx]:.6e} P={pre_stats.max_pressure[idx]:.6e} "
                f"g={pre_stats.max_g[idx]:.6e}"
            )

    summary_path = out_dir / "summary.json"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "base": {"index": args.base, "time": base.time, "step": base.step, "nblocks": int(base.levels.size), "levels": level_hist(base.levels)},
        "pre": {"index": args.pre, "time": pre.time, "step": pre.step, "nblocks": int(pre.levels.size), "levels": level_hist(pre.levels)},
        "post": {"index": args.post, "time": post.time, "step": post.step, "nblocks": int(post.levels.size), "levels": level_hist(post.levels)},
        "threshold": threshold,
        "gravity_level_used": gravity_level,
        "refined_parent_count": nref,
        "base_max_indicator": quantiles(base_stats.max_i),
        "pre_max_indicator": quantiles(pre_stats.max_i),
        "post_max_indicator": quantiles(post_stats.max_i),
        "base_boundary_max_indicator": quantiles(base_stats.max_boundary_i),
        "pre_boundary_max_indicator": quantiles(pre_stats.max_boundary_i),
        "post_boundary_max_indicator": quantiles(post_stats.max_boundary_i),
        "all_same_mesh_summary": all_summary,
        "refined_parent_same_mesh_summary": ref_summary,
    }
    with summary_path.open("w") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
    print(f"Wrote summary JSON: {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
